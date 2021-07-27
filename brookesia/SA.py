"""
    Brookesia
    Reduction and optimization of kinetic mechanisms

    Copyright (C) 2019  Matynia, Delaroque, Chakravarty
    contact : alexis.matynia@sorbonne-universite.fr

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""


import cantera as ct
import numpy as np
import brookesia.Class_def as cdef
from  brookesia.Class_def import print_
import pandas as pd
import copy
import sys
import os
import time as timer



def sensitivities_computation_SA(red_data, mech_data,red_results):

    conditions = red_results.conditions
    mp = conditions.main_path
    print_('\n' + conditions.config + " sensitivity analysis", mp)

#    gas_ref = red_data.gas_ref
    n_sp_ref  = red_data.gas_ref.n_species
    n_r_ref   = red_data.gas_ref.n_reactions
#    n_sp_red  = red_data.gas_red.n_species
    tsp_idx   = red_data.targetSpeciesIdx
    n_tsp     = len(tsp_idx)
    tsp_name  = red_data.tspc

    gas_ref = conditions.composition.gas_ref
    gas_red = red_results.gas

    X_red = conditions.composition.X
    gas_red.TPX = conditions.state_var.T, conditions.state_var.P,X_red

    pts_scatter = red_results.pts_scatter
    n_points = len(pts_scatter)
    n_points_SA = int(red_data.red_op.n_points)

    nu_f = gas_ref.reactant_stoich_coeffs()
    nu_r = gas_ref.product_stoich_coeffs()
    nu   = nu_f-nu_r

    # Sensitivity coeff matrix
    if 'SARGEP' in red_data.reduction_operator:
	    S_react        = np.zeros((n_sp_ref,n_r_ref))
	    S_react_x      = np.zeros((n_points_SA,n_sp_ref,n_r_ref))
	    norm_S_react_x = np.zeros((n_points_SA,n_sp_ref,n_r_ref))
	    S_AB_tsp       = np.zeros((n_sp_ref,n_sp_ref))
	    S_AB_tsp_z     = np.zeros((n_points_SA,n_sp_ref))
	    S_AB_tsp_z_r   = np.zeros((n_points_SA,n_r_ref,n_sp_ref,n_sp_ref))
	    sensi_T_z      = np.zeros((n_points_SA,n_r_ref))
	    norm_sensi_T_z = np.zeros((n_points_SA,n_r_ref))
    else:
	    S_react        = np.zeros((n_tsp,n_r_ref))
	    S_react_x      = np.zeros((n_points_SA,n_tsp,n_r_ref))
	    norm_S_react_x = np.zeros((n_points_SA,n_tsp,n_r_ref))
	    S_AB_tsp       = np.zeros((n_tsp,n_sp_ref))
	    S_AB_tsp_z     = np.zeros((n_points_SA,n_sp_ref))
	    S_AB_tsp_z_r   = np.zeros((n_points_SA,n_r_ref,n_tsp,n_sp_ref))
	    sensi_T_z      = np.zeros((n_points_SA,n_r_ref))
	    norm_sensi_T_z = np.zeros((n_points_SA,n_r_ref))
    sensi_T        = np.zeros((n_r_ref))
    sensi_scatter  = []

    time_start = timer.time()

    if 'flame' in conditions.config:

        f = red_results.f

        # Sensitivity calculation
        if 'free' in conditions.config:
            sensi_Sl = f.get_flame_speed_reaction_sensitivities()
            red_data.red_op.sensi_Sl = sensi_Sl

        if red_data.sens_method == 'adjoint':
            try:
                time_1 = timer.time() ; T = f.T

                def get_species_reaction_sensitivities(f, tsp, grid_point):
                    r"""
                    Compute the normalized sensitivities of the species production
                    :math:`s_{i, spec}` with respect to the reaction rate constants :math:`k_i`:
                    .. math::
                    s_{i, spec} = \frac{k_i}{[X]} \frac{d[X]}{dk_i}
                    """
                    def g(sim):
                        if tsp!='T':
                            return sim.X[f.gas.species_index(tsp), grid_point]
                        else:
                            return sim.T[grid_point]

                    Nvars = sum(D.n_components * D.n_points for D in f.domains)

                    # Index of spec in the global solution vector
                    try:
                        i_spec = f.inlet.n_components + f.flame.component_index(tsp) + f.domains[1].n_components*grid_point
                    except:
                        i_spec = f.fuel_inlet.n_components + f.flame.component_index(tsp) + f.domains[1].n_components*grid_point

                    dgdx = np.zeros(Nvars)
                    dgdx[i_spec] = 1

                    spec_0 = g(f)

                    def perturb(sim, i, dp):
                        sim.gas.set_multiplier(1+dp, i)

                    S_react_x_zi = f.solve_adjoint(perturb, f.gas.n_reactions, dgdx,g,1e-5) / spec_0

                    return S_react_x_zi


                if conditions.simul_param.verbose >=2 :
                    title = 'Sensitivity computation (adjoint)  '
                    bar = cdef.ProgressBar(n_points, title)

                z_i=0
                for z in range(n_points):
                    if z%int(max(n_points/red_data.red_op.n_points,1))==0    \
                    and 0.01*(max(red_results.T)-min(red_results.T)) \
                                        <T[z]-T[0]<           \
                    0.99*(max(red_results.T)-min(red_results.T)) :            # ~ (Ti+1%)<T_pert<(Tf-1%)
                        sensi_scatter.append(z)
                        if conditions.simul_param.verbose >=2: bar.update(z)
                        sensi_T_z_adj = get_species_reaction_sensitivities(f, 'T', z)

                        if 'SARGEP' in red_data.reduction_operator:
                            n_spc = n_sp_ref
                        else:
                            n_spc = len(tsp_name)
                        S_react_x_adj = [[]]*n_spc
                        for spA in range(n_spc):
                            if mech_data.spec.activ_m[spA] :
                                if 'SARGEP' in red_data.reduction_operator:
                                    S_react_x_adj[spA]=get_species_reaction_sensitivities(f, gas_ref.species_name(spA), z)
                                else:
                                    S_react_x_adj[spA]=get_species_reaction_sensitivities(f, tsp_name[spA], z)
                                r_red = 0
                                for r in range(n_r_ref):
                                    if mech_data.react.activ_m[r]:
                                        sensi_T_z[z_i][r] = sensi_T_z_adj[r_red]
                                        S_react_x[z_i][spA][r] = S_react_x_adj[spA][r_red]
                                        r_red+=1
                                    if '_sp' in red_data.reduction_operator:
                                        for r in range(n_r_ref):
                                            if mech_data.react.activ_m[r]:
                                                # Inter-species sensitivity calculation:
                                                for spB in range(n_sp_ref):
                                                    if (mech_data.spec.activ_m[spB] and nu[spB,r]!=0):
                                                        # collecting sensitivities of reactions involving both species
                                                        S_AB_tsp_z_r[z_i][r][spA,spB]=S_react_x[z_i][spA][r]
                        z_i+=1
                if conditions.simul_param.verbose >=2: bar.update(n_points)

                time_2 = timer.time()
                if conditions.simul_param.verbose >=4 :
                    print_("\n      time for adjoint method SA computation: "+str(round(time_2-time_1))+'s',mp)
            except:
                # if problem during the adjoint method, go for the brute-force method
                print_('Warning: error during the sensitivity analysis with adjoint method',mp)
                print_('         brute force method used instead\n',mp)
                red_data.sens_method = 'brute_force'
                # reinitialisation of the sensitivity coeff matrix
                if 'SARGEP' in red_data.reduction_operator:
            	    S_react        = np.zeros((n_sp_ref,n_r_ref))
            	    S_react_x      = np.zeros((n_points_SA,n_sp_ref,n_r_ref))
            	    norm_S_react_x = np.zeros((n_points_SA,n_sp_ref,n_r_ref))
            	    S_AB_tsp       = np.zeros((n_sp_ref,n_sp_ref))
            	    S_AB_tsp_z     = np.zeros((n_points_SA,n_sp_ref))
            	    S_AB_tsp_z_r   = np.zeros((n_points_SA,n_r_ref,n_sp_ref,n_sp_ref))
            	    sensi_T_z      = np.zeros((n_points_SA,n_r_ref))
            	    norm_sensi_T_z = np.zeros((n_points_SA,n_r_ref))
                else:
            	    S_react        = np.zeros((n_tsp,n_r_ref))
            	    S_react_x      = np.zeros((n_points_SA,n_tsp,n_r_ref))
            	    norm_S_react_x = np.zeros((n_points_SA,n_tsp,n_r_ref))
            	    S_AB_tsp       = np.zeros((n_tsp,n_sp_ref))
            	    S_AB_tsp_z     = np.zeros((n_points_SA,n_sp_ref))
            	    S_AB_tsp_z_r   = np.zeros((n_points_SA,n_r_ref,n_tsp,n_sp_ref))
            	    sensi_T_z      = np.zeros((n_points_SA,n_r_ref))
            	    norm_sensi_T_z = np.zeros((n_points_SA,n_r_ref))
                sensi_T        = np.zeros((n_r_ref))
                sensi_scatter  = []

#        red_data.sens_method = 'brute_force'
        if red_data.sens_method == 'brute_force':

            time_1 = timer.time()
            r_red=-1
            title = 'Sensitivity computation (brute force) '
            bar = cdef.ProgressBar(n_r_ref,title)

            for r in range(n_r_ref):
                if mech_data.react.activ_m[r]:
                    gas_red.set_multiplier(1) # initialisation
                    r_red+=1
                    dk = 5e-2 #perturbation
                    gas_red.set_multiplier(1+dk, r_red) # set the multiplier of kf for r reaction


                    os.chdir(conditions.main_path)
                    os.chdir('Flame_ref_results')

                    fn = str(conditions.num)
                    if 'free' in conditions.config or 'burner' in conditions.config:
                        fn += 'ff_'
                    elif 'tp_' in conditions.config:
                        fn =+ 'tp_'
                    else:
                        fn += 'cf_'

                    if 'diff' in conditions.config: phi = 'diff'
                    else: phi = '%.2f' %conditions.composition.phi

                    if conditions.state_var.P>10000:
                        fn+=conditions.composition.fuel.replace('/','').split('(')[0]\
                        +'_'+phi\
                        +'_'+'%.0f'%conditions.state_var.T+'_'+'%.2f'%(conditions.state_var.P/1e5)\
                        +'.xml'
                    else:
                        fn+=conditions.composition.fuel.replace('/','').split('(')[0]\
                        +'_'+phi\
                        +'_'+'%.0f'%conditions.state_var.T+'_'+'%.0f'%(conditions.state_var.P)\
                        +'.xml'

                    #=================================================
                    # supress console output during the simulation
                    if conditions.simul_param.verbose < 9:
                        old_stdout = sys.stdout ; old_stderr = sys.stderr
                        with open(os.devnull, "w") as devnull: sys.stdout = devnull ; sys.stderr = devnull
                    f.restore(fn, 'ref_solution')
                    # Initialize and solve
                    try:
                        f.solve(auto = False, loglevel = 0, refine_grid = False)
                    except:
                        if conditions.simul_param.verbose >= 3:
                            print_("\n     WARNING: No solution found\n",mp)
                    if conditions.simul_param.verbose < 9:
                        sys.stdout = old_stdout ; sys.stderr = old_stderr
                    #=================================================

                    T = f.T
                    z_i=0
                    for z in range(n_points):
                        if z%int(max(n_points/red_data.red_op.n_points,1))==0    \
                        and 0.01*(max(red_results.T)-min(red_results.T)) \
                                                <T[z]-T[0]<           \
                            0.99*(max(red_results.T)-min(red_results.T)) :            # ~ (Ti+1%)<T_pert<(Tf-1%)
                            if sensi_scatter==[]:
                                sensi_scatter.append(z)
                            elif z<np.max(sensi_scatter):
                                sensi_scatter.append(z)
                            f.set_gas_state(z)
                            conc_red = red_results.conc[z]
                            kf_ref = red_results.kf[z][r]
                            kf_pert = gas_red.forward_rate_constants[r_red]
                            if kf_pert!=kf_ref:
                                conc_pert = gas_red.concentrations
                                #Temperature and species sensibilities to reactions
                                sensi_T_z[z_i][r] = (kf_ref/red_results.T[z])\
                                        *((T[z]-red_results.T[z])/(kf_pert-kf_ref))
                                tsp=-1
                                for spA in range(n_sp_ref):
                                    if mech_data.spec.activ_m[spA] :
                                        if 'SARGEP' in red_data.reduction_operator:
            		                            ired = gas_red.species_index(gas_ref.species_name(spA))
            		                            if conc_pert[ired]>conditions.simul_param.tol_ss[1]: # tol_ss_flame[1] : atol for steady-state problem
            		                                S_react_x[z_i][spA][r]=kf_ref/conc_red[spA]\
            		                                  *((np.max([conc_pert[ired],0])-conc_red[spA])\
            		                                     /(kf_pert-kf_ref))
            		                                # Inter-species sensitivity calculation:
            		                                for spB in range(n_sp_ref):
            		                                    if (mech_data.spec.activ_m[spB] and nu[spB,r]!=0):
            		                                        # collecting sensitivities of reactions involving both species
            		                                        S_AB_tsp_z_r[z_i][r][spA,spB]=S_react_x[z_i][spA][r]
                                        else:
                                            if spA in tsp_idx:
                                                tsp+=1
                                                ired = gas_red.species_index(tsp_name[tsp])
                                                if conc_pert[ired]>conditions.simul_param.tol_ss[1]: # tol_ss_flame[1] : atol for steady-state problem
                                                    S_react_x[z_i][tsp][r]=kf_ref/conc_red[spA]\
                                                      *((np.max([conc_pert[ired],0])-conc_red[spA])\
                                                         /(kf_pert-kf_ref))
                                                    # Inter-species sensitivity calculation:
                                                    for spB in range(n_sp_ref):
                                                        if (mech_data.spec.activ_m[spB] and nu[spB,r]!=0):
                                                            # collecting sensitivities of reactions involving both species
                                                            S_AB_tsp_z_r[z_i][r][tsp,spB]=S_react_x[z_i][tsp][r]
                            z_i+=1
                bar.update(r)

            time_2 = timer.time()
            if conditions.simul_param.verbose >=4 :
                print_("\n      time for brute-force SA computation: "+str(round(time_2-time_1))+'s',mp)



#       Sensitivity analysis
        title = '2- Analysis of the sensitivity coeff   '


        if 'SARGEP' in red_data.reduction_operator:
            n_spc = n_sp_ref
            if conditions.simul_param.verbose>=2:
                bar = cdef.ProgressBar(mech_data.spec.activ_m.count(True),title)
        else :
            n_spc = n_tsp
            if conditions.simul_param.verbose>=2:
                bar = cdef.ProgressBar(n_tsp*n_sp_ref,title)
        sensi_T_done = False ; sp_c = 0
        for spA in range(n_spc):
            if mech_data.spec.activ_m[spA]:
                sensi_tsp_done = False ; sp_c+=1
                S_AB_tsp_z = np.zeros((n_points_SA,n_sp_ref)) # initialisation of the matrix of interactions between A and B at z
                # 1 - keep only the max sens coeffs of all x points
                for spB in range(n_sp_ref):
                    if mech_data.spec.activ_m[spB] :
                        z_i = 0
                        for z in range(n_points):
                            if z%int(n_points/red_data.red_op.n_points)==0\
                            and 0.01*abs(max(red_results.T)-min(red_results.T)) \
                                                    <abs(T[z]-T[0])<           \
                                0.99*abs(max(red_results.T)-min(red_results.T)) :        # ~ (Ti+1%)<T_pert<(Tf-1%)

                                max_sensi_T_z   = max(abs(sensi_T_z[z_i]))
                                max_sensi_tsp_z = max(abs(S_react_x[z_i][spA]))
                                for r in range(n_r_ref):
                                    S_AB_tsp_z[z_i][spB] += S_AB_tsp_z_r[z_i][r][spA,spB]
                                    if not sensi_T_done: # to do only at the first loop
                                        # Sensi_T normalisation on z
                                        norm_sensi_T_z[z_i][r] = sensi_T_z[z_i][r]/max_sensi_T_z
                                        # keep the max sens coeffs of all z points
                                        if abs(norm_sensi_T_z[z_i][r])>abs(sensi_T[r]):
                                            sensi_T[r] = norm_sensi_T_z[z_i][r]

                                    if not sensi_tsp_done and max_sensi_tsp_z>1e-10: # to do only at the first spB loop
                                        # Sensi_react normalisation on z
                                        norm_S_react_x[z_i][spA][r] = S_react_x[z_i][spA][r]/max_sensi_tsp_z
                                        # keep the max sens coeffs of all z points
                                        if abs(norm_S_react_x[z_i][spA][r])>abs(S_react[spA][r]):
                                            S_react[spA][r] = norm_S_react_x[z_i][spA][r]
                                z_i += 1

                        sensi_T_done   = True
                        sensi_tsp_done = True
                    if conditions.simul_param.verbose>=2:
                        if 'SARGEP' in red_data.reduction_operator: bar.update(sp_c)
                        else: bar.update(spA*n_sp_ref+spB)

                # AB sensitivity normalisation on z
                z_i = 0
                for z in range(n_points):
                    if z%int(n_points/red_data.red_op.n_points)==0\
                    and 0.01*abs(max(red_results.T)-min(red_results.T)) \
                                            <abs(T[z]-T[0])<           \
                        0.99*abs(max(red_results.T)-min(red_results.T)) :        # ~ (Ti+1%)<T_pert<(Tf-1%)
                        S_AB_tsp_z_max = max(abs(S_AB_tsp_z[z_i]))
                        if S_AB_tsp_z_max>1e-10:
                            for spB in range(n_sp_ref):
                                S_AB_tsp_z[z_i][spB] = S_AB_tsp_z[z_i][spB]/S_AB_tsp_z_max
                                # keep only the max sens coeffs of all z points
                                if S_AB_tsp_z[z_i][spB] > S_AB_tsp[spA,spB]:
                                    S_AB_tsp[spA,spB] = S_AB_tsp_z[z_i][spB]
                        z_i += 1

        if conditions.simul_param.verbose>=2 and 'SARGEP' not in red_data.reduction_operator:
            bar.update(n_tsp*n_sp_ref,title)
        red_data.red_op.sensi_T = sensi_T
        del norm_S_react_x ; del norm_sensi_T_z ; del S_AB_tsp_z


    elif 'reactor' in conditions.config:

        if conditions.config == "reactor_UV":
            reactor = ct.IdealGasReactor(gas_red)
        elif conditions.config == "reactor_HP":
            reactor = ct.IdealGasConstPressureReactor(gas_red)
        sim = ct.ReactorNet([reactor])

        # set the tolerances for the solution and for the sensitivity coefficients
        sim.rtol = conditions.simul_param.tol_ts[0]
        sim.atol = conditions.simul_param.tol_ts[1]
        if red_data.red_op.tol_ts[0]:
            sim.rtol_sensitivity = red_data.red_op.tol_ts[0]
            sim.atol_sensitivity = red_data.red_op.tol_ts[1]
        else:
            sim.rtol_sensitivity = conditions.simul_param.tol_ts[0]
            sim.atol_sensitivity = conditions.simul_param.tol_ts[1]

        for r in range(gas_red.n_reactions):
            reactor.add_sensitivity_reaction(r)

        # Sensitivity coeff matrix
        if 'SARGEP' in red_data.reduction_operator:
            S_react    = np.zeros((n_sp_ref,n_r_ref))
            S_AB_tsp   = np.zeros((n_sp_ref,n_sp_ref))
            S_AB_tsp_t = np.zeros((n_points_SA,n_sp_ref,n_sp_ref))
        else:
            S_react    = np.zeros((n_tsp,n_r_ref))
            S_AB_tsp   = np.zeros((n_tsp,n_sp_ref))
            S_AB_tsp_t = np.zeros((n_points_SA,n_tsp,n_sp_ref))


        # fuel1 index
        ind_fuel  = gas_red.species_index(conditions.composition.fuel.split('/')[0].split('(')[0])
        fuel_conc = [gas_red.concentrations[ind_fuel]]

        bar = cdef.ProgressBar(n_points-int(n_points/red_data.red_op.n_points-1))
        t_i = -1
        for t in range(n_points-int(n_points/red_data.red_op.n_points+1)):
#            try:
            sim.advance(pts_scatter[t])
            #normalized sensitivity for reactions
            if t!=0 and t%int(max(n_points/red_data.red_op.n_points,1))==0 :
                t_i+=1
                sensi_scatter.append(t)
                fuel_conc.append(gas_red.concentrations[ind_fuel])
                sensi_r_t  = sim.sensitivities()
    #             IdealGasReactor
    #                    0 - mass
    #                    1 - volume
    #                    2 - internal energy or temperature
    #                    3+ - mass fractions of the species
    #             CstPressureReactor
    #                    0 - mass
    #                    1 - enthalpy or temperature
    #                    2+ - mass fractions of the species

                #normalized sensitivity for species
    #            S_AB_t    = np.zeros((n_sp,n_sp))

                sp_red=-1 ; tsp=-1
                if 'SARGEP' in red_data.reduction_operator:
                    n_spc = n_sp_ref
                else:
                    n_spc = n_tsp
                for spA in range(n_spc):
                    if 'SARGEP' not in red_data.reduction_operator or mech_data.spec.activ_m[spA]:
                        if 'SARGEP' not in red_data.reduction_operator:
                            sp_red = gas_red.species_index(tsp_name[spA])
                        else:
                            sp_red+=1
                        r_red=-1
                        for r in range(n_r_ref):
                            if mech_data.react.activ_m[r]:
                                r_red+=1
                                # Reaction sensitivity calculation:
                                if conditions.config == "reactor_UV":
                                    try:
                                        S_react_x[t_i][spA][r]=sensi_r_t[sp_red+3,r_red]
                                    except:
                                        print('toto')
                                elif conditions.config == "reactor_HP":
                                    S_react_x[t_i][spA][r]=sensi_r_t[sp_red+2,r_red]
                                for spB in range(n_sp_ref):
                                    if mech_data.spec.activ_m[spB] and nu[spB, r] != 0:
                                        # Inter-species sensitivity calculation:
                                        # (sum of the sensitivities of reactions involving both species)
                                        if conditions.config == "reactor_UV":
                                            S_AB_tsp_t[t_i][spA,spB]+=abs(sensi_r_t[sp_red+3, r_red])
                                        elif conditions.config == "reactor_HP":
                                            S_AB_tsp_t[t_i][spA,spB]+=abs(sensi_r_t[sp_red+2, r_red])

                        # S_react_x sensitivity normalisation on t
                        max_S_react_x = max(abs(S_react_x[t_i][spA]))
                        if max_S_react_x>0:
                            for r in range(n_r_ref):
                                S_react_x[t_i][spA][r]=S_react_x[t_i][spA][r]/max_S_react_x
                        # S_AB_t sensitivity normalisation on t
                        max_S_AB_tsp_t = max(abs(S_AB_tsp_t[t_i][spA]))
                        if max_S_AB_tsp_t>0:
                            for spB in range(n_sp_ref):
                                if mech_data.spec.activ_m[spB] and nu[spB, r] != 0:
                                    S_AB_tsp_t[t_i][spA,spB]=S_AB_tsp_t[t_i][spA,spB]/max_S_AB_tsp_t
                                # keep only the max sens coeffs in all time steps
                                if S_AB_tsp_t[t_i][spA,spB] > S_AB_tsp[spA,spB]:
                                    S_AB_tsp[spA,spB] = S_AB_tsp_t[t_i][spA,spB]

                if conditions.error_param.T_check:
                    # temperature sensitivity
                    r_red=-1
                    for r in range(n_r_ref):
                        if mech_data.react.activ_m[r]:
                            r_red+=1
                            # Reaction sensitivity calculation:
                            if conditions.config == "reactor_UV":
                                sensi_T_z[t_i][r]=sensi_r_t[2,r_red]
                            elif conditions.config == "reactor_HP":
                                sensi_T_z[t_i][r]=sensi_r_t[1,r_red]

                    # S_react_x sensitivity normalisation on t
                    max_sensi_T = max(abs(sensi_T_z[t_i]))
                    if max_sensi_T>0:
                        for r in range(n_r_ref):
                            sensi_T_z[t_i][r]=sensi_T_z[t_i][r]/max_sensi_T
#            except:
#                print_('Warning: error during in sensitivities calculation at t = '+'%.3f'%(pts_scatter[t]*1000)+' ms',mp)
            bar.update(t)


        # Sensitivity analysis
        # keep the maximal sensitivities normalized among each timestep
        S_react_filled = False # flag to avoid S_react matrix unfilled when time iterations and fuel reactivity are decorelated

        for spA in range(n_spc):
            for r in range(n_r_ref):
                t_i=-1
                for t in range(n_points-int(n_points/red_data.red_op.n_points-1)):
                    #normalized sensitivity for reactions
                    if t!=0 and t%int(n_points/red_data.red_op.n_points)==0:
                        t_i+=1
                        try: # ... if error during in sensitivities calculation
                            if 0.01*abs(fuel_conc[0]-fuel_conc[-1])     \
                                   <abs(fuel_conc[0]-fuel_conc[t_i])<    \
                                   0.99*abs(fuel_conc[0]-fuel_conc[-1]):
                                if abs(S_react_x[t_i][spA][r])>abs(S_react[spA][r]):    # ~ (C_fuel_ti+1%)<C_fuel<(C_fuel_tf-1%)
                                        S_react[spA][r]=S_react_x[t_i][spA][r]
                                        S_react_filled = True
                                if conditions.error_param.T_check:
                                    if abs(sensi_T_z[t_i][r])>abs(sensi_T[r]):    # ~ (C_fuel_ti+1%)<C_fuel<(C_fuel_tf-1%)
                                            sensi_T[spA][r]=sensi_T_z[t_i][r]
                        except:
                            a=2
        if S_react_filled == False:
            for spA in range(n_spc):
                for r in range(n_r_ref):
                    t_i=-1
                    for t in range(n_points-int(n_points/red_data.red_op.n_points-1)):
                        #normalized sensitivity for reactions
                        if t!=0 and t%int(n_points/red_data.red_op.n_points)==0:
                            t_i+=1
                            try: # ... if error during in sensitivities calculation
                                if abs(S_react_x[t_i][spA][r])>abs(S_react[spA][r]):    # ~ (C_fuel_ti+1%)<C_fuel<(C_fuel_tf-1%)
                                        S_react[spA][r]=S_react_x[t_i][spA][r]
                            except:
                                a=2

        red_data.red_op.sensi_T = sensi_T


    time_end = timer.time()
    if conditions.simul_param.verbose >=4 :
        print_("\n      time for SA computation+analysis: "+str(round(time_end-time_start))+'s',mp)

    print_("\n",mp)



    red_data.red_op.sensi_sp = S_AB_tsp
    red_data.red_op.sensi_r  = S_react

    if red_data.write_results:
        if 0 in sensi_scatter:
            sensi_scatter.remove(0)
        write_sensitivities(red_data,sensi_scatter,S_react_x,sensi_T_z,S_AB_tsp_z_r,red_results)
        return sensi_scatter,S_react_x,sensi_T_z,red_data
    else:
        return red_data



def speciesWithdrawal(conditions, red_data, red_method, mech_data, eps):
        #gas, fuel, diluant, sensi_species, eps, target_species, active_species_main, verbose = 3) :

    mp = conditions.main_path

    # main variables
    gas  = red_data.gas_ref


    #sensi_AB = red_data.red_op.sensi_sp
    sensi_AB = copy.deepcopy(red_data.red_op.sensi_sp)
    for sp in range(len(sensi_AB)):
        sensi_AB[sp]=sensi_AB[sp]/max(sensi_AB[sp])



    tsp_idx = red_data.targetSpeciesIdx
    verbose=conditions.simul_param.verbose

    ns = gas.n_species
#    nr = gas.n_reactions
    active_species = list(mech_data.spec.activ_p)


    # Add fuel / ox /diluent to the conserved species
    init_spec = conditions.composition.X
    init_spec = init_spec.split(',')
    for spec in init_spec:
        ind_spec = gas.species_index(spec.split(":")[0].replace(' ',''))
        if not active_species[ind_spec]: active_species[ind_spec]=True
    if conditions.composition.X2:
        init_spec = conditions.composition.X2
        init_spec = init_spec.split(',')
        for spec in init_spec:
            ind_spec = gas.species_index(spec.split(":")[0].replace(' ',''))
            if not active_species[ind_spec]: active_species[ind_spec]=True


    if red_method == 'SARGEP_sp':

    # Dijkstra algorithm
        tspi =-1
        for tsp in tsp_idx:
                tspi +=1
                mpq_list = []

                # max-priority queue (mpq) construction
                OIC = []
                for spA in range(ns):
                    OIC.append(sensi_AB[tsp, spA])
                    if OIC[spA] > eps[tspi]:
                        mpq_list.append(spA)

                for mpq in range(len(mpq_list)):
    #                R = np.array((ns,ns))
                    R_max=[] ; visited_sp=[]
                    for sp in range(ns):
                        R_max.append(OIC[mpq]*sensi_AB[mpq,sp])
                        visited_sp.append([tspi,sp])
                        if OIC[mpq]*sensi_AB[mpq,sp] > eps[tspi]:
                            active_species[sp] = True
                    R_max[tspi]=0 ; R_max[mpq]=0


                    # graph search
                    idx_node = R_max.index(max(R_max))
                    while max(R_max)>eps[tspi]:
                        # node informations
                        visited_sp_node = visited_sp[idx_node]

                        # interaction coefficients of the step
                        Ri = []
                        for sp in range(ns):
                            if sp not in visited_sp_node:
                                Ri.append(R_max[idx_node]*sensi_AB[idx_node,sp])
                                if Ri[-1]>eps[tspi]:
                                    active_species[sp] = True
                            else:
                                Ri.append(0)

                        # record informations on the present step of graph analysis
                        R_max[idx_node]=max(Ri)
                        visited_sp[idx_node].append(Ri.index(max(Ri)))
                        if max(Ri)>max(R_max):
                            idx_node = Ri.index(max(Ri))
                        else:
                            idx_node = R_max.index(max(R_max))

        ## Species sensitivities based withdrawal
        #for t in range(len(tsp_idx)):
            #spA=tsp_idx[t]
            #if not active_species[spA]: active_species[spA]=True
            #for spB in range(ns):
                #if not mech_data.spec.activ_m[spB]:
                    #active_species[spB] = False
                #else:
                    #if sensi_AB[t][spB]>eps[t] and not active_species[spB]:
                        #active_species[spB]=True

    if '_sp' in red_method:
        # Species sensitivities based withdrawal
        for t in range(len(tsp_idx)):
            spA=tsp_idx[t]
            if not active_species[spA]: active_species[spA]=True
            for spB in range(ns):
                if not mech_data.spec.activ_m[spB]:
                    active_species[spB] = False
                else:
                    if sensi_AB[t][spB]>eps[t] and not active_species[spB]:
                        active_species[spB]=True

    return active_species







def reactionWithdrawal(conditions, mech_data,active_species,red_data,red_method,eps_r):
        #gas, sensi_reactions, active_species , eps, target_species, active_reactions_main, verbose=3) :

    # main variables
    gas_ref  = red_data.gas_ref

    sensi_r = red_data.red_op.sensi_r
#    target_species = red_data.targetSpeciesIdx

    ns = gas_ref.n_species
    nr = gas_ref.n_reactions

    tsp_idx = red_data.targetSpeciesIdx

    nu_f =gas_ref.reactant_stoich_coeffs()
    nu_r =gas_ref.product_stoich_coeffs()

    active_reactions = list(mech_data.react.activ_p)

    # Add fuel / ox /diluent to the conserved species
    init_spec = conditions.composition.X
    init_spec = init_spec.split(',')
    for spec in init_spec:
        ind_spec = gas_ref.species_index(spec.split(":")[0].replace(' ',''))
        if not active_species[ind_spec]: active_species[ind_spec]=True
    if conditions.composition.X2:
        init_spec = conditions.composition.X2
        init_spec = init_spec.split(',')
        for spec in init_spec:
            ind_spec = gas_ref.species_index(spec.split(":")[0].replace(' ',''))
            if not active_species[ind_spec]: active_species[ind_spec]=True

    # define new activated species / reactions based on species sensitivities
    if '_sp' in red_method:
        active_reactions = [True]*len(mech_data.react.activ_p)
        for r in range(nr):
            if not mech_data.react.activ_m[r]:
                active_reactions[r] = False
            else:
                for sp in range(ns):
                    # remove reactions involving non active species
                    if not active_species[sp]:
                        if nu_f[sp][r]!=0 or nu_r[sp][r]!=0:
                            active_reactions[r] = False
                            break


    # define new activated species / reactions based on reactions sensitivities
    if '_r' in red_method:
        for tsp in range(len(tsp_idx)):
            spA=tsp_idx[tsp]
            if not active_species[spA]: active_species[spA]=True
            sensi_sorted = list(sensi_r[tsp])
            sensi_sorted.sort()
            while 0 in sensi_sorted: sensi_sorted.remove(0)
            lim_val = sensi_sorted[int((len(sensi_sorted)-1)
                                           *min(abs(eps_r[tsp]),1))]
            for r in range(nr):
                if not mech_data.react.activ_m[r]:
                    active_reactions[r] = False
                else:
                    if sensi_r[tsp][r]>lim_val and not active_reactions[r]:
                        active_reactions[r]=True
        try:
            sensi_Sl_sorted = list(red_data.red_op.sensi_Sl)
            sensi_Sl_sorted.sort()
            while 0 in sensi_Sl_sorted: sensi_Sl_sorted.remove(0)
            lim_val = sensi_Sl_sorted[int((len(sensi_Sl_sorted)-1)\
                                           *min(abs(eps_r[tsp]),1))]
            for r in range(nr):
                if not mech_data.react.activ_m[r]:
                    active_reactions[r] = False
                else:
                    if red_data.red_op.sensi_Sl[r]>lim_val and not active_reactions[r]:
                        active_reactions[r]=True
        except:
            A=2
        try:
            sensi_T_sorted = list(red_data.red_op.sensi_T)
            sensi_T_sorted.sort()
            while 0 in sensi_T_sorted: sensi_T_sorted.remove(0)
            lim_val = sensi_T_sorted[int((len(sensi_T_sorted)-1)\
                                           *min(abs(eps_r[tsp]),1))]
            for r in range(nr):
                if not mech_data.react.activ_m[r]:
                    active_reactions[r] = False
                else:
                    if red_data.red_op.sensi_T[r]>lim_val and not active_reactions[r]:
                        active_reactions[r]=True
        except:
            A=2

        for sp in range(ns):
            if not active_species[sp]:
                if mech_data.spec.activ_m[sp]:
                    for r in range(nr):
                        if active_reactions[r]:
                            if nu_f[sp][r] != 0 or nu_r[sp][r] != 0:
                                active_species[sp]=True
                                break

    # check threebody exception (+AR) (+HE) etc.
    for r in range(len(mech_data.react.formula)):
        if mech_data.react.type[r]=="falloff_reaction" \
        and type(mech_data.react.tbe[r]) is str:
            for sp in range(len(mech_data.spec.name)):
                if mech_data.react.tbe[r]==mech_data.spec.name[sp]\
                and not active_species[sp]:
                    active_reactions[r]=False



    return active_reactions, active_species




def write_sensitivities(red_data,sensi_scatter,S_react_x,T_react_z,S_AB_tsp_z_r,red_results):

    import os

    #main variables
    conditions = red_results.conditions
    config  = conditions.config
    fuel    = conditions.composition.fuel
    phi     = conditions.composition.phi
    T       = conditions.state_var.T
    P       = conditions.state_var.P
    tsp     = red_data.tspc

    r_path = os.getcwd()
    main_path = r_path + '/SA_results'

    if not os.path.exists(main_path):  os.mkdir(main_path)
    os.chdir(main_path)

    s_i=0
    for s in sensi_scatter:
        s_i+=1
        csv_file = 'Sensitivities_'+config+'_'+fuel+'_'+'%0.1f' %phi\
                    +'_T'+'%0.0f' %T+'_P'+'%0.2f' %(P/10e5) \
                    +'_('+str(s_i)+').csv'

        fichier_data = open(csv_file, 'w')

        # Write conditions:
        l1_conditions = "config;fuel;oxidant;diluent;phi;diluent_ratio;P(Pa)"
        if "flame" in conditions.config:
            l1_conditions += ";z(m);Ti(K);rtol_ts;atol_ts;rtol_ss;atol_ss;transport_model;Sl0(m/s)"
        elif "reactor" in conditions.config:
            l1_conditions += ";t(s);Ti(K);rtol_ts;atol_ts;ig_time(s)"
        elif "JSR" in conditions.config:
            l1_conditions += ";time(s);rtol_ts;atol_ts"
        l2_conditions = "\n"+conditions.config+";"\
                +conditions.composition.fuel+";"\
                +conditions.composition.oxidant+";"\
                +conditions.composition.diluent+";"\
                +str(conditions.composition.phi)+";"\
                +str(conditions.composition.diluent_ratio)+";"\
                +str(conditions.state_var.P)+";"
        if "flame "in conditions.config:
            l2_conditions+=str(red_results.pts_scatter[s])+";"\
                    +str(conditions.state_var.T)+";"\
                    +str(conditions.simul_param.tol_ts[0])+";"\
                    +str(conditions.simul_param.tol_ts[1])+";"\
                    +str(conditions.simul_param.tol_ss[0])+";"\
                    +str(conditions.simul_param.tol_ss[1])+";"\
                    +conditions.simul_param.transport_model+";"\
                    +str(red_results.Sl)
        elif "reactor "in conditions.config:
            l2_conditions+=str(red_results.pts_scatter[s])+";"\
                    +str(conditions.state_var.T)+";"\
                    +str(conditions.simul_param.tol_ts[0])+";"\
                    +str(conditions.simul_param.tol_ts[1])+";"\
                    +str(red_results.ign_time)
        elif "JSR":
            l2_conditions+=str(red_results.pts_scatter[s])+";"\
                    +str(conditions.simul_param.end_sim)+";"\
                    +str(conditions.simul_param.tol_ts[0])+";"\
                    +str(conditions.simul_param.tol_ts[1])
        fichier_data.write(l1_conditions)
        fichier_data.write(l2_conditions)
        fichier_data.write("\n")

        # l1 headers
        if "flame" in config:
            fichier_data.write("n_react;T;Sl")
        else:
            fichier_data.write("n_react")
        for tsp_i in tsp:
            fichier_data.write(";" + str(tsp_i))
        fichier_data.write("\n")


        # Sensitivities
        for r_i in range(len(S_react_x[0][0])):
            txt = str(r_i+1)
            if "flame" in config:
                if 'burner' not in config:
                    txt += ';'+ '%0.5e' %T_react_z[s_i-1][r_i] + ';' + '%0.5e' %red_data.red_op.sensi_Sl[r_i]
                else:
                    txt += ';'+ '%0.5e' %T_react_z[s_i-1][r_i]
            for tsp_i in range(len(S_react_x[0])):
                txt += ';'+ '%0.5e' %S_react_x[s_i-1][tsp_i][r_i]
            fichier_data.write(txt+'\n')

    os.chdir(r_path)



def plot_sentitivities(sp2plot,red_results,mech_data,sensi_scatter,S_react_x,T_react_z,\
                 red_data,reaction_2plot="all",reaction_label="number"):

    import os
    import matplotlib.pyplot as plt; plt.rcdefaults()

    #main variables
    conditions = red_results.conditions
    config  = conditions.config
    fuel    = conditions.composition.fuel
    phi     = conditions.composition.phi
    T       = conditions.state_var.T
    P       = conditions.state_var.P
    verbose = conditions.simul_param.verbose


    r_path = os.getcwd()
    main_path = r_path + '/SA_results'

    if not os.path.exists(main_path):  os.mkdir(main_path)
    os.chdir(main_path)



    # global sensitivities
    n_sp=-1
    for sp in sp2plot:
        sens_pos=[]; sens_neg=[]
        if sp!='T' and sp!='Sl':
            n_sp+=1
            for r in range(len(red_data.red_op.sensi_r[n_sp])):
                sens_pos.append((r+1,red_data.red_op.sensi_r[n_sp][r]))
                sens_neg.append((r+1,red_data.red_op.sensi_r[n_sp][r]))
        elif sp=='T' and 'flame' in config:
            for r in range(len(red_data.red_op.sensi_T)):
                sens_pos.append((r+1,red_data.red_op.sensi_T[r]))
                sens_neg.append((r+1,red_data.red_op.sensi_T[r]))
        elif sp=='Sl' and 'flame' in config:
            for r in range(len(red_data.red_op.sensi_Sl)):
                sens_pos.append((r+1,red_data.red_op.sensi_Sl[r]))
                sens_neg.append((r+1,red_data.red_op.sensi_Sl[r]))

        if len(sens_pos)!=0:
            # sort negative and positive sensitivities
            sens_pos.sort(reverse=True, key=lambda col: col[1])
            sens_neg.sort(reverse=False, key=lambda col: col[1])

            if reaction_2plot == "all": reaction_2plot = len(sens_pos)
            sens_vect = [] ; reactions = [] ; n_neg = 0 ; n_pos = 0
            for r in range(reaction_2plot):
                # get maximal negative and positive sensitivities
                if sens_pos[n_pos][1]>abs(sens_neg[n_neg][1]):
                    sens_vect.append(sens_pos[n_pos][1]) ; n_pos+=1
                    if reaction_label == "number":
                        reactions.append(str(mech_data.react.number[sens_pos[n_pos][0]]))
                    elif reaction_label == "formula":
                        reactions.append(str(mech_data.react.formula[sens_pos[n_pos][0]]))
                    elif reaction_label == "reaction+formula":
                        reactions.append("("+str(mech_data.react.number[sens_pos[n_pos][0]])\
                                         +") "+mech_data.react.formula[sens_pos[n_pos][0]])
                else:
                    sens_vect.append(sens_neg[n_neg][1]) ; n_neg+=1
                    if reaction_label == "number":
                        reactions.append(str(mech_data.react.number[sens_neg[n_neg][0]]))
                    elif reaction_label == "formula":
                        reactions.append(str(mech_data.react.formula[sens_neg[n_neg][0]]))
                    elif reaction_label == "reaction+formula":
                        reactions.append("("+str(mech_data.react.number[sens_neg[n_neg][0]])\
                                         +") "+mech_data.react.formula[sens_neg[n_neg][0]])
            y_pos = np.arange(len(reactions))

            fig, ax = plt.subplots()
            fig.set_size_inches(7, 2+len(reactions)/6)
            fig.subplots_adjust(left = 0.5, bottom = 0.15,
                           right = 0.9, top = 0.95, wspace = 0, hspace = 0.5)

            ax.barh(y_pos, sens_vect, align='center', alpha=1)
            ax.set_yticks(y_pos)
            ax.set_yticklabels(reactions)
            xlabel = sp + ' global reaction sensitivity'
            ax.set_xlabel(xlabel)

            if verbose>4: plt.show()

            fig_title = 'Global_Sensitivities_'+config+'_'+fuel+'_'+'%0.1f' %phi   \
                    +'_T'+'%0.0f' %T+'_P'+'%0.2f' %(P/10e5)+'_'+sp+'.png'

            fig.savefig(fig_title)






    # local sensitivities
    for s in sensi_scatter:
        n_sp=-1
        for sp in sp2plot:
            sens_pos=[]; sens_neg=[]
            if sp!='T' and sp!='Sl':
                n_sp+=1
                for r in range(len(S_react_x[s][n_sp])):
                    sens_pos.append((r+1,S_react_x[s][n_sp][r]))
                    sens_neg.append((r+1,S_react_x[s][n_sp][r]))
            elif sp=='T' and 'flame' in config:
                for r in range(len(T_react_z[s])):
                    sens_pos.append((r+1,T_react_z[s][r]))
                    sens_neg.append((r+1,T_react_z[s][r]))

            if len(sens_pos)!=0 and sp!='Sl':
                # sort negative and positive sensitivities
                sens_pos.sort(reverse=True, key=lambda col: col[1])
                sens_neg.sort(reverse=False, key=lambda col: col[1])

                if reaction_2plot == "all": reaction_2plot = len(sens_pos)
                sens_vect = [] ; reactions = [] ; n_neg = 0 ; n_pos = 0
                for r in range(reaction_2plot):
                    # get maximal negative and positive sensitivities
                    if sens_pos[n_pos][1]>abs(sens_neg[n_neg][1]):
                        sens_vect.append(sens_pos[n_pos][1]) ; n_pos+=1
                        if reaction_label == "number":
                            reactions.append(str(mech_data.react.number[sens_pos[n_pos][0]]))
                        elif reaction_label == "formula":
                            reactions.append(str(mech_data.react.formula[sens_pos[n_pos][0]]))
                        elif reaction_label == "reaction+formula":
                            reactions.append("("+str(mech_data.react.number[sens_pos[n_pos][0]])\
                                             +") "+mech_data.react.formula[sens_pos[n_pos][0]])
                    else:
                        sens_vect.append(sens_neg[n_neg][1]) ; n_neg+=1
                        if reaction_label == "number":
                            reactions.append(str(mech_data.react.number[sens_neg[n_neg][0]]))
                        elif reaction_label == "formula":
                            reactions.append(str(mech_data.react.formula[sens_neg[n_neg][0]]))
                        elif reaction_label == "reaction+formula":
                            reactions.append("("+str(mech_data.react.number[sens_neg[n_neg][0]])\
                                             +") "+mech_data.react.formula[sens_neg[n_neg][0]])
                y_pos = np.arange(len(reactions))

                fig, ax = plt.subplots()
                fig.set_size_inches(7, 2+len(reactions)/6)
                fig.subplots_adjust(left = 0.5, bottom = 0.15,
                               right = 0.9, top = 0.95, wspace = 0, hspace = 0.5)

                ax.barh(y_pos, sens_vect, align='center', alpha=1)
                ax.set_yticks(y_pos)
                ax.set_yticklabels(reactions)
                xlabel = sp + ' reaction sensitivity at '
                if 'reactor' in red_results.conditions.config :
                    xlabel+='t = '+'%0.2f'%(red_results.pts_scatter[s]*1000)+' ms'
                elif 'flame' in red_results.conditions.config :
                    xlabel+='z = '+'%0.1f'%(red_results.pts_scatter[s]*1000)+' mm'
                elif 'JSR' in red_results.conditions.config :
                    xlabel+='T = '+'%.0f'%(red_results.pts_scatter[s])+'K'
                ax.set_xlabel(xlabel)

                if verbose>4: plt.show()

                fig_title = 'Sensitivities_'+config+'_'+fuel+'_'+'%0.1f' %phi   \
                        +'_T'+'%0.0f' %T+'_P'+'%0.2f' %(P/10e5)+'_'+sp \
                        +'_('+str(s)+').png'

                fig.savefig(fig_title)

    os.chdir(r_path)




