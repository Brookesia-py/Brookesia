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

import numpy as np
import time as timer
import brookesia.Class_def as cdef
from  brookesia.Class_def import print_
import multiprocessing
import os



def dic_par(dic_par_arg_i):

    dic_par_arg     = dic_par_arg_i[0]

#    k_f             = dic_par_arg_i[1]
#    k_r             = dic_par_arg_i[2]
#    conc            = dic_par_arg_i[3]
    reactionRate    = dic_par_arg_i[4]
    i               = dic_par_arg_i[5]
    ns              = dic_par_arg[0]
    nr              = dic_par_arg[1]
    react_activ     = dic_par_arg[2]
    spec_activ            = dic_par_arg[3]
    nu              = dic_par_arg[4]
#    nu_f            = dic_par_arg[5]
#    nu_r            = dic_par_arg[6]
    kronecker       = dic_par_arg[7]
#    bar             = dic_par_arg[10]

    dic=False

    # reaction rates computation
#    reactionRate = np.zeros(nr)
#    fRate = np.ones(nr)
#    rRate = np.ones(nr)
    num = np.zeros((ns,ns))
    den = np.zeros(ns)
    dic = np.zeros((ns, ns))
#    for r in range(nr):
#        if react_activ[r]:
#            for sp in range(ns):
#                if nu_f[sp, r]!=0 or nu_r[sp, r]!=0:
#                    fRate[r] *= (1e3*conc[sp])**nu_f[sp, r]
#                    rRate[r] *= (1e3*conc[sp])**nu_r[sp, r]
#            fRate[r] *= k_f[r]
#            rRate[r] *= k_r[r]
#            reactionRate[r] = fRate[r] - rRate[r]



    # Computation of Direct Interaction Coefficients at time n
    for r in range(nr):
        if react_activ[r]:
            for spA in range(ns):
                for spB in range(ns):
                    num[spA, spB] += (nu[spA, r]*reactionRate[r] \
                                     * kronecker[spB, r])

    for spA in range(ns):
        if spec_activ[spA]:
            PA = 0
            CA = 0
            for r in range(nr):
                PA += max(0, nu[spA, r]*reactionRate[r])
                CA += max(0, -nu[spA, r]*reactionRate[r])
            den[spA] = max(PA, CA)
    for spA in range(ns):
        if spec_activ[spA]:
            if den[spA] == 0:
                dic[spA, :] =0
            else:
                for spB in range(ns):
                    dic[spA, spB] = abs(num[spA, spB])/den[spA]

    return (dic,i)




def dic(red_data,mech_data,results):
    """r_AB =sum(nu_iA rate_i delta_Bi) /max(prod A, conso A) """

    mp = results.conditions.main_path

    gas_ref = red_data.gas_ref
    gas_red = red_data.gas_red

    verbose = red_data.verbose


    if verbose >=2 :
        print_("\n  direct interaction coefficients computation ...",mp)

    timeDRG_1 = timer.time()


    ns = gas_ref.n_species                      # number of species
    nr = gas_ref.n_reactions                    # number of elementary reactions
    # stoechiometric coefficients : direct, reverse, net
    nu_f = gas_ref.reactant_stoich_coeffs()
    nu_r = gas_ref.product_stoich_coeffs()
    nu = nu_f - nu_r

    kronecker = np.zeros((ns,nr))
    for r in range(nr):
        for sp in range(ns):
            if nu_f[sp,r]!=0 or nu_r[sp,r]!=0:
                kronecker[sp,r] = 1


    if gas_red.n_species > 25 :
        bar = cdef.ProgressBar(len(results.pts_scatter))
        bar.update(0)
    drg_i=[]
    n_points = len(results.r_rate)

    ib = []; dic_par_arg_i=[]
    dic_par_arg = [ns, nr,
                   mech_data.react.activ_m, mech_data.spec.activ_m,
                   nu,nu_f,nu_r,
                   kronecker]#,

    for i in range(len(results.r_rate)):
        if i%int(max(n_points/red_data.red_op.n_points,1))==0 and i!=0:
            drg_i.append(i)
            ib.append(i)
            dic_par_arg_i.append([dic_par_arg,results.kf[i],results.kr[i],\
                                  results.conc[i],results.r_rate[i],i])


    num_cores = multiprocessing.cpu_count()

    dic=[]
    if os.name == 'nt': multiprocessing.get_context('spawn')
    with multiprocessing.Pool(num_cores) as p:
        dic=p.map(dic_par, dic_par_arg_i)

    dic.sort(reverse=False, key=lambda col: col[1])

    interactionCoefficients = []

    for i in range(len(dic)):
        interactionCoefficients.append(dic[i][0])

    if gas_red.n_species > 25 :
        bar.update(len(results.pts_scatter));print_("\n",mp)


#        bar.update(n)
    timeDRG_2 = timer.time()
    #  Display options
    if verbose >=4:
        print_("      time for the DIC computation: "+str(round(timeDRG_2-timeDRG_1))+'s',mp)
    elif verbose >=3 :
        print_("   completed",mp)


    red_data.red_op.interaction_coeffs = list(interactionCoefficients)

    del interactionCoefficients

    if red_data.write_results:
        write_sensitivities(red_data,results,drg_i)

    return red_data



def ric(red_data, mech_data, red_results):
    """
    Reaction Interaction Coefficient
    """

    mp = red_results.conditions.main_path

    # main variables
    gas_ref  = red_data.gas_ref
    verbose=red_data.verbose
    n_points = len(red_results.r_rate)

#    n_sp_ref  = gas_ref.n_species
    n_r_ref   = gas_ref.n_reactions
    tspc   = red_data.red_op.new_targets_4_DRG_r
    tsp_idx=[]
    while '' in tsp_idx: tspc.remove('')
    for i in range(len(tspc)):
        	tsp_idx.append(gas_ref.species_index(tspc[i]))


    n_tsp     = len(tsp_idx)

    nu_f = gas_ref.reactant_stoich_coeffs()
    nu_r = gas_ref.product_stoich_coeffs()
    nu = nu_f - nu_r

    r_interCoeff =  np.zeros((n_tsp,n_r_ref))

    div_DRG_points = round(n_points/(red_data.red_op.n_points-1))

    if verbose >=8 :
        print_("   reactions direct interaction coefficients computation ...",mp)

    bar = cdef.ProgressBar(n_points)
    bar.update(0)
    for i in range(n_points):
        if i%max(div_DRG_points,1)==0:
            reactionRate = red_results.r_rate[i]

#            k_f = red_results.kf[i]
#            k_r = red_results.kr[i]
            # reaction rates computation
#            reactionRate = np.zeros(n_r_ref)
#            fRate = np.ones(n_r_ref)
#            rRate = np.ones(n_r_ref)
#            # Reaction rate calculation
#            r_red=-1
#            for r in range(n_r_ref):
#                if mech_data.react.activ_p[r]\
#                or True not in mech_data.react.activ_p: # (1st reduction)
#                    r_red+=1
#                    for spA in range(n_sp_ref):
#                        if mech_data.spec.activ_p[spA] \
#                        or True not in mech_data.spec.activ_p:
#                            if nu_f[spA, r]!=0 or nu_r[spA, r]!=0:
#                                fRate[r]*=(1e3*red_results.conc[i][spA])**nu_f[spA, r]
#                                rRate[r]*=(1e3*red_results.conc[i][spA])**nu_r[spA, r]
#                    fRate[r] *= k_f[r]
#                    rRate[r] *= k_r[r]
#                    reactionRate[r] = fRate[r] - rRate[r]
#
            # Computation of Direct Interaction Coefficients at time n
            for t_sp in range(n_tsp):
                PA = 0
                CA = 0
                for r in range(n_r_ref):
                    if mech_data.react.activ_p[r]\
                    or True not in mech_data.react.activ_p: # (1st reduction)
                        PA += max(0, nu[tsp_idx[t_sp], r]*reactionRate[r])
                        CA += max(0, -nu[tsp_idx[t_sp], r]*reactionRate[r])
                den = max(PA, CA)
                for r in range(n_r_ref):
                    if mech_data.react.activ_p[r]\
                    or True not in mech_data.react.activ_p: # (1st reduction)
                        num = abs(nu[tsp_idx[t_sp], r]*reactionRate[r])
                        if den > 0:
                            r_interCoeff[t_sp,r]=max(num/den,r_interCoeff[t_sp,r])
        bar.update(i)

    bar.update(n_points);print_("\n",mp)


    red_data.red_op.r_interaction_coeffs = list(r_interCoeff)
    del r_interCoeff
    return red_data





def graphSearch(conditions,red_data,mech_data,eps):

    # main variables
    gas     = red_data.gas_ref
    mp = conditions.main_path

    new_targets_4_DRG_r = red_data.red_op.new_targets_4_DRG_r

    IC = red_data.red_op.interaction_coeffs

    target_species = red_data.targetSpeciesIdx
    verbose=conditions.simul_param.verbose

    if verbose >=7 :
            print_("\n  performing"" DRG reduction with eps = "+\
                  str(np.around(eps, decimals=3))+" ...",mp)

    timeDRG_1 = timer.time()

    active_species = list(mech_data.spec.activ_p)

    # Graph construction
    points = len(IC)

    ns = len(IC[0][1,:])

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

    new_target_tsp = False

    if red_data.red_op.graph_search == 'DFS':

        # DFS graph analysis
        OIC_sp = np.zeros((len(target_species),ns))

        for s in range(len(target_species)):
            IC_tsp = 0
            for n in range(points):
                graph = []
                OIC = []
                for spA in range(ns):
                    graph.append([])
                    OIC.append(IC[n][s, spA])
                    if IC[n][s, spA]>OIC_sp[s][spA]:
                        OIC_sp[s][spA]=IC[n][s, spA]
                    # new targets for DRG_r:
                    if spA not in new_targets_4_DRG_r:
                        if IC[n][s, spA] > IC_tsp:
                            IC_tsp = IC[n][s, spA]
                            new_target_tsp = spA
                    for spB in range(ns):
                        if IC[n][spA, spB] > eps[s]:
                            graph[spA].append(spB)
                start_nodes = []
                start_nodes.append(target_species[s])
                visited = []
                while start_nodes!=[]:
                    node = start_nodes.pop()
                    if node not in visited:
                        visited.append(node)
                        adj = graph[node]
                        for k in range(len(adj)):
                            if adj[k] not in visited:
                                start_nodes.append(graph[node][k])

                for ind in range(len(visited)):
                    k = visited[ind]
                    if not active_species[k] and mech_data.spec.activ_m[k]:
                        active_species[k]=True
            if new_target_tsp:
                if new_target_tsp not in new_targets_4_DRG_r:
                    new_targets_4_DRG_r.append(new_target_tsp)


    elif red_data.red_op.graph_search == 'Dijkstra':
        # Dijkstra algorithm
        OIC_sp = np.zeros((len(target_species),ns))

        for tsp in range(len(target_species)):
#            print(target_species)
            tsp_idx = target_species[tsp]
            IC_tsp = 0
            for p in range(points):
                mpq_list = []
                # max-priority queue (mpq) construction
                OIC = []
                for spA in range(ns):
                    OIC.append(IC[p][tsp_idx, spA])
                    if IC[p][tsp_idx, spA]>OIC_sp[tsp][spA]:
                        OIC_sp[tsp][spA]=IC[p][tsp_idx, spA]
                    if OIC[spA] > eps[tsp]:
                        mpq_list.append(spA)
                    # new targets for DRG_r:
                    if spA not in new_targets_4_DRG_r:
                        if IC[p][tsp_idx, spA] > IC_tsp:
                            IC_tsp = IC[p][tsp_idx, spA]
                            new_target_tsp = spA

                for mpq in range(len(mpq_list)):
                    mpq_idx = mpq_list[mpq]
                    R_max=[] ; visited_sp=[]
                    for sp in range(ns):
                        R_max.append(IC[p][mpq_idx,sp])
                        visited_sp.append([tsp,sp])
                        if OIC[mpq]*IC[p][mpq_idx,sp] > eps[tsp]:
                            active_species[sp] = True
                    R_max[tsp]=0 ; R_max[mpq]=0


                    # graph search
                    idx_node = R_max.index(max(R_max))
                    while max(R_max)>eps[tsp]:
#                        print(max(R_max))
#                        print(eps[tsp])
                        # node informations
                        visited_sp_node = visited_sp[idx_node]
                        # interaction coefficients of the step
                        Ri = []
                        for sp in range(ns):
                            if sp not in visited_sp_node:
                                Ri.append(IC[p][idx_node,sp])
                                if Ri[-1]>eps[tsp]:
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
            if new_target_tsp:
                if new_target_tsp not in new_targets_4_DRG_r:
                    new_targets_4_DRG_r.append(new_target_tsp)

    timeDRG_2 = timer.time()
    #  Display options
    if verbose >= 5:
        print_("  time for graph search : "+'%3.0f' %(timeDRG_2-timeDRG_1)+'s',mp)
    if verbose >= 6:
        print_("     "+"kept species:",mp)
        kept_species="     "
        for k in range(len(active_species)):
            if active_species[k]:
                kept_species+=gas.species_name(k)+", "
        print_(kept_species,mp)

    red_data.red_op.OIC_sp = list(OIC_sp)
    red_data.red_op.new_targets_4_DRG_r = new_targets_4_DRG_r

    return active_species




def graphSearch_DRGEP(conditions,red_data,mech_data,eps):

    # main variables
    gas     = red_data.gas_ref

    new_targets_4_DRG_r = red_data.red_op.new_targets_4_DRG_r


    IC = red_data.red_op.interaction_coeffs
    target_species = red_data.targetSpeciesIdx
    verbose=conditions.simul_param.verbose


    if verbose >=7 :
            print(  "\n  performing"" DRG reduction with eps = ", \
                  np.around(eps, decimals=3) , " ...")

    timeDRG_1 = timer.time()

    active_species = list(mech_data.spec.activ_p)

    points = len(IC)

    ns = len(IC[0][1,:])

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

    # Dijkstra algorithm
    OIC_sp = np.zeros((len(target_species),ns))
    new_target_tsp = False
    for tsp in range(len(target_species)):
        tsp_idx = target_species[tsp]
        IC_tsp = 0
        for p in range(points):
            mpq_list = []

            # max-priority queue (mpq) construction
            OIC = []
            for spA in range(ns):
                OIC.append(IC[p][tsp_idx, spA])
                if IC[p][tsp_idx, spA]>OIC_sp[tsp][spA]:
                    OIC_sp[tsp][spA]=IC[p][tsp_idx, spA]
                if OIC[spA] > eps[tsp]:
                    mpq_list.append(spA)
                # new targets for DRG_r:
                if spA not in new_targets_4_DRG_r:
                    if IC[p][tsp_idx, spA] > IC_tsp:
                        IC_tsp = IC[p][tsp_idx, spA]
                        new_target_tsp = spA

            for mpq in range(len(mpq_list)):
                mpq_idx = mpq_list[mpq]
#                R = np.array((ns,ns))
                R_max=[] ; visited_sp=[]
                for sp in range(ns):
                    R_max.append(OIC[mpq]*IC[p][mpq_idx,sp])
                    visited_sp.append([tsp,sp])
                    if OIC[mpq]*IC[p][mpq_idx,sp] > eps[tsp]:
                        active_species[sp] = True
                R_max[tsp]=0 ; R_max[mpq]=0

                # graph search
                idx_node = R_max.index(max(R_max))
                while max(R_max)>eps[tsp]:
                    # node informations
                    visited_sp_node = visited_sp[idx_node]

                    # interaction coefficients of the step
                    Ri = []
                    for sp in range(ns):
                        if sp not in visited_sp_node:
                            Ri.append(R_max[idx_node]*IC[p][idx_node,sp])
                            if Ri[-1]>eps[tsp]:
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
        if new_target_tsp:
            if new_target_tsp not in new_targets_4_DRG_r:
                new_targets_4_DRG_r.append(new_target_tsp)


    timeDRG_2 = timer.time()
    #  Display options
    if verbose >= 5:
        print("  time for graph search : ", '%3.0f' %(timeDRG_2-timeDRG_1), 's')
    if verbose >= 6:
        print("     ", "kept species:")
        kept_species="     "
        for k in range(len(active_species)):
            if active_species[k]:
                kept_species+=gas.species_name(k)+", "
        print(kept_species)

    red_data.red_op.OIC_sp = list(OIC_sp)
    red_data.red_op.new_targets_4_DRG_r = new_targets_4_DRG_r


    return active_species




def reactionWithdrawal(mech_data,red_data,red_method,eps_r,conditions,active_species=False,red_data_list=False,step_n=False):
    """
    Cas 1 DRG_reactionWithdraw = True : retrait des reactions avec un DIC faible + retrait des reactions impliquant des espèces deja retirees
    Cas 2 DRG_reactionWithdraw = False : retrait des reactions impliquant des espèces deja retirees
    """

    # main variables
    gas_ref         = red_data.gas_ref
    verbose         = red_data.verbose

    time_1 = timer.time()
    n_sp_ref  = gas_ref.n_species
    n_r_ref   = gas_ref.n_reactions

    nu_f =gas_ref.reactant_stoich_coeffs()
    nu_r =gas_ref.product_stoich_coeffs()
#    nu = nu_f-nu_r


    if verbose >=8 :
        print(  "  reactions removal ...")


    if not red_data.red_op.OIC_sp:
        IC              = red_data.red_op.interaction_coeffs
        points          = len(IC)
        ns              = len(IC[0][1,:])
        target_species  = red_data.targetSpeciesIdx
        OIC_sp          = np.zeros((len(target_species),ns))
        for tsp in range(len(target_species)):
            tsp_i = target_species[tsp]
            for p in range(points):
                OIC = []
                for spA in range(ns):
                    OIC.append(IC[p][tsp_i, spA])
                    if IC[p][tsp_i, spA]>OIC_sp[tsp][spA]:
                        OIC_sp[tsp][spA]=IC[p][tsp_i, spA]
        red_data.red_op.OIC_sp = list(OIC_sp)

    # define new activated species / reactions based on species sensitivities
    if '_sp' in red_method:
        active_reactions = [True]*len(mech_data.react.activ_p)
        for r in range(n_r_ref):
            if not mech_data.react.activ_m[r]:
                active_reactions[r] = False
            else:
                for sp in range(n_sp_ref):
                    # remove reactions involving non active species
                    if not active_species[sp]:
                        if nu_f[sp][r]!=0 or nu_r[sp][r]!=0:
                            active_reactions[r] = False
                            break


    # define new activated species / reactions based on reaction interaction coeffs
    if '_r' in red_method:
        tsp_idx         = red_data.red_op.new_targets_4_DRG_r
#        print(tsp_idx)
#        print(red_data.targetSpeciesIdx)
        r_inter_coeffs  = red_data.red_op.r_interaction_coeffs

        if red_data.red_op.first_step_DRG_r:
            active_reactions = list(mech_data.react.activ_m)
        else:
            active_reactions = list(mech_data.react.activ_p)
        active_species = list(mech_data.spec.activ_p)

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



        for t in range(len(tsp_idx)):
            spA=tsp_idx[t]
            if t<len(eps_r): eps = eps_r[t]
            else:            eps = np.mean(eps_r)

            if not active_species[spA]: active_species[spA]=True
            r_int_c_sorted = list(r_inter_coeffs[t])
            r_int_c_sorted.sort()
            while 0 in r_int_c_sorted: r_int_c_sorted.remove(0)
            if len(r_int_c_sorted)>0:  # if sens are not equal to zero
                idx_lim_val = int((len(r_int_c_sorted)-1)*min(abs(eps),1))
                lim_val = r_int_c_sorted[idx_lim_val]
                for r in range(n_r_ref):
                    if not mech_data.react.activ_m[r]:
                        active_reactions[r] = False
                    else:
                        if r_inter_coeffs[t][r]!=0:
                            if red_data.red_op.first_step_DRG_r:
                                if r_inter_coeffs[t][r]<lim_val:
                                    active_reactions[r]=False
                            else:
                                if r_inter_coeffs[t][r]>=lim_val and not active_reactions[r]:
                                    active_reactions[r]=True

        # if the DRG_r cannot integrate enougth reactions :
        # previous DRG_r assessment for additional reactions integration
        if step_n>1\
        and conditions.number>1\
        and idx_lim_val==0\
        and active_reactions.count(True) == mech_data.react.activ_p.count(True):


            for _i in range(conditions.number):
                r_inter_coeffs_i = red_data_list[_i].red_op.r_interaction_coeffs
                for t in range(len(tsp_idx)):
                    for r in range(n_r_ref):
                        if r_inter_coeffs_i[t][r]>r_inter_coeffs[t][r]:
                            r_inter_coeffs[t][r]=r_inter_coeffs_i[t][r]

            for t in range(len(tsp_idx)):
                spA=tsp_idx[t]
                if t<len(eps_r): eps = eps_r[t]
                else:            eps = np.mean(eps_r)

                if not active_species[spA]: active_species[spA]=True
                r_int_c_sorted = list(r_inter_coeffs[t])
                r_int_c_sorted.sort()
                while 0 in r_int_c_sorted: r_int_c_sorted.remove(0)
                if len(r_int_c_sorted)>0:  # if sens are not equal to zero
                    idx_lim_val = int((len(r_int_c_sorted)-1)*min(abs(eps),1))
                    lim_val = r_int_c_sorted[idx_lim_val]
                    for r in range(n_r_ref):
                        if not mech_data.react.activ_m[r]:
                            active_reactions[r] = False
                        else:
                            if r_inter_coeffs[t][r]!=0:
                                if red_data.red_op.first_step_DRG_r:
                                    if r_inter_coeffs[t][r]<lim_val:
                                        active_reactions[r]=False
                                else:
                                    if r_inter_coeffs[t][r]>=lim_val and not active_reactions[r]:
                                        active_reactions[r]=True


        for sp in range(n_sp_ref):
            if not active_species[sp]:
                if mech_data.spec.activ_m[sp]:
                    for r in range(n_r_ref):
                        if active_reactions[r]:
                            if nu_f[sp][r]!=0 or nu_r[sp][r]!=0:
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


    time_2 = timer.time()

    #  Display options
    if verbose >=6:
        print("   time for graph search: "+'%3.0f' %(time_2-time_1)+'s')
    if verbose >=8 :
        print("   done")


    return active_reactions, active_species



def write_sensitivities(red_data,red_results,drg_i):

    import os

    #main variables
    conditions = red_results.conditions
    config  = conditions.config
    fuel    = conditions.composition.fuel
    phi     = conditions.composition.phi
    T       = conditions.state_var.T
    P       = conditions.state_var.P
    tsp     = red_data.tspc
    tsp_idx = red_data.targetSpeciesIdx
    gas     = red_data.gas_ref
    drg_i.remove(0)
    interaction_coeffs = red_data.red_op.interaction_coeffs

    r_path = os.getcwd()
    main_path = r_path + '/DRG_results'

    if not os.path.exists(main_path):  os.mkdir(main_path)
    os.chdir(main_path)

    for s in range(len(drg_i)):
        node_csv_file = 'DRG_dic_nodes_'+config+'_'+fuel+'_'+'%0.1f' %phi\
                    +'_T'+'%0.0f' %T+'_P'+'%0.2f' %(P/10e5) \
                    +'_('+str(drg_i[s])+').csv'
        node_file = open(node_csv_file, 'w')
        edge_csv_file = 'DRG_dic_edges_'+config+'_'+fuel+'_'+'%0.1f' %phi\
            +'_T'+'%0.0f' %T+'_P'+'%0.2f' %(P/10e5) \
            +'_('+str(drg_i[s])+').csv'
        edge_file = open(edge_csv_file, 'w')

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

        # nodes
        #    Id,Label,timeset,modularity_class
        #    11,Valjean,,1
        node_file.write("Id;Label;timeset;modularity_class\n")
        for sp in range(len(interaction_coeffs[s])):
            txt = str(sp+1) +";"+gas.species_name(sp)+";"+";"
            if sp in tsp_idx:
                txt += "10\n"
            else:
                txt += "2\n"
            node_file.write(txt)

        # edge
        #    Source,Target,Type,Id,Label,timeset,Weight
        #    1,0,Undirected,0,,,1
        #    2,0,Undirected,1,,,8
        Id = 0
        edge_file.write("Source;Target;Type;Id;Label;timeset;Weight\n")
        for spA in range(len(interaction_coeffs[s])):
            spB = spA+1
            while spB < len(interaction_coeffs[s]):
                txt = str(spA+1) +";"+str(spB+1)+";Undirected;"+str(Id)+";;;"\
                      +str(interaction_coeffs[s][spA][spB])+"\n"
                spB+=1 ; Id+=1
                edge_file.write(txt)

    os.chdir(r_path)

