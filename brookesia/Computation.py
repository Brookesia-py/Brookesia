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
import time as timer
import os
import sys
import pandas as pd

def ref_computation(conditions, verbose=0,act_sp=False,act_r=False):
                                          # act_sp and act_r used if previous red mech given
    mp = conditions.main_path

    # set language parameters (for xml cantera files)
    import locale as lc
    try:
        lc.setlocale(lc.LC_ALL, 'en_US.utf8')
    except:
        try:
            lc.setlocale(lc.LC_ALL, 'en_US.UTF-8')
        except:
            a=False


    if 'free_flame' in conditions.config or 'burner_flame' in conditions.config:

        # Main variables
        gas            = conditions.composition.gas
        tol_ss         = conditions.simul_param.tol_ss
        tol_ts         = conditions.simul_param.tol_ts
        xmax           = conditions.simul_param.end_sim

        if not conditions.exp_data:
            initial_points = list(conditions.simul_param.pts_scatter)
            flam_initial_grid = np.array(initial_points)*xmax
        else:
            initial_points = list(conditions.simul_param.pts_scatter)
            flam_initial_grid = np.array(initial_points)

        if conditions.import_data and not conditions.exp_data: refine_grid = False
        else:                                                  refine_grid = True

        # Flame computation
        print_("   Reference mecanism flame computation",mp)


        loglevel = 0
        gas.TPX = conditions.state_var.T,conditions.state_var.P,conditions.composition.X
        if 'free' in conditions.config:
            f = ct.FreeFlame(gas, flam_initial_grid)
            f.inlet.X = conditions.composition.X
            f.inlet.T = conditions.state_var.T
            f.flame.P = conditions.state_var.P
        elif 'burner' in conditions.config:
            f = ct.BurnerFlame(gas, flam_initial_grid)
            f.burner.T    = conditions.simul_param.T_profile[0]
            f.burner.X    = conditions.composition.X
            f.burner.mdot = conditions.simul_param.mdot    #mass flow rate per unit area [kg/m^2/s]
            z      = np.array(initial_points)/max(initial_points)
            T_prof = np.array(conditions.simul_param.T_profile)
            f.flame.set_fixed_temp_profile(z,T_prof)
        try:
            gas.transport_model = conditions.simul_param.transport_model
            f.transport_model = conditions.simul_param.transport_model
        except :
            print_("Unrecognized transport model. Proceeding with Mix model : ",mp)
            f.transport_model = "Mix"
        f.flame.set_steady_tolerances(default=tol_ss)
        f.flame.set_transient_tolerances(default=tol_ts)


        simul_success = False
        # -------- flame computation
        if conditions.simul_param.restore_flame:
            # -------- restore flame results
            try:
                os.chdir(conditions.r_path + '/_results_input/_flame_results/' + conditions.simul_param.flame_res_folder)
            except:
                try:
                    os.chdir(conditions.simul_param.flame_res_folder)
                except:
                    print_('Warning: restore flame folder not found',mp)

            flame_file = False
            for _file in os.listdir('.'):
                if conditions.num in _file[0:3]:
                    flame_file = _file
                    break
            if flame_file:
                if verbose >2:
                    print_('Restored file: ' + flame_file, mp)
                    # -----
                    # supress console output during the simulation
                if verbose<9:
                    old_stdout = sys.stdout ; old_stderr = sys.stderr
                    with open(os.devnull, "w") as devnull:
                        sys.stdout = devnull ; sys.stderr = devnull
                # -----

                try:
                    f.restore(flame_file, 'ref_solution')

                    os.chdir(conditions.main_path)

                    f.solve(auto = False, loglevel = 0, refine_grid = False)
                    simul_success = True

                    # ---- restore console output
                    if verbose<9: sys.stdout = old_stdout ; sys.stderr = old_stderr
                    if verbose >= 6:
                        print_("     Problem solved on ["+str(f.flame.n_points)+"] point grid",mp)
                except:
                    simul_success = False
                    # ---- restore console output
                    if verbose<9: sys.stdout = old_stdout ; sys.stderr = old_stderr
                    if verbose >= 3:
                        print("\n     WARNING: No solution found starting from the restored file\n",mp)

        if not conditions.simul_param.restore_flame or not flame_file or not simul_success:
            # Grid refinement
            ##  energy disabled
            f.energy_enabled = False
            f.set_refine_criteria(ratio = 7.0, slope = 1, curve = 1)
    #        f.set_time_step(5.0e-6,[10,20,50,80,120,150])
            f.solve(loglevel, refine_grid)
            if verbose >=5 : print_("Problem solved on ["+ str(f.flame.n_points)+ "] point grid",mp)

            ##  energy enabled
            if 'burner' not in conditions.config:
                f.energy_enabled = True
                _auto            = True
            else: _auto = False

            f.set_refine_criteria(ratio = 7.0, slope = 1, curve = 1)
            f.solve(loglevel, auto = _auto)
            if verbose >=5 : print_("Problem solved on ["+ str(f.flame.n_points)+ "] point grid",mp)

            f.set_refine_criteria(ratio = 5.0, slope = 0.5, curve = 0.5)
            f.solve(loglevel, auto = _auto)
            if verbose >=5 : print_("Problem solved on ["+ str(f.flame.n_points)+ "] point grid",mp)

            f.set_refine_criteria(ratio=2.0, slope=0.1, curve=0.1, prune=0.01)
            f.solve(loglevel, auto = _auto)
            if verbose >=5 : print_("Problem solved on ["+ str(f.flame.n_points)+ "] point grid",mp)

            ratio_ff = conditions.simul_param.ratio_ff
            slope_ff = conditions.simul_param.slope_ff
            curve_ff = conditions.simul_param.curve_ff
            prune_ff = conditions.simul_param.prune_ff
            f.set_refine_criteria(ratio=ratio_ff, slope=slope_ff, curve=curve_ff, prune=prune_ff)
            f.solve(loglevel, refine_grid)
            if verbose >=2 : print_("Problem solved on ["+ str(f.flame.n_points)+ "] point grid",mp)


        npoints = f.flame.n_points
        conc = [] ; kf = [] ; kr = [] ; r_rate = []

        for n in range(npoints):
            f.set_gas_state(n)
            if act_sp:
                conc_n = []; sa_i = 0
                for s_i in range(len(act_sp)):
                    if act_sp[s_i]:
                        conc_n.append(gas.concentrations[sa_i])
                        sa_i += 1
                    else:
                        conc_n.append(0)
                conc.append(conc_n)
            else:
                conc.append(gas.concentrations)

            if act_r:
                kf_r = []; kr_r = []; r_rate_r = []; ra_i = 0
                for r_i in range(len(act_r)):
                    if act_r[r_i]:
                        kf_r.append(gas.forward_rate_constants[ra_i])
                        kr_r.append(gas.reverse_rate_constants[ra_i])
                        r_rate_r.append(gas.net_rates_of_progress[ra_i])
                        ra_i += 1
                    else:
                        kf_r.append(0)
                        kr_r.append(0)
                        r_rate_r.append(0)
                kf.append(kf_r)
                kr.append(kr_r)
                r_rate.append(r_rate_r)
            else:
                kf.append(gas.forward_rate_constants)
                kr.append(gas.reverse_rate_constants)
                r_rate.append(gas.net_rates_of_progress)

        results = cdef.Sim_Results(conditions, gas, np.array(f.flame.grid), \
                          list(f.T), f.P, list(conc), list(kf), list(kr))
        results.Sl      = f.u[0]
        results.f       = f
        results.r_rate  = r_rate


        conditions.simul_param.pts_scatter = np.array(f.flame.grid)

        # Save results
        os.chdir(conditions.main_path)
        if 'Flame_ref_results' not in os.listdir():
            os.mkdir('Flame_ref_results')
        os.chdir('Flame_ref_results')
        if conditions.state_var.P>10000:
            fn = conditions.num+'ff_'+conditions.composition.fuel.replace('/','').split('(')[0]\
            +'_'+'%.2f' %conditions.composition.phi\
            +'_'+'%.0f'%conditions.state_var.T+'_'+'%.2f'%(conditions.state_var.P/1e5)\
            +'.xml'
        else:
            fn = conditions.num+'ff_'+conditions.composition.fuel.replace('/','').split('(')[0]\
            +'_'+'%.2f' %conditions.composition.phi\
            +'_'+'%.0f'%conditions.state_var.T+'_'+'%.0f'%(conditions.state_var.P)\
            +'.xml'
        f.save(fn,'ref_solution')
        os.chdir(conditions.main_path)

    if 'pp_flame' in conditions.config :
        if conditions.composition.X == conditions.composition.X2:
            print_('Warning: same mixture composition in burner 1 and 2',mp)
            print_('         -> condition config is modified to "tp_flame"',mp)
            conditions.config = 'tp_flame'

    if 'diff_flame' in conditions.config or 'pp_flame' in conditions.config :

        # Main variables
        gas   = conditions.composition.gas
        tol_ss         = conditions.simul_param.tol_ss
        tol_ts         = conditions.simul_param.tol_ts

        width = conditions.simul_param.end_sim

        if conditions.exp_data and 'cflow' in conditions.config:
            initial_points = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
        else:
            initial_points = np.array(conditions.simul_param.pts_scatter)\
                            /np.max(conditions.simul_param.pts_scatter)

        # =============================================================================
        # # PART 1: INITIALIZATION
        # =============================================================================

        # Set up an initial hydrogen-oxygen counterflow flame at 1 bar and low strain
        # rate (maximum axial velocity gradient = 2414 1/s)

        f = ct.CounterflowDiffusionFlame(gas, width=width)
#        f = ct.CounterflowPremixedFlame(gas, width=width)


        # Define the operating pressure and boundary conditions
        f.P = conditions.state_var.P  # 1 bar
        if 'diff_' in conditions.config:
            f.fuel_inlet.mdot = conditions.simul_param.mdot  # kg/m^2/s
            f.fuel_inlet.X =  conditions.composition.X
            f.fuel_inlet.T = conditions.state_var.T  # K
            f.oxidizer_inlet.mdot = conditions.simul_param.mdot2  # kg/m^2/s
            f.oxidizer_inlet.X =  conditions.composition.X2
            f.oxidizer_inlet.T = conditions.state_var.T2  # K
        else:
            if conditions.composition.phi > conditions.composition.phi2:
                f.fuel_inlet.mdot = conditions.simul_param.mdot  # kg/m^2/s
                f.fuel_inlet.X =  conditions.composition.X
                f.fuel_inlet.T = conditions.state_var.T  # K
                f.oxidizer_inlet.mdot = conditions.simul_param.mdot2  # kg/m^2/s
                f.oxidizer_inlet.X =  conditions.composition.X2
                f.oxidizer_inlet.T = conditions.state_var.T2  # K
            else:
                f.fuel_inlet.mdot = conditions.simul_param.mdot2  # kg/m^2/s
                f.fuel_inlet.X =  conditions.composition.X2
                f.fuel_inlet.T = conditions.state_var.T2  # K
                f.oxidizer_inlet.mdot = conditions.simul_param.mdot  # kg/m^2/s
                f.oxidizer_inlet.X =  conditions.composition.X
                f.oxidizer_inlet.T = conditions.state_var.T  # K
        f.flame.grid = np.array(initial_points)*width

        f.flame.set_steady_tolerances(default=tol_ss)
        f.flame.set_transient_tolerances(default=tol_ts)

        try:
            gas.transport_model = conditions.simul_param.transport_model
            f.transport_model = conditions.simul_param.transport_model
        except :
            print_("Unrecognized transport model. Proceeding with Mix model : ",mp)
            f.transport_model = "Mix"


        # Define a limit for the maximum temperature below which the flame is
        # considered as extinguished and the computation is aborted
        # This increases the speed of refinement is enabled
        temperature_limit_extinction = max(conditions.state_var.T,conditions.state_var.T2)+500  # K
        def interrupt_extinction(t):
            temperature_limit_extinction = max(conditions.state_var.T,conditions.state_var.T2)+500  # K
            if np.max(f.T) < temperature_limit_extinction:
                raise Exception('Flame extinguished')
            return 0.
        f.set_interrupt(interrupt_extinction)


        energy_enabled = True
        refine_grid = True
        loglevel = 0

        f.set_refine_criteria(ratio = 7.0, slope = 1, curve = 1)
        try:
            f.solve(loglevel, refine_grid, auto = energy_enabled)
        except:
            print_('Warning: inital computation failed',mp)
            print_('try to recompute with an alternative approach',mp)
            print_('     - first step:  pure undiluted diffusion flame (fuel/oxygen) ....',mp)
            f.oxidizer_inlet.X = 'O2:1'
            f.fuel_inlet.X     = conditions.composition.fuel + ':1'
            f.solve(loglevel, refine_grid, auto = energy_enabled)
            f.set_refine_criteria(ratio = 5.0, slope = 0.5, curve = 0.5)
            f.solve(loglevel, refine_grid)
            print_('                    flame simulation successful',mp)
            print_('     - second step: simulation of the initial counterflow flame',mp)
            f.fuel_inlet.X =  conditions.composition.X2
            f.oxidizer_inlet.X =  conditions.composition.X
            f.solve(loglevel, refine_grid, auto = energy_enabled)

#        f.set_time_step(5.0e-6,[10,20,50,80,120,150])
#        f.solve(loglevel, refine_grid, auto = energy_enabled)
        if verbose >=5 : print_("Problem solved on ["+ str(f.flame.n_points)+ "] point grid",mp)

        f.set_refine_criteria(ratio = 5.0, slope = 0.5, curve = 0.5)
        f.solve(loglevel, refine_grid)
        if verbose >=5 : print_("Problem solved on ["+ str(f.flame.n_points)+ "] point grid",mp)

        f.set_refine_criteria(ratio=3.0, slope=0.1, curve=0.2, prune=0.03)
        f.solve(loglevel, refine_grid)
        if verbose >=5 : print_("Problem solved on ["+str(f.flame.n_points)+ "] point grid",mp)

        ratio_ff = conditions.simul_param.ratio_ff
        slope_ff = conditions.simul_param.slope_ff
        curve_ff = conditions.simul_param.curve_ff
        prune_ff = conditions.simul_param.prune_ff
        f.set_refine_criteria(ratio=ratio_ff, slope=slope_ff, curve=curve_ff, prune=prune_ff)
        f.solve(loglevel, refine_grid)
        if verbose >=2 : print_("Problem solved on ["+ str(f.flame.n_points)+ "] point grid",mp)


        # Save data
        npoints = f.flame.n_points
        conc = [] ; kf = [] ; kr = [] ; r_rate = []
        for n in range(npoints):
            f.set_gas_state(n)

            if act_sp:
                conc_n = []; sa_i = 0
                for s_i in range(len(act_sp)):
                    if act_sp[s_i]:
                        conc_n.append(gas.concentrations[sa_i])
                        sa_i += 1
                    else:
                        conc_n.append(0)
                conc.append(conc_n)
            else:
                conc.append(gas.concentrations)

            if act_r:
                kf_r = []; kr_r = []; r_rate_r = []; ra_i = 0
                for r_i in range(len(act_r)):
                    if act_r[r_i]:
                        kf_r.append(gas.forward_rate_constants[ra_i])
                        kr_r.append(gas.reverse_rate_constants[ra_i])
                        r_rate_r.append(gas.net_rates_of_progress[ra_i])
                        ra_i += 1
                    else:
                        kf_r.append(0)
                        kr_r.append(0)
                        r_rate_r.append(0)
                kf.append(kf_r)
                kr.append(kr_r)
                r_rate.append(r_rate_r)
            else:
                kf.append(gas.forward_rate_constants)
                kr.append(gas.reverse_rate_constants)
                r_rate.append(gas.net_rates_of_progress)

        results = cdef.Sim_Results(conditions, gas, np.array(f.flame.grid), \
                          list(f.T), f.P, list(conc), list(kf), list(kr))
        results.K_max   = f.strain_rate('max')
        results.f       = f
        results.r_rate  = r_rate


        conditions.simul_param.pts_scatter = np.array(f.flame.grid)


        # Save results
        os.chdir(conditions.main_path)
        if 'Flame_ref_results' not in os.listdir():
            os.mkdir('Flame_ref_results')
        os.chdir('Flame_ref_results')
        if 'diff' in conditions.config: phi = 'diff'
        else: phi = '%.2f' %conditions.composition.phi
        if conditions.state_var.P>10000:
            fn = conditions.num+'cf_'+conditions.composition.fuel.replace('/','').split('(')[0]\
            +'_'+phi\
            +'_'+'%.0f'%conditions.state_var.T+'_'+'%.2f'%(conditions.state_var.P/1e5)\
            +'.xml'
        else:
            fn = conditions.num+'cf_'+conditions.composition.fuel.replace('/','').split('(')[0]\
            +'_'+phi\
            +'_'+'%.0f'%conditions.state_var.T+'_'+'%.0f'%(conditions.state_var.P)\
            +'.xml'
        f.save(fn,'ref_solution')



        if conditions.error_param.K_check:
            if verbose>3:
                print_('Extinction strain rate computation:',mp)

            # =============================================================================
            # # PART 3: STRAIN RATE LOOP
            # =============================================================================

            # Compute counterflow diffusion flames at increasing strain rates at 1 bar
            # The strain rate is assumed to increase by 25% in each step until the flame is
            # extinguished
            strain_factor = 10

            # Exponents for the initial solution variation with changes in strain rate
            # Taken from Fiala and Sattelmayer (2014) (doi:10.1155/2014/484372)
            exp_d_a = - 1. / 2.
            exp_u_a = 1. / 2.
            exp_V_a = 1.
            exp_lam_a = 2.
            exp_mdot_a = 1. / 2.

            # Restore initial solution
            if verbose<8:
                old_stdout = sys.stdout ; old_stderr = sys.stderr
                with open(os.devnull,"w") as devnull: sys.stdout=devnull;sys.stderr=devnull
            f.restore(fn, 'ref_solution')
            if verbose<8: sys.stdout = old_stdout ; sys.stderr = old_stderr      # restore console output



            strain_accuracy = conditions.error_param.strain_accuracy
            max_it = 30 ; it_n = 0 ; n = 0 ; fn = 'Kref_'+fn
            pct_var_strain=strain_accuracy+1 ;
            strain_rate = 1 ; strain_rate_prev = 0
            it_strain = 1 ; restart_sim=False   #009 46  025 41 05 43  2 46   20 49

            # Do the strain rate loop
            clock = cdef.Clock('Stretch') ; clock.start()

            while pct_var_strain>strain_accuracy:
                it_n += 1
                if it_n>max_it: break
                while np.max(f.T) > temperature_limit_extinction or restart_sim:
#                    while np.max(f.T) > temperature_limit_extinction or restart_sim:
                    # supress console output during the simulation
                    if verbose<9:
                        old_stdout = sys.stdout ; old_stderr = sys.stderr
                        with open(os.devnull,"w") as devnull: sys.stdout=devnull;sys.stderr=devnull
                    if n!=0: f.restore(fn, 'solution')
                    if verbose<9: sys.stdout = old_stdout ; sys.stderr = old_stderr      # restore console output
                    n += 1

                    # Create an initial guess based on the previous solution
                    # Update grid
                    f.flame.grid *= strain_factor ** exp_d_a
                    normalized_grid = f.grid / (f.grid[-1] - f.grid[0])
                    # Update mass fluxes
                    f.fuel_inlet.mdot *= strain_factor ** exp_mdot_a
                    f.oxidizer_inlet.mdot *= strain_factor ** exp_mdot_a
                    # Update velocities
                    f.set_profile('u', normalized_grid, f.u * strain_factor ** exp_u_a)
                    f.set_profile('V', normalized_grid, f.V * strain_factor ** exp_V_a)
                    # Update pressure curvature
                    f.set_profile('lambda', normalized_grid, f.L * strain_factor ** exp_lam_a)
                    try:
                        # Try solving the flame
                        f.solve(loglevel=0)
                        if verbose<9:
                            old_stdout = sys.stdout ; old_stderr = sys.stderr
                            with open(os.devnull,"w") as devnull: sys.stdout=devnull;sys.stderr=devnull
                        f.save(fn,'solution')
                        if verbose<9: sys.stdout = old_stdout ; sys.stderr = old_stderr      # restore console output
                        strain_factor=1+it_strain ;
                        if not restart_sim:
                            strain_rate = f.strain_rate('max') # the maximum axial strain rate
                            if verbose>3:
                                print('\rStrain rate: ' + format(strain_rate, '.2e') + ' 1/s',end='')
                            strain_rate_prev2 = strain_rate_prev
                            strain_rate_prev  = strain_rate
                            temperature_limit_extinction=.75*np.max(f.T)
                        else:
                            restart_sim=False
                    except Exception as e:
                        if e.args[0] == 'Flame extinguished':
                            if verbose>3: print('\nFlame extinguished')
                            strain_factor=1-it_strain/(it_strain+1)
                            restart_sim=True
                        else:
                            print('Error occurred while solving:', e)
                        break
                it_strain /=2
                pct_var_strain = abs((strain_rate_prev2-strain_rate_prev)/strain_rate_prev)
                if verbose>3: print_('Accuracy: '+  format(pct_var_strain*100, '.1f') + '%',mp)
            if verbose>3: print_('Extinction strain rate: ' + format(strain_rate, '.2e') + ' 1/s',mp)
            results.K_ext = strain_rate_prev
            clock.stop()
            if verbose>2: clock.display()
        os.chdir(conditions.main_path)


    if 'tp_flame' in conditions.config:

        gas         = conditions.composition.gas
        tol_ss      = conditions.simul_param.tol_ss
        tol_ts      = conditions.simul_param.tol_ts
        width       = conditions.simul_param.end_sim

        ratio = conditions.simul_param.ratio_ff
        slope = conditions.simul_param.slope_ff
        curve = conditions.simul_param.curve_ff
        prune = conditions.simul_param.prune_ff
        gas.TPX=(conditions.state_var.T,conditions.state_var.P,conditions.composition.X)

        f = ct.CounterflowTwinPremixedFlame(gas, width=width)
#        f.inlet.X = conditions.composition.X
#        f.inlet.T = conditions.state_var.T
#        f.flame.P = conditions.state_var.P
        f.reactants.mdot = conditions.simul_param.mdot  # kg/m^2/s
        f.flame.set_steady_tolerances(default=tol_ss)
        f.flame.set_transient_tolerances(default=tol_ts)
        try:
            gas.transport_model = conditions.simul_param.transport_model
            f.transport_model = conditions.simul_param.transport_model
        except :
            print_("Unrecognized transport model. Proceeding with Mix model : ",mp)
            f.transport_model = "Mix"
        f.set_refine_criteria(ratio=ratio, slope=slope, curve=curve, prune=prune)


#        f.show_solution()
        f.solve(loglevel=0, auto=True)

        # Save data
        npoints = f.flame.n_points
        conc = [] ; kf = [] ; kr = [] ; r_rate = []

        for n in range(npoints):
            f.set_gas_state(n)

            f.set_gas_state(n)

            if act_sp:
                conc_n = []; sa_i = 0
                for s_i in range(len(act_sp)):
                    if act_sp[s_i]:
                        conc_n.append(gas.concentrations[sa_i])
                        sa_i += 1
                    else:
                        conc_n.append(0)
                conc.append(conc_n)
            else:
                conc.append(gas.concentrations)

            if act_r:
                kf_r = []; kr_r = []; r_rate_r = []; ra_i = 0
                for r_i in range(len(act_r)):
                    if act_r[r_i]:
                        kf_r.append(gas.forward_rate_constants[ra_i])
                        kr_r.append(gas.reverse_rate_constants[ra_i])
                        r_rate_r.append(gas.net_rates_of_progress[ra_i])
                        ra_i += 1
                    else:
                        kf_r.append(0)
                        kr_r.append(0)
                        r_rate_r.append(0)
                kf.append(kf_r)
                kr.append(kr_r)
                r_rate.append(r_rate_r)
            else:
                kf.append(gas.forward_rate_constants)
                kr.append(gas.reverse_rate_constants)
                r_rate.append(gas.net_rates_of_progress)


        results = cdef.Sim_Results(conditions, gas, np.array(f.flame.grid), \
                          list(f.T), f.P, list(conc), list(kf), list(kr))
        results.f       = f
        results.r_rate  = r_rate

        conditions.simul_param.pts_scatter = np.array(f.flame.grid)


        # Save results
        os.chdir(conditions.main_path)
        if 'Flame_ref_results' not in os.listdir():
            os.mkdir('Flame_ref_results')
        os.chdir('Flame_ref_results')
        if conditions.state_var.P>10000:
            fn = conditions.num+'tp_'+conditions.composition.fuel.replace('/','').split('(')[0]\
            +'_'+'%.2f' %conditions.composition.phi\
            +'_'+'%.0f'%conditions.state_var.T+'_'+'%.2f'%(conditions.state_var.P/1e5)\
            +'.xml'
        else:
            fn = conditions.num+'tp_'+conditions.composition.fuel.replace('/','').split('(')[0]\
            +'_'+'%.2f' %conditions.composition.phi\
            +'_'+'%.0f'%conditions.state_var.T+'_'+'%.0f'%(conditions.state_var.P)\
            +'.xml'
        f.save(fn,'ref_solution')





    elif 'reactor' in conditions.config or 'PFR' in conditions.config:

        # import the gas model and set the initial conditions
        gas = conditions.composition.gas
        init = conditions.state_var.T, conditions.state_var.P,conditions.composition.X
        gas.TPX = init

        if 'PFR' in conditions.config:
            length      = conditions.simul_param.end_sim    # *approximate* PFR length [m]
            u_0         = conditions.simul_param.u_0        # inflow velocity [m/s]
            area        = conditions.simul_param.area       # cross-sectional area [m**2]
            n_steps     = conditions.simul_param.n_pts      # number of time steps considered for the simulation
            mass_flow_rate1 = u_0 * gas.density * area

        if 'PFR' in conditions.config and not conditions.simul_param.PFR_auto_time:
            # Plug Flow reactor simulation:
            # https://cantera.org/examples/python/reactors/pfr.py.html

            #####################################################################
            # Method 1: Lagrangian Particle Simulation
            #####################################################################
            # A Lagrangian particle is considered which travels through the PFR. Its
            # state change is computed by upwind time stepping. The PFR result is produced
            # by transforming the temporal resolution into spatial locations.
            # The spatial discretization is therefore not provided a priori but is instead
            # a result of the transformation.

            tic_0 = timer.time()

            T, P, conc, kf, kr, r_rate = [],[],[],[],[],[]


            # create a new reactor
            r1 = ct.IdealGasConstPressureReactor(gas)
            # create a reactor network for performing time integration
            sim = ct.ReactorNet([r1])

            sim.rtol = conditions.simul_param.tol_ts[0]
            sim.atol = conditions.simul_param.tol_ts[1]

            # approximate a time step to achieve a similar resolution as in the next method
            t_total = length / u_0
            dt = t_total / n_steps
            # define time, space, and other information vectors
            timeVec = (np.arange(n_steps) + 1) * dt
            z1 = np.zeros_like(timeVec)
            u1 = np.zeros_like(timeVec)

            # save data at t0
            P.append(r1.thermo.P)
            T.append(r1.T)
            conc.append(r1.thermo.concentrations.tolist())
            kf.append(gas.forward_rate_constants)
            kr.append(gas.reverse_rate_constants)

            for n1, t_i in enumerate(timeVec):
                # perform time integration
                sim.advance(t_i)
                # compute velocity and transform into space
                u1[n1] = mass_flow_rate1 / area / r1.thermo.density
                z1[n1] = z1[n1 - 1] + u1[n1] * dt
                P.append(r1.thermo.P)
                T.append(r1.T)


                if act_sp:
                    conc_n = []; sa_i = 0
                    for s_i in range(len(act_sp)):
                        if act_sp[s_i]:
                            conc_n.append(r1.thermo.concentrations.tolist()[sa_i])
                            sa_i += 1
                        else:
                            conc_n.append(0)
                    conc.append(conc_n)
                else:
                    conc.append(r1.thermo.concentrations.tolist())
                if act_r:
                    kf_r = []; kr_r = []; r_rate_r = []; ra_i = 0
                    for r_i in range(len(act_r)):
                        if act_r[r_i]:
                            kf_r.append(gas.forward_rate_constants[ra_i])
                            kr_r.append(gas.reverse_rate_constants[ra_i])
                            r_rate_r.append(gas.net_rates_of_progress[ra_i])
                            ra_i += 1
                        else:
                            kf_r.append(0)
                            kr_r.append(0)
                            r_rate_r.append(0)
                    kf.append(kf_r)
                    kr.append(kr_r)
                    r_rate.append(r_rate_r)
                else:
                    kf.append(gas.forward_rate_constants)
                    kr.append(gas.reverse_rate_constants)
                    r_rate.append(gas.net_rates_of_progress)


            tic_1 = timer.time()
            #print('Simulation time: '+ "%5.1f" %(tic_1-tic_0) + ' s')
            #####################################################################
            timeVec = np.insert(timeVec,0,0)
            z1      = np.insert(z1,0,0)


        else:

            init_points = int(conditions.simul_param.tign_nPoints)     # max number of time step
            n_pts       = conditions.simul_param.n_pts
            delta_npts  = conditions.simul_param.delta_npts
            t_max_coeff = conditions.simul_param.t_max_coeff
            Scal_ref    = conditions.simul_param.Scal_ref  # Scalar selected for the time vector computation
            grad_curv_ratio = conditions.simul_param.grad_curv_ratio
            tmax_react  = conditions.simul_param.t_max_react

        #   Auto ignition time Detection
                 # Record the scalar evolution and check ignition
                 # with H2O variation (considered if variation >20%)
            time1 = timer.time()
            print_("  Computation...",mp)
            gas.TPX = init
            if conditions.config == 'reactor_UV':
                reactor = ct.IdealGasReactor(gas)
            elif conditions.config == 'reactor_HP' or conditions.config == 'PFR':
                reactor = ct.IdealGasConstPressureReactor(gas)
            sim = ct.ReactorNet([reactor])
            sim.rtol = conditions.simul_param.tol_ts[0]
            sim.atol = conditions.simul_param.tol_ts[1]
            time_init = tmax_react
            timeVec_ignit = []
            for i in range(init_points):
                timeVec_ignit.append(time_init)
                time_init/=1.05
            timeVec_ignit.reverse()
            H2O = np.zeros(init_points-1, 'd')
            X_Scal = np.zeros(init_points-1, 'd')
            dH2O = np.zeros(init_points-2, 'd')


            if Scal_ref == "T" or Scal_ref == "T(K)" or Scal_ref == "Temp":
                X_Scal[0] = reactor.T
            else:
                X_Scal[0] = gas.X[gas.species_index(Scal_ref)]
            try:
                H2O[0] = gas.X[gas.species_index("H2O")]
            except:
                H2O[0] = 0

            for n in range(1, init_points-1):
    #            timeVec_ignit[n] = timeVec_ignit[n-1] + dt
                sim.advance(timeVec_ignit[n])

                if Scal_ref == "T" or Scal_ref == "T(K)" or Scal_ref == "Temp":
                    X_Scal[n] = reactor.T
                else:
                    X_Scal[n] = gas.X[gas.species_index(Scal_ref)]

                try:
                    H2O[n] = gas.X[gas.species_index("H2O")]
                    dH2O[n-1] = (H2O[n]-H2O[n-1]) *0.5
                except:
                    H2O[n] = 0
                    dH2O[n-1] = np.abs(X_Scal[n]-X_Scal[n-1]) *0.5

            idx_maxgrad = np.where(dH2O==max(dH2O)) #find max grad index
            tign = timeVec_ignit[idx_maxgrad[0][0]]
            if verbose >=4:
                if H2O[0]*1.2>H2O[-1] or tign>0.97*timeVec_ignit[-1]:
                    print_("    WARNING : No ignition was detected in the first "+
                          str(tmax_react)+"s" +".",mp)
                elif verbose>7:
                    print_("    - first ignition time estimation:  "+"%5.3f" %(tign*1e6)+'µs',mp)


        # Computing the time vector according to the selected scalar variation
            tmax = min(tmax_react,t_max_coeff*tign)


            # Grad and curve vectors calculation
            grad_vect_0 = [] ; curv_vect = []
            # grad calculation
            for i in range(len(X_Scal)-1):
                grad_vect_0.append((X_Scal[i+1]-X_Scal[i])/(timeVec_ignit[i+1]-timeVec_ignit[i]))
            # curv calculation
            for i in range(len(grad_vect_0)-1):
                curv_vect.append((grad_vect_0[i+1]-grad_vect_0[i])/(2*(timeVec_ignit[i+1]-timeVec_ignit[i])))
            # grad recomputed, for corresponding to the curv values
            grad_vect = []
            for i in range(len(X_Scal)-2):
                grad_vect.append((X_Scal[i+2]-X_Scal[i])/(2*(timeVec_ignit[i+1]-timeVec_ignit[i])))

            if verbose>7:
                print_("grad and curve calculation done",mp)

            time_stepping_calc=True; pt_coeff=2e5
            coeff_dtmax = n_pts/10 ; coeff_dtmin = n_pts*2
            time_stepping_calc_try=0
            while time_stepping_calc:
                # time stepping
                timeVec = [0]
                # first point => considering only grad_vect_0 (curve data not available)
                timeVec.append(tmax/((grad_vect_0[0]/np.sum([grad_vect_0[0]]))*n_pts))
                # Next points => considering only grad_vect_0 (curve data not available)
                while timeVec[-1]<tmax:
                    idx = min((np.abs(timeVec_ignit-timeVec[-1])).argmin()-2,len(grad_vect)-1)
                    dt = (tmax/(n_pts*pt_coeff))/(abs((grad_vect[idx]/np.sum(np.abs([grad_vect])))*
                         grad_curv_ratio+curv_vect[idx]/np.sum(np.abs([curv_vect]))*(1-grad_curv_ratio)))
                    if dt<(tign/coeff_dtmin):
                        timeVec.append(timeVec[-1]+tign/200)
                    elif dt>(tign/coeff_dtmax):
                        timeVec.append(timeVec[-1]+tign/5)
                    else:
                        timeVec.append(timeVec[-1]+dt)

                if len(timeVec)>n_pts-delta_npts and len(timeVec)<n_pts+delta_npts:
                    time_stepping_calc=False
                else:
                    if verbose>4:
                        print_("Bad number of time steps:"+str(len(timeVec))+". Recomputing time vector",mp)
                    pt_coeff=pt_coeff*(n_pts/len(timeVec))**2
                    time_stepping_calc=True
                    time_stepping_calc_try+=1
                    if time_stepping_calc_try>50:  time_stepping_calc=False
                    if pt_coeff<1e-5:
                        print_("Warning : larger time steps imposed for the calculation",mp)
                        pt_coeff = 2e5 ; coeff_dtmax = coeff_dtmax/4
                    if pt_coeff>1e20:
                        print_("Warning : smaller time steps imposed for the calculation",mp)
                        pt_coeff = 2e5 ; coeff_dtmin = coeff_dtmin*4
            if verbose>3:
                print_("    - time vector contains: "+str(len(timeVec))+" points",mp)



        # Data computation
            gas.TPX = init
            if conditions.config == 'reactor_UV' :
                reactor = ct.IdealGasReactor(gas)
            elif conditions.config == 'reactor_HP' or conditions.config == 'PFR':
                reactor = ct.IdealGasConstPressureReactor(gas)
            sim = ct.ReactorNet([reactor])
            T = []
            P = []
            conc = [] ; kf = [] ; kr = [] ; r_rate = []
            heat_release = [0]
            T.append(reactor.T)
            P.append(reactor.thermo.P)
            if act_sp:
                conc_n = []; sa_i = 0
                for s_i in range(len(act_sp)):
                    if act_sp[s_i]:
                        conc_n.append(reactor.thermo.concentrations.tolist()[sa_i])
                        sa_i += 1
                    else:
                        conc_n.append(0)
                conc.append(conc_n)
            else:
                conc.append(reactor.thermo.concentrations.tolist())
            if act_r:
                kf_r = []; kr_r = []; r_rate_r = []; ra_i = 0
                for r_i in range(len(act_r)):
                    if act_r[r_i]:
                        kf_r.append(gas.forward_rate_constants[ra_i])
                        kr_r.append(gas.reverse_rate_constants[ra_i])
                        r_rate_r.append(gas.net_rates_of_progress[ra_i])
                        ra_i += 1
                    else:
                        kf_r.append(0)
                        kr_r.append(0)
                        r_rate_r.append(0)
                kf.append(kf_r)
                kr.append(kr_r)
                r_rate.append(r_rate_r)
            else:
                kf.append(gas.forward_rate_constants)
                kr.append(gas.reverse_rate_constants)
                r_rate.append(gas.net_rates_of_progress)

#            conc.append(reactor.thermo.concentrations.tolist())
#            kf.append(gas.forward_rate_constants)
#            kr.append(gas.reverse_rate_constants)
#            r_rate.append(gas.net_rates_of_progress)
            fuel = conditions.composition.fuel.split('/')[0].split('(')[0]
            target_ign_idx = gas.species_index(fuel)
            target_ign_ = []
            grad_fuel = []
            z1 = np.zeros_like(timeVec)
            u1 = np.zeros_like(timeVec)

            for n in range(1,len(timeVec)):
                sim.advance(timeVec[n])
                T.append(reactor.T)
                P.append(reactor.thermo.P)

                if act_sp:
                    conc_n = []; sa_i = 0
                    for s_i in range(len(act_sp)):
                        if act_sp[s_i]:
                            conc_n.append(reactor.thermo.concentrations.tolist()[sa_i])
                            sa_i += 1
                        else:
                            conc_n.append(0)
                    conc.append(conc_n)
                else:
                    conc.append(reactor.thermo.concentrations.tolist())
                if act_r:
                    kf_r = []; kr_r = []; r_rate_r = []; ra_i = 0
                    for r_i in range(len(act_r)):
                        if act_r[r_i]:
                            kf_r.append(gas.forward_rate_constants[ra_i])
                            kr_r.append(gas.reverse_rate_constants[ra_i])
                            r_rate_r.append(gas.net_rates_of_progress[ra_i])
                            ra_i += 1
                        else:
                            kf_r.append(0)
                            kr_r.append(0)
                            r_rate_r.append(0)
                    kf.append(kf_r)
                    kr.append(kr_r)
                    r_rate.append(r_rate_r)
                else:
                    kf.append(gas.forward_rate_constants)
                    kr.append(gas.reverse_rate_constants)
                    r_rate.append(gas.net_rates_of_progress)

#                conc.append(reactor.thermo.concentrations.tolist())
#                kf.append(gas.forward_rate_constants)
#                kr.append(gas.reverse_rate_constants)
#                r_rate.append(gas.net_rates_of_progress)

                target_ign_.append(gas.X[target_ign_idx])

                if 'PFR' in conditions.config:
                    u1[n] = mass_flow_rate1 / area / reactor.thermo.density
                    dt = timeVec[n]-timeVec[n-1]
                    z1[n] = z1[n - 1] + u1[n] * dt

                # heat release calculation
                hr="ok"
                try:
                    heat_release.append(-np.dot(gas.net_rates_of_progress,\
                                               gas.delta_enthalpy))
                except:
                    if n==1:
                        print_("warning: heat release calculation issues",mp)
                        hr="no_heat_release"
            # ignition time calculation
               # 1- based on heat release (default)
            if hr!= "no_heat_release":
                ign_time_hr = timeVec[heat_release.index(max(heat_release))]
            else: ign_time_hr=False
            ign_time_hr = timeVec[heat_release.index(max(heat_release))]
               # 2- based on fuel gradients (if heat relase calculation troubles)
            for t in range(len(timeVec)-3):
                grad_fuel.append((target_ign_[t+2]-target_ign_[t])/(timeVec[t+2]-timeVec[t]))
            ign_time_sp = timeVec[grad_fuel.index(max(grad_fuel))+2]
            if 'PFR' in conditions.config:
                ign_dist_sp = z1[grad_fuel.index(max(grad_fuel))+2]

            time2 = timer.time()
            if verbose>3:
                if 'PFR' in conditions.config:
                    print_("    - ignition distance :  "+"%5.3f" %(ign_dist_sp*1e3)+'mm',mp)
                else:
                    print_("    - ignition time :  "+"%5.3f" %(ign_time_sp*1e6)+'µs',mp)
                print_("      Time to compute reference data : "+"%5.2f" %(time2-time1)+'s',mp)


        results = cdef.Sim_Results(conditions, gas, np.array(timeVec), list(T), \
                             list(P), list(conc), list(kf), list(kr))
        results.r_rate = list(r_rate)

        conditions.simul_param.pts_scatter=np.array(timeVec)

        if 'PFR' in conditions.config:
            results.z1 = list(z1)
        else:
            results.ign_time_hr = ign_time_hr
            results.ign_time_sp = ign_time_sp



    elif 'JSR' in conditions.config:
        # adapted from https://www.cantera.org/examples/jupyter/reactors/stirred_reactor.ipynb.html

        # =============================================================================
        # Initial conditions
        # =============================================================================
        gas  = conditions.composition.gas
        T_list = list(conditions.simul_param.pts_scatter)
        gas.TPX = T_list[0], conditions.state_var.P,conditions.composition.X
        time1 = timer.time()
        # Reactor parameters
        residenceTime = conditions.simul_param.end_sim #s
        reactorVolume = 30.5*(1e-2)**3 #m3
        # Instrument parameters
        pressureValveCoefficient = .01
        # Simulation termination criterion
        maxSimulationTime = 50 # seconds

        # =============================================================================
        # Reactor model
        # =============================================================================
        #Initialize the stirred reactor and connect all peripherals
        fuelAirMixtureTank = ct.Reservoir(gas)
        exhaust = ct.Reservoir(gas)
        stirredReactor = ct.IdealGasReactor(gas, energy='off', volume=reactorVolume)
        massFlowController = ct.MassFlowController(upstream=fuelAirMixtureTank,
                                                   downstream=stirredReactor,
                                                   mdot=stirredReactor.mass/residenceTime)
        pressureRegulator = ct.Valve(upstream=stirredReactor,
                                     downstream=exhaust,
                                     K=pressureValveCoefficient)
        reactorNetwork = ct.ReactorNet([stirredReactor])
        # now compile a list of all variables for which we will store data
        columnNames = [stirredReactor.component_name(item) for item in range(stirredReactor.n_vars)]
        columnNames = ['pressure'] + columnNames
        # use the above list to create a DataFrame
        timeHistory = pd.DataFrame(columns=columnNames)

        # =============================================================================
        # Modeling data
        # =============================================================================
        T = [] ; P = [] ; conc = []; kf = [] ; kr = [] ; r_rate = []

        # Create a data frame to store values for the above points
        tempDependence = pd.DataFrame(columns=timeHistory.columns)
        tempDependence.index.name = 'Temperature'
        inletConcentrations = conditions.composition.X
        concentrations      = inletConcentrations
        tic_sim = timer.time()
        for temperature in T_list:
            #Re-initialize the gas
            reactorTemperature = temperature #Kelvin
            gas.TPX = reactorTemperature,conditions.state_var.P,inletConcentrations
            # Re-initialize the dataframe used to hold values
            timeHistory = pd.DataFrame(columns=columnNames)
            # Re-initialize all the reactors, reservoirs, etc
            fuelAirMixtureTank = ct.Reservoir(gas)
            exhaust            = ct.Reservoir(gas)
            # use concentrations from the previous iteration to speed up convergence
            gas.TPX = reactorTemperature,conditions.state_var.P,concentrations

            stirredReactor = ct.IdealGasReactor(gas, energy='off', volume=reactorVolume)
            massFlowController = ct.MassFlowController(upstream=fuelAirMixtureTank,
                                                       downstream=stirredReactor,
                                                       mdot=stirredReactor.mass/residenceTime)
            pressureRegulator  = ct.Valve(upstream=stirredReactor,
                                         downstream=exhaust,
                                         K=pressureValveCoefficient)
            reactorNetwork     = ct.ReactorNet([stirredReactor])
            tic = timer.time()

            # Re-run the isothermal simulations
            t = 0
            while t < maxSimulationTime:
                t = reactorNetwork.step()
                if verbose>7:
                    write_t = '\r'+str(t) + '                '
                    print(write_t ,end="")
            if verbose>7:print(' ')
            state = np.hstack([stirredReactor.thermo.P,
                               stirredReactor.mass,
                               stirredReactor.volume,
                               stirredReactor.T,
                               stirredReactor.thermo.X])
            diff_P = conditions.state_var.P - stirredReactor.thermo.P
            if(abs(diff_P/conditions.state_var.P) > 0.01):
                print_("WARNING: Non-trivial pressure rise in the reactor. You may adjust pressureValveCoefficient value in ./packages/Computation.py",mp)
            concentrations = stirredReactor.thermo.X
            # Store the result in the dataframe that indexes by temperature
            tempDependence.loc[temperature] = state
            toc = timer.time()

            if verbose>4:
                print_('Simulation at T={}K took {:3.2f}s to compute'.format(temperature, toc-tic),mp)
            else:
                print("\r  Computation at T="+str(stirredReactor.T)+'K ...  ',end='')

            # results
            T.append(stirredReactor.T)
            P.append(stirredReactor.thermo.P)
            if act_sp:
                conc_n = []; sa_i = 0
                for s_i in range(len(act_sp)):
                    if act_sp[s_i]:
                        conc_n.append(stirredReactor.thermo.concentrations.tolist()[sa_i])
                        sa_i += 1
                    else:
                        conc_n.append(0)
                conc.append(conc_n)
            else:
                conc.append(stirredReactor.thermo.concentrations.tolist())
            if act_r:
                kf_r = []; kr_r = []; r_rate_r = []; ra_i = 0
                for r_i in range(len(act_r)):
                    if act_r[r_i]:
                        kf_r.append(gas.forward_rate_constants[ra_i])
                        kr_r.append(gas.reverse_rate_constants[ra_i])
                        r_rate_r.append(gas.net_rates_of_progress[ra_i])
                        ra_i += 1
                    else:
                        kf_r.append(0)
                        kr_r.append(0)
                        r_rate_r.append(0)
                kf.append(kf_r)
                kr.append(kr_r)
                r_rate.append(r_rate_r)
            else:
                kf.append(gas.forward_rate_constants)
                kr.append(gas.reverse_rate_constants)
                r_rate.append(gas.net_rates_of_progress)

#            conc.append(stirredReactor.thermo.concentrations.tolist())
#            kf.append(gas.forward_rate_constants)
#            kr.append(gas.reverse_rate_constants)
#            r_rate.append(gas.net_rates_of_progress)

        toc_sim = timer.time()


        results = cdef.Sim_Results(conditions, gas, np.array(T_list), list(T), \
                             list(P), list(conc), list(kf), list(kr))
        results.r_rate = r_rate
        conditions.simul_param.time_jsr_ref_lim = toc_sim-tic_sim



    return results,conditions




def red_computation(conditions, gas_red, act_sp,act_r):

    verbose    = conditions.simul_param.verbose
    mp = conditions.main_path

    if 'free_flame' in conditions.config or 'burner_flame' in conditions.config:

        # Main variables
        gas_ref    = conditions.composition.gas_ref
        tol_ss     = conditions.simul_param.tol_ss
        tol_ts     = conditions.simul_param.tol_ts

        if 'free' in conditions.config:
            grid      = list(conditions.simul_param.pts_scatter)
            f         = ct.FreeFlame(gas_red, grid)
        elif 'burner' in conditions.config:
            grid             = list(conditions.simul_param.pts_scatter_i)
            f                = ct.BurnerFlame(gas_red,grid)
            f.burner.mdot    = conditions.simul_param.mdot    #mass flow rate per unit area [kg/m^2/s]
            f.energy_enabled = False
            f.flame.set_fixed_temp_profile(np.array(grid)/max(grid),conditions.simul_param.T_profile)

        try:
            gas_red.transport_model = conditions.simul_param.transport_model
            f.transport_model = conditions.simul_param.transport_model
        except :
            print_("Unrecognized transport model. Proceeding with Mix model : ",mp)
            f.transport_model = "Mix"
        f.flame.set_steady_tolerances(default=tol_ss)
        f.flame.set_transient_tolerances(default=tol_ts)


        # Get saved results
        os.chdir(conditions.main_path)
        os.chdir('Flame_ref_results')
        if conditions.state_var.P>10000:
            fn = conditions.num+'ff_'+conditions.composition.fuel.replace('/','').split('(')[0]\
            +'_'+'%.2f' %conditions.composition.phi\
            +'_'+'%.0f'%conditions.state_var.T+'_'+'%.2f'%(conditions.state_var.P/1e5)\
            +'.xml'
        else:
            fn = conditions.num+'ff_'+conditions.composition.fuel.replace('/','').split('(')[0]\
            +'_'+'%.2f' %conditions.composition.phi\
            +'_'+'%.0f'%conditions.state_var.T+'_'+'%.0f'%(conditions.state_var.P)\
            +'.xml'


        # --------------------     Simulation     ----------------------

        # supress console output during the simulation
        if verbose<9:
            old_stdout = sys.stdout ; old_stderr = sys.stderr
            with open(os.devnull, "w") as devnull:
                sys.stdout = devnull ; sys.stderr = devnull
        f.restore(fn, 'ref_solution')

        os.chdir(conditions.main_path)

        try:
            f.solve(auto = False, loglevel = 0, refine_grid = False)
            simul_success = True

            # restore console output
            if verbose<9: sys.stdout = old_stdout ; sys.stderr = old_stderr
            if verbose >= 6:
                print_("     Problem solved on ["+str(f.flame.n_points)+"] point grid",mp)
        except:
            simul_success = False
            # restore console output
            if verbose<9: sys.stdout = old_stdout ; sys.stderr = old_stderr
            if verbose >= 3:
                print("\n     WARNING: No solution found\n",mp)
        grid = f.flame.grid

        # ------------------     end of simulation     --------------------

        # saving spec concentrations at
        conc=[] ; kf=[] ; kr=[] ; r_rate=[] ; npoints = f.flame.n_points
        for n in range(npoints):
            f.set_gas_state(n)
            # saving spec concentrations
            conc_x = []
            for sp in range(len(act_sp)):
                if act_sp[sp] and simul_success:
                    sp_ind=gas_red.species_index(gas_ref.species_name(sp))
                    conc_x.append(gas_red.concentrations[sp_ind])
                else:
                    conc_x.append(0)
            conc.append(conc_x)
            # saving kf and kr values
            kf_x=[];kr_x=[];r_rate_x=[];r_red=0
            for r in range(len(act_r)):
                if act_r[r]:
                        kf_x.append(gas_red.forward_rate_constants[r_red])
                        kr_x.append(gas_red.reverse_rate_constants[r_red])
                        r_rate_x.append(gas_red.net_rates_of_progress[r_red])
                        r_red+=1
                else:
                    kf_x.append(0); kr_x.append(0); r_rate_x.append(0)
            kf.append(kf_x); kr.append(kr_x) ; r_rate.append(r_rate_x)

        results = cdef.Sim_Results(conditions, gas_red, list(f.flame.grid), \
                          list(f.T), f.P, list(conc), list(kf), list(kr))
        results.Sl = f.u[0]
        results.f  = f
        results.r_rate = r_rate

    if 'tp_flame' in conditions.config:

        # Main variables
        gas_ref   = conditions.composition.gas_ref


        f = ct.CounterflowTwinPremixedFlame(gas_red)

        # Get saved results
        os.chdir(conditions.main_path)
        os.chdir('Flame_ref_results')
        if conditions.state_var.P>10000:
            fn = conditions.num+'tp_'+conditions.composition.fuel.replace('/','').split('(')[0]\
            +'_'+'%.2f' %conditions.composition.phi\
            +'_'+'%.0f'%conditions.state_var.T+'_'+'%.2f'%(conditions.state_var.P/1e5)\
            +'.xml'
        else:
            fn = conditions.num+'tp_'+conditions.composition.fuel.replace('/','').split('(')[0]\
            +'_'+'%.2f' %conditions.composition.phi\
            +'_'+'%.0f'%conditions.state_var.T+'_'+'%.0f'%(conditions.state_var.P)\
            +'.xml'


        # supress console output during the simulation
        if verbose<9:
            old_stdout = sys.stdout ; old_stderr = sys.stderr
            with open(os.devnull, "w") as devnull:
                sys.stdout = devnull ; sys.stderr = devnull

        f.restore(fn, 'ref_solution')

        # Initialize and solve
        try:
            f.solve(auto = False, loglevel = 0, refine_grid = False)
            simul_success = True

            # restore console output
            if verbose<9: sys.stdout = old_stdout ; sys.stderr = old_stderr
            if verbose >= 6:
                print_("     Problem solved on ["+str(f.flame.n_points)+"] point grid",mp)
        except:
            simul_success = False
            # restore console output
            if verbose<9: sys.stdout = old_stdout ; sys.stderr = old_stderr
            if verbose >= 3:
                print("\n     WARNING: No solution found\n",mp)

        grid = f.flame.grid


        # saving spec concentrations at
        conc=[] ; kf=[] ; kr=[] ; r_rate=[] ; npoints = f.flame.n_points
        for n in range(npoints):
            f.set_gas_state(n)
            # saving spec concentrations
            conc_x = []
            for sp in range(len(act_sp)):
                if act_sp[sp] and simul_success:
                    sp_ind=gas_red.species_index(gas_ref.species_name(sp))
                    conc_x.append(gas_red.concentrations[sp_ind])
                else:
                    conc_x.append(0)
            conc.append(conc_x)
            # saving kf and kr values
            kf_x=[];kr_x=[];r_rate_x=[];r_red=0
            for r in range(len(act_r)):
                if act_r[r]:
                        kf_x.append(gas_red.forward_rate_constants[r_red])
                        kr_x.append(gas_red.reverse_rate_constants[r_red])
                        r_rate_x.append(gas_red.net_rates_of_progress[r_red])
                        r_red+=1
                else:
                    kf_x.append(0); kr_x.append(0); r_rate_x.append(0)
            kf.append(kf_x); kr.append(kr_x) ; r_rate.append(r_rate_x)

        results = cdef.Sim_Results(conditions, gas_red, list(f.flame.grid), \
                          list(f.T), f.P, list(conc), list(kf), list(kr))
        results.f       = f
        results.r_rate  = r_rate


        os.chdir(conditions.main_path)




    if 'diff_flame' in conditions.config or 'pp_flame' in conditions.config :

        # Main variables
        gas_ref  = conditions.composition.gas_ref
        width    = conditions.simul_param.end_sim


        # =============================================================================
        # # PART 1: INITIALIZATION
        # =============================================================================

        f = ct.CounterflowDiffusionFlame(gas_red, width=width)

        # Get saved results
        os.chdir(conditions.main_path)
        os.chdir('Flame_ref_results')
        if 'diff' in conditions.config: phi = 'diff'
        else: phi = '%.2f' %conditions.composition.phi
        if conditions.state_var.P>10000:
            fn = conditions.num+'cf_'+conditions.composition.fuel.replace('/','').split('(')[0]\
            +'_'+phi\
            +'_'+'%.0f'%conditions.state_var.T+'_'+'%.2f'%(conditions.state_var.P/1e5)\
            +'.xml'
        else:
            fn = conditions.num+'cf_'+conditions.composition.fuel.replace('/','').split('(')[0]\
            +'_'+phi\
            +'_'+'%.0f'%conditions.state_var.T+'_'+'%.0f'%(conditions.state_var.P)\
            +'.xml'

#        # Define a limit for the maximum temperature below which the flame is
#        # considered as extinguished and the computation is aborted
#        # This increases the speed of refinement is enabled
        temperature_limit_extinction = max(conditions.state_var.T,conditions.state_var.T2)+500  # K
        def interrupt_extinction(t):
            temperature_limit_extinction = max(conditions.state_var.T,conditions.state_var.T2)+500  # K
            if np.max(f.T) < temperature_limit_extinction:
                raise Exception('Flame extinguished')
            return 0.
        f.set_interrupt(interrupt_extinction)

        # supress console output during the simulation
        if verbose<9:
            old_stdout = sys.stdout ; old_stderr = sys.stderr
            with open(os.devnull, "w") as devnull:
                sys.stdout = devnull ; sys.stderr = devnull

        f.restore(fn, 'ref_solution')

        # Initialize and solve
        try:
            f.solve(auto = False, loglevel = 0, refine_grid = False)
            simul_success = True

            # restore console output
            if verbose<9: sys.stdout = old_stdout ; sys.stderr = old_stderr
            if verbose >= 6:
                print_("     Problem solved on ["+str(f.flame.n_points)+"] point grid",mp)
        except:
            simul_success = False
            # restore console output
            if verbose<9: sys.stdout = old_stdout ; sys.stderr = old_stderr
            if verbose >= 3:
                print("\n     WARNING: No solution found\n",mp)

        grid = f.flame.grid


        # saving spec concentrations at
        conc = [] ; kf = [] ; kr = [] ; r_rate=[] ; npoints = f.flame.n_points
        for n in range(npoints):
            f.set_gas_state(n)
            # saving spec concentrations
            conc_x = []
            for sp in range(len(act_sp)):
                if act_sp[sp] and simul_success:
                    sp_ind=gas_red.species_index(gas_ref.species_name(sp))
                    conc_x.append(gas_red.concentrations[sp_ind])
                else:
                    conc_x.append(0)
            conc.append(conc_x)

            # saving kf and kr values
            kf_x=[];kr_x=[];r_rate_x=[];r_red=0
            for r in range(len(act_r)):
                if act_r[r]:
                        kf_x.append(gas_red.forward_rate_constants[r_red])
                        kr_x.append(gas_red.reverse_rate_constants[r_red])
                        r_rate_x.append(gas_red.net_rates_of_progress[r_red])
                        r_red+=1
                else:
                    kf_x.append(0); kr_x.append(0); r_rate_x.append(0)
            kf.append(kf_x); kr.append(kr_x) ; r_rate.append(r_rate_x)
#
#            # saving kf and kr values
#            kf_x=[];kr_x=[];r_red=0
#            for r in range(len(act_r)):
#                if act_r[r]:
#                        kf_x.append(gas_red.forward_rate_constants[r_red])
#                        kr_x.append(gas_red.reverse_rate_constants[r_red])
#                        r_red+=1
#                else:
#                    kf_x.append(0); kr_x.append(0)
#            kf.append(kf_x); kr.append(kr_x)

        results = cdef.Sim_Results(conditions, gas_red, list(f.flame.grid), \
                          list(f.T), f.P, list(conc), list(kf), list(kr))
        results.K_max   = f.strain_rate('max')
        results.f       = f
        results.r_rate  = r_rate


        if conditions.error_param.K_check and simul_success:
            if verbose>4:
                print_('Extinction strain rate computation:',mp)
            # =============================================================================
            # # PART 3: STRAIN RATE LOOP
            # =============================================================================

            # Compute counterflow diffusion flames at increasing strain rates at 1 bar
            # The strain rate is assumed to increase by 25% in each step until the flame is
            # extinguished
            strain_factor = 0.95

            # Exponents for the initial solution variation with changes in strain rate
            # Taken from Fiala and Sattelmayer (2014) (doi:10.1155/2014/484372)
            exp_d_a = - 1. / 2.
            exp_u_a = 1. / 2.
            exp_V_a = 1.
            exp_lam_a = 2.
            exp_mdot_a = 1. / 2.

            # Counter to identify the loop
            strain_accuracy = conditions.error_param.strain_accuracy
            max_it = 30
            if conditions.simul_param.par_ind: Ki = conditions.simul_param.par_ind
            else:                              Ki = ''
            it_n = 0 ; fnK = 'Kref_'+fn ; fnKi = 'K'+Ki+'_'+fn
            pct_var_strain=strain_accuracy+1 ;
            strain_rate = 1 ; strain_rate_prev = 0
            it_strain = .05 ; restart_sim = False   #009 46  025 41 05 43  2 46   20 49


            if simul_success:
                first_it = True
                while pct_var_strain>strain_accuracy:
                    it_n += 1
                    if it_n>max_it: break
                    it_sn = 0
                    while np.max(f.T) > temperature_limit_extinction or restart_sim:
                        it_sn +=1
                        if it_sn>40:
                            it_n = max_it*2
                            break
                        # supress console output during the simulation
                        if verbose<9:
                            old_stdout = sys.stdout ; old_stderr = sys.stderr
                            with open(os.devnull,"w") as devnull: sys.stdout=devnull;sys.stderr=devnull
                        if first_it: f.restore(fnK, 'solution')
                        else:        f.restore(fnKi, 'solution')
                        if verbose<9: sys.stdout = old_stdout ; sys.stderr = old_stderr      # restore console output

                        # Create an initial guess based on the previous solution
                        # Update grid
                        f.flame.grid    *= strain_factor ** exp_d_a
                        normalized_grid  = f.grid / (f.grid[-1] - f.grid[0])
                        # Update mass fluxes
                        f.fuel_inlet.mdot     *= strain_factor ** exp_mdot_a
                        f.oxidizer_inlet.mdot *= strain_factor ** exp_mdot_a
                        # Update velocities
                        f.set_profile('u', normalized_grid, f.u * strain_factor ** exp_u_a)
                        f.set_profile('V', normalized_grid, f.V * strain_factor ** exp_V_a)
                        # Update pressure curvature
                        f.set_profile('lambda', normalized_grid, f.L * strain_factor ** exp_lam_a)
                        try:
                            # Try solving the flame
                            f.solve(loglevel=0)
                            if verbose<9:
                                old_stdout = sys.stdout ; old_stderr = sys.stderr
                                with open(os.devnull,"w") as devnull: sys.stdout=devnull;sys.stderr=devnull
                            f.save(fnKi,'solution')
                            if verbose<9: sys.stdout = old_stdout ; sys.stderr = old_stderr      # restore console output
                            strain_factor=1+it_strain ;
                            first_it = False
                            if not restart_sim:
                                strain_rate = f.strain_rate('max') # the maximum axial strain rate
                                if verbose>5:
                                    print('\rStrain rate:' + format(strain_rate, '.2e') + ' 1/s',end='')
                                strain_rate_prev2 = strain_rate_prev
                                strain_rate_prev  = strain_rate
                                temperature_limit_extinction=.75*np.max(f.T)
                            else:
                                restart_sim=False

                        except Exception as e:
                            if e.args[0] == 'Flame extinguished':
                                if verbose>5: print('\nFlame extinguished')
                                strain_factor=1-it_strain/(it_strain+1)
                                restart_sim=True
                                # if extinction at first iteration
                                if first_it:
                                    it_strain *= 4
                                    strain_rate_prev2 = 1 ; strain_rate_prev=0.1
                                else:
                                    it_strain = abs(it_strain)
                            else:
                                print('Error occurred while solving:', e)
                            break
                    it_strain /=2
                    pct_var_strain = abs((strain_rate_prev2-strain_rate_prev)/strain_rate_prev)
                    if verbose>6: print_('Accuracy: '+ format(pct_var_strain*100, '.1f') + '%',mp)
                if verbose>5: print_('Extinction strain rate: ' + format(strain_rate, '.2e') + ' 1/s',mp)
                results.K_ext = strain_rate_prev

        os.chdir(conditions.main_path)


    if 'reactor' in conditions.config or 'PFR' in conditions.config:

        X_red = conditions.composition.X
        gas_red.TPX = conditions.state_var.T, conditions.state_var.P,X_red
        gas_ref     = conditions.composition.gas_ref

        timeVec = conditions.simul_param.pts_scatter

        if 'PFR' in conditions.config:
        # Plug Flow reactor simulation:
        # https://cantera.org/examples/python/reactors/pfr.py.html
        #
        # A Lagrangian particle is considered which travels through the PFR. Its
        # state change is computed by upwind time stepping. The PFR result is produced
        # by transforming the temporal resolution into spatial locations.
        # The spatial discretization is therefore not provided a priori but is instead
        # a result of the transformation.
            length      = conditions.simul_param.end_sim    # *approximate* PFR length [m]
            u_0         = conditions.simul_param.u_0        # inflow velocity [m/s]
            area        = conditions.simul_param.area       # cross-sectional area [m**2]
            n_steps     = conditions.simul_param.n_pts      # number of time steps considered for the simulation
            mass_flow_rate1 = u_0 * gas_red.density * area


        if conditions.config == 'reactor_UV' :
            reactor = ct.IdealGasReactor(gas_red)
        elif conditions.config == 'reactor_HP' or conditions.config == 'PFR':
            reactor = ct.IdealGasConstPressureReactor(gas_red)

        sim = ct.ReactorNet([reactor])
        sim.rtol = conditions.simul_param.tol_ts[0]
        sim.atol = conditions.simul_param.tol_ts[1]


        T = [] ; P = [] ; conc = [] ; kf = [] ; kr = [] ; r_rate = []

        # saving data at t=0
        heat_release = [0]
        T.append(reactor.T)
        P.append(reactor.thermo.P)

        # saving spec concentrations at t=0
        conc_t=[]
        for sp in range(len(act_sp)):
            if act_sp[sp]:
                sp_ind=gas_red.species_index(gas_ref.species_name(sp))
                conc_t.append(reactor.thermo.concentrations[sp_ind])
            else:
                conc_t.append(0)
        conc.append(conc_t)

        # saving kf and kr values at t=0
        kf_t=[];kr_t=[];r_rate_t=[];r_red=0
        for r in range(len(act_r)):
            if act_r[r]:
                    kf_t.append(gas_red.forward_rate_constants[r_red])
                    kr_t.append(gas_red.reverse_rate_constants[r_red])
                    r_rate_t.append(gas_red.net_rates_of_progress[r_red])
                    r_red+=1
            else:
                kf_t.append(0); kr_t.append(0); r_rate_t.append(0)
        kf.append(kf_t); kr.append(kr_t) ; r_rate.append(r_rate_t)

        # definition of fuel for ignition delay detection
        fuel = X_red.split(':')[0]
        target_ign_idx = gas_red.species_index(fuel)
        target_ign_ = []
        grad_fuel = []
        if 'PFR' in conditions.config:
            # define space and velocity vectors
            z1 = np.zeros_like(timeVec)
            u1 = np.zeros_like(timeVec)

        for n in range(1,len(timeVec)):
            time = timeVec[n]

            # supress console output during the simulation
            if verbose<9:
                old_stdout = sys.stdout ; old_stderr = sys.stderr
                with open(os.devnull, "w") as devnull:
                    sys.stdout = devnull ; sys.stderr = devnull
            try:
                sim.advance(time)
                simul_success = True
                # restore console output
                if verbose<9: sys.stdout = old_stdout ; sys.stderr = old_stderr
            except:
                simul_success = False
                # restore console output
                if verbose<9: sys.stdout = old_stdout ; sys.stderr = old_stderr
                if verbose >= 3:
                    print("\n     WARNING: No solution found\n",mp)

            # compute velocity and transform into space
            if 'PFR' in conditions.config:
                u1[n] = mass_flow_rate1 / area / reactor.thermo.density
                dt = timeVec[n]-timeVec[n-1]
                z1[n] = z1[n - 1] + u1[n] * dt


            T.append(reactor.T)
            P.append(reactor.thermo.P)
            # saving spec concentrations
            conc_t=[]
            for sp in range(len(act_sp)):
                if act_sp[sp] and simul_success:
                    sp_ind=gas_red.species_index(gas_ref.species_name(sp))
                    conc_t.append(reactor.thermo.concentrations[sp_ind])
                else:
                    conc_t.append(0)
            conc.append(conc_t)
            target_ign_.append(gas_red.X[target_ign_idx])

            # saving kf and kr values
            kf_t=[];kr_t=[];r_rate_t=[];r_red=0
            for r in range(len(act_r)):
                if act_r[r]:
                        kf_t.append(gas_red.forward_rate_constants[r_red])
                        kr_t.append(gas_red.reverse_rate_constants[r_red])
                        r_rate_t.append(gas_red.net_rates_of_progress[r_red])
                        r_red+=1
                else:
                    kf_t.append(0); kr_t.append(0); r_rate_t.append(0)
            kf.append(kf_t); kr.append(kr_t) ; r_rate.append(r_rate_t)

#            # saving kf and kr values
#            kf_t=[];kr_t=[];r_red=0
#            for r in range(len(act_r)):
#                if act_r[r]:
#                        kf_t.append(gas_red.forward_rate_constants[r_red])
#                        kr_t.append(gas_red.reverse_rate_constants[r_red])
#                        r_red+=1
#                else:
#                    kf_t.append(0); kr_t.append(0)
#            kf.append(kf_t); kr.append(kr_t)
            # heat release calculation
            hr="ok"
            try:
                heat_release.append(-np.dot(gas_red.net_rates_of_progress,\
                                           gas_red.delta_enthalpy))
            except:
                if n==1:
                    print_("warning: heat release calculation issues",mp)
                    hr="no_heat_release"

        # define time and position at the beginning (-> 0):
        if 'PFR' in conditions.config and not conditions.simul_param.PFR_auto_time:
            timeVec = np.insert(timeVec,0,0)
            z1      = np.insert(z1,0,0)


        # ignition time calculation
           # 1- based on heat release (default)
        if hr!= "no_heat_release":
            ign_time_hr = timeVec[heat_release.index(max(heat_release))]
        else: ign_time_hr=False
           # 2- based on fuel (if heat relase calculation troubles)
        for t in range(len(timeVec)-3):
            try:
                grad_fuel.append((target_ign_[t+2]-target_ign_[t])/(timeVec[t+2]-timeVec[t]))
            except:
                print('t:'+str(t))
                print('len timeVec:'+str(len(timeVec)))
                print('target_ing:'+str(len(target_ign_)))

        ign_time_sp = timeVec[grad_fuel.index(max(grad_fuel))+2]


        results = cdef.Sim_Results(conditions, gas_red, list(timeVec), list(T),\
                             list(P), list(conc), list(kf), list(kr))

        if 'PFR' in conditions.config:
            results.z1 = list(z1)
        else:
            results.ign_time_hr = ign_time_hr
            results.ign_time_sp = ign_time_sp

        results.r_rate      = list(r_rate)


    elif 'JSR' in conditions.config:
        # adapted from https://www.cantera.org/examples/jupyter/reactors/stirred_reactor.ipynb.html

        # =============================================================================
        # Initial conditions
        # =============================================================================
        X_red = conditions.composition.X

        gas_ref  = conditions.composition.gas_ref
        T_list = list(conditions.simul_param.pts_scatter)
        gas_red.TPX = T_list[0], conditions.state_var.P,X_red
        # Reactor parameters
        residenceTime = conditions.simul_param.end_sim #s
        reactorVolume = 30.5*(1e-2)**3 #m3
        # Instrument parameters
        pressureValveCoefficient = .01
        # Simulation termination criterion
        maxSimulationTime = 50 # seconds

        # =============================================================================
        # Reactor model
        # =============================================================================
        #Initialize the stirred reactor and connect all peripherals
        fuelAirMixtureTank = ct.Reservoir(gas_red)
        exhaust = ct.Reservoir(gas_red)
        stirredReactor = ct.IdealGasReactor(gas_red, energy='off', volume=reactorVolume)
        massFlowController = ct.MassFlowController(upstream=fuelAirMixtureTank,
                                                   downstream=stirredReactor,
                                                   mdot=stirredReactor.mass/residenceTime)
        pressureRegulator = ct.Valve(upstream=stirredReactor,
                                     downstream=exhaust,
                                     K=pressureValveCoefficient)
        reactorNetwork = ct.ReactorNet([stirredReactor])
        # now compile a list of all variables for which we will store data
        columnNames = [stirredReactor.component_name(item) for item in range(stirredReactor.n_vars)]
        columnNames = ['pressure'] + columnNames
        # use the above list to create a DataFrame
        timeHistory = pd.DataFrame(columns=columnNames)

        # =============================================================================
        # Modeling data
        # =============================================================================
        T = [] ; P = [] ; conc = []; kf = [] ; kr = [] ; r_rate = []

        # Create a data frame to store values for the above points
        tempDependence = pd.DataFrame(columns=timeHistory.columns)
        tempDependence.index.name = 'Temperature'
        inletConcentrations = X_red
        concentrations      = inletConcentrations
        tic_sim = timer.time()
        for temperature in T_list:
            conc_T = []
            #Re-initialize the gas
            reactorTemperature = temperature #Kelvin
            reactorVolume = 30.5*(1e-2)**3 #m3
            gas_red.TPX = reactorTemperature,conditions.state_var.P,inletConcentrations
            # Re-initialize the dataframe used to hold values
            timeHistory = pd.DataFrame(columns=columnNames)
            # Re-initialize all the reactors, reservoirs, etc
            fuelAirMixtureTank = ct.Reservoir(gas_red)
            exhaust = ct.Reservoir(gas_red)
            # use concentrations from the previous iteration to speed up convergence
            gas_red.TPX = reactorTemperature,conditions.state_var.P,concentrations

            stirredReactor = ct.IdealGasReactor(gas_red, energy='off', volume=reactorVolume)
            massFlowController = ct.MassFlowController(upstream=fuelAirMixtureTank,
                                                       downstream=stirredReactor,
                                                       mdot=stirredReactor.mass/residenceTime)
            pressureRegulator  = ct.Valve(upstream=stirredReactor,
                                         downstream=exhaust,
                                         K=pressureValveCoefficient)
            reactorNetwork     = ct.ReactorNet([stirredReactor])

            # supress console output during the simulation
            if verbose<9:
                old_stdout = sys.stdout ; old_stderr = sys.stderr
                with open(os.devnull, "w") as devnull:
                    sys.stdout = devnull ; sys.stderr = devnull

            try:
                # Re-run the isothermal simulations
                t = 0
                simul_success = True
                while t < maxSimulationTime:
                    toc_sim = timer.time()
                    if (toc_sim-tic_sim) > (100*conditions.simul_param.time_jsr_ref_lim):  # if calculation time becomes too long
                        simul_success = False
                        t = maxSimulationTime + 1
                        if verbose >= 3:
                            print("\n     WARNING: No solution found\n")
                            print(str(toc_sim-tic_sim))
                    else:
                        t = reactorNetwork.step()

                # restore console output
                if verbose<9: sys.stdout = old_stdout ; sys.stderr = old_stderr
            except:
                simul_success = False
                # restore console output
                if verbose<9: sys.stdout = old_stdout ; sys.stderr = old_stderr
                if verbose >= 3:
                    print("\n     WARNING: No solution found\n",mp)

            state = np.hstack([stirredReactor.thermo.P,
                               stirredReactor.mass,
                               stirredReactor.volume,
                               stirredReactor.T,
                               stirredReactor.thermo.X])
            if simul_success:
                concentrations = stirredReactor.thermo.X
            # Store the result in the dataframe that indexes by temperature
            tempDependence.loc[temperature] = state

            # results
            T.append(stirredReactor.T)
            P.append(stirredReactor.thermo.P)
            # saving spec concentrations
            for sp in range(len(act_sp)):
                if act_sp[sp] and simul_success:
                    sp_ind=gas_red.species_index(gas_ref.species_name(sp))
                    conc_T.append(stirredReactor.thermo.concentrations[sp_ind])
                else:
                    conc_T.append(0)
            conc.append(conc_T)

            # saving kf and kr values
            kf_t=[];kr_t=[];r_rate_t=[];r_red=0
            for r in range(len(act_r)):
                if act_r[r]:
                        kf_t.append(gas_red.forward_rate_constants[r_red])
                        kr_t.append(gas_red.reverse_rate_constants[r_red])
                        r_rate_t.append(gas_red.net_rates_of_progress[r_red])
                        r_red+=1
                else:
                    kf_t.append(0); kr_t.append(0); r_rate_t.append(0)
            kf.append(kf_t); kr.append(kr_t) ; r_rate.append(r_rate_t)

#            # saving kf and kr values
#            kf_t=[];kr_t=[];r_red=0
#            for r in range(len(act_r)):
#                if act_r[r]:
#                        kf_t.append(gas_red.forward_rate_constants[r_red])
#                        kr_t.append(gas_red.reverse_rate_constants[r_red])
#                        r_red+=1
#                else:
#                    kf_t.append(0); kr_t.append(0)
#            kf.append(kf_t); kr.append(kr_t)

        results = cdef.Sim_Results(conditions, gas_red, list(T_list), list(T), \
                             list(P), list(conc), list(kf), list(kr))
        results.r_rate      = list(r_rate)


#    elif 'PFR' in conditions.config:
#
#        # Plug Flow reactor simulation:
#        # https://cantera.org/examples/python/reactors/pfr.py.html
#
#        #####################################################################
#        # Method 1: Lagrangian Particle Simulation
#        #####################################################################
#        # A Lagrangian particle is considered which travels through the PFR. Its
#        # state change is computed by upwind time stepping. The PFR result is produced
#        # by transforming the temporal resolution into spatial locations.
#        # The spatial discretization is therefore not provided a priori but is instead
#        # a result of the transformation.
#
#        tic_0 = timer.time()
#
#        T = [] ; P = [] ; conc = []; kf = [] ; kr = [] ; r_rate = []
#
#        length      = conditions.simul_param.end_sim    # 1.5e-7  # *approximate* PFR length [m]
#        u_0         = conditions.simul_param.u_0        # .006  # inflow velocity [m/s]
#        area        = conditions.simul_param.area       # 1.e-4  # cross-sectional area [m**2]
#        n_steps     = conditions.simul_param.n_pts      #2000 # number of time steps considered for the simulation
#
#        X_red = conditions.composition.X
#        gas_red.TPX = conditions.state_var.T, conditions.state_var.P,X_red
#        gas_ref     = conditions.composition.gas
#
#        timeVec = conditions.simul_param.pts_scatter
#
#
#        mass_flow_rate1 = u_0 * gas_red.density * area
#
#        # create a new reactor
#        r1 = ct.IdealGasConstPressureReactor(gas_red)
#        # create a reactor network for performing time integration
#        sim = ct.ReactorNet([r1])
#
#        sim.rtol = conditions.simul_param.tol_ts[0]
#        sim.atol = conditions.simul_param.tol_ts[1]
#
#
#        # approximate a time step to achieve a similar resolution as in the next method
#        t_total = length / u_0
#        dt = t_total / n_steps
#        # define time, space, and other information vectors
#        z1 = np.zeros_like(timeVec)
#        u1 = np.zeros_like(timeVec)
#
#
#
#        # saving data at t=0
#        T.append(r1.T)
#        P.append(r1.thermo.P)
#        # conc
#        conc_t=[]
#        for sp in range(len(act_sp)):
#            if act_sp[sp]:
#                sp_ind=gas_red.species_index(gas_ref.species_name(sp))
#                conc_t.append(r1.thermo.concentrations[sp_ind])
#            else:
#                conc_t.append(0)
#        conc.append(conc_t)
#
#
#        # saving kf and kr values
#        kf_t=[];kr_t=[];r_rate_t=[];r_red=0
#        for r in range(len(act_r)):
#            if act_r[r]:
#                    kf_t.append(gas_red.forward_rate_constants[r_red])
#                    kr_t.append(gas_red.reverse_rate_constants[r_red])
#                    r_rate_t.append(gas_red.net_rates_of_progress[r_red])
#                    r_red+=1
#            else:
#                kf_t.append(0); kr_t.append(0); r_rate_t.append(0)
#        kf.append(kf_t); kr.append(kr_t) ; r_rate.append(r_rate_t)
##        # reaction rate coefficients
##        kf_t=[];kr_t=[];r_red=0
##        for r in range(len(act_r)):
##            if act_r[r]:
##                    kf_t.append(gas_red.forward_rate_constants[r_red])
##                    kr_t.append(gas_red.reverse_rate_constants[r_red])
##                    r_red+=1
##            else:
##                kf_t.append(0); kr_t.append(0)
##        kf.append(kf_t); kr.append(kr_t)
#
#
#        for n1, t_i in enumerate(timeVec):
#
#            # perform time integration
#
#            # supress console output during the simulation
#            if verbose<9:
#                old_stdout = sys.stdout ; old_stderr = sys.stderr
#                with open(os.devnull, "w") as devnull:
#                    sys.stdout = devnull ; sys.stderr = devnull
#            try:
#                sim.advance(t_i)
#                simul_success = True
#                # restore console output
#                if verbose<9: sys.stdout = old_stdout ; sys.stderr = old_stderr
#            except:
#                simul_success = False
#                # restore console output
#                if verbose<9: sys.stdout = old_stdout ; sys.stderr = old_stderr
#                if verbose >= 3:
#                    print("\n     WARNING: No solution found\n",mp)
#
#            # compute velocity and transform into space
#            u1[n1] = mass_flow_rate1 / area / r1.thermo.density
#            z1[n1] = z1[n1 - 1] + u1[n1] * dt
##            states1.append(r1.thermo.state)
#
#            # save simulation data
#            T.append(r1.T)
#            P.append(r1.thermo.P)
#            #conc
#            conc_t=[]
#            for sp in range(len(act_sp)):
#                if act_sp[sp] and simul_success:
#                    sp_ind=gas_red.species_index(gas_ref.species_name(sp))
#                    conc_t.append(r1.thermo.concentrations[sp_ind])
#                else:
#                    conc_t.append(0)
#            conc.append(conc_t)
#
#
#            # reaction rate coefficients
#            kf_t=[];kr_t=[];r_red=0
#            for r in range(len(act_r)):
#                if act_r[r]:
#                        kf_t.append(gas_red.forward_rate_constants[r_red])
#                        kr_t.append(gas_red.reverse_rate_constants[r_red])
#                        r_red+=1
#                else:
#                    kf_t.append(0); kr_t.append(0)
#            kf.append(kf_t); kr.append(kr_t)
#
#
#        tic_1 = timer.time()
#        timeVec = np.insert(timeVec,0,0)
#        z1      = np.insert(z1,0,0)
#
#
#        results = cdef.Sim_Results(conditions, gas_red, timeVec, list(T),\
#                             list(P), list(conc), list(kf), list(kr))
#        results.z1 = list(z1)
#        conditions.simul_param.pts_scatter=np.array(timeVec)


    return results


def transportModel(gas_ref, transport_model):
    try:
        gas_ref.transport_model = transport_model
    except :
        print("Warning, unrecognized or unsupported transport model. Proceeding with the default model : ",gas_ref.transport_model)
    return gas_ref

