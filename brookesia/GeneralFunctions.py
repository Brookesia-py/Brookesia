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

#==============================================================================
#   Packages
#==============================================================================
import numpy                       as np
import os                          as os
import sys                         as sys
import cantera                     as ct
import time                        as timer
import brookesia.Computation      as comp
import brookesia.SA               as sa
import brookesia.DRG              as drg
import brookesia.GeneticAlgorithm as ga
import brookesia.PSO              as pso
import brookesia.Class_def        as cdef
import gc
import copy
import datetime
import traceback

from  brookesia.Class_def   import print_
from scipy.interpolate       import interp1d
from shutil import copyfile



#==============================================================================
#   Functions
#==============================================================================

def main_redopt_algo(filename,WD_path):
#    WD_path
    print('\n'*100)
    os.chdir(WD_path)

    #==============================================================================
    #%%      Read reduction conditions
    #==============================================================================

    conditions_list, red_data_list, ref_results_list = get_reduction_parameters(filename)


    #==============================================================================
    #%%      New folder creation
    #==============================================================================

    mf = conditions_list[0].main_folder
    mp = conditions_list[0].main_path
    mech      = conditions_list[0].mech
    mech_prev_red = conditions_list[0].mech_prev_red
    if mf:
        if not os.path.exists(mf):  os.mkdir(mf)
    if not os.path.exists(mp):  os.mkdir(mp)
    try: copyfile('_kinetic_mech/'+mech,mp+'/'+mech)
    except:  copyfile(mech,mp+'/'+mech.split('/')[-1])
    if mech_prev_red:
        try:    copyfile('_kinetic_mech/'+mech_prev_red,mp+'/'+mech_prev_red)
        except: copyfile(mech_prev_red,mp+'/'+mech_prev_red.split('/')[-1])
    try:    copyfile('_conditions_input/'+filename,mp+'/Conditions_redopt.inp')
    except: copyfile(filename,mp+'/'+filename.split('/')[-1])
    if os.path.isfile('_uncertainties/uncertainties.csv'):
        copyfile('_uncertainties/uncertainties.csv',mp+'/uncertainties.csv')
    os.chdir(conditions_list[0].main_path)

    print_('Computed with :\n * Cantera  '+ct.__version__+'\n * Brookesia 1.5\n\n',mp)

    try:
        #==============================================================================
        #%%      Mechanism informations
        #==============================================================================

        verbose = conditions_list[0].simul_param.verbose
        mech_data = cdef.Mech_data(mech,verbose)

        if mech_prev_red:
            mech_data_prev = cdef.Mech_data(mech_prev_red,verbose)
            mech_data.compare_red(mech_data_prev,mp)
#            mech_data.get_kin_data(mech_prev_red)

        #==============================================================================
        #%%      Information displaying
        #==============================================================================

        gas_ref = conditions_list[0].composition.gas
        ns_ref = gas_ref.n_species
        nr_ref = gas_ref.n_reactions

        if verbose >=1:
            print_('\n\n\nmechanism: ' + mech,mp)
            print_(str(ns_ref)+' species',mp)
            print_(str(nr_ref)+' reactions',mp)


        if verbose >=1:
            print_('Number of evaluated conditions: '+ '%.0f' %len(conditions_list),mp)

        #==============================================================================
        #%%     Reference data computation
        #==============================================================================

        if not ref_results_list:  # if no external data provided
            ref_results_list = []
            for i in range(len(conditions_list)):
                results,conditions = computation_reference(conditions_list[i],verbose)
                conditions_list[i] = conditions
                ref_results_list.append(results)

    #    ref_results_list[0].plotData_opt(['CH4'])

        #==============================================================================
        #%%    Reduction procedure
        ##=============================================================================
        os.chdir(mp)

        ns_red,nr_red,red_results_list,red_data_list\
            =reduction(conditions_list,ref_results_list,red_data_list,mech_data)


        #============================================================================
        #%%    Display results
        ##===========================================================================
        _reduction = False
        for op in range(len(red_data_list)):
            if red_data_list[op][0].reduction_operator != 'NULL': _reduction = True
        if _reduction:
            print_('\n\n Species reduction:  '+ '%4.0f' %ns_ref+ '  -> '+ '%5.0f' %ns_red,mp)
            print_('\n\n Reaction reduction: '+ '%4.0f' %nr_ref+ '  -> '+ '%5.0f' %nr_red,mp)


        #============================================================================
        #%%    Saving results in csv file if no reduction
        ##===========================================================================
        if len(red_data_list)==0:
            ref_results_list[0].write_mech_info('')
            for l in range(len(ref_results_list)):
                ref_results_list[l].write_case_data('Reference','',l+1)
        os.chdir(conditions_list[0].main_path)

    except:
        print_(traceback.format_exc(),mp)
    os.chdir(conditions_list[0].main_path)


def computation_reference(conditions, verbose=1, act_sp=False, act_r=False):
    mp = conditions.main_path
    if verbose >=1:
        print_('\n\n============================ ',mp)
        print_('Configuration: '+ conditions.config,mp)
        if conditions.config=='JSR':
            print_('T   = '+ str(conditions.simul_param.pts_scatter[0])+ '-'+ \
                            str(conditions.simul_param.pts_scatter[-1])+ ' K',mp)
        else:
            print_('T   = '+str(conditions.state_var.T)+ ' K',mp)
        print_('P   = '+str(conditions.state_var.P)+ ' Pa',mp)
        if conditions.config=='pp_flame':
            print_('phis = '+str(conditions.composition.phi)+' / '+str(conditions.composition.phi2),mp)
        else:
            print_('phi = '+str(conditions.composition.phi),mp)
        print_('============================ \n\n',mp)

        results,conditions = comp.ref_computation(conditions,verbose,act_sp,act_r)

        return results,conditions


def reduction(conditions_list,ref_results_list,red_data_list,mech_data):

    verbose = conditions_list[0].simul_param.verbose
    mp = conditions_list[0].main_path

    if conditions_list[0].import_data:
        ref_results_list_ext_data = copy_ref_results(ref_results_list)
        conditions_list, mech_results_list, ref_results_list \
                  = input_results_treatment(conditions_list,ref_results_list,mech_data)

    if 'mech_results_list' in locals():
        red_results_list = copy_ref_results(mech_results_list)
    else:
        red_results_list = copy_ref_results(ref_results_list)

    if conditions_list[0].mech_prev_red:
        gas_prev_mech = ct.Solution(conditions_list[0].mech_prev_red)
        for i in range(len(red_results_list)):
            red_results_list[i].gas = gas_prev_mech


    red_errors_list = []


    for op in range(len(red_data_list)):

        red_results_list[0].write_mech_info(op)
        red_data_list[op][0].op = op

        red_method = red_data_list[op][0].reduction_operator
        print_ ('\n\n ----------------------------------------------------',mp)
        if red_method == 'DRG_sp':
            print_('            Species DRG reduction \n',mp)
            folderCreation(conditions_list[0].mech,'','','_','__DRG_sp',verbose)
        if red_method == 'DRGEP_sp':
            print_('            Species DRGEP reduction \n',mp)
            folderCreation(conditions_list[0].mech,'','','_','__DRGEP_sp',verbose)
        elif red_method == 'DRG_r':
            print_('            Reaction DRG reduction \n',mp)
            folderCreation(conditions_list[0].mech,'','','_','__DRG_r',verbose)
        elif red_method == 'SAR_sp':
            print_('            Species SA reduction \n',mp)
            folderCreation(conditions_list[0].mech,'','','_','__SAR_sp',verbose)
        elif red_method == 'SARGEP_sp':
            print_('            Species SARGEP reduction \n',mp)
            folderCreation(conditions_list[0].mech,'','','_','__SARGEP_sp',verbose)
        elif red_method == 'SAR_r':
            print_('            Reaction SA reduction \n',mp)
            folderCreation(conditions_list[0].mech,'','','_','__SAR_r',verbose)
        elif red_method == 'CSP':
            folderCreation(conditions_list[0].mech,'','','_','__CSP',verbose)
            if len(mech_data.spec.CSP_radicals)>0: mech_data.spec.CSP_radicals = []
        elif red_method == 'NULL' and not conditions_list[0].import_data:
            mech_results_list = ref_results_list


        #Clock
        clock = cdef.Clock(red_method) ; clock.start()

        simulation = 0 ;

        # main variables
        n_tspecies = len(red_data_list[op][0].tspc)
        tsp_idx    = red_data_list[op][0].targetSpeciesIdx
        tspc       = red_data_list[op][0].tspc
        mech_data.spec.activ_p  = [False]*len(mech_data.spec.name)
        mech_data.react.activ_p = [False]*len(mech_data.react.number)


        if red_method != 'NULL':
            for i in range(len(conditions_list)):

                # -------------------------------
                # Reduction loop

                # Simulation condition
                conditions   = conditions_list[i]
                # Reference mech results
                ref_results = ref_results_list[i]
                # Reduced mech results
                red_results = red_results_list[i]
                # Reduction parameters
                red_data    = red_data_list[op][i]

                # loop variables
                try_acc = [5]*n_tspecies ; try_n=1 # Reduction acceleration coeff
                sp_try = np.zeros(n_tspecies)
                sp_inter = np.zeros(n_tspecies)
                T_try=0;ig_try=0;Sl_try=0;K_try=0
                active_sp_pm = mech_data.spec.activ_m
                active_r_pm  = mech_data.react.activ_m
                red_results_loop = red_results
                cross_red_error = False    # interactions between species reduction
                stop_reduction=False

                if verbose >=1:
                    print_('\n\n============================ ',mp)
                    print_('Configuration: '+conditions.config,mp)
                    if conditions.config=='JSR':
                        print_('T   = '+ str(conditions.simul_param.pts_scatter[0])+ '-'+ \
                                        str(conditions.simul_param.pts_scatter[-1])+ ' K',mp)
                    else:
                        print_('T   = '+ str(conditions.state_var.T)+ ' K',mp)
                    print_('P   = '+ str(conditions.state_var.P)+ ' Pa',mp)
                    print_('phi = '+ str(conditions.composition.phi),mp)
                    print_('============================ \n\n',mp)
                    timer.sleep(0.1)


                if 'DRG' in red_method:
                    # Interaction coefficients calculation  (useful also for DRG_r to assess the interactions between target species (TSI loop))
                    red_data = drg.dic(red_data,mech_data,red_results)
                    # Reaction interaction coefficients calculation
                    if '_r' in red_method:
                        if i==0: red_data.red_op.first_step_DRG_r = True
                        else:    red_data.red_op.first_step_DRG_r = False
                        # add main coupled pecies coupled to targets for important reaction analysis
                        if red_data.red_op.new_targets_4_DRG_r == tsp_idx:
                            _eps = np.array(red_data.red_op.eps_init)*1e5
                            _a = drg.graphSearch(conditions,red_data,mech_data,_eps)
                            del _a
                        red_data = drg.ric(red_data, mech_data, red_results)
                    if red_data.optim!='False':
                        if red_data.optim_param.optim_on_meth:
                            red_data = drg.ric(red_data, mech_data, red_results)
                elif 'SA' in red_method:
                    if conditions.config=='JSR':
                        print_('Warning: non sensitivity analysis for JSR configuration',mp)
                        break
                    # Sensitivity coefficients calculation
                    red_data = sa.sensitivities_computation_SA(red_data,mech_data,\
                                                            red_results)
                    if red_data.optim!='False':
                        red_data.optim_param.target_r.append(red_data.red_op.sensi_r)
                elif 'CSP' in red_method:
                    if 'flame' in conditions.config:
                        print_('Warning: no CSP analysis for flame configurations',mp)
                        break
                    red_data, CSP_radicals = csp.CSP_analysis(red_data,mech_data,red_results)
                    mech_data.spec.CSP_radicals.append(CSP_radicals)

                # eps definition
                eps_init   = copy.deepcopy(red_data.red_op.eps_init)
                delta_eps      = copy.deepcopy(red_data.red_op.delta_eps_init)
                max_eps_config = copy.deepcopy(red_data.red_op.eps_max)
                eps = copy.deepcopy(max_eps_config)


                # define eps
                if conditions.error_param.T_check and conditions.config!='JSR':
                    for sp in conditions.error_param.sp_T:
                        idx = tspc.index(sp)
                        if sp not in tspc[0:red_data.n_tspc]:
                            eps[idx] = eps_init[idx]
                if conditions.error_param.Sl_check and 'free_flame' in conditions.config:
                    for sp in conditions.error_param.sp_Sl:
                        idx = tspc.index(sp)
                        if sp not in tspc[0:red_data.n_tspc]:
                            eps[idx] = eps_init[idx]
                if conditions.error_param.K_check \
                   and ('diff_flame' in conditions.config or 'pp_flame' in conditions.config):
                    for sp in conditions.error_param.sp_K:
                        idx = tspc.index(sp)
                        if sp not in tspc[0:red_data.n_tspc]:
                            eps[idx] = eps_init[idx]
                if conditions.error_param.ig_check and 'reactor' in conditions.config:
                    for sp in conditions.error_param.sp_ig:
                        idx = tspc.index(sp)
                        if sp not in tspc[0:red_data.n_tspc]:
                            eps[idx] = eps_init[idx]
                for idx in range(red_data.n_tspc):
                    eps[idx] = eps_init[idx]


    #            under_tol_red = True
                act_sp_prev = active_sp_pm
                act_r_prev  = active_r_pm
                sp_prev = active_sp_pm
                r_prev = active_r_pm
                First_try_T_error  = False
                First_try_Sl_error = False
                First_try_ig_error = False
                First_try_K_error  = False
                First_try_sp_error = [False]*red_data.n_tspc
                one_simul_success  = True
                T_error     = False
                Sl_error    = False
                ig_error    = False
                K_error     = False
                sp_error    = [False]*red_data.n_tspc
                stop_T_red  = False
                stop_Sl_red = False
                stop_ig_red = False
                stop_K_red  = False
                cross_red_error = False
                eps_stop = [True]*n_tspecies   # stop reduction if errors>tolerance
                for sp in range(red_data.n_tspc): eps_stop[sp]=False
                for sp in conditions.error_param.sp_T:
                    if conditions.error_param.T_check and not 'JSR' in conditions.config:
                        if sp in tspc: idx = tspc.index(sp) ; eps_stop[idx] = False
                for sp in conditions.error_param.sp_Sl:
                    if conditions.error_param.Sl_check and 'free_flame'  in conditions.config:
                        if sp in tspc: idx = tspc.index(sp) ; eps_stop[idx] = False
                for sp in conditions.error_param.sp_K:
                    if conditions.error_param.K_check and ('diff_flame' in conditions.config or 'pp_flame' in conditions.config):
                        if sp in tspc: idx = tspc.index(sp) ; eps_stop[idx] = False
                for sp in conditions.error_param.sp_ig:
                    if conditions.error_param.ig_check and 'reactor' in conditions.config:
                        if sp in tspc: idx = tspc.index(sp) ; eps_stop[idx] = False

                eps_cre  = [False]*n_tspecies   # stop reduction if cross reduction errors
                eps_prev = list(eps)

                while not stop_reduction:
                    print_('\nStep:'+str(try_n)+' ('+red_method+')',mp); try_n+=1
                    sp_try[:] += 1
                    T_try+=1;ig_try+=1;Sl_try+=1;K_try+=1
                    sp_inter_flag = [False]*len(tspc)

                    # Active species and reactions definition
                    if 'DRG' in red_method:
                        # species selection
                        if red_method == 'DRGEP_sp':
                            active_sp_pm   = drg.graphSearch_DRGEP(conditions,red_data,mech_data,eps)
                        elif red_method == 'DRG_sp':
                            active_sp_pm   = drg.graphSearch(conditions,red_data,mech_data,eps)
                        # reaction selection
                        if '_sp' in red_method:
                            active_r_pm, active_sp_pm = drg.reactionWithdrawal\
                                (mech_data,red_data,red_method,eps,conditions,active_sp_pm)
                        else:
                            active_r_pm, active_sp_pm = drg.reactionWithdrawal\
                                (mech_data,red_data,red_method,eps,conditions)
                    elif 'SA' in red_method:
                        active_sp_pm = sa.speciesWithdrawal\
                            (conditions,red_data,red_method,mech_data,eps)
                        active_r_pm,active_sp_pm = sa.reactionWithdrawal\
                            (conditions,mech_data,active_sp_pm,red_data,red_method,eps)

                    # testing new mech
                    if active_sp_pm   != sp_prev \
                    or active_r_pm    != r_prev  \
                    or try_n==2:
                        if verbose>3:
                            print_("  "+str(active_sp_pm.count(True))+\
                                  " species, "+ str(active_r_pm.count(True))+\
                                  " reactions remaining",mp)
                        os.chdir(conditions_list[0].main_path+'/__'+red_method)
                        mech_data.write_new_mech("temp.cti",active_sp_pm,active_r_pm)
                        # --------------------------------------------------------------------------------
                        # interpretation of the new mech
                        if verbose<8:
                        # supress console output during the simulation
                            old_stdout = sys.stdout ; old_stderr = sys.stderr
                            with open(os.devnull, "w") as devnull: sys.stdout = devnull ; sys.stderr = devnull

                        red_data.red_op.gas = ct.Solution('temp.cti')

                        if verbose<8:
                        # restore console output
                            sys.stdout = old_stdout ; sys.stderr = old_stderr
                        # --------------------------------------------------------------------------------

                        red_results_loop = comp.red_computation(\
                                    conditions, red_data.red_op.gas,\
                                    active_sp_pm,active_r_pm)
                        errors = cdef.Errors(conditions,ref_results,\
                                      red_results_loop,red_data,red_data.red_op)
                        # error calculation
#                        if errors.under_tol:
#                            os.system('cp temp.cti previous_mech.cti')
                        if conditions.simul_param.show_plots:
                            plotData(tspc[0:red_data.n_tspc],ref_results,red_results_loop)

                    elif  active_sp_pm == mech_data.spec.activ_p       \
                      and active_r_pm  == mech_data.react.activ_p      \
                      and errors.under_tol:
                        if verbose > 2:
                            print_('     All remaining species and reactions are to be kept',mp)
                        break

                    else:
                        if verbose > 5:
                            print_("     Mechanism already evaluated. No further reduction was achieved",mp)

                    # =============================================================
                    #                  Check Errors

                    # Error threshold on T
                    if conditions.error_param.T_check and conditions.config!='JSR':
                        if errors.under_tol_T : # under tol -> continuing reduction
                            T_error = False ; First_try_T_error = False
                        else:                   # above tol
                            if T_try<=1:
                                First_try_T_error=True ; T_error=True ; T_try=0
                            else:
                                First_try_T_error = False ; T_error = True
                                if stop_T_red:
                                    cross_red_error = True    # check cross reduction induced error
                                    if red_data.red_op.inter_sp_inter:
                                        # Define the most interactive spec among the target species
    #                                    OIC_sp = copy.deepcopy(red_data.red_op.OIC_sp)
                                        if 'SA' in red_method:
                                            OIC_sp = copy.deepcopy(red_data.red_op.sensi_T)
                                        elif 'DRG' in red_method:
                                            OIC_sp = copy.deepcopy(red_data.red_op.OIC_sp)
                                        OIC_sp_tsp = np.zeros(len(tsp_idx))
                                        id_sp=0
                                        for sp_tsp in range(len(mech_data.spec.activ_m)):
                                            if sp_tsp in tsp_idx:
                                                if 'SA' in red_method:
                                                    OIC_sp_tsp[id_sp] = max(OIC_sp_tsp[id_sp],OIC_sp[sp_tsp])+1e-6
                                                elif 'DRG' in red_method:
                                                    OIC_sp_tsp[id_sp] = max(OIC_sp_tsp[id_sp],OIC_sp[idx][sp_tsp])+1e-6
                                                id_sp+=1
                                        for sp_i in range(int(sp_inter[idx])):
                                            OIC_sp_tsp[list(OIC_sp_tsp).index(max(OIC_sp_tsp))]=0
                                        best_ICsp = list(OIC_sp_tsp).index(max(OIC_sp_tsp))
                                        eps_cre[best_ICsp]=True
                                        sp_inter[idx] += 1
                                stop_T_red = True

                    # Error threshold on Sl
                    if conditions.error_param.Sl_check and 'free_flame' in conditions.config:
                        if errors.under_tol_Sl: # under tol -> continuing reduction
                            Sl_error = False ; First_try_Sl_error = False
                        else:                   # above tol
                            if Sl_try<=1:
                                First_try_Sl_error=True ; Sl_error=True ; Sl_try=0
                            else:
                                First_try_Sl_error = False ; Sl_error = True
                                if stop_Sl_red:
                                    cross_red_error = True    # check cross reduction induced error
                                    if red_data.red_op.inter_sp_inter:
                                        # Define the most interactive spec among the target species
                                        if 'SA' in red_method:
                                            OIC_sp = copy.deepcopy(red_data.red_op.sensi_Sl)
                                        elif 'DRG' in red_method:
                                            OIC_sp = copy.deepcopy(red_data.red_op.OIC_sp)
                                        OIC_sp_tsp = np.zeros(len(tsp_idx))
                                        id_sp=0
                                        idx = tspc.index(conditions.error_param.sp_Sl[0])
                                        for sp_tsp in range(len(mech_data.spec.activ_m)):
                                            if sp_tsp in tsp_idx:
                                                if 'SA' in red_method:
                                                    OIC_sp_tsp[id_sp] = max(OIC_sp_tsp[id_sp],OIC_sp[sp_tsp])+1e-6
                                                elif 'DRG' in red_method:
                                                    OIC_sp_tsp[id_sp] = max(OIC_sp_tsp[id_sp],OIC_sp[idx][sp_tsp])+1e-6
                                                id_sp+=1
                                        for sp_i in range(int(sp_inter[idx])):
                                            OIC_sp_tsp[list(OIC_sp_tsp).index(max(OIC_sp_tsp))]=0
                                        best_ICsp = list(OIC_sp_tsp).index(max(OIC_sp_tsp))
                                        eps_cre[best_ICsp]=True
                                        if not sp_inter_flag[idx]: sp_inter[idx]+=1; sp_inter_flag[idx]=True
                                stop_Sl_red = True

                    # Error threshold on K
                    if conditions.error_param.K_check \
                    and ('diff_flame' in conditions.config or 'pp_flame' in conditions.config):
                        if errors.under_tol_K: # under tol -> continuing reduction
                            K_error = False ; First_try_K_error = False
                        else:                   # above tol
                            if K_try<=1:
                                First_try_K_error=True ; K_error=True ; K_try=0
                            else:
                                First_try_K_error = False ; K_error = True
                                if stop_K_red:
                                    cross_red_error = True    # check cross reduction induced error
                                    if red_data.red_op.inter_sp_inter:
                                        # Define the most interactive spec among the target species
                                        if 'SA' in red_method:
                                            OIC_sp = copy.deepcopy(red_data.red_op.sensi_sp)
                                        elif 'DRG' in red_method:
                                            OIC_sp = copy.deepcopy(red_data.red_op.OIC_sp)
                                        OIC_sp_tsp = np.zeros(len(tsp_idx))
                                        id_sp=0
                                        idx = tspc.index(conditions.error_param.sp_K[0])
                                        for sp_tsp in range(len(mech_data.spec.activ_m)):
                                            if sp_tsp in tsp_idx:
                                                OIC_sp_tsp[id_sp] = max(OIC_sp_tsp[id_sp],OIC_sp[idx][sp_tsp])+1e-6
                                                id_sp+=1
                                        for sp_i in range(int(sp_inter[idx])):
                                            OIC_sp_tsp[list(OIC_sp_tsp).index(max(OIC_sp_tsp))]=0
                                        best_ICsp = list(OIC_sp_tsp).index(max(OIC_sp_tsp))
                                        eps_cre[best_ICsp]=True
                                        if not sp_inter_flag[idx]: sp_inter[idx]+=1; sp_inter_flag[idx]=True
                                stop_K_red = True


                    # Error threshold on Ignition time
                    if conditions.error_param.ig_check and 'reactor' in conditions.config:
                        if errors.under_tol_ig: # under tol -> continuing reduction
                            ig_error = False ; First_try_ig_error = False
                        else:                   # above tol
                            if ig_try<=1:
                                First_try_ig_error=True ; ig_error=True ; ig_try=0
                            else :
                                First_try_ig_error = False ; ig_error = True
                                if stop_ig_red:
                                    cross_red_error = True    # check cross reduction induced error
                                    if red_data.red_op.inter_sp_inter:
                                        # Define the most interactive spec among the target species
                                        if 'SARGEP' in red_method:
                                            OIC_sp = []
                                            for sp_oic in range(len(active_r_pm)):
                                                if sp_oic in tsp_idx:
                                                    OIC_sp.append(copy.deepcopy(red_data.red_op.sensi_sp[sp_oic]))
                                        elif 'SA' in red_method:
                                            OIC_sp = copy.deepcopy(red_data.red_op.sensi_sp)
                                        elif 'DRG' in red_method:
                                            OIC_sp = copy.deepcopy(red_data.red_op.OIC_sp)
    #                                    OIC_sp = copy.deepcopy(red_data.drg.OIC_sp)
                                        OIC_sp_tsp = np.zeros(len(tsp_idx))
                                        id_sp=0
                                        idx = tspc.index(conditions.error_param.sp_ig[0])
                                        for sp_tsp in range(len(mech_data.spec.activ_m)):
                                            if sp_tsp in tsp_idx:
                                                OIC_sp_tsp[id_sp] = max(OIC_sp_tsp[id_sp],OIC_sp[idx][sp_tsp])+1e-6
                                                id_sp+=1
                                        for sp_i in range(int(sp_inter[idx])):
                                            OIC_sp_tsp[list(OIC_sp_tsp).index(max(OIC_sp_tsp))]=0
                                        best_ICsp = list(OIC_sp_tsp).index(max(OIC_sp_tsp))
                                        eps_cre[best_ICsp]=True
                                        if not sp_inter_flag[idx]: sp_inter[idx]+=1; sp_inter_flag[idx]=True
                                stop_ig_red = True

                    # Error threshold on Species
                    for idx in range(red_data.n_tspc):
    #                    sp = tspc[idx]
                        if errors.under_tol_s[idx]: # under tol -> continuing reduction
                            sp_error[idx] = False ; First_try_sp_error[idx]=False
                        else:                       # above tol
                            if sp_try[idx]<=1:
                                First_try_sp_error[idx]=True ; sp_error[idx]=True
                                sp_try[idx] = 0
                            else:
                                First_try_sp_error[idx]=False ; sp_error[idx]=True
                                if eps_stop[idx]:
                                    cross_red_error = True    # check cross reduction induced error
                                    if red_data.red_op.inter_sp_inter:
                                        # Define the most interactive spec among the target species
                                        if 'SA' in red_method:
                                            OIC_sp = copy.deepcopy(red_data.red_op.sensi_sp)
                                        elif 'DRG' in red_method:
                                            OIC_sp = copy.deepcopy(red_data.red_op.OIC_sp)
                                        OIC_sp_tsp = np.zeros(len(tsp_idx))
                                        id_sp=0
                                        for sp_tsp in range(len(mech_data.spec.activ_m)):
                                            if sp_tsp in tsp_idx:#[0:red_data.n_tspc]:
                                                OIC_sp_tsp[id_sp] = max(OIC_sp_tsp[id_sp],OIC_sp[idx][sp_tsp])+1e-6
                                                id_sp+=1
                                        for sp_i in range(int(sp_inter[idx])):
                                            OIC_sp_tsp[list(OIC_sp_tsp).index(max(OIC_sp_tsp))]=0
                                        best_ICsp = list(OIC_sp_tsp).index(max(OIC_sp_tsp))
                                        eps_cre[best_ICsp]=True
                                        if not sp_inter_flag[idx]: sp_inter[idx]+=1; sp_inter_flag[idx]=True

                                eps_stop[idx] = True; eps[idx]=eps_prev[idx]


                    # =============================================================
                    #                   eps evolution

                    eps_evol = ['equal']*n_tspecies
                    first_try =   [First_try_T_error]  + [First_try_ig_error]\
                                + [First_try_Sl_error] +  [First_try_K_error]\
                                +  First_try_sp_error
                    if not cross_red_error:
                        if True in first_try or not errors.under_tol: # if reduced mech fails to respect the accuracy requirements
                            if (conditions.error_param.T_check and conditions.config!='JSR')\
                            and (First_try_T_error or T_error):#errors.under_tol):
                                for sp in conditions.error_param.sp_T:
                                    idx = tspc.index(sp)
                                    eps_evol[idx] = 'decrease'
                            if (conditions.error_param.Sl_check and 'free_flame' in conditions.config)\
                            and (First_try_Sl_error or Sl_error):#errors.under_tol):
                                for sp in conditions.error_param.sp_Sl:
                                    idx = tspc.index(sp)
                                    eps_evol[idx] = 'decrease'
                            if (conditions.error_param.K_check \
                                and ('diff_flame' in conditions.config or 'pp_flame' in conditions.config))\
                            and (First_try_K_error or K_error):#errors.under_tol):
                                for sp in conditions.error_param.sp_K:
                                    idx = tspc.index(sp)
                                    eps_evol[idx] = 'decrease'
                            if (conditions.error_param.ig_check and 'reactor' in conditions.config)\
                            and (First_try_ig_error or ig_error):
                                for sp in conditions.error_param.sp_ig:
                                    idx = tspc.index(sp)
                                    eps_evol[idx] = 'decrease'
                            for sp in range(red_data.n_tspc):
                                if First_try_sp_error[sp] or sp_error[sp]:#errors.under_tol:
                                    eps_evol[sp] = 'decrease'
                        else: # if reduced mech respects the accuracy requirements
                            one_simul_success = True
                            if conditions.error_param.T_check and conditions.config!='JSR':
                                for sp in conditions.error_param.sp_T:
                                    idx = tspc.index(sp)
                                    if not eps_stop[idx] :
                                        if T_error:
                                            eps_evol[idx] = 'decrease'
                                            eps_stop[idx] = True;eps[idx]=eps_prev[idx]
                                        else:
                                            eps_evol[idx] = 'increase'
                            if conditions.error_param.Sl_check and 'free_flame' in conditions.config:
                                for sp in conditions.error_param.sp_Sl:
                                    idx = tspc.index(sp)
                                    if not eps_stop[idx] and eps_evol[idx]!='decrease':
                                        if Sl_error:
                                            eps_evol[idx] = 'decrease'
                                            eps_stop[idx] = True;eps[idx]=eps_prev[idx]
                                        else:
                                            eps_evol[idx] = 'increase'
                            if conditions.error_param.K_check \
                            and ('diff_flame' in conditions.config or 'pp_flame' in conditions.config):
                                for sp in conditions.error_param.sp_K:
                                    idx = tspc.index(sp)
                                    if not eps_stop[idx] and eps_evol[idx]!='decrease':
                                        if K_error:
                                            eps_evol[idx] = 'decrease'
                                            eps_stop[idx] = True;eps[idx]=eps_prev[idx]
                                        else:
                                            eps_evol[idx] = 'increase'
                            if conditions.error_param.ig_check and 'reactor' in conditions.config:
                                for sp in conditions.error_param.sp_ig:
                                    idx = tspc.index(sp)
                                    if not eps_stop[idx] and eps_evol[idx]!='decrease':
                                        if ig_error:
                                            eps_evol[idx] = 'decrease'
                                            eps_stop[idx] = True;eps[idx]=eps_prev[idx]
                                        else:
                                            eps_evol[idx] = 'increase'
                            for idx in range(red_data.n_tspc):
                                if not eps_stop[idx] and eps_evol[idx]!='decrease':
                                    if sp_error[idx] or eps_cre[idx]:
                                        eps_evol[idx] = 'decrease'
                                        eps_stop[idx] = True; eps[idx]=eps_prev[idx]
                                    else:
                                        eps_evol[idx] = 'increase'
                    elif red_data.red_op.inter_sp_inter:
                        eps[best_ICsp]      = eps_prev[best_ICsp]
                        if not eps_stop[best_ICsp]:
                            eps_evol[best_ICsp] = 'decrease'
                            if eps[best_ICsp]-delta_eps[best_ICsp]>1e-2:
                           	    eps[best_ICsp] -= delta_eps[best_ICsp]
                            else:
                                eps[best_ICsp] /= 2
                            eps_stop[best_ICsp] = True


                    eps_prev = list(eps)

                    # =============================================================
                    #                   New eps calculation

                    for idx in range(len(eps)):
                        if not eps_stop[idx]:
                            if eps_evol[idx] == 'increase':
                                if (eps[idx]+delta_eps[idx])<max_eps_config[idx]:
                                    eps[idx] += delta_eps[idx]
                                else:
                                    eps[idx]=max_eps_config[idx]
                                    eps_stop[idx] = True
                            elif eps_evol[idx] == 'decrease':
                                if eps[idx]-delta_eps[idx]>1e-2:
                                    eps[idx] -= delta_eps[idx]
                                else:
                                    eps[idx] /= 2

                        #Acceleration of the iterative process
                        if sp_try[idx] > try_acc[idx] and not eps_stop[idx]:
                            delta_eps[idx] = delta_eps[idx]*2
                            try_acc[idx]+=5

                    # =============================================================
                    #                   Display informations

                    # conditions / errors / eps_stop / eps_evol / tspc / First_try_(...)_error
                    if verbose>=3:
                        if conditions.error_param.T_check and conditions.config!='JSR':
                            txt_T='  Temperature error:   '+'%0.1f' %(errors.qoi_T*100)
                            if errors.under_tol_T: txt_T+='% < '
                            else:                  txt_T+='% > '
                            txt_T+=str(red_data.red_op.max_error_T)+'%    '
                            if First_try_T_error: txt_T+='(!) 1st try '
                            txt_T += '-> '
                            for sp in conditions.error_param.sp_T:
                                idx = tspc.index(sp)
                                txt_T += 'eps('+tspc[idx]+')'
                                if eps_stop[idx]:                 txt_T += '(#):'
                                elif eps_evol[idx] == 'decrease': txt_T += '(-):'
                                elif eps_evol[idx] == 'equal':    txt_T += '(=):'
                                elif eps_evol[idx] == 'increase': txt_T += '(+):'
                                txt_T += '%2.3f' %eps[idx]+'  '
                            print_(txt_T,mp)
                        if conditions.error_param.Sl_check and 'free_flame' in conditions.config:
                            txt_Sl='  Flame speed error:   '+'%0.1f' %(errors.qoi_Sl*100)
                            if errors.under_tol_Sl: txt_Sl+='% < '
                            else:                  txt_Sl+='% > '
                            txt_Sl+=str(red_data.red_op.max_error_Sl)+'%    '
                            if First_try_Sl_error: txt_Sl+='(!) 1st try '
                            txt_Sl += '-> '
                            for sp in conditions.error_param.sp_Sl:
                                idx = tspc.index(sp)
                                txt_Sl += 'eps('+tspc[idx]+')'
                                if eps_stop[idx]:                 txt_Sl += '(#):'
                                elif eps_evol[idx] == 'decrease': txt_Sl += '(-):'
                                elif eps_evol[idx] == 'equal':    txt_Sl += '(=):'
                                elif eps_evol[idx] == 'increase': txt_Sl += '(+):'
                                txt_Sl += '%2.3f' %eps[idx]+'  '
                            print_(txt_Sl,mp)
                        if conditions.error_param.K_check \
                        and ('diff_flame' in conditions.config or 'diff_flame' in conditions.config):
                            txt_K='  Extinction rate error:   '+'%0.1f' %(errors.qoi_K*100)
                            if errors.under_tol_K: txt_K+='% < '
                            else:                  txt_K+='% > '
                            txt_K+=str(red_data.red_op.max_error_K)+'%    '
                            if First_try_K_error: txt_K+='(!) 1st try '
                            txt_K += '-> '
                            for sp in conditions.error_param.sp_K:
                                idx = tspc.index(sp)
                                txt_K += 'eps('+tspc[idx]+')'
                                if eps_stop[idx]:                 txt_K += '(#):'
                                elif eps_evol[idx] == 'decrease': txt_K += '(-):'
                                elif eps_evol[idx] == 'equal':    txt_K += '(=):'
                                elif eps_evol[idx] == 'increase': txt_K += '(+):'
                                txt_K += '%2.3f' %eps[idx]+'  '
                            print_(txt_K,mp)

                        if conditions.error_param.ig_check and 'reactor' in conditions.config:
                            txt_ig='  Ignition delay error:   '+'%0.1f' %(errors.qoi_ig*100)
                            if errors.under_tol_ig: txt_ig+='% < '
                            else:                  txt_ig+='% > '
                            txt_ig+=str(red_data.red_op.max_error_ig)+'%    '
                            if First_try_ig_error: txt_ig+='(!) 1st try '
                            txt_ig += '-> '
                            for sp in conditions.error_param.sp_ig:
                                idx = tspc.index(sp)
                                txt_ig += 'eps('+tspc[idx]+')'
                                if eps_stop[idx]:                 txt_ig += '(#):'
                                elif eps_evol[idx] == 'decrease': txt_ig += '(-):'
                                elif eps_evol[idx] == 'equal':    txt_ig += '(=):'
                                elif eps_evol[idx] == 'increase': txt_ig += '(+):'
                                txt_ig += '%2.3f' %eps[idx]+'  '
                            print_(txt_ig,mp)
                        for idx in range(red_data.n_tspc):
                            sp = tspc[idx]
                            txt_sp='  ' + sp +' '*(5-len(sp))+'error:   '\
                            +'%0.1f' %(errors.qoi_s[idx]*100)
                            if errors.under_tol_s[idx]: txt_sp+='% < '
                            else:                       txt_sp+='% > '
                            txt_sp+=str(red_data.red_op.max_error_sp[idx])+'%    '
                            if First_try_sp_error[idx]: txt_sp+='(!) 1st try '
                            txt_sp += '-> '
                            if eps_stop[idx]:                 txt_sp += 'eps(#):'
                            elif eps_evol[idx] == 'decrease': txt_sp += 'eps(-):'
                            elif eps_evol[idx] == 'equal':    txt_sp += 'eps(=):'
                            elif eps_evol[idx] == 'increase': txt_sp += 'eps(+):'
                            txt_sp += '%2.3f' %eps[idx]+'  '
                            print_(txt_sp,mp)


                    # Check if the reduced mechanism can be recorded
                    if errors.under_tol:
                        act_sp_prev = active_sp_pm
                        act_r_prev  = active_r_pm
                        sp_prev     = active_sp_pm
                        r_prev      = active_r_pm
                    else :
                        sp_prev     = active_sp_pm
                        r_prev      = active_r_pm

                    n_it_max = 100
                    if max_eps_config == eps \
                    or False not in eps_stop\
                    or errors.above_tol\
                    or max(eps)<=min(eps_init)/n_it_max\
                    or max(sp_try)>n_it_max or T_try>n_it_max or ig_try>n_it_max or Sl_try>n_it_max or K_try>n_it_max:
                        stop_reduction=True

                    if cross_red_error:
                        print_('\n    ERROR FROM AN OTHER SPECIES REDUCTION',mp)
                        if red_data.red_op.inter_sp_inter and not eps_stop[idx]:
                            print_('\n    Main coupled species:'+tspc[best_ICsp]+\
                                  ' -> eps(-): '+'%0.3f' %(eps[best_ICsp]),mp)
                            cross_red_error = False
                        elif red_data.red_op.inter_sp_inter and eps_stop[idx]:
                            print_('\n    Main coupled species:'+tspc[best_ICsp]+\
                                  ' -> eps(#): '+'%0.3f' %(eps[best_ICsp]),mp)
                            cross_red_error = False



                # =========================
                # Saving results in objects
                if verbose>8: print_("\n\n# Saving results",mp)
                red_data_list[op][i]    = red_data
                if True in first_try or not one_simul_success:
                    print_('\n\nWarning: no reduction passed the accuracy threshold',mp)
                    print_('Previous reduced mech conserved',mp)
                    mech_data.spec.activ_p  = list(mech_data.spec.activ_m)
                    mech_data.react.activ_p = list(mech_data.react.activ_m)
                else:
                    mech_data.spec.activ_p = list(act_sp_prev)
                    mech_data.react.activ_p = list(act_r_prev)

                # Display informations
                if verbose>3:
                    print_("\n=>  "+str(mech_data.spec.activ_p.count(True))+\
                          " species, "+ str(mech_data.react.activ_p.count(True))+\
                          " reactions remaining",mp)
                if verbose>8:
                    species_txt='Remaining species: '
                    for sp in range(len(mech_data.spec.name)):
                        if mech_data.spec.activ_p[sp]:
                            species_txt+=mech_data.spec.name[sp]+" "
                    print_(species_txt,mp)
                    reaction_txt='Remaining reaction: '
                    for r in range(len(mech_data.react.number)):
                        if mech_data.react.activ_p[r]:
                            reaction_txt+=str(mech_data.react.number[r])+" "
                    print_(reaction_txt,mp)

                simulation += 1

            # clock
            clock.stop()

            # Save method's reduced mech
            mech_data.spec.activ_m = list(mech_data.spec.activ_p)
            mech_data.react.activ_m = list(mech_data.react.activ_p)

            # Save the non optimised mech
            filename = conditions_list[0].mech_ext
            new_filename = str(op+1) + '_' + filename
            os.chdir(conditions_list[0].main_path)
            if op==0:  os.mkdir('Red_mech')
            os.chdir('Red_mech')
            mech_data.write_new_mech(new_filename)
            if conditions_list[0].simul_param.write_ck:
                mech_data.write_chemkin_mech(new_filename)
            os.chdir(conditions_list[0].main_path)


            # Calculation of the all data with last reduced mech
            for i in range(len(conditions_list)):
                # --------------------------------------------------------------------------------
                # interpretation of the new mech

                # supress console output during the simulation
                old_stdout = sys.stdout ; old_stderr = sys.stderr
                with open(os.devnull, "w") as devnull: sys.stdout = devnull ; sys.stderr = devnull

                red_results_list[i].gas = ct.Solution('Red_mech/'+new_filename)

                # restore console output
                sys.stdout = old_stdout ; sys.stderr = old_stderr
                # --------------------------------------------------------------------------------

                red_results_list[i] = comp.red_computation(conditions_list[i],\
                                    red_results_list[i].gas,\
                                    mech_data.spec.activ_p,mech_data.react.activ_p)
                errors = cdef.Errors(conditions_list[i],ref_results_list[i],\
                                  red_results_list[i],red_data_list[op][i],red_data_list[op][i].red_op)
                red_errors_list.append(errors)



        # flush memory
        gc.collect()

        # =====================================================================
        # Optimization
        if red_data_list[op][0].optim:
            if red_method == 'NULL':
                if 'DRG' in red_data_list[op][0].optim_param.optim_on_meth:
                    for i in range(len(red_data_list[op])):
                        red_data_list[op][i]=drg.dic(red_data_list[op][i],mech_data,mech_results_list[i])
                        red_data_list[op][i]=drg.ric(red_data_list[op][i],mech_data,mech_results_list[i])
                if 'SA'  in red_data_list[op][0].optim_param.optim_on_meth:
                    for i in range(len(red_data_list[op])):
                        red_data_list[op][i]=sa.sensitivities_computation_SA(red_data_list[op][i],\
                                     mech_data,mech_results_list[i])


            clock_opt = cdef.Clock(red_data_list[op][0].optim) ; clock_opt.start()
            if red_data_list[op][0].optim == 'GA':
                mech_data.react.kin,opt_results_list,opt_errors_list,fitness\
                = ga.geneticAlgorithm(conditions_list,mech_data,\
                                  ref_results_list,red_data_list[op])
            elif red_data_list[op][0].optim == 'PSO':
                mech_data.react.kin,opt_results_list,opt_errors_list,fitness\
                = pso.PSO(conditions_list,mech_data,\
                                  ref_results_list,red_data_list[op])

            clock_opt.stop()

        # Directory
        os.chdir(conditions_list[0].main_path)

        # =====================================================================
        # Saving results in csv file
        if 'GA' in red_data_list[op][0].optim or 'PSO' in red_data_list[op][0].optim:
            opt_results_list[-1].write_method(clock,red_method+'_opt',op,\
                            fitness,clock_opt,red_data_list[op][0].optim_param)
        else:
            red_results_list[-1].write_method(clock,red_method,op)
        for l in range(len(red_results_list)):
            if conditions_list[0].exp_data:
                ref_results_list_ext_data[l].gas = ref_results_list[l].gas
                ref_results_list_ext_data[l].write_case_data('Reference',op,l+1)
            else: ref_results_list[l].write_case_data('Reference',op,l+1)
            if red_method != 'NULL':
                red_results_list[l].write_case_data('Reduction',op,False,red_errors_list[l])
            if 'GA' in red_data_list[op][0].optim or 'PSO' in red_data_list[op][0].optim:
                opt_results_list[l].write_case_data('Optimization',op,\
                                False,opt_errors_list[l])


    return  mech_data.spec.activ_p.count(True), \
            mech_data.react.activ_p.count(True), \
            red_results_list,red_data_list



def input_results_treatment(conditions_list,ref_results_list,mech_data):

    ''' return:
        - condition_list            which contain the vecto pts_scatter for the
                                    scattering definition of future simulation
        - mech_results_list         which contains the kinetic informations of
                                    the evaluated cases for DRG/SA computation
        - ref_results_list_smooth   wich contains the reference results
                                    interpolated for comparison and error
                                    computation

     computation of the different cases:
     -> to get the kinetic informations for DRG/SA computation
     -> to define the point scattering for future simulations (reduction/optimization)
     -> to interpolate the experimental data on this scattering for future comparison
     -> to shift the experimental flame position
    '''

    mech_results_list = []
    if conditions_list[-1].mech_prev_red:
        gas = ct.Solution(conditions_list[-1].mech_prev_red)
        for i in range(len(conditions_list)):
            conditions_list[i].composition.gas = gas
    for i in range(len(conditions_list)):
        results,conditions = computation_reference(conditions_list[i],conditions_list[i].simul_param.verbose,\
                                mech_data.spec.activ_m,mech_data.react.activ_m)
        mech_results_list.append(results)

    if conditions_list[0].exp_data :
    # shift free flame_data:
        for i in range(len(conditions_list)):
            conditions   = conditions_list[i]
            if conditions.config == 'free_flame' and len(ref_results_list[i].X)>4:
                # find z(1/2CH4)
                fuel_idx = conditions_list[i].composition.gas_ref.species_index(\
                        ref_results_list[i].conditions.composition.fuel.split('(')[0].split('/')[0])
                c_i = ref_results_list[i].X[0][fuel_idx]
                c_f = ref_results_list[i].X[-1][fuel_idx]
                c_12 = (c_i-c_f)/2
                for pt in range(len(ref_results_list[i].pts_scatter)):
                    if ref_results_list[i].X[pt][fuel_idx]<=c_12:
                        pt_12_exp=ref_results_list[i].pts_scatter[pt]
                        c_12 = ref_results_list[i].X[pt][fuel_idx]
                        break
                for pt in range(len(mech_results_list[i].pts_scatter)):
                    if mech_results_list[i].X[pt][fuel_idx]<=c_12:
                        pt_12_sim=mech_results_list[i].pts_scatter[pt]
                        break
                # shift the experimental flame position:
                shift = pt_12_sim-pt_12_exp
                ref_results_list[i].pts_scatter = ref_results_list[i].pts_scatter+shift
                # add points at the origin:
                ref_results_list[i].pts_scatter = np.insert\
                                                  (ref_results_list[i].pts_scatter,0,0)
                ref_results_list[i].X=[ref_results_list[i].X[0]]+ref_results_list[i].X
                ref_results_list[i].T=[ref_results_list[i].T[0]]+ref_results_list[i].T

                # remove potential extra experimental points:
                for pt in range(len(ref_results_list[i].pts_scatter)):
                    pts_scatter_end = mech_results_list[i].pts_scatter[-1]
                    conc_end        = ref_results_list[i].X[-(pt+1)]
                    T_end           = ref_results_list[i].T[-(pt+1)]
                    if ref_results_list[i].pts_scatter[-(pt+1)]>=mech_results_list[i].pts_scatter[-1]:
                        ref_results_list[i].pts_scatter = np.delete(ref_results_list[i].pts_scatter,-1)
                        del ref_results_list[i].X[-1]
                        del ref_results_list[i].T[-1]
                    else:
                        ref_results_list[i].pts_scatter=\
                            np.concatenate((ref_results_list[i].pts_scatter,[pts_scatter_end]))
                        ref_results_list[i].X+=list([conc_end])
                        ref_results_list[i].T.append(T_end)
                        ref_results_list[i].X2conc()
                        break

        # Smoothed interpolation of experimental results
        ref_results_list_smooth=[]
        for i in range(len(conditions_list)):
            ref_results_list_smooth.append(copy.copy(ref_results_list[i]))
            pts_scatter_mech = list(mech_results_list[i].pts_scatter)
            ns = conditions_list[i].composition.gas_ref.n_species
            conc_smooth = np.zeros((len(pts_scatter_mech),ns))
            for sp in range(len(mech_results_list[i].X[0])):
                y=[]
                for pt in range(len(ref_results_list[i].X)):
                    y.append(ref_results_list[i].X[pt][sp])
                x=np.array(ref_results_list[i].pts_scatter)
                y=np.array(y)
                if len(x)==1: x=np.append(x,0.2) ; y=np.append(y,y)
                conc_smooth_sp = interp1d(x, y, kind='linear')
                for pt in range(len(pts_scatter_mech)):
                    if pts_scatter_mech[pt]<=ref_results_list[i].pts_scatter[-1]:
                        conc_smooth[pt][sp]=conc_smooth_sp(pts_scatter_mech[pt])
                conc_smooth[-1] = conc_smooth[-2]

            if 'False' not in ref_results_list[i].T:
                T_v = np.array(ref_results_list[i].T)
                if len(T_v)==1: T_v=np.append(T_v,T_v)
                T = interp1d(x, T_v, kind='linear')
                T_smooth = np.zeros((len(pts_scatter_mech)))
                for pt in range(len(pts_scatter_mech)):
                    if pts_scatter_mech[pt]<=ref_results_list[i].pts_scatter[-1]:
                        T_smooth[pt]=T(pts_scatter_mech[pt])
                T_smooth[-1]=T_smooth[-2]
                ref_results_list_smooth[i].T = list(T_smooth)
            else: ref_results_list_smooth[i].T = ref_results_list[i].T

            ref_results_list_smooth[i].X = list(conc_smooth)
            ref_results_list_smooth[i].X2conc()

            ref_results_list_smooth[i].pts_scatter=list(pts_scatter_mech)
            ref_results_list_smooth[i].kf=list(mech_results_list[i].kf)
            ref_results_list_smooth[i].kr=list(mech_results_list[i].kr)
            conditions_list[i].simul_param.pts_scatter=list(pts_scatter_mech)
#        else:
#            ref_results_list_smooth = copy.copy(ref_results_list)
#            for i in range(len(conditions_list)):
#                ref_results_list_smooth[i].kf=list(mech_results_list[i].kf)
#                ref_results_list_smooth[i].kr=list(mech_results_list[i].kr)
#                ref_results_list_smooth[i].X2conc()
#                mech_results_list[i].X2conc()
    else:
        # saving and suppression of unpickable variables
        gas = ref_results_list[0].gas
        f=[]
        for res in range(len(ref_results_list)):
            f.append(ref_results_list[res].f)
            del ref_results_list[res].gas
            del ref_results_list[res].f
            del ref_results_list[res].conditions

        ref_results_list_smooth = copy.deepcopy(ref_results_list)

        for res in range(len(ref_results_list)):
            ref_results_list[res].conditions        = conditions_list[res]
            ref_results_list_smooth[res].conditions = conditions_list[res]
            ref_results_list[res].gas               = gas
            ref_results_list_smooth[res].gas        = gas
            ref_results_list[res].f                 = f[res]
            ref_results_list_smooth[res].f          = f[res]

    return conditions_list, mech_results_list, ref_results_list_smooth


def folderCreation(mech, T, P, phi, name='', verbose=3):
    """ Creation of the folder to store the reduced mechanisms"""

    if name == '':
        folder = "T"+str(T)+"p"+str(P)+"phi"+str(phi)
    else:
        folder = name
#        fs = open(mech, 'a')
#        fs.write('\n#End')
#        fs.close()
    if verbose >=8:
        print("creation of folder: ", folder)
    if not os.path.exists(folder):
        os.mkdir(folder)
    os.chdir(folder)
    os.system('cp ../'+ mech + " .")

def spec2check(conditions, tspc):

    spc2rem=0
    tspc
    if 'reactor' in conditions.config and conditions.error_param.Sl_check:
        for sp in conditions.error_param.sp_Sl:
            if sp not in tspc[0:conditions.error_param.n_tspc]:
                spc2rem += 1
    if 'reactor' in conditions.config and conditions.error_param.K_check:
        for sp in conditions.error_param.sp_K:
            if sp not in tspc[0:conditions.error_param.n_tspc]:
                spc2rem += 1
    if 'flame' in conditions.config and conditions.error_param.ig_check:
        for sp in conditions.error_param.sp_ig:
            if sp not in tspc[0:conditions.error_param.n_tspc]:
                spc2rem += 1
    if 'free_flame' in conditions.config and conditions.error_param.K_check:
        for sp in conditions.error_param.sp_K:
            if sp not in tspc[0:conditions.error_param.n_tspc]:
                spc2rem += 1
    if ('diff_flame' in conditions.config or 'pp_flame' in conditions.config) \
    and conditions.error_param.Sl_check:
        for sp in conditions.error_param.sp_Sl:
            if sp not in tspc[0:conditions.error_param.n_tspc]:
                spc2rem += 1
    if 'JSR' in conditions.config and conditions.error_param.ig_check:
        for sp in conditions.error_param.sp_ig:
            if sp not in tspc[0:conditions.error_param.n_tspc]:
                spc2rem += 1
    if 'JSR' in conditions.config and conditions.error_param.Sl_check:
        for sp in conditions.error_param.sp_Sl:
            if sp not in tspc[0:conditions.error_param.n_tspc]:
                spc2rem += 1
    if 'JSR' in conditions.config and conditions.error_param.T_check:
        for sp in conditions.error_param.sp_T:
            if sp not in tspc[0:conditions.error_param.n_tspc]:
                spc2rem += 1
    if 'JSR' in conditions.config and conditions.error_param.K_check:
        for sp in conditions.error_param.sp_K:
            if sp not in tspc[0:conditions.error_param.n_tspc]:
                spc2rem += 1
    if 'PFR' in conditions.config and conditions.error_param.ig_check:
        for sp in conditions.error_param.sp_ig:
            if sp not in tspc[0:conditions.error_param.n_tspc]:
                spc2rem += 1
    if 'PFR' in conditions.config and conditions.error_param.Sl_check:
        for sp in conditions.error_param.sp_Sl:
            if sp not in tspc[0:conditions.error_param.n_tspc]:
                spc2rem += 1
    if 'PFR' in conditions.config and conditions.error_param.K_check:
        for sp in conditions.error_param.sp_K:
            if sp not in tspc[0:conditions.error_param.n_tspc]:
                spc2rem += 1
    tspc_2_check = len(tspc)-spc2rem

    return tspc_2_check




def read_ref_data(data_file_name,gas,conc_unit,ext_data_type,tspc,mech,verbose=4):

    import csv


    #==========================================================================
    ###    Get data
    #==========================================================================

    # flags
    case_nb = -1
    step_nb = -1

    # initiating lists
    data = []
    data_save = []
    steps_save = []
    headers_steps=[]
    headers_steps_save = []
    case_titles = []
    case_data = []
    case_steps = []
    K_ext = []

    try:    open(data_file_name)
    except:
        data_file_name = '_results_input/'+data_file_name

    with open(data_file_name) as csvfile:
        readcsv = csv.reader(csvfile, delimiter=';')

        read_data=0
        for line in readcsv:
            if line!=[]:
                line_com = line
                while ' ' in line_com: line_com=line_com.replace(' ','')
                if len(line_com[0])>0:
                    if '#' in line_com[0][0]:
                        a=2
                    else:
                        if "Case" in line[0]:
                            # if applicable, save previous config data
                            if case_nb>-1:
                                data_save.append(data) ; data = [] ;
                                steps_save.append(step_title) ;
                                case_data.append(data_save); data_save = []
                                case_steps.append(steps_save); steps_save = []
                                headers_steps_save.append(headers_steps); headers_steps=[]

                            if len(K_ext)<len(case_data):
                                K_ext.append(False)
                            case_nb += 1
                            step_nb = -1
                            read_data=0   # if applicable, stop recording data in data list

                        if "reactor" in line[0] or "flame" in line[0] or "JSR" in line[0] or "PFR" in line[0]:
                            case_titles.append(line)


                        if "Step" in line[0]:
                            # if applicable, save previous step data
                            if step_nb>-1:
                                data_save.append(data) ; data = []
                                steps_save.append(step_title)

                            step_title = line[0].split(": ")[1]
                            step_nb +=1
                            read_data=0      # if applicable, stop recording data in data list

                        if "K_ext(1/s):" in line[0] and step_nb==0: #get K_ext ref
                            K_ext.append(float(line[1]))

                        if line[0] == "":
                            read_data=0      # if applicable, stop recording data in data list

                        if "T(K)" in line or "Tf(K)" in line:   # means the headers line is reached
                            headers_steps.append(line)
                            read_data=1

                        elif read_data == 1:
                            data.append(line) # Construct the data list


    # Save previous config data
    data_save.append(data)
    case_data.append(data_save)
    steps_save.append(step_title)
    case_steps.append(steps_save)
    headers_steps_save.append(headers_steps)
    if len(K_ext)<len(case_data): K_ext.append(False)

    #==========================================================================
    ###    Store data
    #==========================================================================

    conditions_list=[] ; ref_results_list=[]

    for c in range(len(case_titles)):
        config        = case_titles[c][0]
        fuel          = case_titles[c][1]
        oxidant       = case_titles[c][2]
        diluent       = case_titles[c][3]
        if "diff_" not in config:
            phi           = float(case_titles[c][4])
        else: phi = False
        if "/" in case_titles[c][5]:
            diluent_ratio = [case_titles[c][5].split('[')[1].split(',')[0],\
            float(case_titles[c][5].split('[')[1].split(',')[1].split(']')[0])]
        else:
            try:
                diluent_ratio = float(case_titles[c][5])
            except:
                diluent_ratio = ''
        if case_titles[c][6]=='':
            mixt = False
        else:
            mixt   = True
            mixt_X = case_titles[c][6]
        P             = float(case_titles[c][7])
        if "diff_" not in config and "pp_" not in config:
            Ti            = float(case_titles[c][8])
            rtol_ts       = float(case_titles[c][9])
            atol_ts       = float(case_titles[c][10])
        else:
            Ti            = float(case_titles[c][16])
            rtol_ts       = float(case_titles[c][17])
            atol_ts       = float(case_titles[c][18])
        if "flame" in config and "diff_" not in config and "pp_" not in config:
            rtol_ss   = float(case_titles[c][11])
            atol_ss   = float(case_titles[c][12])
            transport_model = case_titles[c][13]
        elif "JSR" in config:
            end_sim   = float(case_titles[c][8])
        elif "PFR" in config:
            end_sim   = float(case_titles[c][11])
            u_0       = float(case_titles[c][12])
            area      = float(case_titles[c][13])
        elif "diff_" in config or "pp_" in config:
            mdot      = float(case_titles[c][8])
            fuel_2    = case_titles[c][9]
            oxidant_2 = case_titles[c][10]
            diluent_2 = case_titles[c][11]
            if "diff_" not in config:
                phi_2 = float(case_titles[c][12])
            else: phi_2 = False
            diluent_r_2 = case_titles[c][13]
            mixt2     = case_titles[c][14]
            mdot2     = float(case_titles[c][15])
            Ti        = float(case_titles[c][16])
            rtol_ts   = float(case_titles[c][17])
            atol_ts   = float(case_titles[c][18])
            rtol_ss   = float(case_titles[c][19])
            atol_ss   = float(case_titles[c][20])
            transport_model = case_titles[c][21]
        pts_scatter=[] ; T=[] ; conc=[] ; z1=[]
        for pt in range(len(case_data[c][0])):
            pts_scatter.append(float(case_data[c][0][pt][0]))
            if "PFR" in config:
                try:
                    T.append(float(case_data[c][0][pt][2]))
                    z1.append(float(case_data[c][0][pt][1]))
                except:
                    T.append(270.5)
                    z1.append(0)
            else:
                if case_data[c][0][pt][1] != 'False':
                    try:     T.append(float(case_data[c][0][pt][1]))
                    except:
                        T.append(270.5)
                else: T=['False']

        # conc matric construction
        conc=[]
        for pt in range(len(case_data[c][0])):
            conc_l=[]
            for sp_mech_idx in range(gas.n_species):
                sp_name = gas.species_name(sp_mech_idx)
                spec_found=False
                for sp_data_idx in range(len(headers_steps_save[c][0])):
                    if sp_name==headers_steps_save[c][0][sp_data_idx]:
                        try:
                            conc_l.append(float(case_data[c][0][pt][sp_data_idx]))
                            spec_found=True ; break
                        except: a=True
                if not spec_found: conc_l.append(0)
            conc.append(list(conc_l))

        # -----------------------------
        #   conditions
        conditions_list.append(cdef.conditions(mech,config,fuel,oxidant,diluent,\
                                    phi,diluent_ratio,Ti,P,tspc))
        conditions_list[-1].import_data = True
        conditions_list[-1].mixt = mixt
        if mixt: conditions_list[-1].composition.X = mixt_X
        else: conditions_list[-1].composition.X = conditions_list[-1].composition.molarFraction(False)

        conditions_list[-1].simul_param.verbose = verbose
        if "JSR" in config:
            conditions_list[-1].simul_param.end_sim = end_sim
            conditions_list[-1].simul_param.pts_scatter = np.array(pts_scatter)
        elif "PFR" in config:
            conditions_list[-1].simul_param.end_sim = end_sim
            conditions_list[-1].simul_param.u_0     = u_0
            conditions_list[-1].simul_param.area    = area
        elif pts_scatter[-1]!=0:  conditions_list[-1].simul_param.end_sim = pts_scatter[-1]
        if "burner" in config:
            conditions_list[-1].simul_param.mdot = float(case_titles[c][14])
            conditions_list[-1].simul_param.T_profile = list(T)
            conditions_list[-1].simul_param.pts_scatter = np.array(pts_scatter)
            conditions_list[-1].simul_param.pts_scatter_i = np.array(pts_scatter)
        if ext_data_type!="Experimental_data":
            conditions_list[-1].simul_param.pts_scatter   = np.array(pts_scatter)
            conditions_list[-1].simul_param.pts_scatter_i = np.array(pts_scatter)
            conditions_list[-1].simul_param.rtol_ts       = rtol_ts
            conditions_list[-1].simul_param.atol_ts       = atol_ts
        else:
            conditions_list[-1].exp_data = True
        if "flame" in config:
            conditions_list[-1].simul_param.rtol_ss = rtol_ss
            conditions_list[-1].simul_param.atol_ss = atol_ss
            conditions_list[-1].simul_param.transport_model=transport_model
        if "diff_"  in config or "pp_" in config:
            conditions_list[-1].simul_param.mdot           = mdot
            conditions_list[-1].composition.fuel2          = fuel_2
            conditions_list[-1].composition.oxidant2       = oxidant_2
            conditions_list[-1].composition.diluent2       = diluent_2
            conditions_list[-1].composition.phi2           = phi_2
            conditions_list[-1].composition.diluent_ratio2 = diluent_r_2
            conditions_list[-1].composition.X2             = mixt2
            conditions_list[-1].simul_param.mdot2          = mdot2  # kg/m^2/s
            conditions_list[-1].state_var.T2               = Ti  # K

        conditions_list[-1].ext_data  = True
        conditions_list[-1].conc_unit = conc_unit

        # -----------------------------
        #    results
        ref_results_list.append(cdef.Sim_Results(conditions_list[-1]))
        ref_results_list[-1].pts_scatter = np.array(pts_scatter)
        ref_results_list[-1].T = list(T)
        ref_results_list[-1].P = P
        if conc_unit=="mol_m3":
            ref_results_list[-1].conc = list(conc) ; ref_results_list[-1].conc2X()
        if conc_unit=="Molar_fraction":
            ref_results_list[-1].X = list(conc) ;
            if 'False' not in ref_results_list[-1].T:
                ref_results_list[-1].X2conc()
        ref_results_list[-1].kf = False
        ref_results_list[-1].kf = False
        if "reactor" in config:
            try:
                ref_results_list[-1].ign_time_hr=float(case_titles[c][11])
                ref_results_list[-1].ign_time_sp=float(case_titles[c][11])
            except:
                ref_results_list[-1].ign_time_hr=False
                ref_results_list[-1].ign_time_sp=False
        elif "free_flame" in config:
            try:     ref_results_list[-1].Sl=float(case_titles[c][14])
            except:  ref_results_list[-1].Sl=0
        elif "diff_flame" in config or "pp_flame" in config :
            ref_results_list[-1].K_ext = K_ext[c]
            try:     ref_results_list[-1].K_max=float(case_titles[c][22])
            except:  ref_results_list[-1].K_max=0
        elif "PFR" in config:
            ref_results_list[-1].z1 = list(z1)


    return list(conditions_list), list(ref_results_list)



def plotData(spec2plot,ref_results,red_results=False,opt_results=False):

    import matplotlib
    matplotlib.use('Agg')

    import matplotlib.pyplot as plt
#        import matplotlib.cm as cm      #color management

    linestyles = [(0, ()),       #'solid'                                   0
         (0, (5, 1)),            #('densely dashed',      ),                1
         (0, (3, 1, 1, 1)),      # ('densely dashdotted',  ),               2
         (0, (5, 5)),            #('dashed',              ),                3
         (0, (1, 1)),            #('densely dotted',      ),                4
         (0, (3, 5, 1, 5, 1, 5)),       #('dashdotdotted',         ),       5
         (0, (3, 1, 1, 1, 1, 1)),       # ('densely dashdotdotted', )])     6
         (0, (5, 10)),                  #'loosely dashed',      ),          7
         (0, (3, 5, 1, 5)),             # 'dashdotted',          ),         8
         (0, (3, 10, 1, 10)),           # ('loosely dashdotted',  ),        9
         (0, (1, 5)),                   #'dotted'                           10
         (0, (3, 10, 1, 10, 1, 10)),    #('loosely dashdotdotted', ),       11
         (0, (1, 10))]*2                  #'loosely dotted'                   12
    scatter_styles = ["o", "^", "*",  "+", "s", "d", "v", "x", "p", "."]*3
    colors_styles = ["blue","red","green","grey","purple","black", "darkturkoise",\
                     "sandybrown","saddlebrown","magenta"]*3



    fig, ax = plt.subplots()
    if "PFR" in ref_results.conditions.config:
        ref_abs = np.array(ref_results.z1)*1000
        if red_results: red_abs = np.array(red_results.z1)*1000
        if opt_results: opt_abs = np.array(opt_results.z1)*1000
        ax.set_xlabel('Z (mm)')
    else:
        ref_abs = ref_results.pts_scatter
        if red_results: red_abs = red_results.pts_scatter
        if opt_results: opt_abs = opt_results.pts_scatter
        if "reactor" in ref_results.conditions.config:
            ax.set_xlabel('t (s)')
        elif "flame" in ref_results.conditions.config:
            ax.set_xlabel('Z (mm)')
        elif "JSR" in ref_results.conditions.config:
            ax.set_xlabel('T (K)')

    sp2plt_idx=[]
    for sp in spec2plot:
        sp2plt_idx.append(ref_results.conditions.composition.gas.species_index(sp))

    # plot results
    for i in range(len(sp2plt_idx)):
        Xi=[];idx=sp2plt_idx[i]
        for step_X in ref_results.X:
            Xi.append(step_X[idx])

        ax.plot(ref_abs[0:len(Xi)], Xi, \
                       linestyle=linestyles[i], color=colors_styles[0], \
                       linewidth=2, label=spec2plot[i])

    # plot red_results
    if red_results:
        for i in range(len(sp2plt_idx)):
            Xi=[];idx=sp2plt_idx[i]
            for step_X in red_results.X:
                Xi.append(step_X[idx])
            ax.plot(red_abs[0:len(Xi)], Xi, \
                           linestyle=linestyles[i], color=colors_styles[1],\
                           linewidth=2, label=spec2plot[i])
    # plot opt_results
    if opt_results:
        for i in range(len(sp2plt_idx)):
            Xi=[];idx=sp2plt_idx[i]
            for step_X in opt_results.X:
                Xi.append(step_X[idx])
            ax.plot(opt_abs[0:len(Xi)], Xi, \
                           linestyle=linestyles[i], color=colors_styles[0], \
                           linewidth=2, label=spec2plot[i])



    lines, labels = ax.get_legend_handles_labels()
    ax.legend(lines, labels, loc=0)

    plt.show()



def get_reduction_parameters(filename):

    try :   fs = open('_conditions_input/'+filename, 'r')
    except: fs = open(filename, 'r')

    caution_opt_jsr      = True
    caution_opt_fflame   = True
    caution_opt_cf_flame = True
    first_it             = True

    txt = fs.readline()

# =============================================================================
#     # Main parameters
# =============================================================================
    sp_T=[];sp_Sl=[];sp_K=[];sp_ig=[];
    T_check=False;Sl_check=False;K_check=False;ig_check=False;external_results=False
    concentration_units='Molar_fraction';ext_data_type='Experimental_data'
    ref_results_list=False;
    while 'Case' not in txt[-1] and 'Operators' not in txt[-1] and txt[0] != '':
        txt = fs.readline().split('=')
        txt[0]=txt[0].replace(' ','')

        if txt[0] == 'main_path':           main_path_ext       = clean_txt(txt[1])
        if txt[0] == 'mech':                mech                = clean_txt(txt[1])
        if txt[0] == 'mech_prev_red':       mech_prev_red       = clean_txt2(txt[1])
        if txt[0] == 'ext_results_file':    external_results    = clean_txt2(txt[1])
        if txt[0] == 'conc_units':          concentration_units = clean_txt2(txt[1])
        if txt[0] == 'ext_data_type':       ext_data_type       = clean_txt2(txt[1])
        if txt[0] == 'verbose':             verbose             = int(txt[1])
        if txt[0] == 'show_plots':
            show_plots  = str2bool(txt[1])
        if txt[0] == 'write_ck':
            write_ck    = str2bool(txt[1])
        if txt[0] == 'tspc':
            tspc          = txt2list_string(txt[1])
            n_tspc = len(tspc)
        if txt[0] == 'T_check':           T_check       = str2bool(txt[1])
        if txt[0] == 'sp_T':              sp_T          = txt2list_string(txt[1])
        if txt[0] == 'Sl_check':          Sl_check      = str2bool(txt[1])
        if txt[0] == 'sp_Sl':             sp_Sl         = txt2list_string(txt[1])
        if txt[0] == 'ig_check':          ig_check      = str2bool(txt[1])
        if txt[0] == 'sp_ig':             sp_ig         = txt2list_string(txt[1])
        if txt[0] == 'K_check':           K_check       = str2bool(txt[1])
        if txt[0] == 'sp_K':              sp_K          = txt2list_string(txt[1])
        if txt[0] == 'error_calculation': error_calculation = clean_txt(txt[1])
        if txt[0] == 'error_coupling':    error_coupling    = clean_txt(txt[1])
        if txt[0] == 'restore_flame_folder': restore_flame_folder = clean_txt(txt[1]);

    # gas
    try:
        gas_ref = ct.Solution('_kinetic_mech/'+mech)
    except:
        gas_ref = ct.Solution(mech)


# =============================================================================
#     # Simulation cases
# =============================================================================
    conditions_list = [] ; save_conditions=False ; case_n=0
    while 'Operators' not in txt[-1] and txt[0] != '':
        txt = fs.readline().split('=')
        # get data
        while 'Case' not in txt[-1] and 'Operators' not in txt[-1] and txt[0] != '':
            txt[0]=txt[0].replace(' ','')
            if txt[0]== 'config':             config        = clean_txt(txt[1])

            if txt[0] == 'fuel'          or txt[0] == 'fuel_1':              fuel          = clean_txt(txt[1])
            if txt[0] == 'oxidant'       or txt[0] == 'oxidant_1':           oxidant       = clean_txt(txt[1])
            if txt[0] == 'diluent'       or txt[0] == 'diluent_1':           diluent       = clean_txt(txt[1])
            if txt[0] == 'diluent_ratio' or txt[0] == 'diluent_ratio_1':     diluent_ratio = diluent_fct(txt[1])
            if txt[0] == 'Ts'            or txt[0] == 'Ts_1':                Ts            = txt2list_float(txt[1])
            if txt[0] == 'T_min'         or txt[0] == 'T_min_1':             T_min         = float(txt[1])
            if txt[0] == 'T_max'         or txt[0] == 'T_max_1':             T_max         = float(txt[1])
            if txt[0] == 'T_incr'        or txt[0] == 'T_incr_1':
                T_incr        = float(txt[1])
                Ts = np.arange(T_min, T_max+T_incr/2, T_incr)

            if txt[0] == 'phis'          or txt[0] == 'phis_1':              phis          = txt2list_float(txt[1])
            if txt[0] == 'phi_min'       or txt[0] == 'phi_min_1':           phi_min       = float(txt[1])
            if txt[0] == 'phi_max'       or txt[0] == 'phi_max_1':           phi_max       = float(txt[1])
            if txt[0] == 'phi_incr'      or txt[0] == 'phi_incr_1':
                phi_incr      = float(txt[1])
                phis = np.arange(phi_min, phi_max+phi_incr/2, phi_incr)
            if txt[0] == 'mixt'          or txt[0] == 'mixt_1':              mixt          = txt2mixt(txt[1])

            if txt[0] == 'P_min':             P_min         = float(txt[1])
            if txt[0] == 'P_max':             P_max         = float(txt[1])
            if txt[0] == 'P_incr':            P_incr        = float(txt[1])
            if txt[0] == 'P_incr':  Ps = np.arange(P_min, P_max+P_incr/2, P_incr)
            if txt[0] == 'Ps':                Ps            = txt2list_float(txt[1])
            if txt[0] == 'tol_ts':            tol_ts        = txt2list_float(txt[1])
            # options for reactor and PFR
            if txt[0] == 'n_pts':             n_pts           = float(txt[1])
            if txt[0] == 'delta_npts':        delta_npts      = float(txt[1])
            if txt[0] == 't_max_coeff':       t_max_coeff     = float(txt[1])
            if txt[0] == 'Scal_ref':          Scal_ref        = clean_txt(txt[1])
            if txt[0] == 'grad_curv_ratio':   grad_curv_ratio = float(txt[1])
            if txt[0] == 'tign_nPoints':      tign_nPoints    = float(txt[1])
            if txt[0] == 'tign_det_tmax_sec': t_max_react     = float(txt[1])
            # options for jsr
            if txt[0] == 't_max':             t_max           = float(txt[1]); caution_opt_jsr=False
            # options for flame
            if txt[0] == 'tol_ss':            tol_ss          = txt2list_float(txt[1])
            if txt[0] == 'transport_model':   transport_model = clean_txt(txt[1])
            if txt[0] == 'pts_scatter':       pts_scatter     = txt2list_float(txt[1])
            if txt[0] == 'slope':             slope_ff        = float(txt[1])
            if txt[0] == 'curve':             curve_ff        = float(txt[1])
            if txt[0] == 'ratio':             ratio_ff        = float(txt[1])
            if txt[0] == 'prune':             prune_ff        = float(txt[1])
            # option for free_flame / burner_flame
            if txt[0] == 'restore_flame_folder': restore_flame_folder = clean_txt(txt[1]);

            # option for burner_flame
            if txt[0] == 'T_profile':         T_profile       = txt2list_float(txt[1]);
            # option for cflow_flame
            if txt[0] == 'width':             width           = float(txt[1]);  caution_opt_cf_flame=False
            if txt[0] == 'mdot' or txt[0]=='mdots' or txt[0]=='mdots_1': mdot1 = txt2list_float(txt[1])
            if txt[0] == 'mdots_2':           mdot2           = txt2list_float(txt[1])
            if txt[0] == 'fuel_2':            fuel2           = clean_txt(txt[1])
            if txt[0] == 'oxidant_2':         oxidant2        = clean_txt(txt[1])
            if txt[0] == 'diluent_2':         diluent2        = clean_txt(txt[1])
            if txt[0] == 'diluent_ratio_2':   diluent_ratio2  = diluent_fct(txt[1])
            if txt[0] == 'T2_min':            T2_min          = float(txt[1])
            if txt[0] == 'T2_max':            T2_max          = float(txt[1])
            if txt[0] == 'T2_incr':           T2_incr         = float(txt[1])
            if txt[0] == 'T2_incr':  Ts2 = np.arange(T2_min, T2_max+T2_incr/2, T_incr)
            if txt[0] == 'Ts2':               Ts2             = txt2list_float(txt[1])
            if txt[0] == 'P2_min':            P2_min          = float(txt[1])
            if txt[0] == 'P2_max':            P2_max          = float(txt[1])
            if txt[0] == 'P2_incr':           P2_incr         = float(txt[1])
            if txt[0] == 'P2_incr':  Ps2 = np.arange(P2_min, P2_max+P2_incr/2, P_incr)
            if txt[0] == 'Ps2':               Ps2             = txt2list_float(txt[1])
            if txt[0] == 'phi2_min':          phi2_min        = float(txt[1])
            if txt[0] == 'phi2_max':          phi2_max        = float(txt[1])
            if txt[0] == 'phi2_incr':         phi2_incr       = float(txt[1])
            if txt[0] == 'phi2_incr':  phis2 = np.arange(phi2_min, phi2_max+phi2_incr/2, phi_incr)
            if txt[0] == 'phis_2':            phis2           = txt2list_float(txt[1])
            if txt[0] == 'mixt_2':            mixt2           = txt2mixt(txt[1])
            if txt[0] == 'T_lim':             T_lim           = float(txt[1])
            # options for PFR
            if txt[0] == 'pfr_autom':         pfr_autom       = str2bool(txt[1])
            if txt[0] == 'length':            xmax            = float(txt[1])
            if txt[0] == 'u_0':               u_0             = float(txt[1])
            if txt[0] == 'area':              area            = float(txt[1])


            txt = fs.readline().split('=')
            save_conditions = True

        if save_conditions and not external_results:
            case_n+=1

            for sp in sp_T:
                if sp not in tspc and T_check:   tspc.append(sp)
            for sp in sp_Sl:
                if sp not in tspc and Sl_check:  tspc.append(sp)
            for sp in sp_K:
                if sp not in tspc and K_check:   tspc.append(sp)
            for sp in sp_ig:
                if sp not in tspc and ig_check:  tspc.append(sp)

            if config == 'JSR':
                T_scatter = copy.deepcopy(Ts)
                Ts = [Ts[0]]

            if 'diff_' in config:
                phis=mdot1
                oxidant=False ; fuel2 = False

            # store data
            it_T=-1
            for T in Ts:
                it_T+=1; it_P=-1 ;
                for P in Ps:
                    it_P+=1; it_phi=-1
                    for phi in phis:
                        it_phi+=1
                        if 'fuel' not in locals():
                            fuel = mixt[0].split(',')[0].split(":")[0].replace(' ','')
                        if 'oxidant' not in locals():
                            try:    oxidant = mixt[0].split(',')[1].split(":")[0].replace(' ','')
                            except: oxidant = fuel
                        if 'diluent' not in locals():
                            try:
                                diluent_ratio = 1
                                diluent = mixt[0].split(',')[1].split(":")[0].replace(' ','')
                            except:
                                diluent = oxidant
                                diluent_ratio = 1


                        if 'diff_' in config:
                            phi=[False]; phis2=[False]*len(mdot1)
                        conditions_list.append(cdef.conditions(mech,config,\
                                        fuel,oxidant,diluent,\
                                        phi,diluent_ratio,\
                                        T,P,tspc))
                        if 'tol_ts' in locals():
                            conditions_list[-1].simul_param.tol_ts = tol_ts
                        if 'mech_prev_red' in locals():
                            if mech_prev_red != 'False':
                                conditions_list[-1].mech_prev_red = mech_prev_red
                            else: conditions_list[-1].mech_prev_red = False
                        # options for reactor or PFR
                        if 'n_pts' in locals():
                            conditions_list[-1].simul_param.n_pts = n_pts
                        # options for reactor
                        if 'delta_npts' in locals():
                            conditions_list[-1].simul_param.delta_npts = delta_npts
                        if 't_max_coeff' in locals():
                            conditions_list[-1].simul_param.t_max_coeff = t_max_coeff
                        if 'Scal_ref' in locals():
                            conditions_list[-1].simul_param.Scal_ref = Scal_ref
                        if 'grad_curv_ratio' in locals():
                            conditions_list[-1].simul_param.grad_curv_ratio = grad_curv_ratio
                        if 'tign_nPoints' in locals():
                            conditions_list[-1].simul_param.tign_nPoints = int(tign_nPoints)
                        if 'tign_dt' in locals():
                            conditions_list[-1].simul_param.tign_dt = tign_dt
                        # options for jsr
                        if 't_max' in locals():
                            caution_opt_jsr=False
                            conditions_list[-1].simul_param.end_sim = t_max
                        if config == 'JSR':
                            conditions_list[-1].simul_param.pts_scatter = T_scatter
                        # options for flame
                        if 'tol_ss' in locals():
                            conditions_list[-1].simul_param.tol_ss = tol_ss
                        if 'transport_model' in locals():
                            conditions_list[-1].simul_param.transport_model  = "Mix"
                        if 'pts_scatter' in locals():
                            conditions_list[-1].simul_param.pts_scatter = pts_scatter
                        if 'slope_ff' in locals():
                            conditions_list[-1].simul_param.slope_ff = slope_ff
                        if 'curve_ff' in locals():
                            conditions_list[-1].simul_param.curve_ff = curve_ff
                        if 'ratio_ff' in locals():
                            conditions_list[-1].simul_param.ratio_ff = ratio_ff
                        if 'prune_ff' in locals():
                            conditions_list[-1].simul_param.prune_ff = prune_ff
                        if 'restore_flame_folder' in locals():
                            if restore_flame_folder!= 'False':
                                conditions_list[-1].simul_param.restore_flame    = True
                                conditions_list[-1].simul_param.flame_res_folder = restore_flame_folder
#                            del restore_flame_folder
                        # options for free_flame or PFR
                        if 'xmax' in locals():
                            caution_opt_fflame=False
                            conditions_list[-1].simul_param.end_sim = xmax
                        # options for burner_flame
                        if 'T_profile' in locals():
                            conditions_list[-1].simul_param.T_profile = T_profile
                        # options for cflow_flame
                        if 'width' in locals():
                            caution_opt_cf_flame=False
                            conditions_list[-1].simul_param.end_sim = width
                        if 'mdot1' in locals():
                            try:    conditions_list[-1].simul_param.mdot = mdot1[it_phi]
                            except:
                                if type(mdot1) is list: conditions_list[-1].simul_param.mdot = mdot1[0]
                                else: conditions_list[-1].simul_param.mdot = mdot1
                        if 'mdot2' in locals():
                            try:    conditions_list[-1].simul_param.mdot2 = mdot2[it_phi]
                            except:
                                if type(mdot2) is list: conditions_list[-1].simul_param.mdot2 = mdot2[0]
                                else: conditions_list[-1].simul_param.mdot2 = mdot2
                        if 'fuel2' in locals():
                            conditions_list[-1].composition.fuel2 = fuel2
                        if 'oxidant2' in locals():
                            conditions_list[-1].composition.oxidant2 = oxidant2
                        if 'diluent2' in locals():
                            conditions_list[-1].composition.diluent2 = diluent2
                        if 'phis2' in locals():
                            conditions_list[-1].composition.phi2 = phis2[it_phi]
                        if 'mixt' in locals():
                            conditions_list[-1].composition.mixt = mixt[it_phi].replace(' ','').replace('\n','')
                            conditions_list[-1].composition.X = mixt[it_phi].replace(' ','').replace('\n','')
                            conditions_list[-1].composition.fuel = fuel.split('/')[0].split(',')[0].split(' ')[0].replace('(','').replace(')','')
                        else:
                            conditions_list[-1].composition.X = conditions_list[-1].composition.molarFraction(False)
                        if 'mixt2' in locals():
                            conditions_list[-1].composition.mixt2 = mixt2[it_phi].replace(' ','').replace('\n','')
                        if 'diluent_ratio2' in locals():
                            conditions_list[-1].composition.diluent_ratio2 = diluent_ratio2
                        if 'Ts2' in locals():
                            conditions_list[-1].state_var.T2 = Ts2[it_T]
                        if 'Ps2' in locals():
                            conditions_list[-1].simul_param.P2 = Ps2[it_T]
                        if 'phis2' in locals():
                            conditions_list[-1].simul_param.phi2 = phis2[it_phi]
                        if 'T_lim' in locals():
                            conditions_list[-1].simul_param.T_lim = T_lim
                        if 'mixt2' in locals():
                            conditions_list[-1].composition.X2 = mixt2[it_phi].replace(' ','').replace('\n','')
                        elif 'fuel2' in locals():
                            if "(" in str(conditions_list[-1].composition.fuel2)+str(conditions_list[-1].composition.oxidant2) and gas_ref:
                                X2 = conditions_list[-1].composition.molarFraction_fromvalues(gas_ref,True)
                                conditions_list[-1].composition.X2 = X2
                            elif gas_ref:
                                X2 = conditions_list[-1].composition.molarFraction(gas_ref,True)
                                conditions_list[-1].composition.X2 = X2
                        # options for PFR
                        if 'pfr_autom' in locals():
                            conditions_list[-1].simul_param.PFR_auto_time = pfr_autom
                        if 'u_0' in locals():
                            conditions_list[-1].simul_param.u_0 = u_0
                        if 'area' in locals():
                            conditions_list[-1].simul_param.area = area

                        # main parameters
                        if 'verbose' in locals():
                            conditions_list[-1].simul_param.verbose = verbose
                        if 'show_plots' in locals():
                            conditions_list[-1].simul_param.show_plots = show_plots
                        if 'write_ck' in locals():
                            conditions_list[-1].simul_param.write_ck = write_ck
                        if 'error_calculation' in locals():
                            conditions_list[-1].error_param.error_calculation = error_calculation
                        if 'error_coupling' in locals():
                            conditions_list[-1].error_param.error_coupling = error_coupling
                        if 'T_check' in locals():
                            conditions_list[-1].error_param.T_check = T_check
                        if 'ig_check' in locals():
                            conditions_list[-1].error_param.ig_check = ig_check
                        if 'Sl_check' in locals():
                            conditions_list[-1].error_param.Sl_check = Sl_check
                        if 'K_check' in locals():
                            conditions_list[-1].error_param.K_check = K_check
                        if 'sp_T' in locals():
                            conditions_list[-1].error_param.sp_T = sp_T
                        if 'sp_ig' in locals():
                            conditions_list[-1].error_param.sp_ig = sp_ig
                        if 'sp_Sl' in locals():
                            conditions_list[-1].error_param.sp_Sl = sp_Sl
                        if 'sp_K' in locals():
                            conditions_list[-1].error_param.sp_K = sp_K

                        conditions_list[-1].error_param.tspc   = tspc
                        conditions_list[-1].error_param.n_tspc = n_tspc
                        conditions_list[-1].main_path = main_path_ext

            if caution_opt_jsr and config=='JSR':
                print('WARNING: resident time for JSR not given (t_max) !')
                print('default value: t_max = 0.02 s')
                conditions_list[-1].simul_param.end_sim = 0.02
            if caution_opt_fflame and config=='free_flame':
                print('WARNING: xmax for free flame not given (xmax) !')
                print('default value: t_max = 0.02 m')
                conditions_list[-1].simul_param.end_sim = 0.02

            if 'config'           in locals(): del config
            if 'fuel'             in locals(): del fuel
            if 'oxidant'          in locals(): del oxidant
            if 'diluent'          in locals(): del diluent
            if 'diluent_ratio'    in locals(): del diluent_ratio
            if 'Ts'               in locals(): del Ts
            if 'T_min'            in locals(): del T_min
            if 'T_max'            in locals(): del T_max
            if 'T_incr'           in locals(): del T_incr
            if 'Ts'               in locals(): del Ts
            if 'phis'             in locals(): del phis
            if 'phi_min'          in locals(): del phi_min
            if 'phi_max'          in locals(): del phi_max
            if 'phi_incr'         in locals(): del phi_incr
            if 'phis2'            in locals(): del phis2
            if 'mixt'             in locals(): del mixt
            if 'P_min'            in locals(): del P_min
            if 'P_max'            in locals(): del P_max
            if 'P_incr'           in locals(): del P_incr
            if 'Ps'               in locals(): del Ps
            if 'tol_ts'           in locals(): del tol_ts
            # options for reactor
            if 'n_pts'            in locals(): del n_pts
            if 'delta_npts'       in locals(): del delta_npts
            if 't_max_coeff'      in locals(): del t_max_coeff
            if 'Scal_ref'         in locals(): del Scal_ref
            if 'grad_curv_ratio'  in locals(): del grad_curv_ratio
            if 'tign_nPoints'     in locals(): del tign_nPoints
            if 'tign_dt'          in locals(): del tign_dt
            # options for jsr
            if 't_max'            in locals(): del t_max ; caution_opt_jsr=True
            # options for flame
            if 'tol_ss'           in locals(): del tol_ss
            if 'transport_model'  in locals(): del transport_model
            if 'pts_scatter'      in locals(): del pts_scatter
            if 'slope_ff'         in locals(): del slope_ff
            if 'curve_ff'         in locals(): del curve_ff
            if 'ratio_ff'         in locals(): del ratio_ff
            if 'prune_ff'         in locals(): del prune_ff
            # option for free_flame
            if 'xmax'             in locals(): del xmax  ; caution_opt_fflame=True
            # option for burner_flame
            if 'T_profile'        in locals(): del T_profile
            # option for cflow_flame
            if 'width'            in locals(): del width ; caution_opt_cf_flame=True
            if 'mdot1'            in locals(): del mdot1
            if 'mdot2'            in locals(): del mdot2
            if 'fuel2'            in locals(): del fuel2
            if 'oxidant2'         in locals(): del oxidant2
            if 'diluent2'         in locals(): del diluent2
            if 'diluent_ratio2'   in locals(): del diluent_ratio2
            if 'T2_min'           in locals(): del T2_min
            if 'T2_max'           in locals(): del T2_max
            if 'T2_incr'          in locals(): del T2_incr
            if 'Ts2'              in locals(): del Ts2
            if 'Ts2'              in locals(): del Ts2
            if 'P2_min'           in locals(): del P2_min
            if 'P2_max'           in locals(): del P2_max
            if 'P2_incr'          in locals(): del P2_incr
            if 'Ps2'              in locals(): del Ps2
            if 'phi2_min'         in locals(): del phi2_min
            if 'phi2_max'         in locals(): del phi2_max
            if 'phi2_incr'        in locals(): del phi2_incr
            if 'phis2'            in locals(): del phis2
            if 'mixt2'            in locals(): del mixt2
            if 'T_lim'            in locals(): del T_lim
            if 'restore_flame_folder' in locals(): del restore_flame_folder
            save_conditions = False

    if external_results:
        if 'verbose' in locals():
            conditions_list, ref_results_list = read_ref_data\
            (external_results,gas_ref,concentration_units,ext_data_type,tspc,mech,verbose)
        else:
            conditions_list, ref_results_list = read_ref_data\
            (external_results,gas_ref,concentration_units,ext_data_type,tspc,mech)

        for sp in sp_T:
            if sp not in tspc and T_check:  tspc.append(sp)
        for sp in sp_Sl:
            if sp not in tspc and Sl_check: tspc.append(sp)
        for sp in sp_K:
            if sp not in tspc and K_check:  tspc.append(sp)
        for sp in sp_ig:
            if sp not in tspc and ig_check: tspc.append(sp)
        if 'mech_prev_red' in locals():
            if mech_prev_red != 'False':
                for _ci in range(len(conditions_list)):
                    conditions_list[_ci].mech_prev_red = mech_prev_red
            else: conditions_list[-1].mech_prev_red = False
        if 'restore_flame_folder' in locals():
            if restore_flame_folder!= 'False':
                for _ci in range(len(conditions_list)):
                    conditions_list[_ci].simul_param.restore_flame    = True
                    conditions_list[_ci].simul_param.flame_res_folder = restore_flame_folder
            del restore_flame_folder



    # main path
    r_path = os.getcwd()
    now    = datetime.datetime.now()
    date   = '/'+now.strftime("%Y%m%d_%H%M_")
    if 'main_path_ext' not in locals(): main_path_ext = ''
    if '/' in main_path_ext:
        main_folder = r_path + '/' + main_path_ext.split('/')[0]
        main_path = r_path + '/' + main_path_ext.split('/')[0] + date + main_path_ext.split('/')[1]
    else:
        main_path = r_path + date + main_path_ext
        main_folder = False

    for cond in conditions_list:
        # main parameters
        if 'verbose' in locals():
            cond.simul_param.verbose = verbose
        if 'show_plots' in locals():
            cond.simul_param.show_plots = show_plots
        if 'error_calculation' in locals():
            cond.error_param.error_calculation = error_calculation
        if 'error_coupling' in locals():
            cond.error_param.error_coupling = error_coupling
        if 'T_check' in locals():
            cond.error_param.T_check = T_check
        if 'ig_check' in locals():
            cond.error_param.ig_check = ig_check
        if 'Sl_check' in locals():
            cond.error_param.Sl_check = Sl_check
        if 'K_check' in locals():
            cond.error_param.K_check = K_check
        if 'sp_T' in locals():
            cond.error_param.sp_T = sp_T
        if 'sp_ig' in locals():
            cond.error_param.sp_ig = sp_ig
        if 'sp_Sl' in locals():
            cond.error_param.sp_Sl = sp_Sl
        if 'sp_K' in locals():
            cond.error_param.sp_K = sp_K
        cond.error_param.tspc   = tspc
        cond.error_param.n_tspc = n_tspc
        cond.main_path          = main_path
        cond.r_path             = r_path



    while '> Op:' not in txt[-1] and txt[0] != '':
        txt = fs.readline().split('=')

# =============================================================================
#     # Reduction / optimization operators
# =============================================================================
    red_data_list = [] ; save_op = False ; op_n=0
    if txt[0] == '':
        print('Error: operators are missing')
        print(' you must add, to your input file, the text: ')
        print('\n\n\n#=============================================\n')
        print('#           Operators\n')
        print('#=============================================\n\n')
        print('for simulation only, write: ')
        print('operator        = NULL \n')

    while txt[0] != '':
        while '> Op:' not in txt[-1] and txt[0] != '':
            txt[0]=txt[0].replace(' ','')

            if txt[0]== 'operator':           reduction_operator = clean_txt(txt[1])
            if txt[0] == 'eps':
                eps = txt2list_float(txt[1])
                while len(eps)<len(tspc): eps.append(eps[-1])
            if txt[0] == 'delta_eps':
                delta_eps = txt2list_float(txt[1])
                while len(delta_eps)<len(tspc): delta_eps.append(delta_eps[-1])
            if txt[0] == 'n_points':
                n_points           = float(txt[1])
            if txt[0] == 'max_error_sp':
                max_error_sp = txt2list_float(txt[1])
                while len(max_error_sp)<n_tspc: max_error_sp.append(max_error_sp[-1])
            if txt[0] == 'max_error_T':       max_error_T        = float(txt[1])
            if txt[0] == 'max_error_ig':      max_error_ig       = float(txt[1])
            if txt[0] == 'max_error_Sl':      max_error_Sl       = float(txt[1])
            if txt[0] == 'max_error_K':       max_error_K        = float(txt[1])
            if txt[0] == 'inter_sp_inter':    inter_sp_inter     = str2bool(txt[1])
            if txt[0] == 'optim':             optim              = clean_txt(txt[1])

            if txt[0] == 'csp_method':        csp_method         = clean_txt(txt[1])
            if txt[0] == 'tol_csp':           tol_csp            = float(txt[1])
            if txt[0] == 'time_resolution':   time_resolution    = float(txt[1])
            if txt[0] == 'epsilon_rel':
                epsilon_rel        = float(txt[1])
            if txt[0] == 'epsilon_abs':       epsilon_abs        = float(txt[1])
            if txt[0] == 'csp_refin_iter':    csp_refin_iter     = int(txt[1])

            if txt[0] == 'ttol_sensi':
                try:    ttol_sensi = txt2list_float(txt[1])
                except: ttol_sensi = txt2list_bool(txt[1])
            if txt[0] == 'ttol_sensi':
                try:    ttol_sensi = txt2list_float(txt[1])
                except: ttol_sensi = txt2list_bool(txt[1])

            if 'optim' in locals():
                if 'GA' in optim or 'PSO' in optim:
                    if txt[0] == 'n_gen':               n_gen              = int(txt[1])
                    if txt[0] == 'n_it':                n_gen              = int(txt[1])
                    if txt[0] == 'n_indiv':             n_indiv            = int(txt[1])
                    if txt[0] == 'error_fitness':       error_fitness      = clean_txt(txt[1])
                    if txt[0] == 'Arrh_max_variation':  Arrh_max_variation = txt2list_float(txt[1])
                    if txt[0] == 'optim_on_meth':       optim_on_meth      = clean_txt(txt[1])
                    if txt[0] == 'nb_r2opt':            nb_r2opt           = float(txt[1])
                    if txt[0] == 'sub_mech_sel':        sub_mech_sel       = txt2list_string(txt[1])

                    # genetic algorithm options
                    if txt[0] == 'selection_operator':  selection_operator = clean_txt(txt[1])
                    if txt[0] == 'selection_options':   selection_options  = txt2list_float(txt[1])
                    if txt[0] == 'Xover_operator':      Xover_operator     = txt2list_string(txt[1])
                    if txt[0] == 'Xover_pct':           Xover_pct          = txt2list_int(txt[1])
                    if txt[0] == 'Xover_opt':           Xover_opt          = txt2list_float(txt[1])
                    if txt[0] == 'mut_operator':        mut_operator       = txt2list_string(txt[1])
                    if txt[0] == 'mut_pct':             mut_pct            = txt2list_int(txt[1])
                    if txt[0] == 'mut_opt':             mut_opt            = txt2list_float(txt[1])
                    if txt[0] == 'mut_intensity':       mut_intensity      = float(txt[1])
                    # PSO options
                    if txt[0] == 'inertia_score':       inertia_score      = str2bool(txt[1])
                    if txt[0] == 'inertia_min':         inertia_min        = float(txt[1])
                    if txt[0] == 'inertia_max_i'\
                    or txt[0] == 'inertia_i':           inertia_i          = float(txt[1])
                    if txt[0] == 'inertia_max_endi'\
                    or txt[0] == 'inertia_end':         inertia_end        = float(txt[1])
                    if txt[0] == 'cognitive_accel_i':   cognitive_accel_i   = float(txt[1])
                    if txt[0] == 'score':               score_integ         = str2bool(txt[1])
                    if txt[0] == 'cognitive_accel_end': cognitive_accel_end = float(txt[1])
                    if txt[0] == 'social_accel_i':      social_accel_i      = float(txt[1])
                    if txt[0] == 'social_accel_end':    social_accel_end    = float(txt[1])

            else:
                optim = 'False'
            txt = fs.readline().split('=')
            save_op = True


        # store data
        if save_op:
            red_data_list.append([])
            for i in range(len(conditions_list)):
                red_data = cdef.Red_data(gas_ref,mech,tspc,n_tspc,\
                                         reduction_operator,optim,verbose)
                red_data_list[-1].append(red_data)
                if 'eps'            in locals(): red_data.red_op.eps_init       = eps
                if 'delta_eps'      in locals(): red_data.red_op.delta_eps_init = delta_eps
                if 'n_points'       in locals(): red_data.red_op.n_points       = n_points
                if 'max_error_sp'   in locals(): red_data.red_op.max_error_sp   = max_error_sp
                if 'max_error_T'    in locals(): red_data.red_op.max_error_T    = max_error_T
                if 'max_error_ig'   in locals(): red_data.red_op.max_error_ig   = max_error_ig
                if 'max_error_Sl'   in locals(): red_data.red_op.max_error_Sl   = max_error_Sl
                if 'max_error_K'    in locals(): red_data.red_op.max_error_K    = max_error_K
                if 'inter_sp_inter' in locals(): red_data.red_op.inter_sp_inter = inter_sp_inter
                if 'optim'          in locals(): red_data.red_op.optim          = optim
                if 'ttol_sensi'     in locals(): red_data.red_op.tol_ts         = ttol_sensi

                if 'csp_method'     in locals(): red_data.red_op.csp_method     = csp_method
                if 'tol_csp'        in locals(): red_data.red_op.tol_csp        = tol_csp
                if 'time_resolution'in locals(): red_data.red_op.TimeResolution = time_resolution
                if 'epsilon_rel'    in locals(): red_data.red_op.epsilon_rel    = epsilon_rel
                if 'epsilon_abs'    in locals(): red_data.red_op.epsilon_abs    = epsilon_abs
                if 'csp_refin_iter' in locals(): red_data.red_op.csp_refin_iter = csp_refin_iter


                if 'optim' in locals():
                    if 'GA' in optim or 'PSO' in optim:
                        red_data.optim_param = cdef.Optim_param(tspc, n_tspc)
                        if 'n_gen'              in locals(): red_data.optim_param.n_gen              = n_gen
                        if 'n_indiv'            in locals(): red_data.optim_param.n_ind              = n_indiv
                        if 'error_fitness'      in locals(): red_data.optim_param.error_fitness      = error_fitness
                        if 'Arrh_max_variation' in locals(): red_data.optim_param.Arrh_max_variation = Arrh_max_variation
                        if 'optim_on_meth'      in locals(): red_data.optim_param.optim_on_meth      = optim_on_meth
                        if 'nb_r2opt'           in locals(): red_data.optim_param.nb_r2opt           = nb_r2opt
                        if 'sub_mech_sel' in locals():
                            # C0 submech
                            if 'H2' not in sub_mech_sel: red_data.optim_param.opt_subm_C[0] = False
                            # C1-C30 submech
                            for x in range(29):
                                if 'C'+str(x+1) not in sub_mech_sel: red_data.optim_param.opt_subm_C[x+1] = False
                            # CO submech
                            if 'CO' not in sub_mech_sel: red_data.optim_param.opt_subm_CO = False
                            # N submech
                            for x in range(30):
                                if 'N' not in sub_mech_sel: red_data.optim_param.opt_subm_N[x] = False
                            # S / Si submech
                            if 'S' not in sub_mech_sel: red_data.optim_param.opt_subm_S = False
                            if 'Si' not in sub_mech_sel: red_data.optim_param.opt_subm_Si = False
                            del sub_mech_sel
                        if red_data.optim_param.optim_on_meth:
                            red_data.optim_param.nr = round((red_data.optim_param.nb_r2opt/100)*gas_ref.n_reactions)
                        red_data.optim_param.main_path = main_path
                        if external_results:
                            red_data.optim_param.exp_data = True

                        # genetic algorithm options
                        if 'selection_operator' in locals(): red_data.optim_param.selection_operator = selection_operator
                        if 'selection_options'  in locals(): red_data.optim_param.selection_options  = selection_options
                        if 'Xover_operator'     in locals(): red_data.optim_param.Xover_operator     = Xover_operator
                        if 'Xover_pct'          in locals(): red_data.optim_param.Xover_pct          = Xover_pct
                        if 'Xover_opt'          in locals(): red_data.optim_param.Xover_option       = Xover_opt
                        if 'mut_operator'       in locals(): red_data.optim_param.mut_operator       = mut_operator
                        if 'mut_pct'            in locals(): red_data.optim_param.mut_pct            = mut_pct
                        if 'mut_opt'            in locals(): red_data.optim_param.mut_option         = mut_opt
                        if 'mut_intensity'      in locals(): red_data.optim_param.mut_intensity      = mut_intensity
                        # PSO options
                        if 'inertia_score'       in locals(): red_data.optim_param.inertia_score      = inertia_score
                        if 'inertia_min'         in locals(): red_data.optim_param.inertia_min           = inertia_i
                        if 'inertia_i'           in locals(): red_data.optim_param.inertia_i           = inertia_i
                        if 'inertia_end'         in locals(): red_data.optim_param.inertia_end         = inertia_end
                        if 'cognitive_accel_i'   in locals(): red_data.optim_param.cognitive_accel_i   = cognitive_accel_i
                        if 'cognitive_accel_end' in locals(): red_data.optim_param.cognitive_accel_end = cognitive_accel_end
                        if 'social_accel_i'      in locals(): red_data.optim_param.social_accel_i      = social_accel_i
                        if 'social_accel_end'    in locals(): red_data.optim_param.social_accel_end    = social_accel_end


            if 'eps'                in locals(): del eps
            if 'delta_eps'          in locals(): del delta_eps
            if 'n_points'           in locals(): del n_points
            if 'max_error_sp'       in locals(): del max_error_sp
            if 'max_error_T'        in locals(): del max_error_T
            if 'max_error_ig'       in locals(): del max_error_ig
            if 'max_error_Sl'       in locals(): del max_error_Sl
            if 'max_error_K'        in locals(): del max_error_K
            if 'inter_sp_inter'     in locals(): del inter_sp_inter
            if 'optim'              in locals(): del optim
            if 'ttol_sensi'         in locals(): del ttol_sensi
            if 'n_gen'              in locals(): del n_gen
            if 'n_indiv'            in locals(): del n_indiv
            if 'error_fitness'      in locals(): del error_fitness
            if 'Arrh_max_variation' in locals(): del Arrh_max_variation
            if 'optim_on_meth'      in locals(): del optim_on_meth
            if 'nb_r2opt'           in locals(): del nb_r2opt
            if 'selection_operator' in locals(): del selection_operator
            if 'selection_options'  in locals(): del selection_options
            if 'Xover_operator'     in locals(): del Xover_operator
            if 'Xover_pct'          in locals(): del Xover_pct
            if 'Xover_opt'          in locals(): del Xover_opt
            if 'mut_operator'       in locals(): del mut_operator
            if 'mut_pct'            in locals(): del mut_pct
            if 'mut_opt'            in locals(): del mut_opt
            if 'mut_intensity'      in locals(): del mut_intensity
            if 'sub_mech_sel'       in locals(): del sub_mech_sel
            if 'inertia_score'       in locals(): del inertia_score
            if 'inertia_min'         in locals(): del inertia_min
            if 'inertia_i'           in locals(): del inertia_i
            if 'inertia_end'         in locals(): del inertia_end
            if 'cognitive_accel_i'   in locals(): del cognitive_accel_i
            if 'cognitive_accel_end' in locals(): del cognitive_accel_end
            if 'social_accel_i'      in locals(): del social_accel_i
            if 'social_accel_end'    in locals(): del social_accel_end

            save_op = False

        txt = fs.readline().split('=')

    for num in range(len(conditions_list)):
        conditions_list[num].num = str(num+1)+'_'

    return conditions_list, red_data_list, ref_results_list



def clean_txt(_txt):
    _txt  = _txt.replace(' ','')
    _txt  = _txt.replace('\n','')
    _txt  = _txt.replace('[','')
    _txt  = _txt.replace(']','')
    return _txt

def clean_txt2(_txt):
    while _txt[0] == ' ':
        _txt = _txt[1:-1]
    return _txt

def txt2list_string(_txt):
    _txt  = _txt.replace(' ','')
    _txt  = _txt.replace('[','')
    _txt  = _txt.replace(']','')
    _txt  = _txt.replace('\n','')
    _txt  = _txt.replace("'","")
    _txt  = _txt.replace('"','')
    _list = _txt.split(',')
    return _list

def txt2list_float(_txt):
    _txt  = _txt.replace(' ','')
    _txt  = _txt.replace('[','')
    _txt  = _txt.replace(']','')
    _txt  = _txt.replace('\n','')
    _list = _txt.split(',')
    for l in range(len(_list)):
        if _list[l]=='': _list[l]=0
    _list = [float(num) for num in _list]
    return _list

def txt2list_int(_txt):
    _txt  = _txt.replace(' ','')
    _txt  = _txt.replace('[','')
    _txt  = _txt.replace(']','')
    _txt  = _txt.replace('\n','')
    _list = _txt.split(',')
    for l in range(len(_list)):
        if _list[l]=='': _list[l]=0
    _list = [int(num) for num in _list]
    return _list

def txt2list_bool(_txt):
    _txt  = _txt.replace(' ','')
    _txt  = _txt.replace('[','')
    _txt  = _txt.replace('\n','')
    _txt  = _txt.replace(']','')

    _list = _txt.split(',')
    _list = [str2bool(_bool) for _bool in _list]
    return _list

def txt2mixt(_txt):
    _list = _txt.split(';')
    return _list

def diluent_fct(_txt):
    _list = _txt.split(' ')
    while '' in _list : _list.remove('')
    _list[0].replace(',','')
    try:
        diluent_ratio = float(_list[0])
    except:
        _list[1].replace(',','')
        diluent_ratio = [_list[0], float(_list[1])]
    return diluent_ratio

def str2bool(v):
  v = v.replace(' ','')
  v = v.replace('\n','')

  return v.lower() in ("yes", "true", "t", "1")

def copy_ref_results(ref_results_list):

    # 1- save unpickable variables:
    gas     = ref_results_list[0].conditions.composition.gas
    gas_ref = ref_results_list[0].conditions.composition.gas_ref
    f_save = []
    for ref_results in ref_results_list:
        f_save.append(ref_results.f)

    # 2- remove unpickable variables from ref_results_list:
    for ref_results in ref_results_list:
        del ref_results.gas
        del ref_results.conditions.composition.gas
        del ref_results.conditions.composition.gas_ref
        del ref_results.f

    # 3- create red_results_list as a copy ref_results_list:
    ref_results_copy = copy.deepcopy(ref_results_list)

    # 4- create red_results_list as a copy ref_results_list:
    for i in range(len(ref_results_copy)):
        ref_results_copy[i].gas                             = gas_ref
        ref_results_copy[i].conditions.composition.gas      = gas
        ref_results_copy[i].conditions.composition.gas_ref  = gas_ref
        ref_results_copy[i].f                               = f_save[i]
        ref_results_list[i].gas                             = gas_ref
        ref_results_list[i].conditions.composition.gas      = gas
        ref_results_list[i].conditions.composition.gas_ref  = gas_ref
        ref_results_list[i].f                               = f_save[i]


    return ref_results_copy


def copytree(src, dst, symlinks=False, ignore=None):
    import shutil
    if not os.path.exists(dst):
        os.makedirs(dst)
    for item in os.listdir(src):
        s = os.path.join(src, item)
        d = os.path.join(dst, item)
        if os.path.isdir(s):
            copytree(s, d, symlinks, ignore)
        else:
            if not os.path.exists(d) or os.stat(s).st_mtime - os.stat(d).st_mtime > 1:
                shutil.copy2(s, d)
