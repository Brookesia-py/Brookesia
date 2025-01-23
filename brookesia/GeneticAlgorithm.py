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

import brookesia.Computation as comp
import brookesia.Class_def as cdef
from  brookesia.Class_def import print_
import os
import sys
import multiprocessing
from multiprocessing import Pool
from pebble import ProcessPool
from concurrent.futures import TimeoutError


from shutil import copyfile
import time as timer
import cantera as ct

import random
import copy
#import gc
import operator
import csv

import warnings
#warnings.filterwarnings("ignore") # suppress warnings with exponential calculation
import pdb  # pdb.set_trace()

def geneticAlgorithm(conditions_list,mech_data,ref_results_list,red_data_list):


    verbose     = conditions_list[0].simul_param.verbose
    optim_param = red_data_list[0].optim_param
    mp          = optim_param.main_path

    print_ ('\n\n ----------------------------------------------------',mp)
    print_('           Genetic Algorithm  \n',mp)
    print_(str(optim_param.n_gen)+ ' generations  /  '+\
                        str(optim_param.n_ind)+' individuals\n',mp)

    os.chdir(mp)
    if not os.path.exists("GA"): os.mkdir("GA")
    os.chdir("GA")


    # Population creation, composed of:
    # 1- selected individuals
    # 2- childs from Xover
    # 3- childs from mutations
    size_ind   = optim_param.n_ind
    optim_param.count_Xover()
    size_Xover = optim_param.total_Xover
    optim_param.count_mut()
    size_mut   = optim_param.total_mut
    size_tot   = int(size_ind + size_Xover + size_mut)
    pop = Population(conditions_list,mech_data,red_data_list,ref_results_list,\
                     size_tot)

    pop2keep = pop.def_pop2keep()

    # Reference ind
    mech_data = pop.get_uncertainty(mech_data,optim_param,conditions_list)
    ref_ind = Chromosome(conditions_list,mech_data,\
                         ref_results_list,red_data_list,False)
    ref_ind.fitness = ref_ind.fitness_eval(conditions_list,optim_param,ref_results_list,0,True)
    if verbose >= 1 :
        print_("Non optimized reduced mechanism fitness: "+"%.3f"%(ref_ind.fitness),mp)
    best_ind = copy.deepcopy(ref_ind)


    # Population evaluation:
    pop.fitness_eval_newpop(optim_param,conditions_list,ref_results_list)


    # Find new best ind and save the mech,
    # if not, replace worst ind of the current pop by the previous best ind
    best_ind,new_best_ind = pop.compare_best_ind(best_ind,optim_param,verbose)


    gen=0
    # save and display convergence informations
    warnings.filterwarnings("ignore", category=ResourceWarning)

    if verbose >= 1 : print_('Initial population:',mp)
    pop.convergence_information(gen,optim_param,verbose)

    pop.selection(optim_param,verbose)

    time_1 = timer.time()

    for gen in range(1, optim_param.n_gen+1):
        print_("\n\nGeneration:" + str(gen),mp)
        pop.Xover(optim_param,conditions_list,ref_results_list,verbose)
        pop.mutation(optim_param,conditions_list,ref_results_list,gen,verbose)
        pop.fitness_eval_newchilds(optim_param,conditions_list,ref_results_list)
        pop.selection(optim_param,verbose)
        pop.check_early_conv(optim_param,gen,verbose)
        best_ind,new_best_ind = pop.compare_best_ind(best_ind,optim_param,verbose)
        if new_best_ind and optim_param.exp_data:                              # optimization of the time iteration for reactor models
            best_ind.fitness = best_ind.fitness_eval(conditions_list,optim_param,ref_results_list)
        if gen<optim_param.keep_until_gen:
            pop.insert_pop2keep(optim_param, pop2keep)
        if verbose > 5: pop.display(optim_param)
        pop.convergence_information(gen,optim_param,verbose)

    time_2 = timer.time()
    if verbose >= 5 :
        print_("\n      time for optimization: "+str(round(time_2-time_1))+'s',mp)


    new_filename = str(red_data_list[0].op+1) + 'opt_' + conditions_list[0].mech_ext

    os.chdir(conditions_list[0].main_path)
    if not os.path.exists('Red_mech'):  os.mkdir('Red_mech')
    os.chdir('Red_mech')
    if '.cti' in new_filename:
        best_ind.mech.write_new_mech(new_filename)
    else:
        best_ind.mech.write_yaml_mech(new_filename)

    if conditions_list[0].simul_param.write_ck:
        best_ind.mech.write_chemkin_mech(new_filename,conditions_list[0].version)
    os.chdir(conditions_list[0].main_path)

    Opt_results_list, errors_list, fitness = \
        best_ind.export_data(conditions_list,optim_param,ref_results_list)

    del pop ; del ref_ind

    os.chdir(conditions_list[0].main_path)

    return list(best_ind.mech.react.kin),Opt_results_list,errors_list,fitness



def plotConvergence(optim_param):
    import matplotlib.pyplot as plt

    fig = plt.figure()
    fig_title = 'genetic algorithm convergence'
    fig.text(0.5,0.95,fig_title,fontsize=22,horizontalalignment='center')
    plt.plot(optim_param.genVec, optim_param.best_fitness, label='best fitness')#, marker='o')
    plt.plot(optim_param.genVec, optim_param.mean_fitness, label='average fitness')#,marker='o')
    plt.xlabel('generation')
    plt.ylabel("fitness")
    plt.legend()
    plt.show()




class Chromosome:
    def __init__(self,conditions_list,mech_data,ref_results_list,red_data_list,\
                 rand_kin=True):
        optim_param = red_data_list[0].optim_param
        if type(rand_kin) is bool:
            self.mech    = copy.deepcopy(mech_data)
        else: # if the option to import kinetic mechanism have been selected
            self.mech    = copy.deepcopy(rand_kin)
        self.r2opt   = []
        self.find_r2opt(red_data_list,optim_param,conditions_list)

        if rand_kin == True:
            self.randomize_kin(optim_param)

        self.fitness = 0


    def find_r2opt(self,red_data_list,optim_param,conditions_list):


        if optim_param.optim_on_meth!='False':
            n_tspc     = red_data_list[0].n_tspc
            n_r2opt    = optim_param.nb_r2opt                    # total number of react to opt
            self.mech.react.modif = [False]*len(self.mech.react.modif)
            mp         = optim_param.main_path

            # Get list of reactions to optimize
            if 'SA' in red_data_list[0].reduction_operator \
            or optim_param.optim_on_meth=='SA':
                # keep maximal sensitivities
                if len(red_data_list[0].red_op.sensi_r)!=0:
                    # ------------ 10/07/2023
                    add_col = 0
                    for l in range(len(red_data_list)):
                        if red_data_list[l].red_op.sensi_Sl is not False\
                        and conditions_list[l].error_param.Sl_check:
                            if red_data_list[l].red_op.sensi_T is not False\
                            and conditions_list[l].error_param.T_check:
                                if red_data_list[l].red_op.sensi_igt is not False\
                                and conditions_list[l].error_param.ig_check:
                                    # sensitivity analysis on: Sl, T, igt
                                    add_col = 3
                                    idx_Sl = -3 ; idx_T = -2 ; idx_igt = -1
                                else:
                                    # sensitivity analysis on: Sl, T
                                    add_col = max(add_col,2)
                                    idx_Sl = -2 ; idx_T = -1
                            elif red_data_list[l].red_op.sensi_igt is not False\
                            and conditions_list[l].error_param.ig_check:
                                # sensitivity analysis on: Sl, igt
                                add_col = max(add_col,2)
                                idx_Sl = -2 ; idx_igt = -1
                            else:
                                # sensitivity analysis on: Sl
                                add_col = max(add_col,1)
                                idx_Sl = -1
                        elif red_data_list[l].red_op.sensi_T is not False\
                        and conditions_list[l].error_param.T_check:
                            if red_data_list[l].red_op.sensi_igt is not False\
                            and conditions_list[l].error_param.ig_check:
                                # sensitivity analysis on: T, igt
                                add_col = max(add_col,2)
                                idx_T = -2 ; idx_igt = -1
                            else:
                                # sensitivity analysis on: T
                                add_col = max(add_col,1)
                                idx_T = -1
                        elif red_data_list[l].red_op.sensi_igt is not False\
                        and conditions_list[l].error_param.ig_check:
                            # sensitivity analysis on: igt
                            add_col = max(add_col,1)
                            idx_igt = -1

                    # create max_sens_list
                    max_sens_list=np.zeros((len(red_data_list[0].red_op.sensi_r)+add_col,len(red_data_list[0].red_op.sensi_r[0])))


                    # fill max_sens_list   - species sensitivities
                    for l in range(len(red_data_list)):
                        for idx in range(len(red_data_list[l].red_op.sensi_r)):
                            for r in range(len(red_data_list[l].red_op.sensi_r[idx])):
                                if abs(red_data_list[l].red_op.sensi_r[idx][r])>max_sens_list[idx][r]:
                                    max_sens_list[idx][r]=abs(red_data_list[l].red_op.sensi_r[idx][r])
                                # if l == len(red_data_list[l].red_op.sensi_r-1):
                                #     try:
                                #         if abs(red_data_list[l].red_op.sensi_Sl[r])>max_sens_list[idx][r]\
                                #         and conditions_list[l].error_param.Sl_check:
                                #             max_sens_list[idx][r]=abs(red_data_list[l].red_op.sensi_Sl[r])
                                #     except:
                                #         A=2
                                #     try:
                                #         if abs(red_data_list[l].red_op.sensi_T[r])>max_sens_list[idx][r]\
                                #         and conditions_list[l].error_param.T_check:
                                #             max_sens_list[idx][r]=abs(red_data_list[l].red_op.sensi_T[r])
                                #     except:
                                #         A=2
                                #     try:
                                #         if abs(red_data_list[l].red_op.sensi_igt[r])>max_sens_list[idx][r]\
                                #         and conditions_list[l].error_param.ig_check:
                                #             max_sens_list[idx][r]=abs(red_data_list[l].red_op.sensi_igt[r])
                                #     except:
                                #         A=2

                    # fill max_sens_list...
                    for l in range(len(red_data_list)):
                        #   - Sl sensitivities
                        if red_data_list[l].red_op.sensi_Sl is not False\
                        and conditions_list[l].error_param.Sl_check:
                            if len(max_sens_list[idx_Sl])==1:
                                max_sens_list[idx_Sl] = red_data_list[l].red_op.sensi_Sl
                            else:
                                for r in range(len(red_data_list[l].red_op.sensi_Sl)):
                                    if abs(red_data_list[l].red_op.sensi_Sl[r])>max_sens_list[idx_Sl][r]\
                                    and conditions_list[l].error_param.Sl_check:
                                        max_sens_list[idx_Sl][r]=abs(red_data_list[l].red_op.sensi_Sl[r])
                        #   - T sensitivities
                        if red_data_list[l].red_op.sensi_T is not False\
                        and conditions_list[l].error_param.T_check:
                            if len(max_sens_list[idx_T])==1:
                                max_sens_list[idx_T] = red_data_list[l].red_op.sensi_T
                            else:
                                for r in range(len(red_data_list[l].red_op.sensi_T)):
                                    if abs(red_data_list[l].red_op.sensi_T[r])>max_sens_list[idx_T][r]\
                                    and conditions_list[l].error_param.T_check:
                                        max_sens_list[idx_T][r]=abs(red_data_list[l].red_op.sensi_T[r])
                        #   - igt sensitivities
                        if red_data_list[l].red_op.sensi_igt is not False\
                        and conditions_list[l].error_param.ig_check:
                            if len(max_sens_list[idx_igt])==1:
                                max_sens_list[idx_igt] = red_data_list[l].red_op.sensi_igt
                            else:
                                for r in range(len(red_data_list[l].red_op.sensi_igt)):
                                    if abs(red_data_list[l].red_op.sensi_igt[r])>max_sens_list[idx_igt][r]\
                                    and conditions_list[l].error_param.ig_check:
                                        max_sens_list[idx_igt][r]=abs(red_data_list[l].red_op.sensi_igt[r])
                    if len(max_sens_list)==1:
                        print_('Warning, no sensitivity data',mp)

                    # ------------------------


# =============================================================================
#                     # --------------   OLD   --------------
#                     max_sens_list=np.zeros((len(red_data_list[0].red_op.sensi_r),len(red_data_list[0].red_op.sensi_r[0])))
#                     n_r2opt_sp_max = round(n_r2opt/len(red_data_list[0].red_op.sensi_r))    # number of react to opt per target data (spec / Sl / ...)
#                     for l in range(len(red_data_list)):
#                         for idx in range(len(red_data_list[l].red_op.sensi_r)):
#                             for r in range(len(red_data_list[l].red_op.sensi_r[idx])):
#                                 if abs(red_data_list[l].red_op.sensi_r[idx][r])>max_sens_list[idx][r]:
#                                     max_sens_list[idx][r]=abs(red_data_list[l].red_op.sensi_r[idx][r])
#                                 if l == len(red_data_list[l].red_op.sensi_r-1):
#                                     try:
#                                         if abs(red_data_list[l].red_op.sensi_Sl[r])>max_sens_list[idx][r]\
#                                         and conditions_list[l].error_param.Sl_check:
#                                             max_sens_list[idx][r]=abs(red_data_list[l].red_op.sensi_Sl[r])
#                                     except:
#                                         A=2
#                                     try:
#                                         if abs(red_data_list[l].red_op.sensi_T[r])>max_sens_list[idx][r]\
#                                         and conditions_list[l].error_param.T_check:
#                                             max_sens_list[idx][r]=abs(red_data_list[l].red_op.sensi_T[r])
#                                     except:
#                                         A=2
#                 else:
#                     max_sens_list = [0]
#                     for l in range(len(red_data_list)): #
#                         if red_data_list[l].red_op.sensi_Sl is not False\
#                         and conditions_list[l].error_param.Sl_check:
#                             if len(max_sens_list)==1:
#                                 max_sens_list = red_data_list[l].red_op.sensi_Sl
#                             else:
#                                 for r in range(len(red_data_list[l].red_op.sensi_Sl)):
#                                     if abs(red_data_list[l].red_op.sensi_Sl[r])>max_sens_list[r]\
#                                     and conditions_list[l].error_param.Sl_check:
#                                         max_sens_list[r]=abs(red_data_list[l].red_op.sensi_Sl[r])
#                         if red_data_list[l].red_op.sensi_T is not False\
#                         and conditions_list[l].error_param.T_check:
#                             if len(max_sens_list)==1:
#                                 max_sens_list = red_data_list[l].red_op.sensi_T
#                             else:
#                                 for r in range(len(red_data_list[l].red_op.sensi_T)):
#                                     if abs(red_data_list[l].red_op.sensi_T[r])>max_sens_list[r]\
#                                     and conditions_list[l].error_param.T_check:
#                                         max_sens_list[r]=abs(red_data_list[l].red_op.sensi_T[r])
#                     if len(max_sens_list)==1:
#                         print_('Warning, no sensitivity data',mp)
#
#                     # --------------   OLD   --------------
# =============================================================================


                n_r2opt_sp_max = round(n_r2opt/len(red_data_list[0].red_op.sensi_r))    # number of react to opt per target data (spec / Sl / ...)
                react2mod = []
                for idx in range(len(max_sens_list)):
                    react2mod.append([])
                    sensi_r=[]
                    for r in range(len(max_sens_list[idx])):
                        sensi_r.append((max_sens_list[idx][r],r+1))
                    sensi_r.sort()

                    n_r2opt_sp=0
                    for r1 in range(len(sensi_r)):
                        x = -(r1+1)
                        for r2 in range(len(max_sens_list[idx])):
                            if max_sens_list[idx][r2]==sensi_r[x][0]:  # find most sensitive reactions
                                if not self.mech.react.modif[r2]:
                                    # check if sub mechs can be modified:
                                    n_C = self.mech.react.subm_C[r2]
                                    sub_H = self.mech.react.subm_C[r2]==0 and not self.mech.react.subm_CO[r2]
                                    n_N = self.mech.react.subm_N[r2]
                                    n_S = self.mech.react.subm_S[r2]
                                    n_Si = self.mech.react.subm_Si[r2]

                                    if    n_C>0 or sub_H: sub_C = optim_param.opt_subm_C[n_C]
                                    else:                 sub_C = False
                                    if self.mech.react.subm_CO[r2]:
                                        sub_CO = optim_param.opt_subm_CO
                                    else:
                                        sub_CO = False
                                    if    n_N>0: sub_N = optim_param.opt_subm_N[n_N]
                                    else:        sub_N = False
                                    if    n_S>0: sub_S = optim_param.opt_subm_N[n_S]
                                    else:        sub_S = False
                                    if    n_Si>0: sub_Si = optim_param.opt_subm_Si[n_Si]
                                    else:         sub_Si = False
                                    # if sub mechs can be modified:
                                    if sub_C or sub_CO or sub_N or sub_S or sub_Si:
                                        self.mech.react.modif[r2]=True ; n_r2opt_sp+=1
                                        react2mod[-1].append(r2)
                            if n_r2opt_sp==n_r2opt_sp_max: break
                        if n_r2opt_sp==n_r2opt_sp_max: break

            if 'DRG' in red_data_list[0].reduction_operator \
            or optim_param.optim_on_meth=='DRG':
                if n_tspc != 0: n_r2opt_sp_max = round(n_r2opt/n_tspc)    # number of react to opt per target data (spec / Sl / ...)
                max_coeffs_list=np.zeros((n_tspc,len(red_data_list[0].red_op.r_interaction_coeffs[0])))
                for l in range(len(red_data_list)):
                    for idx in range(n_tspc):
                        for r in range(len(red_data_list[l].red_op.r_interaction_coeffs[idx])):
                            if abs(red_data_list[l].red_op.r_interaction_coeffs[idx][r])>max_coeffs_list[idx][r]:
                                max_coeffs_list[idx][r]=abs(red_data_list[l].red_op.r_interaction_coeffs[idx][r])


                react2mod = []
                for idx in range(n_tspc):
                    react2mod.append([])
                    drg_coeffs=[]
                    for j in range(len(max_coeffs_list[idx])):
                        drg_coeffs.append((max_coeffs_list[idx][j],j+1))
                    drg_coeffs.sort()

                    #DRG_sorted.sort()
                    n_r2opt_sp=0
                    for r1 in range(len(drg_coeffs)):
                        x = -(r1+1)
                        for r2 in range(len(max_coeffs_list[idx])):
                            if max_coeffs_list[idx][r2]==drg_coeffs[x][0]:
                                if not self.mech.react.modif[r2]:
                                    # check if sub mechs can be modified:
                                    n_C = self.mech.react.subm_C[r2]
                                    sub_H = self.mech.react.subm_C[r2]==0 and not self.mech.react.subm_CO[r2]
                                    n_N = self.mech.react.subm_N[r2]
                                    n_S = self.mech.react.subm_S[r2]
                                    n_Si = self.mech.react.subm_Si[r2]

                                    if    n_C>0 or sub_H: sub_C = optim_param.opt_subm_C[n_C]
                                    else:                 sub_C = False
                                    if self.mech.react.subm_CO[r2]:
                                        sub_CO = optim_param.opt_subm_CO
                                    else:
                                        sub_CO = False
                                    if    n_N>0: sub_N = optim_param.opt_subm_N[n_N]
                                    else:        sub_N = False
                                    if    n_S>0: sub_S = optim_param.opt_subm_N[n_S]
                                    else:        sub_S = False
                                    if    n_Si>0: sub_Si = optim_param.opt_subm_Si[n_Si]
                                    else:         sub_Si = False
                                    # if sub mechs can be modified:
                                    if sub_C or sub_CO or sub_N or sub_S or sub_Si:
                                        self.mech.react.modif[r2]=True ; n_r2opt_sp+=1
                                        react2mod[-1].append(r2)

                            if n_r2opt_sp==n_r2opt_sp_max: break
                        if n_r2opt_sp==n_r2opt_sp_max: break

            if optim_param.display_react2opt:
                n_r_opt=0
                for _r in range(len(self.mech.react.modif)):
                    if self.mech.react.modif[_r] == True:n_r_opt+=1


                print_('-------------------------------- \n'+str(n_r_opt)+' reactions to optimize:',mp)
                for _tn in range(len(red_data_list[0].tspc)):
                    print_(' * for target: '+red_data_list[0].tspc[_tn]+':',mp)
                    for _r in react2mod[_tn]:
                        if   int(self.mech.react.number[_r])<10:    spaces='    -   '
                        elif int(self.mech.react.number[_r])<100:   spaces='   -   '
                        elif int(self.mech.react.number[_r])<1000:  spaces='  -   '
                        elif int(self.mech.react.number[_r])<10000: spaces=' -   '
                        print_(str(self.mech.react.number[_r]) + spaces + self.mech.react.equation[_r],mp)

                # ------------ 10/07/2023
                if 'idx_Sl' in locals():
                    print_(' * for target Flame speed'+':',mp)
                    for _r in react2mod[idx_Sl]:
                        if   int(self.mech.react.number[_r])<10:    spaces='    -   '
                        elif int(self.mech.react.number[_r])<100:   spaces='   -   '
                        elif int(self.mech.react.number[_r])<1000:  spaces='  -   '
                        elif int(self.mech.react.number[_r])<10000: spaces=' -   '
                        print_(str(self.mech.react.number[_r]) + spaces + self.mech.react.equation[_r],mp)

                if 'idx_T' in locals():
                    print_(' * for target Temperature'+':',mp)
                    for _r in react2mod[idx_T]:
                        if   int(self.mech.react.number[_r])<10:    spaces='    -   '
                        elif int(self.mech.react.number[_r])<100:   spaces='   -   '
                        elif int(self.mech.react.number[_r])<1000:  spaces='  -   '
                        elif int(self.mech.react.number[_r])<10000: spaces=' -   '
                        print_(str(self.mech.react.number[_r]) + spaces + self.mech.react.equation[_r],mp)
                # ------------------------




#                for _r in range(len(self.mech.react.modif)):
#                    if self.mech.react.modif[_r] == True:
#                        if   int(self.mech.react.number[_r])<10:    spaces='    -   '
#                        elif int(self.mech.react.number[_r])<100:   spaces='   -   '
#                        elif int(self.mech.react.number[_r])<1000:  spaces='  -   '
#                        elif int(self.mech.react.number[_r])<10000: spaces=' -   '
#
#                        print_(str(self.mech.react.number[_r]) + spaces + self.mech.react.equation[_r],mp)
                print_('--------------------------------\n\n',mp)
                optim_param.display_react2opt=False


        else:   #selection of submech in gui


            if False in optim_param.opt_subm_C \
            or False in optim_param.opt_subm_N \
            or False in optim_param.opt_subm_S \
            or False in optim_param.opt_subm_Si\
            or optim_param.opt_subm_CO == False: # opt_subm_C : selection of submech in gui
                self.mech.react.modif = [True]*len(self.mech.react.modif)
                # CxHyOz sub-mechanisms (H2 submech is considered as a C0 submech)
                for n_C in range(len(optim_param.opt_subm_C)):
                    if optim_param.opt_subm_C[n_C]==False:
                        for r in range(len(self.mech.react.equation)):
                            if self.mech.react.subm_C[r] == n_C \
                            and not self.mech.react.subm_CO[r]  \
                            and self.mech.react.subm_N[r]  == 0 \
                            and self.mech.react.subm_S[r]  == 0 \
                            and self.mech.react.subm_Si[r] == 0 :
                                self.mech.react.modif[r] = False
                # CO sub-mechanism
                if optim_param.opt_subm_CO==False:
                    for r in range(len(self.mech.react.equation)):
                        if self.mech.react.subm_CO[r] == True:
                            self.mech.react.modif[r] = False
                # N sub-mechanisms
                for n_N_ in range(len(optim_param.opt_subm_N)-1):
                    n_N = n_N_ + 1
                    if optim_param.opt_subm_N[n_N]==False:
                        for r in range(len(self.mech.react.equation)):
                            if self.mech.react.subm_N[r] == n_N:
                                self.mech.react.modif[r] = False
                # S sub-mechanism
                if optim_param.opt_subm_S==False:
                    for r in range(len(self.mech.react.equation)):
                        if self.mech.react.subm_S[r] == True:
                            self.mech.react.modif[r] = False
                # Si sub-mechanism
                if optim_param.opt_subm_Si==False:
                    for r in range(len(self.mech.react.equation)):
                        if self.mech.react.subm_Si[r] == True:
                            self.mech.react.modif[r] = False
            else:
                self.mech.react.modif = [True]*len(self.mech.react.modif)

        if type(optim_param.reactions2opt) is not bool:
            if False not in self.mech.react.modif: # if no restriction is defined, prevent the modification of not specified reactions
                self.mech.react.modif = [False]*len(self.mech.react.modif)
            # allow the modification of specified reactions
            for r2mod in optim_param.reactions2opt:
                for r in range(len(self.mech.react.modif)):
                    if r2mod == r+1:
                        self.mech.react.modif[r] = True






#    def get_values(self,optim_param,n_ind):
#
#        n_ind+=1
#
#        # Directory
#        os.chdir(optim_param.main_path)
#        table = [] ; new_values =[] ; num_react=[]
#
#        try :
#            f=open('new_values.csv')
#            csv_f=csv.reader(f, delimiter=';')
#            for row in csv_f:
#                table.append(row)
#            f.close()
#
#            for i in range(len(table)-1):
#                num_react.append(table[i+1][0])      # [0] reaction number column
#                new_values.append(table[i+1][n_ind].split(','))    # [4] uncertainties column
#        except:
#            a=2
#            # print('Uncertainties.csv absent: all uncertainties are fixed to ',str(optim_param.Arrh_max_variation))
#
##        self.mech.react.incert = []
##        for r in range(len(self.mech.react.equation)):
##            if not react_found:
##                self.mech.react.incert.append(optim_param.Arrh_max_variation)
#
#        for r in range(len(self.mech.react.type)):
#            for r_inc in range(len(new_values)):
#                if int(num_react[r_inc])==self.mech.react.number[r]\
#                and new_values[r_inc]!='':
#                    if self.mech.react.type[r] == "three_body_reaction"\
#                    or self.mech.react.type[r] == "reaction":
#                        self.mech.react.kin[r]=new_values[r]
#
#                    elif self.mech.react.type[r] == "falloff_reaction"\
#                    or self.mech.react.type[r] == "chemically_activated_reaction"\
#                    or self.mech.react.type[r] == "chebyshev"\
#                    or self.mech.react.type[r] == "pdep_arrhenius":
#                        for k1 in range(len(self.mech.react.kin[r])):
#                            for k2 in range(len(self.mech.react.kin[r][k1])):
#                                self.mech.react.kin[r][k1] = new_values[r][k1*3:len(self.mech.react.kin[r][k1])]
#
#
#
#        os.chdir("GA")



    def shift_flame_data(self,conditions_list,optim_param,ref_results_list):
        verbose = conditions_list[0].simul_param.verbose

        os.chdir(conditions_list[0].main_path+'/GA')

        if '.cti' in conditions_list[0].mech:
            filename = 'temp.cti'
            self.mech.write_new_mech(filename)
        else:
            filename = 'temp.yaml'
            self.mech.write_yaml_mech(filename)


        # --------------------------------------------------------------------------------
        # interpretation of the new mech
        gas = cdef.get_gas_ct(filename)
        # --------------------------------------------------------------------------------

        qoi_tot = [] ; qoi_tot_pond =  [] ; pond=0


        for i in range(len(conditions_list)):
            conditions   = conditions_list[i]
            if conditions.config == 'free_flame':

                T_check  = conditions.error_param.T_check

                # Simulation conditions
                ref_results = ref_results_list[i]

                Opt_results = comp.red_computation(ref_results.conditions, \
                                   gas,self.mech.spec.activ_m,self.mech.react.activ_m)

                end_sim     = ref_results.conditions.simul_param.end_sim
                shifting    = -end_sim
                original_pts_scatter = np.array(ref_results.pts_scatter)
                pts_scatter = np.array(ref_results.pts_scatter)
                fitness     = 0
                for shif_it in range(100):
                    shifting    += end_sim/100
                    ref_results.pts_scatter            = pts_scatter+shifting
                    conditions.simul_param.pts_scatter = pts_scatter+shifting
                    errors = cdef.Errors(conditions,ref_results,Opt_results,\
                                     optim_param)
                    for sp in range(optim_param.n_tspc):
                        qoi_tot.append(errors.qoi_s[sp])
                        qoi_tot_pond.append(errors.qoi_s[sp]*optim_param.coeff_s[sp])
                        pond+=optim_param.coeff_s[sp]
                    if T_check:
                        qoi_tot.append(errors.qoi_T)
                        qoi_tot_pond.append(errors.qoi_T*optim_param.coeff_T)
                        pond+=optim_param.coeff_T
                    if conditions.error_param.error_type_fit == 'mean':
                         fitness_i = 1/(np.sum(qoi_tot)/pond)
                    elif conditions.error_param.error_type_fit == 'max':
                         fitness_i = 1/np.max(qoi_tot)
                    if fitness_i>fitness:
                        best_shift = shifting
                ref_results.pts_scatter            = original_pts_scatter+best_shift
                conditions.simul_param.pts_scatter = original_pts_scatter+best_shift
                ref_results.conditions.simul_param.shift = best_shift

        return conditions_list,ref_results_list

    def time_step_optim(self,conditions_list,ref_results_list):
        mp = conditions_list[0].main_path
        verbose = conditions_list[0].simul_param.verbose

        print_('time step optimization',mp)
        os.chdir(conditions_list[0].main_path+'/GA')
        if '.cti' in conditions_list[0].mech:
            filename = 'temp.cti'
            self.mech.write_new_mech(filename)
        else:
            filename = 'temp.yaml'
            self.mech.write_yaml_mech(filename)

        # --------------------------------------------------------------------------------
        # interpretation of the new mech

        for i in range(len(conditions_list)):
            conditions   = conditions_list[i]
            conditions.composition.gas = cdef.get_gas_ct(filename)
            if 'reactor' in conditions.config:
                opt_results, conditions = comp.ref_computation(conditions)
                ref_results = comp.red_computation(conditions, \
                            conditions.composition.gas_ref,\
                            self.mech.spec.activ_m,self.mech.react.activ_m)
                ref_results_list[i] = ref_results
                conditions_list[i]  = conditions
        # --------------------------------------------------------------------------------


        return conditions_list, ref_results_list


    def randomize_kin(self,optim_param):

        for r in range(len(self.mech.react.type)):
            try_r = 0 ; valid = False ; damping = 0.75 ; max_try = 10
            while not valid and try_r < max_try:
                var_range = 1-(try_r/max_try)
                if self.mech.react.modif[r]:
                    incert_r = [u/100 for u in self.mech.react.incert[r]]
                    if self.mech.react.type[r] == "three_body_reaction"\
                    or self.mech.react.type[r] == "reaction":
                        for k in range(len(self.mech.react.kin[r])):
                            self.mech.react.kin[r][k]=self.mech.react.ref_kin[r][k]\
                            +self.mech.react.ref_kin[r][k]*random.uniform(-var_range,var_range)*incert_r[k]*(damping**try_r)

                    elif self.mech.react.type[r] == "falloff_reaction"\
                    or self.mech.react.type[r] == "pdep_arrhenius"\
                    or self.mech.react.type[r] == "chemically_activated_reaction"\
                    or self.mech.react.type[r] == "chebyshev":
                        for k1 in range(len(self.mech.react.kin[r])):
                            for k2 in range(len(self.mech.react.kin[r][k1])):
                                self.mech.react.kin[r][k1][k2]=self.mech.react.ref_kin[r][k1][k2]\
                                +self.mech.react.ref_kin[r][k1][k2]*random.uniform(-var_range,var_range)*incert_r[k2]*(damping**try_r)

                    valid = check_k(self.mech.react,r)
                else:
                    valid = True

                try_r +=1

            if not valid:
                self.mech.react.kin[r] = copy.deepcopy(self.mech.react.ref_kin[r])





    def fitness_eval(self,conditions_list,optim_param,ref_results_list,n_par=0,ref_ind=False):

        verbose = conditions_list[0].simul_param.verbose
        mp = conditions_list[0].main_path
        os.chdir(conditions_list[0].main_path+'/GA')

        filename = 'temp_'+str(n_par)
        # if self.mech.keep4opt == True:
        #     filename = 'keep_' + filename
        if '.cti' in conditions_list[0].mech:
            filename += '.cti'
            self.mech.write_new_mech(filename)
        else:
            filename += '.yaml'
            self.mech.write_yaml_mech(filename)


        # --------------------------------------------------------------------------------
        # interpretation of the new mech

        # suppress console output during the interpretation
        # if verbose<9:
        #     old_stdout = sys.stdout ; old_stderr = sys.stderr
        #     with open(os.devnull, "w") as devnull: sys.stdout = devnull ; sys.stderr = devnull

        gas = cdef.get_gas_ct(filename)
        
        #restore console output
        # if verbose<9: sys.stdout = old_stdout ; sys.stderr = old_stderr
        # --------------------------------------------------------------------------------
        
        txt_f='\n'
        qoi_tot = [] ; pond=0 ; qoi_tot_mean=0
        for i in range(len(conditions_list)):

            # -------------------------------
            # Reduction loop

            T_check  = conditions_list[i].error_param.T_check
            Sl_check = conditions_list[i].error_param.Sl_check
            ig_check = conditions_list[i].error_param.ig_check
            K_check  = conditions_list[i].error_param.K_check


            # Simulation conditions
            conditions                      = conditions_list[i]
            conditions.simul_param.par_ind  = str(n_par)
            ref_results                     = ref_results_list[i]

            cur_path = os.getcwd()
            Opt_results = comp.red_computation(conditions,gas, \
                               self.mech.spec.activ_m,self.mech.react.activ_m)
            os.chdir(cur_path)

            # ####### DEV 
            # if ref_ind:
            #     ref_results.opt_refind_T            = copy.deepcopy(Opt_results.T)
            #     ref_results.opt_refind_conc         = copy.deepcopy(Opt_results.conc)
            #     ref_results.opt_refind_X            = copy.deepcopy(Opt_results.X)
            #     ref_results.opt_refind_ign_time_hr  = copy.deepcopy(Opt_results.ign_time_hr)
            #     ref_results.opt_refind_ign_time_sp  = copy.deepcopy(Opt_results.ign_time_sp)
            #     ref_results.opt_refind_ign_time     = copy.deepcopy(Opt_results.ign_time)
            #     ref_results.opt_refind_Sl           = copy.deepcopy(Opt_results.Sl)
            #     ref_results.opt_refind_K_ext        = copy.deepcopy(Opt_results.K_ext)

            errors = cdef.Errors(conditions,ref_results,Opt_results,\
                                 optim_param)

            # Condition fitness weighting
            if optim_param.coeff_cond:  coeff_cond = optim_param.coeff_cond[i]
            else:                       coeff_cond = 1

            qoi_case, pond_case = 0,0
            cc = conditions.config
            for sp in range(optim_param.n_tspc):
                if errors.qoi_s[sp]:
                    pond_i = optim_param.coeff_s[sp]*coeff_cond
                    txt_f+='Fit ' + cc + '  sp - ' + optim_param.tspc[sp] + ' : ' + "%.3f" %(1-errors.qoi_s[sp]) + '  pond = ' + "%.1f" %pond_i + '\n'
                    qoi_tot_mean += (1-errors.qoi_s[sp])*pond_i
                    qoi_case     += (1-errors.qoi_s[sp])*pond_i
                    qoi_tot.append((1-errors.qoi_s[sp])*min(1,np.ceil(pond_i)))
                    pond+=pond_i ; pond_case+=pond_i
            if 'JSR' not in conditions.config and T_check  and errors.qoi_T:
                pond_i = optim_param.coeff_T*coeff_cond
                txt_f+='Fit ' + cc + ' - T: ' + "%.3f" %(1-errors.qoi_T) + '  pond = ' + "%.1f" %pond_i + '\n'
                qoi_tot_mean += (1-errors.qoi_T)*pond_i
                qoi_case     += (1-errors.qoi_T)*pond_i
                qoi_tot.append((1-errors.qoi_T)*min(1,np.ceil(pond_i)))
                pond+=pond_i ; pond_case+=pond_i
            if 'reactor' in conditions.config and ig_check and errors.qoi_ig:
                pond_i = optim_param.coeff_ig*coeff_cond
                txt_f+='Fit reactor - igt: ' + "%.3f" %(1-errors.qoi_ig) + '  pond = ' + "%.1f" %pond_i + '\n'
                qoi_tot_mean += (1-errors.qoi_ig)*pond_i
                qoi_case     += (1-errors.qoi_ig)*pond_i
                qoi_tot.append((1-errors.qoi_ig)*min(1,np.ceil(pond_i)))
                pond+=pond_i ; pond_case+=pond_i
            if 'free_flame' in conditions.config and Sl_check and errors.qoi_Sl:
                pond_i = optim_param.coeff_Sl*coeff_cond
                txt_f+='Fit free_flame - Sl: ' + "%.3f" %(1-errors.qoi_Sl) + '  pond = ' + "%.1f" %pond_i + '\n'
                qoi_tot_mean += (1-errors.qoi_Sl)*pond_i
                qoi_case     += (1-errors.qoi_Sl)*pond_i
                qoi_tot.append((1-errors.qoi_Sl)*min(1,np.ceil(pond_i)))
                pond+=pond_i ; pond_case+=pond_i
            if ('diff_flame' in conditions.config or 'pp_flame' in conditions.config) and K_check and errors.qoi_K:
                pond_i = optim_param.coeff_K*coeff_cond
                txt_f+='Fit diff_flame - K: ' + "%.3f" %(1-errors.qoi_K) + '  pond = ' + "%.1f" %pond_i + '\n'
                qoi_tot_mean += (1-errors.qoi_K)*pond_i
                qoi_case     += (1-errors.qoi_K)*pond_i
                qoi_tot.append((1-errors.qoi_K)*min(1,np.ceil(pond_i)))
                pond+=pond_i ; pond_case+=pond_i

            # detailed informations on initial fitness (for each case)
#            print_detailed_fit = True
#            if print_detailed_fit:
#                fit_i =  1/(qoi_case/pond_case)
#                txt_fit = 'Fitness condition ' + str(i+1) + ': ' + str(fit_i) \
#                          + '   pond: ' + str(pond_case)
#                print_(txt_fit, mp)
        if verbose>=5:
             print_(txt_f,mp)       
 
        if 'no data' in qoi_tot:  qoi_tot.remove('no data')
        if conditions.error_param.error_type_fit == 'mean':
#             self.fitness = 1/(np.sum(qoi_tot)/pond)
             # fitness = 1/(qoi_tot_mean/pond)
             fitness = qoi_tot_mean/pond
             if verbose>=5:
                  print_('Fitness of the individual= ' + '%.3f' %fitness,mp)
        elif conditions.error_param.error_type_fit == 'max':
             fitness = 1-np.max(qoi_tot)

        # check nan
        if fitness != fitness:
            fitness = 0

        return max(fitness,0)

    def export_data(self,conditions_list,optim_param,ref_results_list,filename='temp.cti'):
        verbose = conditions_list[0].simul_param.verbose

        errors_list=[] ; Opt_results_list = []

        qoi_tot = [] ; pond=0

        os.chdir(conditions_list[0].main_path+'/GA')
        if '.cti' in conditions_list[0].mech:
            self.mech.write_new_mech(filename)
        else:
            if '.cti' in filename:
                filename = filename[:-4] + '.yaml'
            self.mech.write_yaml_mech(filename)

        # --------------------------------------------------------------------------------
        # interpretation of the new mech

        ct.suppress_thermo_warnings()
        warnings.filterwarnings("ignore", category=DeprecationWarning)
        gas = cdef.get_gas_ct(filename)

        # --------------------------------------------------------------------------------


        for i in range(len(conditions_list)):

            # -------------------------------
            # Reduction loop

            # Simulation condition
            conditions   = conditions_list[i]
            ref_results = ref_results_list[i]
            os.chdir(conditions_list[0].main_path+'/GA')
            if '.cti' in filename:
                self.mech.write_new_mech(filename)
            else:
                self.mech.write_yaml_mech(filename)

            T_check  = conditions_list[i].error_param.T_check
            Sl_check = conditions_list[i].error_param.Sl_check
            ig_check = conditions_list[i].error_param.ig_check
            K_check  = conditions_list[i].error_param.K_check

            Opt_results_list.append(comp.red_computation(conditions, \
                       gas,self.mech.spec.activ_m,self.mech.react.activ_m))
            os.chdir(conditions_list[0].main_path+'/GA')
            errors_list.append(cdef.Errors(conditions,ref_results,\
                                         Opt_results_list[-1],optim_param))
            for sp in range(optim_param.n_tspc):
                if errors_list[-1].qoi_s[sp]:
                    qoi_tot.append(errors_list[-1].qoi_s[sp])
                    pond+=optim_param.coeff_s[sp]
            if 'JSR' not in conditions.config and T_check  and errors_list[-1].qoi_T:
                qoi_tot.append(errors_list[-1].qoi_T)
                pond+=optim_param.coeff_T
            if 'reactor' in conditions.config and ig_check and errors_list[-1].qoi_ig:
                qoi_tot.append(errors_list[-1].qoi_ig)
                pond+=optim_param.coeff_ig
            if 'free_flame' in conditions.config and Sl_check and errors_list[-1].qoi_Sl:
                qoi_tot.append(errors_list[-1].qoi_Sl)
                pond+=optim_param.coeff_Sl
            if ('diff_flame' in conditions.config or 'pp_flame' in conditions.config) and K_check and errors_list[-1].qoi_K:
                qoi_tot.append(errors_list[-1].qoi_K)
                pond+=optim_param.coeff_K

        if 'no data' in qoi_tot:  qoi_tot.remove('no data')
        if conditions.error_param.error_type_fit == 'mean':
             self.fitness = 1/np.mean(qoi_tot)
        elif conditions.error_param.error_type_fit == 'max':
             self.fitness = 1/np.max(qoi_tot)

        return Opt_results_list, errors_list, self.fitness


class Population:
    def __init__(self,conditions_list,mech_data,red_data_list,ref_results_list,\
                 size_pop):

        optim_param = red_data_list[0].optim_param
        mp          = optim_param.main_path


        self.population = []

        mech_data = self.get_uncertainty(mech_data,optim_param,conditions_list)


        # option to import modified mecanisms in the first population
        if optim_param.import_mech:
            rand_kin = optim_param.import_mech
            cur_path = os.getcwd()
            os.chdir(mp+'/../_kinetic_mech/'+optim_param.import_mech)
            print_(mp+'/../_kinetic_mech/'+optim_param.import_mech,mp)
            print_('Import kinetic mechanisms from:',mp)
            list_files = os.listdir()
            list_mech = []
            for _file in list_files:
                if '.cti' in _file or '.yaml' in _file:
                    list_mech.append(cdef.Mech_data(_file))

                    list_mech[-1].react.incert      = copy.deepcopy(mech_data.react.incert)

                    list_mech[-1].react.f_max       = copy.deepcopy(mech_data.react.f_max)
                    list_mech[-1].react.f_T_fit     = copy.deepcopy(mech_data.react.f_T_fit)
                    list_mech[-1].react.k0_fit      = copy.deepcopy(mech_data.react.k0_fit)
                    list_mech[-1].react.f_max_lp    = copy.deepcopy(mech_data.react.f_max_lp)
                    list_mech[-1].react.f_T_fit_lp  = copy.deepcopy(mech_data.react.f_T_fit_lp)
                    list_mech[-1].react.k0_fit_lp   = copy.deepcopy(mech_data.react.k0_fit_lp)
                    list_mech[-1].react.f_max_hp    = copy.deepcopy(mech_data.react.f_max_hp)
                    list_mech[-1].react.f_T_fit_hp  = copy.deepcopy(mech_data.react.f_T_fit_hp)
                    list_mech[-1].react.k0_fit_hp   = copy.deepcopy(mech_data.react.k0_fit_hp)

                    list_mech[-1].react.k_ref       = copy.deepcopy(mech_data.react.k_ref)
                    list_mech[-1].react.k_min       = copy.deepcopy(mech_data.react.k_min)
                    list_mech[-1].react.k_max       = copy.deepcopy(mech_data.react.k_max)
                    list_mech[-1].react.k_ref_lp    = copy.deepcopy(mech_data.react.k_ref_lp)
                    list_mech[-1].react.k_min_lp    = copy.deepcopy(mech_data.react.k_min_lp)
                    list_mech[-1].react.k_max_lp    = copy.deepcopy(mech_data.react.k_max_lp)
                    list_mech[-1].react.k_ref_hp    = copy.deepcopy(mech_data.react.k_ref_hp)
                    list_mech[-1].react.k_min_hp    = copy.deepcopy(mech_data.react.k_min_hp)
                    list_mech[-1].react.k_max_hp    = copy.deepcopy(mech_data.react.k_max_hp)
                    if 'keep' in _file.lower():
                        list_mech[-1].keep4opt = True

                    print_('Imported mechanism: ' + _file,mp)
            os.chdir(cur_path)


        for ind in range(size_pop):
            # if import_mech, get the mecanism to import in the first population
            try:    rand_kin = list_mech[ind]
            except: rand_kin = True

            self.population.append(Chromosome(conditions_list,mech_data,\
                                   ref_results_list,red_data_list,rand_kin))


#        optim_from_values = False
#        if optim_from_values:
#            for ind in range(size_pop):
#                self.get_values(red_data_list[0].optim_param,ind)


    def __getitem__(self, i):
        return self.population[i]



    def get_uncertainty(self,mech_data,optim_param,conditions_list):

        # Directory
        os.chdir(optim_param.main_path)
        table = [] ; uncertainty =[] ; num_react=[]


        if 'uncertainties.csv' in os.listdir():
            f=open('uncertainties.csv')
            csv_f=csv.reader(f, delimiter=';')
            for row in csv_f:
                table.append(row)
            f.close()

            for i in range(len(table)-1):
                num_react.append(table[i+1][0])      # [0] reaction number column
                uncertainty.append(table[i+1][4])    # [4] uncertainties column
        else:
            a=2
            # print('Uncertainties.csv absent: all uncertainties are fixed to ',str(optim_param.Arrh_max_variation))

        mech_data.react.incert = []
        for r in range(len(mech_data.react.equation)):
            react_found=False
            for r_inc in range(len(uncertainty)):
                if int(num_react[r_inc])==mech_data.react.number[r]\
                and uncertainty[r_inc]!='':
                    mech_data.react.incert.append([float(u) for u in uncertainty[r_inc].split(',')])
                    react_found = True
                    break
            if not react_found:
                mech_data.react.incert.append(optim_param.Arrh_max_variation)


        # =============================================================================
        # Incertitude _f
        # =============================================================================


        # 1- lire fichier csv
        if 'uncertainties_f.csv' in os.listdir():
            f=open('uncertainties_f.csv')
            csv_f     = csv.reader(f, delimiter=';')
            react_f   = {'equation':[],
                         'f_max':[],'f_T_fit':[],'k0_fit':[],
                         'f_max_lp':[],'f_T_fit_lp':[],'k0_fit_lp':[],
                         'f_max_hp':[],'f_T_fit_hp':[],'k0_fit_hp':[]}
            # new_react = False
            for row in csv_f:
                if 'T_it:' in row:
                    T_min  = int(row[2])
                    T_max  = int(row[3])
                    T_incr = int(row[4])
                    mech_data.react.f_Tit = [T_min,T_max,T_incr]
#                if "'---------------------------" in row:
#                    new_react = True
#                if new_react:
                if 'equation:' in row:
                    react_f['equation'].append(row[3])
                    react_f['f_max'].append(False)
                    react_f['f_max_lp'].append(False)
                    react_f['f_max_hp'].append(False)                    
                    react_f['f_T_fit'].append(False)
                    react_f['f_T_fit_lp'].append(False)
                    react_f['f_T_fit_hp'].append(False)
                    react_f['k0_fit'].append(False)
                    react_f['k0_fit_lp'].append(False)
                    react_f['k0_fit_hp'].append(False)
                    if 'h + ho2 => h2 + o2' in row[3]:
                        print('toto')
                    
                # --- atmospheric pressure
                if 'f max:' in row:
                    f_max = row[3]      # row.split(';')[3]
                    if f_max == 'False':         react_f['f_max'][-1] = optim_param.f_default
                    else:                        react_f['f_max'][-1] = f_max
                if 'f_T_fit:' in row:
                    f_T_fit = []
                    if row[3] != 'False' :
                        for _t in range(9):
                            f_T_fit.append(np.float32(row[3+_t]))
                        react_f['f_T_fit'][-1] = f_T_fit
                    else:
                        react_f['f_T_fit'][-1] = [optim_param.f_default,0,0,0,0,0,0]
                if 'k0_fit:' in row:
                    k0_fit = []
                    if row[3] != 'False' :
                        for _t in range(9):
                            k0_fit.append(np.float32(row[3+_t]))
                        react_f['k0_fit'][-1] = k0_fit
                    else:
                        react_f['k0_fit'][-1] = False
#                        new_react = False
                # --- low pressure
                if 'f max (low P):' in row:
                    f_max_lp = row[3]       # row.split(';')[3]
                    if f_max_lp == 'False':      react_f['f_max_lp'][-1] = optim_param.f_default
                    else:                        react_f['f_max_lp'][-1] = f_max_lp
                if 'f_T_fit_lp:' in row:
                    f_T_fit_lp = []
                    if row[3] != 'False' :
                        for _t in range(9):
                            f_T_fit_lp.append(np.float32(row[3+_t]))
                        react_f['f_T_fit_lp'][-1] = f_T_fit_lp
                    else:
                        react_f['f_T_fit_lp'][-1] = [optim_param.f_default,0,0,0,0,0,0]
                if 'k0_fit_lp:' in row:
                    k0_fit_lp = []
                    if row[3] != 'False' :
                        for _t in range(9):
                            k0_fit_lp.append(np.float32(row[3+_t]))
                        react_f['k0_fit_lp'][-1] = k0_fit_lp
                    else:
                        react_f['k0_fit_lp'][-1] = False
                # --- high pressure
                if 'f max (high P):' in row:
                    f_max_hp = row[3]     # row.split(';')[3]
                    if f_max_hp == 'False':      react_f['f_max_hp'][-1] = optim_param.f_default
                    else:                        react_f['f_max_hp'][-1] = f_max_hp
                if 'f_T_fit_hp:' in row:
                    f_T_fit_hp = []
                    if row[3] != 'False' :
                        for _t in range(9):
                            f_T_fit_hp.append(np.float32(row[3+_t]))
                        react_f['f_T_fit_hp'][-1] = f_T_fit_hp
                    else:
                        react_f['f_T_fit_hp'][-1] = [optim_param.f_default,0,0,0,0,0,0]
                if 'k0_fit_hp:' in row:
                    k0_fit_hp = []
                    if row[3] != 'False' :
                        for _t in range(9):
                            k0_fit_hp.append(np.float32(row[3+_t]))
                        react_f['k0_fit_hp'][-1] = k0_fit_hp
                    else:
                        react_f['k0_fit_hp'][-1] = False
#                        new_react = False

            f.close()


            # 2b - comparer les réactifs du méca avec ceux du fichier
            for _r1,r1_equation in enumerate(mech_data.react.equation):
                if _r1==15:
                    print('toto')
                if '<=>' in r1_equation:
                    reactants_m = r1_equation.split('<=>')[0].upper()
                    products_m  = r1_equation.split('<=>')[1].upper()
                else:
                    reactants_m = r1_equation.split('=>')[0].upper()
                    products_m  = r1_equation.split('=>')[1].upper()
                reactants_meca = cdef.get_react_species(reactants_m)
                products_meca = cdef.get_react_species(products_m)
                reaction_found   = False
                for _r2,r2_equation in enumerate(react_f['equation']):
                    if '<=>' in r2_equation:
                        reactants_f = r2_equation.split('<=>')[0].upper()
                        products_f  = r2_equation.split('<=>')[1].upper()
                    else:
                        reactants_f = r2_equation.split('=>')[0].upper()
                        products_f  = r2_equation.split('=>')[1].upper()
                    reactants_file = cdef.get_react_species(reactants_f)
                    products_file  = cdef.get_react_species(products_f)
                    # check identical reactions
                    if (set(reactants_file) == set(reactants_meca) and set(products_file) == set(products_meca)) \
                    or (set(reactants_file) == set(products_meca)  and set(products_file) == set(reactants_meca) and '<=>' in r2_equation):
                        mech_data.react.f_max.append(react_f['f_max'][_r2])
                        mech_data.react.f_T_fit.append(react_f['f_T_fit'][_r2])
                        mech_data.react.k0_fit.append(react_f['k0_fit'][_r2])
                        mech_data.react.f_max_lp.append(react_f['f_max_lp'][_r2])
                        mech_data.react.f_T_fit_lp.append(react_f['f_T_fit_lp'][_r2])
                        mech_data.react.k0_fit_lp.append(react_f['k0_fit_lp'][_r2])
                        mech_data.react.f_max_hp.append(react_f['f_max_hp'][_r2])
                        mech_data.react.f_T_fit_hp.append(react_f['f_T_fit_hp'][_r2])
                        mech_data.react.k0_fit_hp.append(react_f['k0_fit_hp'][_r2])
                        reaction_found = True
                        break
                    # else:
                    #     if r1_equation.upper() == r2_equation.upper():
                    #         mech_data.react.f_max.append(react_f['f_max'][_r2])
                    #         mech_data.react.f_T_fit.append(react_f['f_T_fit'][_r2])
                    #         mech_data.react.k0_fit.append(react_f['k0_fit'][_r2])
                    #         mech_data.react.f_max_lp.append(react_f['f_max_lp'][_r2])
                    #         mech_data.react.f_T_fit_lp.append(react_f['f_T_fit_lp'][_r2])
                    #         mech_data.react.k0_fit_lp.append(react_f['k0_fit_lp'][_r2])
                    #         mech_data.react.f_max_hp.append(react_f['f_max_hp'][_r2])
                    #         mech_data.react.f_T_fit_hp.append(react_f['f_T_fit_hp'][_r2])
                    #         mech_data.react.k0_fit_hp.append(react_f['k0_fit_hp'][_r2])
                    #         reaction_found = True
                    #         break

                if not reaction_found:
                    mech_data.react.f_max.append(optim_param.f_default)
                    mech_data.react.f_T_fit.append(optim_param.f_default)
                    mech_data.react.k0_fit.append(False)
                    mech_data.react.f_max_lp.append(optim_param.f_default)
                    mech_data.react.f_T_fit_lp.append(optim_param.f_default)
                    mech_data.react.k0_fit_lp.append(False)
                    mech_data.react.f_max_hp.append(optim_param.f_default)
                    mech_data.react.f_T_fit_hp.append(optim_param.f_default)
                    mech_data.react.k0_fit_hp.append(False)

        else:
            for r in range(len(mech_data.react.equation)):
                mech_data.react.f_max.append(optim_param.f_default)
                mech_data.react.f_T_fit.append(optim_param.f_default)
                mech_data.react.k0_fit.append(False)
                mech_data.react.f_max_lp.append(optim_param.f_default)
                mech_data.react.f_T_fit_lp.append(optim_param.f_default)
                mech_data.react.k0_fit_lp.append(False)
                mech_data.react.f_max_hp.append(optim_param.f_default)
                mech_data.react.f_T_fit_hp.append(optim_param.f_default)
                mech_data.react.k0_fit_hp.append(False)

        p_min  = 5e4    # Pa
        p_max  = 5e6    # Pa
        for i in range(len(conditions_list)):
            if conditions_list[i].state_var.P < p_min:
                p_min = conditions_list[i].state_var.P
            if conditions_list[i].state_var.P > p_max:
                p_max = conditions_list[i].state_var.P
        mech_data = compute_kref_and_klim(mech_data,p_min,p_max)

        os.chdir("GA")

        return mech_data



    def check_early_conv(self,optim_param,gen,verbose):

        mp      = optim_param.main_path
        MaxIt   = optim_param.n_gen

        # calculation of the fitness stats
        fitness_list = []
        for p in range(len(self.population)):
            fitness_list.append(self.population[p].fitness)
        best_fitness = np.max(fitness_list) ; std_fit = np.std(fitness_list)

        fitness_list = np.array((fitness_list))

        if gen < .9*MaxIt and std_fit<(.1*best_fitness):
            print_('Early convergence detected, create new random individuals',mp)
            self.sort_fitness()

            for p_i in range(round(len(self.population)/2)):
                 self.population[p_i].randomize_kin(optim_param)



    def display(self, optim_param, ind = -3):
        mp = optim_param.main_path
        """ Affichage d'un ou de tout les individus de la population"""
        if ind == -3:
            print_("Fitness of the new population:\n",mp)
            fit_disp=""
            for i in range(optim_param.n_ind):
                fit_disp = fit_disp + '%5.3f' %(self.population[i].fitness) + ' '
                if (i+1)%10==0 and i!=0:
                    print_(fit_disp,mp)
                    fit_disp=""
            print_(fit_disp,mp)
        else:
            print_("",mp)
            print_('individu'+str(ind)+': '+"%.3f" %(self.population[ind].fitness),mp)
            print_("\n",mp)

    def compare_best_ind(self, best_ind, optim_param, main_path, verbose=0):
        mp = optim_param.main_path

        mp = optim_param.main_path
        n_ind = optim_param.n_ind
        best_idx = self.find_best()

        new_best_ind = False
        # compare the current best population ind to the previous best ind
        if self.population[best_idx].fitness > best_ind.fitness:
            best_ind = copy.deepcopy(self.population[best_idx])
            os.chdir(mp+'/GA')
            if '.cti' in self.population[0].mech.name:
                best_ind.mech.write_new_mech("optim_mech.cti")
            else:
                best_ind.mech.write_yaml_mech("optim_mech.yaml")

            if verbose >= 3:
                print_("New best_ind: "+"%.3f" %(best_ind.fitness),mp)
            new_best_ind = True
        best_idx = self.find_best(n_ind)

        if self.population[best_idx].fitness < best_ind.fitness:
            worst_idx=self.find_worst(best_idx,n_ind)
            self.population[worst_idx]=copy.deepcopy(best_ind)

        return best_ind, new_best_ind


    def def_pop2keep(self):
        pop2keep = []
        for ind in self.population:
            if ind.mech.keep4opt == True:
                pop2keep.append(copy.deepcopy(ind))
        return pop2keep


    def insert_pop2keep(self, optim_param, pop2keep):
        n_ind = optim_param.n_ind

        # Fitness list
        list_fit = []
        for _p in range(n_ind):
            list_fit.append(self.population[_p].fitness)

        # Get sorted element indices
        sort_index = sorted(range(len(list_fit)), key=lambda k: list_fit[k])

        # Create the order list from the sorted indices
        order = [0] * len(list_fit)
        for i, _index in enumerate(sort_index, start=1):
            order[_index] = i

        for _p2k, mech_2keep in enumerate(pop2keep):
            self.population[order[_p2k]]=copy.deepcopy(mech_2keep)


    def find_best(self,n_ind=False):
        if not n_ind: n_ind=len(self.population)
        best_fit=-1
        for p in range(n_ind):
            if self.population[p].fitness > best_fit:
                best_idx = p ; best_fit = self.population[p].fitness

        try:     A = best_idx+1-1
        except:  best_idx = 0

        return best_idx

    def find_worst(self,best_idx,n_ind=False):
        if not n_ind: n_ind=len(self.population)
        best_fit = self.population[best_idx].fitness
        worst_idx=0
        for p in range(n_ind):
            if self.population[p].fitness < best_fit:
                worst_idx = p ; best_fit = self.population[p].fitness
        return worst_idx

    def fitness_eval_par(self,fit_eval_inp):

        optim_param         = fit_eval_inp[0]
        conditions_list     = fit_eval_inp[1]
        ref_results_list    = fit_eval_inp[2]
        bar                 = fit_eval_inp[3]
        ind                 = fit_eval_inp[4]
        is_child            = fit_eval_inp[5]
        if is_child:
            ind               = optim_param.n_ind + ind
            title="New ind evaluation  "
        else:
            title="New pop evaluation  "

        mech = conditions_list[0].mech

        # change dir
        os.chdir(conditions_list[0].main_path)

        gas = cdef.get_gas_ct(mech)

        for _c in range(len(conditions_list)):
            conditions_list[_c].composition.gas     = gas
            conditions_list[_c].composition.gas_ref = gas

        os.chdir("GA")

        try:
            fitness = self.population[ind].fitness_eval(conditions_list,optim_param,ref_results_list,ind)
        except:
            fitness = 0

        bar.update(ind,title)

        return (fitness,ind)

    def fitness_eval_newchilds(self,optim_param,conditions_list,ref_results_list):

        mp          = optim_param.main_path
        num_cores   = multiprocessing.cpu_count()
        child_nb    = optim_param.total_Xover + optim_param.total_mut

        gas         = conditions_list[0].composition.gas
        gas_ref     = conditions_list[0].composition.gas_ref

        # saving and suppression of unpickable variables on fitness eval inputs
        for cond in range(len(conditions_list)):
            del conditions_list[cond].composition.gas
            del conditions_list[cond].composition.gas_ref
        f=[]
        simul_time_limit = 0
        for res in range(len(ref_results_list)):
            f.append(ref_results_list[res].f)
            del ref_results_list[res].gas
            del ref_results_list[res].f
            simul_time_limit += ref_results_list[res].simul_time
#        simul_time_limit = simul_time_limit*(6+np.random.uniform()*16)#*(child_nb/num_cores)
        simul_time_limit = simul_time_limit*2*(child_nb/num_cores) + 30
#        simul_time_limit = simul_time_limit*5 + 30


        bar = cdef.ProgressBar(child_nb, '')
        title="New ind evaluation  "
        bar.update(0,title)

        fit_eval_inp = []
        for ch in range(child_nb):
            fit_eval_inp.append([optim_param,conditions_list,ref_results_list,bar,ch,True])

        # saving of the ref mechanism
        mech = conditions_list[0].mech
        os.chdir(conditions_list[0].main_path)
        copyfile(mech,'GA/'+mech)
        os.chdir("GA")


        # Fitness calculation

        if sys.gettrace() is not None : is_debug_mode = True
        else:                           is_debug_mode = False

        # no parallelisation --------------------------------------------------
        if is_debug_mode:
            fit_list = []
            mech = conditions_list[0].mech

            gas = cdef.get_gas_ct(mech)
            for _c in range(len(conditions_list)):
                conditions_list[_c].composition.gas     = gas
                conditions_list[_c].composition.gas_ref = gas
            for ch in range(child_nb):
                fit_list.append([self.population[ch].fitness_eval(conditions_list,optim_param,ref_results_list,ch),ch])
        else:

            # Parallelisation 1 ---------------------------------------------------
    #        with Pool(processes=num_cores) as pool:
    #            fit_list = [ pool.apply_async(self.fitness_eval_par, args) for args in fit_eval_inp ]

            # Parallelisation 2 ---------------------------------------------------
    #        fit_list = []
    #        def log_result(fit_i):
    #            # This is called whenever foo_pool(i) returns a result.
    #            # result_list is modified only by the main process, not the pool workers.
    #            fit_list.append(fit_i)
    #
    #        pool=Pool(processes=num_cores)
    #        for args in fit_eval_inp:
    #            pool.apply_async(self.fitness_eval_par, args, callback=log_result)
    #        pool.close()

            # Parallelisation 3 ---------------------------------------------------
    #        if os.name == 'nt': multiprocessing.get_context('spawn')
    #        with multiprocessing.Pool(num_cores) as p:
    #            fit_list = p.map(self.fitness_eval_par, fit_eval_inp)


            # Parallelisation 4   (with simulation time check) --------------------
            #  https://pythonhosted.org/Pebble/#pools
            fit_list = []
            with ProcessPool() as pool:
                sim_results = pool.map(self.fitness_eval_par, fit_eval_inp, timeout=simul_time_limit)
                try:
                    for fit in sim_results.result():
                        fit_list.append(fit)
                except TimeoutError:
                    print_('\n\nWarning : simulation time > ' + '%.0f' %simul_time_limit + 's (> 1.5 x ref simulation time)',mp)
                    print_("TimeoutError: aborting remaining computations",mp)
                    fit_list.sort(reverse=False, key=lambda col: col[1])
                    fitness_incomplete = copy.deepcopy(fit_list)
                    list_ind_eval = []
                    for ind_fit_inc in fitness_incomplete:
                        list_ind_eval.append(ind_fit_inc[1])
                    n_sim=len(fitness_incomplete)
                    for _i in range(child_nb):
                        try:
                            if _i+optim_param.n_ind not in list_ind_eval:
                                fit_list.insert(_i,(0,_i))
                        except:
                            fit_list.append((0,_i))
                    print_('Number of individuals evaluated: ' + str(n_sim) +'\n\n',mp)
                    sim_results.cancel()
                except:
                    print_('\n\nWarning : error in simulation',mp)
                    print_("aborting remaining computations",mp)
                    fit_list.sort(reverse=False, key=lambda col: col[1])
                    fitness_incomplete = copy.deepcopy(fit_list)
                    list_ind_eval = []
                    for ind_fit_inc in fitness_incomplete:
                        list_ind_eval.append(ind_fit_inc[1])
                    n_sim=len(fitness_incomplete)
                    for _i in range(child_nb):
                        try:
                            if _i+optim_param.n_ind not in list_ind_eval:
                                fit_list.insert(_i,(0,_i))
                        except:
                            fit_list.append((0,_i))
                    print_('Number of individuals evaluated: ' + str(n_sim) +'\n\n',mp)
                    sim_results.cancel()



#        for _i in range(len(fit_list)):
#            self.population[_i].fitness = fit_list[_i][0]

#        if os.name == 'nt': multiprocessing.get_context('spawn')
#        with multiprocessing.Pool(num_cores) as p:
#            fit_list = p.map(self.fitness_eval_par, fit_eval_inp, timeout=5)
#            p.start()
#            p.join(1)
#
#            # If thread is active
#            if p.is_alive():
#                print_("foo is running... let's kill it...")
#                # Terminate foo
#                p.terminate()
#                # Cleanup
#                p.join()


#        if __name__ == '__main__':
#            # Start foo as a process
#        #    p = multiprocessing.Process(target=simul_flame, name="Foo", args=())
#            p = multiprocessing.Process(target=simul_flame, args=())
#
#            p.start()
#
#            # Wait a maximum of 10 seconds for foo
#            # Usage: join([timeout in seconds])
#            p.join(10)
#
#            # If thread is active
#            if p.is_alive():
#                print("foo is running... let's kill it...")
#
#                # Terminate foo
#                p.terminate()
#                # Cleanup
#                p.join()



        fit_list.sort(reverse=False, key=lambda col: col[1])


        for _i in range(len(fit_list)):
            ind = optim_param.n_ind + _i
            self.population[ind].fitness = fit_list[_i][0]

        bar.update(child_nb,title)
        print('\n')




        for i in range(len(conditions_list)):
            conditions_list[i].composition.gas     = gas
            conditions_list[i].composition.gas_ref = gas_ref
        for i in range(len(ref_results_list)):
            ref_results_list[i].gas = gas
            ref_results_list[i].f   = f[i]

    def fitness_eval_newpop(self,optim_param,conditions_list,ref_results_list):

        mp          = optim_param.main_path

        n_ind = optim_param.n_ind
#        print('total_Xover: '+str(optim_param.total_Xover))
#        print('total_mut: '+str(optim_param.total_mut))
#        print('total_child: '+str(child_nb))

        gas     = conditions_list[0].composition.gas
        gas_ref = conditions_list[0].composition.gas_ref

        # saving and suppression of unpickable variables on fitness eval inputs
        for cond in range(len(conditions_list)):
            del conditions_list[cond].composition.gas
            del conditions_list[cond].composition.gas_ref
        f=[]
        simul_time_limit = 0
        for res in range(len(ref_results_list)):
            f.append(ref_results_list[res].f)
            del ref_results_list[res].gas
            del ref_results_list[res].f
            simul_time_limit += ref_results_list[res].simul_time
        simul_time_limit = simul_time_limit*3

        bar = cdef.ProgressBar(n_ind, '')
        title="New pop evaluation  "
        bar.update(0,title)

        fit_eval_inp = []
        for _i in range(n_ind):
            fit_eval_inp.append([optim_param,conditions_list,ref_results_list,bar,_i,False])

        # saving of the ref mechanism
        mech = conditions_list[0].mech
        os.chdir(conditions_list[0].main_path)
        copyfile(mech,'GA/'+mech)
        os.chdir("GA")


        if sys.gettrace() is not None : is_debug_mode = True
        else:                           is_debug_mode = False

        if is_debug_mode:
            # bypass parallelisation for debugging
            fit_i = []
            mech = conditions_list[0].mech

            gas = cdef.get_gas_ct(mech)
            for _c in range(len(conditions_list)):
                conditions_list[_c].composition.gas     = gas
                conditions_list[_c].composition.gas_ref = gas
            for ind in range(n_ind):
                fit_i.append([self.population[ind].fitness_eval(conditions_list,optim_param,ref_results_list,ind),ind])

        else:
            # Parallelized Fitness calculation
            num_cores   = multiprocessing.cpu_count()
            if os.name == 'nt': multiprocessing.get_context('spawn')
            with multiprocessing.Pool(num_cores) as p:
                fit_i = p.map(self.fitness_eval_par, fit_eval_inp)
    #
            fit_i.sort(reverse=False, key=lambda col: col[1])

            for _i in range(len(fit_i)):
                self.population[_i].fitness = fit_i[_i][0]



        bar.update(n_ind,title)
        print('\n')

        for i in range(len(conditions_list)):
            conditions_list[i].composition.gas     = gas
            conditions_list[i].composition.gas_ref = gas_ref
        for i in range(len(ref_results_list)):
            ref_results_list[i].gas = gas
            ref_results_list[i].f   = f[i]



    def convergence_information(self,gen,optim_param,verbose=0):
        mp = optim_param.main_path

        vec = []
        for p in range(optim_param.n_ind):
            vec.append(self.population[p].fitness)
        min_fit = np.min(vec) ; max_fit = np.max(vec) ; avr_fit = np.mean(vec)

        if verbose>=2:
            print_("worst ind: "+"%.3f" %min_fit+     \
                  "   best ind: "+"%.3f" %max_fit+   \
                  "   mean fitness: "+"%.3f" %avr_fit,mp)

        if gen==0:
            optim_param.genVec.append(0)
            optim_param.best_fitness.append(max_fit)
            optim_param.mean_fitness.append(avr_fit)
            optim_param.worst_fitness.append(min_fit)
        else:
            optim_param.genVec.append(gen)
            optim_param.best_fitness.append(max_fit)
            optim_param.mean_fitness.append(avr_fit)
            optim_param.worst_fitness.append(min_fit)






#%% Selection

    def selection(self,optim_param,verbose):
        mp = optim_param.main_path

        if   optim_param.selection_operator == 'Elitism':     # elitism
            self.select_1_elit()
            if verbose >= 6: print_("Elitism selection",mp)
        elif optim_param.selection_operator == 'Roulette':    # roulette
            self.select_2_roulette(optim_param)
            if verbose >= 6: print_("Roulette selection",mp)
        elif optim_param.selection_operator == 'Rank':    # rank
            self.select_3_rank(optim_param)
            if verbose >= 6: print_("Rank selection",mp)
        elif optim_param.selection_operator == 'Geometric_norm':    # geometric
            self.select_4_geomNorm(optim_param)
            if verbose >= 6: print_("Geometric norm selection",mp)


    def select_1_elit(self):
        self.sort_fitness_rev()

    def select_2_roulette(self,optim_param):

        pop_copy = copy.deepcopy(self)
        pop_copy.sort_fitness()

        random.seed() ; fit = [] ; proba = [] ; rand_sel_vect=[] #; nan_list = []
        size_pop = optim_param.n_ind

        # build fitness vector
        for p in range(len(pop_copy.population)):
            # check if fitness is not nan
#            if pop_copy.population[p].fitness==pop_copy.population[p].fitness:
            fit.append(pop_copy.population[p].fitness)
#            else:
#                nan_list.append(p)
        # calculate selection probability vector
        for p in range(len(fit)):
            proba.append(fit[p]/np.sum(fit))
        proba = np.cumsum(proba)

        # build the random vector for the selection
        for p in range(size_pop):
            rand_sel_vect.append(random.random())
        rand_sel_vect.sort()


        new_ind = 0 ; i_prob =0 ; proba[-1]=1

        while new_ind<size_pop:
            if rand_sel_vect[new_ind]<=proba[i_prob]:
                self.population[new_ind]=copy.deepcopy(pop_copy.population[i_prob])
                new_ind += 1 #; i_prob = 0
            else:
                i_prob += 1
        #self.sort_fitness_rev()

    def select_3_rank(self, optim_param):

        pop_copy = copy.deepcopy(self)
        pop_copy.sort_fitness()

        random.seed() ; rank = [] ; proba = [] ; rand_sel_vect=[]

        # build rank vector
        for p in range(len(pop_copy.population)):
            rank.append(p+1)
        # calculate selection probability vector
        for p in range(len(pop_copy.population)):
            proba.append(rank[p]/np.sum(rank))
        proba = np.cumsum(proba)

        # build the random vector for the selection
        for p in range(len(pop_copy.population)):
            rand_sel_vect.append(random.random())
        rand_sel_vect.sort()

        new_ind = 0 ; i_prob =0 ;  proba[-1]=1
        size_pop = optim_param.n_ind
        while new_ind<size_pop:
            if rand_sel_vect[new_ind]<=proba[i_prob]:
                self.population[new_ind]=copy.deepcopy(pop_copy.population[i_prob])
                new_ind +=1 #; i_prob = 0
            else:
                i_prob+=1

    def select_4_geomNorm(self,optim_param):


        q = optim_param.selection_options[0]
        pop_copy = copy.deepcopy(self)
        pop_copy.sort_fitness

        random.seed() ; proba = [] ; rand_sel_vect=[]

        # calculate selection probability vector
        for r in range(len(pop_copy.population)):
            p = (q/(1-(1-q)**len(pop_copy.population)) )*(1-q)**r
            proba.append(p)
        proba = np.cumsum(proba)/np.sum(proba)
        proba.sort()

        # build the random vector for the selection
        for i in range(len(pop_copy.population)):
            rand_sel_vect.append(random.random())
        rand_sel_vect.sort()

        new_ind = 0 ; i_prob =0 ; proba[-1]=1
        size_pop = optim_param.n_ind
        while new_ind<size_pop:
            if rand_sel_vect[new_ind]<=proba[i_prob]:
                self.population[new_ind]=copy.deepcopy(pop_copy.population[i_prob])
                new_ind +=1 #; i_prob = 0
            else:
                i_prob += 1
        #self.sort_fitness_rev()







#%% CrossOver

    def Xover(self,optim_param,conditions_list,ref_results_list,verbose):
        mp = optim_param.main_path

        created_ind_total = 0

        Xover_num = optim_param.Xover_num


        for op in range(len(optim_param.Xover_operator)):
            created_ind = 0
            operator = optim_param.Xover_operator[op]
            if operator == 'simple_Xover' :     # simple Xover
                while created_ind<sum(Xover_num[0:1]):
                    self.Xover_1_simple(conditions_list,optim_param,\
                                        ref_results_list,created_ind)
                    created_ind_total += 2
                    created_ind += 2
            elif operator == 'multiple_Xover':    # multiple Xover
                while created_ind<sum(Xover_num[0:2]):
                    self.Xover_2_multiple(conditions_list,optim_param,\
                                          ref_results_list,created_ind)
                    created_ind_total += 2
                    created_ind += 2
            elif operator == 'arith_Xover':    # arithmetic Xover
                while created_ind<sum(Xover_num[0:3]):
                    self.Xover_3_arith(conditions_list,optim_param,\
                                       ref_results_list,created_ind)
                    created_ind_total += 2
                    created_ind += 2
            elif operator == 'heuristic_Xover':    # heuristic Xover
                while created_ind<sum(Xover_num[0:4]):
                    self.Xover_4_heuri(conditions_list,optim_param,\
                                       ref_results_list,created_ind)
                    created_ind_total += 2
                    created_ind += 2


    def Xover_1_simple(self,conditions_list,optim_param,ref_results_list,created_ind):

        random.seed()
        size_react = len(self.population[0].mech.react.kin)
        size_pop   = optim_param.n_ind
        parent1    = random.randrange(0, size_pop)
        parent2    = random.randrange(0, size_pop)
        child1     = size_pop+created_ind
        child2     = size_pop+created_ind+1
        while parent1 == parent2:
            parent2 = random.randrange(0, size_pop)

        self.population[child1] = copy.deepcopy(self.population[parent1])
        self.population[child2] = copy.deepcopy(self.population[parent2])

        a = random.randrange(0,size_react)
        for i in range(a,size_react):
            self.population[child1].mech.react.kin[i] = copy.deepcopy(self.population[parent2].mech.react.kin[i])
            self.population[child2].mech.react.kin[i] = copy.deepcopy(self.population[parent1].mech.react.kin[i])
                

    def Xover_2_multiple(self,conditions_list,optim_param,ref_results_list,created_ind):

        random.seed()
        size_react = len(self.population[0].mech.react.kin)
        size_pop   = optim_param.n_ind
        parent1    = random.randrange(0, size_pop)
        parent2    = random.randrange(0, size_pop)
        child1     = size_pop+created_ind
        child2     = size_pop+created_ind+1
        while parent1 == parent2:
            parent2 = random.randrange(0, size_pop)

        self.population[child1] = copy.deepcopy(self.population[parent1])
        self.population[child2] = copy.deepcopy(self.population[parent2])

        for i in range(size_react):
            a = random.randrange(0,2)
            if a ==1:
                self.population[child1].mech.react.kin[i] = copy.deepcopy(self.population[parent2].mech.react.kin[i])
                self.population[child2].mech.react.kin[i] = copy.deepcopy(self.population[parent1].mech.react.kin[i])


    def Xover_3_arith(self,conditions_list,optim_param,ref_results_list,created_ind):

        mp = optim_param.main_path
        random.seed()
        size_react = len(self.population[0].mech.react.kin)
        size_pop   = optim_param.n_ind
        p1    = random.randrange(0, size_pop)
        p2    = random.randrange(0, size_pop)
        child1     = size_pop+created_ind
        child2     = size_pop+created_ind+1
        while p1 == p2:
            p2 = random.randrange(0, size_pop)

        self.population[child1] = copy.deepcopy(self.population[p1])
        self.population[child2] = copy.deepcopy(self.population[p2])


        for r in range(size_react):
            try_r = 0 ; valid1,valid2 = False,False
            while not ((valid1 and valid2) or try_r > 3):
                incert_r = [u/100 for u in self.population[0].mech.react.incert[r]]
                if self.population[0].mech.react.modif[r]:
                    mix = random.random()
                    if self.population[0].mech.react.type[r] =='three_body_reaction' \
                    or self.population[0].mech.react.type[r] =='reaction':
                        nTry = 0; maxRetry = 15
                        while nTry <= maxRetry:
                            val1=[];val2=[]
                            for k in range(3):
                                ref = self.population[0].mech.react.ref_kin[r][k]
                                if ref!= 0.0:
                                    new_val1  = self.population[p1].mech.react.kin[r][k]*mix \
                                               +self.population[p2].mech.react.kin[r][k]*(1-mix)
                                    new_val2  = self.population[p2].mech.react.kin[r][k]*mix \
                                               +self.population[p1].mech.react.kin[r][k]*(1-mix)
                                    if optim_param.Arrh_var:
                                        a = ref*(1-incert_r[k])
                                        b = ref*(1-incert_r[k])
                                        min_val = min(a,b)
                                        max_val = max(a,b)
                                    else:
                                        f_min = self.population[0].mech.react.f_min[r]
                                        T_min = self.population[0].mech.react.f_Tit[0]
                                        T_max = self.population[0].mech.react.f_Tit[1]            
                                        R = 8.314/4.1868 #(cal/mol)                                    
                                        if k==0: # A
                                            if ref>0:
                                                min_val=ref/(10**f_min)
                                                max_val=ref*(10**f_min)
                                            if ref<0:
                                                min_val=ref*(10**f_min)
                                                max_val=ref/(10**f_min)
                                        if k==1: # n
                                            if ref>0:
                                                min_val=max(0,ref-(f_min/np.log10(T_max)))
                                                max_val=ref+(f_min/np.log10(T_max))
                                            if ref<0:
                                                min_val=ref-(f_min/np.log10(T_max))
                                                max_val=min(0,ref+(f_min/np.log10(T_max)))
                                        if k==2: # Ea
                                            if ref>0:
                                                min_val=max(0,ref-(f_min*R*T_min*np.log(10)))
                                                max_val=ref+(f_min*R*T_min*np.log(10))
                                            if ref<0:
                                                min_val=ref-(f_min*R*T_min*np.log(10))
                                                max_val=min(0,ref+(f_min*R*T_min*np.log(10)))                                
                                    if  min_val < new_val1 < max_val\
                                    and min_val < new_val2 < max_val\
                                    and new_val1/ref > 0\
                                    and new_val2/ref > 0:
                                        val1.append(new_val1)
                                        val2.append(new_val2)
                                else :
                                    val1.append(ref)
                                    val2.append(ref)
                            if len(val1)==3:
                                self.population[child1].mech.react.kin[r] = list(val1)
                                self.population[child2].mech.react.kin[r] = list(val2)
                                break
                            else:
                                nTry +=1 ; mix = random.random()

                    if self.population[0].mech.react.type[r] =='falloff_reaction' \
                    or self.population[0].mech.react.type[r] =='chemically_activated_reaction'\
                    or self.population[0].mech.react.type[r] =='chebyshev'\
                    or self.population[0].mech.react.type[r] =='pdep_arrhenius':
                        mix = random.random()                      
                        nTry = 0; maxRetry = 15
                        while nTry <= maxRetry:
                            success = True                        
                            for j in range(len(self.population[0].mech.react.kin[r])):
                                val1=[];val2=[]
                                for k in range(3):
                                    ref=self.population[0].mech.react.ref_kin[r][j][k]
                                    if ref!= 0.0:
                                        new_val1  = self.population[p1].mech.react.kin[r][j][k]*mix\
                                                   +self.population[p2].mech.react.kin[r][j][k]*(1-mix)
                                        new_val2  = self.population[p2].mech.react.kin[r][j][k]*mix \
                                                   +self.population[p1].mech.react.kin[r][j][k]*(1-mix)
                                        if optim_param.Arrh_var or self.population[0].mech.react.type[r] =='chebyshev':
                                            a = ref*(1-incert_r[k])
                                            b = ref*(1-incert_r[k])
                                            min_val = min(a,b)
                                            max_val = max(a,b)
                                        else:
                                            f_min = self.population[0].mech.react.f_min[r]
                                            T_min = self.population[0].mech.react.f_Tit[0]
                                            T_max = self.population[0].mech.react.f_Tit[1]            
                                            R = 8.314/4.1868 #(cal/mol)                                    
                                            if k==0: # A
                                                if ref>0:
                                                    min_val=ref/(10**f_min)
                                                    max_val=ref*(10**f_min)
                                                if ref<0:
                                                    min_val=ref*(10**f_min)
                                                    max_val=ref/(10**f_min)
                                            if k==1: # n
                                                if ref>0:
                                                    min_val=max(0,ref-(f_min/np.log10(T_max)))
                                                    max_val=ref+(f_min/np.log10(T_max))
                                                if ref<0:
                                                    min_val=ref-(f_min/np.log10(T_max))
                                                    max_val=min(0,ref+(f_min/np.log10(T_max)))
                                            if k==2: # Ea
                                                if ref>0:
                                                    min_val=max(0,ref-(f_min*R*T_min*np.log(10)))
                                                    max_val=ref+(f_min*R*T_min*np.log(10))
                                                if ref<0:
                                                    min_val=ref-(f_min*R*T_min*np.log(10))
                                                    max_val=min(0,ref+(f_min*R*T_min*np.log(10)))    
                                        if  min_val < new_val1 < max_val\
                                        and min_val < new_val2 < max_val\
                                        and new_val1/ref > 0\
                                        and new_val2/ref > 0:
                                            val1.append(new_val1)
                                            val2.append(new_val2)
                                        else:
                                            success = False
                                    else :
                                        val1.append(.0)
                                        val2.append(.0)
                            if success:
                                for j in range(len(self.population[0].mech.react.kin[r])):
                                    for k in range(3):
                                        if ref!= 0.0:
                                            new_val1  = self.population[p1].mech.react.kin[r][j][k]*mix\
                                                       +self.population[p2].mech.react.kin[r][j][k]*(1-mix)
                                            new_val2  = self.population[p2].mech.react.kin[r][j][k]*mix \
                                                       +self.population[p1].mech.react.kin[r][j][k]*(1-mix)
                                            self.population[child1].mech.react.kin[r][j][k] = new_val1
                                            self.population[child2].mech.react.kin[r][j][k] = new_val2                                
                                break
                            else:
                                nTry +=1 ; 
                    valid1 = check_k(self.population[child1].mech.react,r)
                    valid2 = check_k(self.population[child2].mech.react,r)
                    try_r += 1
                else:
                    valid1,valid2 = True,True

            if not valid1 or not valid2:
                if not valid1:
                    self.population[child1].mech.react.kin[r] = copy.deepcopy(self.population[p1].mech.react.kin[r])
                if not valid2:
                    self.population[child2].mech.react.kin[r] = copy.deepcopy(self.population[p2].mech.react.kin[r])



    def Xover_4_heuri(self,conditions_list,optim_param,ref_results_list,created_ind):

        random.seed()
        size_react = len(self.population[0].mech.react.kin)
        size_pop   = optim_param.n_ind
        p1    = random.randrange(0, size_pop)
        p2    = random.randrange(0, size_pop)
        damping = 0.75
        while p1 == p2:
            p2 = random.randrange(0, size_pop)
        if self.population[p1].fitness>self.population[p2].fitness:
            best=p1; worse=p2
        else:
            best=p2; worse=p1
        child1     = size_pop+created_ind
        child2     = size_pop+created_ind+1

        self.population[child1] = copy.deepcopy(self.population[worse])
        self.population[child2] = copy.deepcopy(self.population[best])

        for r in range(size_react):
            try_r = 0 ; valid1,valid2 = False,False
            while not ((valid1 and valid2) or try_r > 3):
                incert_r = [u/100 for u in self.population[0].mech.react.incert[r]]
                mix = random.random()
                if self.population[0].mech.react.modif[r]:
                    if self.population[0].mech.react.type[r] =='three_body_reaction' \
                    or self.population[0].mech.react.type[r] =='reaction':
                        nTry = 0; maxRetry = 15
                        while nTry <= maxRetry:
                            success = True
                            val1=[]
                            for k in range(3):
                                ref = self.population[0].mech.react.ref_kin[r][k]
                                if ref!= 0.0:
                                    bestVal  = self.population[best].mech.react.kin[r][k]
                                    worseVal = self.population[worse].mech.react.kin[r][k]
                                    new_val  = mix*(damping**try_r)*(bestVal-worseVal)+bestVal
                                    if optim_param.Arrh_var:
                                        a = ref*(1-incert_r[k])
                                        b = ref*(1-incert_r[k])
                                        min_val = min(a,b)
                                        max_val = max(a,b)
                                    else:
                                        f_min = self.population[0].mech.react.f_min[r]
                                        T_min = self.population[0].mech.react.f_Tit[0]
                                        T_max = self.population[0].mech.react.f_Tit[1]            
                                        R = 8.314/4.1868 #(cal/mol)                                    
                                        if k==0: # A
                                            if ref>0:
                                                min_val=ref/(10**f_min)
                                                max_val=ref*(10**f_min)
                                            if ref<0:
                                                min_val=ref*(10**f_min)
                                                max_val=ref/(10**f_min)
                                        if k==1: # n
                                            if ref>0:
                                                min_val=max(0,ref-(f_min/np.log10(T_max)))
                                                max_val=ref+(f_min/np.log10(T_max))
                                            if ref<0:
                                                min_val=ref-(f_min/np.log10(T_max))
                                                max_val=min(0,ref+(f_min/np.log10(T_max)))
                                        if k==2: # Ea
                                            if ref>0:
                                                min_val=max(0,ref-(f_min*R*T_min*np.log(10)))
                                                max_val=ref+(f_min*R*T_min*np.log(10))
                                            if ref<0:
                                                min_val=ref-(f_min*R*T_min*np.log(10))
                                                max_val=min(0,ref+(f_min*R*T_min*np.log(10)))                                
                                    if  min_val < new_val < max_val\
                                    and new_val/ref > 0:
                                        val1.append(new_val)
                                    else:
                                        success = False
                                else :
                                    val1.append(ref)
                            if success:
                                self.population[child1].mech.react.kin[r] = list(val1)
                                break
                            else:
                                nTry +=1 ; mix = random.random()
                    if self.population[0].mech.react.type[r] =='falloff_reaction' \
                    or self.population[0].mech.react.type[r] =='chemically_activated_reaction'\
                    or self.population[0].mech.react.type[r] =='chebyshev'\
                    or self.population[0].mech.react.type[r] =='pdep_arrhenius':
                        nTry = 0; maxRetry = 15
                        while nTry <= maxRetry:
                            success = True                        
                            for j in range(len(self.population[0].mech.react.kin[r])):
                                val1=[]
                                for k in range(3):
                                    ref=self.population[0].mech.react.ref_kin[r][j][k]
                                    if ref!= 0.0:
                                        bestVal  = self.population[best].mech.react.kin[r][j][k]
                                        worseVal = self.population[worse].mech.react.kin[r][j][k]
                                        new_val  = mix*(damping**try_r)*(bestVal-worseVal)+bestVal
                                        if optim_param.Arrh_var:
                                            a = ref*(1-incert_r[k])
                                            b = ref*(1-incert_r[k])
                                            min_val = min(a,b)
                                            max_val = max(a,b)
                                        else:
                                            f_min = self.population[0].mech.react.f_min[r]
                                            T_min = self.population[0].mech.react.f_Tit[0]
                                            T_max = self.population[0].mech.react.f_Tit[1]            
                                            R = 8.314/4.1868 #(cal/mol)                                    
                                            if k==0: # A
                                                if ref>0:
                                                    min_val=ref/(10**f_min)
                                                    max_val=ref*(10**f_min)
                                                if ref<0:
                                                    min_val=ref*(10**f_min)
                                                    max_val=ref/(10**f_min)
                                            if k==1: # n
                                                if ref>0:
                                                    min_val=max(0,ref-(f_min/np.log10(T_max)))
                                                    max_val=ref+(f_min/np.log10(T_max))
                                                if ref<0:
                                                    min_val=ref-(f_min/np.log10(T_max))
                                                    max_val=min(0,ref+(f_min/np.log10(T_max)))
                                            if k==2: # Ea
                                                if ref>0:
                                                    min_val=max(0,ref-(f_min*R*T_min*np.log(10)))
                                                    max_val=ref+(f_min*R*T_min*np.log(10))
                                                if ref<0:
                                                    min_val=ref-(f_min*R*T_min*np.log(10))
                                                    max_val=min(0,ref+(f_min*R*T_min*np.log(10)))
                                        if  min_val < new_val < max_val\
                                        and new_val/ref > 0:
                                            val1.append(new_val)
                                        else:
                                            success = False
                                    else:
                                        val1.append(.0)
                            if success:
                                for j in range(len(self.population[0].mech.react.kin[r])):
                                    for k in range(3):
                                        bestVal  = self.population[best].mech.react.kin[r][j][k]
                                        worseVal = self.population[worse].mech.react.kin[r][j][k]                                        
                                        new_val  = mix*(damping**try_r)*(bestVal-worseVal)+bestVal
                                        self.population[child1].mech.react.kin[r][j][k]=new_val
                                break
                            else:
                                nTry +=1 ; mix = random.random()
                    valid1 = check_k(self.population[child1].mech.react,r)
                    valid2 = check_k(self.population[child2].mech.react,r)
                    try_r += 1
                else:
                    valid1,valid2 = True,True
            if not valid1 and valid2:
                if not valid1:
                    self.population[child1].mech.react.kin[r] = copy.deepcopy(self.population[p1].mech.react.kin[r])
                if not valid2:
                    self.population[child2].mech.react.kin[r] = copy.deepcopy(self.population[p2].mech.react.kin[r])


#%% Mutation

    def mutation(self,optim_param,conditions_list,ref_results_list,gen,verbose):
        mp = optim_param.main_path

        mut_num = optim_param.mut_num

        created_ind_Xover = optim_param.total_Xover
        created_ind_total = 0 + created_ind_Xover



        for op in range(len(optim_param.mut_operator)):
            operator = optim_param.mut_operator[op]

            if operator == 'uniform_mutation' :     # 1- unifMutation
                while (created_ind_total-created_ind_Xover) < \
                            sum(mut_num[0:1]):
                    self.mut_1_unif(conditions_list,optim_param,ref_results_list,\
                               created_ind_total,optim_param.mut_option[0])
                    created_ind_total += 1

            elif operator == 'non_uniform_mutation':    # 2- nonUnifMutation
                while (created_ind_total-created_ind_Xover) < \
                            sum(mut_num[0:2]):
                    self.mut_2_nonUnif(conditions_list,optim_param,ref_results_list,\
                               created_ind_total,gen,optim_param.mut_option[1])
                    created_ind_total += 1

            elif operator == 'boundary_mutation':    # 3- boundaryMutation
                while (created_ind_total-created_ind_Xover) < \
                            sum(mut_num[0:3]):
                    self.mut_3_bound(conditions_list,optim_param,ref_results_list,\
                               created_ind_total,optim_param.mut_option[2])
                    created_ind_total += 1

    def mut_1_unif(self,conditions_list,optim_param,ref_results_list,created_ind,opt):

        random.seed()
        size_react = len(self.population[0].mech.react.kin)
        size_pop   = int(optim_param.n_ind)
        parent1    = random.randrange(0, size_pop)
        child1     = int(size_pop+created_ind)
        probaMut   = optim_param.mut_intensity
        damping    = 0.75

        self.population[child1] = copy.deepcopy(self.population[parent1])

        for r in range(size_react):
            try_r = 0 ; valid = False
            while not (valid or try_r > 10):
                incert_r = [u/100 for u in self.population[0].mech.react.incert[r]]
                if self.population[0].mech.react.modif[r]:
                    if self.population[0].mech.react.type[r] =='three_body_reaction' \
                    or self.population[0].mech.react.type[r] =='reaction':
                        val=[]
                        for k in range(3):                            
                            ref = self.population[0].mech.react.ref_kin[r][k]
                            rand = random.random()*100
                            if ref==0:
                                rand = probaMut+1
                            if rand < probaMut: # random modif of a kinetic constant
                                if optim_param.Arrh_var:
                                    val.append(random.uniform(ref*(1-incert_r[k]*(damping**try_r))\
                                                            ,ref*(1+incert_r[k]*(damping**try_r))))
                                
                                else:
                                    f_min = self.population[0].mech.react.f_min[r]*(damping**try_r)
                                    T_min = self.population[0].mech.react.f_Tit[0]
                                    T_max = self.population[0].mech.react.f_Tit[1]            
                                    R = 8.314/4.1868 #(cal/mol)                                    
                                    if k==0: # A
                                        val.append(random.uniform(ref/(10**f_min),ref*(10**f_min)))
                                    if k==1: # n
                                        if ref>0:
                                            min_val=max(0,ref-(f_min/np.log10(T_max)))
                                            max_val=ref+(f_min/np.log10(T_max))
                                        if ref<0:
                                            min_val=ref-(f_min/np.log10(T_max))
                                            max_val=min(0,ref+(f_min/np.log10(T_max)))
                                        val.append(random.uniform(min_val,max_val))
                                    if k==2: # Ea
                                        if ref>0:
                                            min_val=max(0,ref-(f_min*R*T_min*np.log(10)))
                                            max_val=ref+(f_min*R*T_min*np.log(10))
                                        if ref<0:
                                            min_val=ref-(f_min*R*T_min*np.log(10))
                                            max_val=min(0,ref+(f_min*R*T_min*np.log(10)))
                                        val.append(random.uniform(min_val,max_val))
                            else:
                                val.append(self.population[child1].mech.react.kin[r][k])
                        self.population[child1].mech.react.kin[r]=list(val)

                    if self.population[0].mech.react.type[r] =='falloff_reaction' \
                    or self.population[0].mech.react.type[r] =='chemically_activated_reaction'\
                    or self.population[0].mech.react.type[r] =='pdep_arrhenius':
                        mod_factor = [] ; j=0
                        # calculate the variation for the three Arrhenius parameters
                        # and then apply the same factors to the other Arrhenius parameters
                        # to keep consistancy on the modifications (see Bertolino et al. Comb and Flame 229 (2021) 111366)
                        for k in range(3):
                            ref = self.population[0].mech.react.ref_kin[r][j][k]
                            rand = random.random()*100
                            if ref==0:
                                rand = probaMut+1
                            if rand < probaMut: # random modif of a kinetic constant
                                if optim_param.Arrh_var:
                                    val = random.uniform(ref*(1-incert_r[k]*(damping**try_r))\
                                                            ,ref*(1+incert_r[k]*(damping**try_r)))
                                    if ref!=0:  mod_factor.append(val/ref)
                                    else:       mod_factor.append(0)
                                else:
                                    f_min = self.population[0].mech.react.f_min[r]*(damping**try_r)
                                    T_min = self.population[0].mech.react.f_Tit[0]
                                    T_max = self.population[0].mech.react.f_Tit[1]            
                                    R = 8.314/4.1868 #(cal/mol)                                    
                                    if k==0: # A
                                        val = random.uniform(ref/(10**f_min),ref*(10**f_min))
                                    if k==1: # n
                                        if ref>0:
                                            min_val=max(0,ref-(f_min/np.log10(T_max)))
                                            max_val=ref+(f_min/np.log10(T_max))
                                        if ref<0:
                                            min_val=ref-(f_min/np.log10(T_max))
                                            max_val=min(0,ref+(f_min/np.log10(T_max)))
                                        val = random.uniform(min_val,max_val)
                                    if k==2: # Ea
                                        if ref>0:
                                            min_val=max(0,ref-(f_min*R*T_min*np.log(10)))
                                            max_val=ref+(f_min*R*T_min*np.log(10))
                                        if ref<0:
                                            min_val=ref-(f_min*R*T_min*np.log(10))
                                            max_val=min(0,ref+(f_min*R*T_min*np.log(10)))
                                        val = random.uniform(min_val,max_val)
                                    if ref!=0:  mod_factor.append(val/ref)
                                    else:       mod_factor.append(0)
                            else:
                                mod_factor.append(1)
                        # Introduction of the new Arrhenius parameters                                
                        for j in range(len(self.population[0].mech.react.kin[r])):
                            for k in range(3):
                                ref = self.population[0].mech.react.ref_kin[r][j][k]
                                self.population[child1].mech.react.kin[r][j][k]=ref*mod_factor[k]
                                
                    if self.population[0].mech.react.type[r] =='chebyshev':
                        for j in range(len(self.population[0].mech.react.kin[r])):
                            val=[]
                            for k in range(3):
                                ref = self.population[0].mech.react.ref_kin[r][j][k]
                                rand = random.random()*100
                                if ref==0:
                                    rand = probaMut+1
                                if rand < probaMut: # random modif of a kinetic constant
                                    val.append(random.uniform(ref*(1-incert_r[k]*(damping**try_r))\
                                                            ,ref*(1+incert_r[k]*(damping**try_r))))
                                else:
                                    val.append(self.population[child1].mech.react.kin[r][j][k])
                            self.population[child1].mech.react.kin[r][j]=list(val)
                                
                # Validation of the new reaction rate                                
                valid = check_k(self.population[child1].mech.react,r)
                try_r += 1

            if not valid:
                self.population[child1].mech.react.kin[r] = copy.deepcopy(self.population[parent1].mech.react.kin[r])



    def mut_2_nonUnif(self,conditions_list,optim_param,ref_results_list,created_ind,\
                      gen,option):

        random.seed()
        size_react = len(self.population[0].mech.react.kin)
        size_pop   = int(optim_param.n_ind)
        parent1    = random.randrange(0, size_pop)
        child1     = int(size_pop+created_ind)
        probaMut   = optim_param.mut_intensity
        ratio      = gen/optim_param.n_gen
        shape      = option
        damping    = 0.75

        self.population[child1] = copy.deepcopy(self.population[parent1])

        for r in range(size_react):
            try_r = 0 ; valid = False
            while not (valid or try_r > 10):
                rand = random.random()*100
                incert_r = [u/100 for u in self.population[0].mech.react.incert[r]]
                if rand < probaMut:
                    if self.population[0].mech.react.modif[r]:
                        if self.population[0].mech.react.type[r] =='three_body_reaction' \
                        or self.population[0].mech.react.type[r] =='reaction':
                            for k in range(3):
                                ref = self.population[0].mech.react.ref_kin[r][k]
                                val = self.population[child1].mech.react.kin[r][k]
                                if optim_param.Arrh_var:
                                    if ref>0:
                                        min_val = ref*(1-incert_r[k]*(damping**try_r))
                                        max_val = ref*(1+incert_r[k]*(damping**try_r))
                                    else:
                                        max_val = ref*(1-incert_r[k]*(damping**try_r))
                                        min_val = ref*(1+incert_r[k]*(damping**try_r))
                                else:
                                    f_min = self.population[0].mech.react.f_min[r]*(damping**try_r)
                                    T_min = self.population[0].mech.react.f_Tit[0]
                                    T_max = self.population[0].mech.react.f_Tit[1]            
                                    R = 8.314/4.1868 #(cal/mol)   
                                    if ref==0:
                                        min_val, max_val = 0,0
                                    else:
                                        if k==0: # A
                                            if ref>0:
                                                min_val=ref/(10**f_min)
                                                max_val=ref*(10**f_min)
                                            if ref<0:
                                                min_val=ref*(10**f_min)
                                                max_val=ref/(10**f_min)
                                        if k==1: # n
                                            if ref>0:
                                                min_val=max(0,ref-(f_min/np.log10(T_max)))
                                                max_val=ref+(f_min/np.log10(T_max))
                                            if ref<0:
                                                min_val=ref-(f_min/np.log10(T_max))
                                                max_val=min(0,ref+(f_min/np.log10(T_max)))
                                        if k==2: # Ea
                                            if ref>0:
                                                min_val=max(0,ref-(f_min*R*T_min*np.log(10)))
                                                max_val=ref+(f_min*R*T_min*np.log(10))
                                            if ref<0:
                                                min_val=ref-(f_min*R*T_min*np.log(10))
                                                max_val=min(0,ref+(f_min*R*T_min*np.log(10)))
                                rand_dir = random.randrange(0,2)
                                if rand_dir == 0: # random modif of a kinetic constant
                                    change=(max_val-val)*(random.random()*(1-ratio))**shape
                                    self.population[child1].mech.react.kin[r][k]+=change
                                elif rand_dir == 1:
                                    change=(val-min_val)*(random.random()*(1-ratio))**shape
                                    self.population[child1].mech.react.kin[r][k]-=change
                                    
                        elif self.population[0].mech.react.type[r] =='falloff_reaction' \
                        or self.population[0].mech.react.type[r] =='chemically_activated_reaction'\
                        or self.population[0].mech.react.type[r] =='pdep_arrhenius':
                            mod_factor = [] ; j=0     
                            # calculate the variation for the three Arrhenius parameters
                            # and then apply the same factors to the other Arrhenius parameters
                            # to keep consistancy on the modifications (see Bertolino et al. Comb and Flame 229 (2021) 111366)
                            for k in range(3):
                                ref = self.population[0].mech.react.ref_kin[r][j][k]
                                val = self.population[child1].mech.react.kin[r][j][k]                                    
                                if optim_param.Arrh_var:
                                    if ref>0:
                                        min_val = ref*(1-incert_r[k]*(damping**try_r))
                                        max_val = ref*(1+incert_r[k]*(damping**try_r))
                                    else:
                                        max_val = ref*(1-incert_r[k]*(damping**try_r))
                                        min_val = ref*(1+incert_r[k]*(damping**try_r))
                                else:
                                    f_min = self.population[0].mech.react.f_min[r]*(damping**try_r)
                                    T_min = self.population[0].mech.react.f_Tit[0]
                                    T_max = self.population[0].mech.react.f_Tit[1]            
                                    R = 8.314/4.1868 #(cal/mol)
                                    if ref==0:
                                        min_val, max_val = 0,0
                                    else:                                        
                                        if k==0: # A
                                            if ref>0:
                                                min_val=ref/(10**f_min)
                                                max_val=ref*(10**f_min)
                                            if ref<0:
                                                min_val=ref*(10**f_min)
                                                max_val=ref/(10**f_min)
                                        if k==1: # n
                                            if ref>0:
                                                min_val=max(0,ref-(f_min/np.log10(T_max)))
                                                max_val=ref+(f_min/np.log10(T_max))
                                            if ref<0:
                                                min_val=ref-(f_min/np.log10(T_max))
                                                max_val=min(0,ref+(f_min/np.log10(T_max)))
                                        if k==2: # Ea
                                            if ref>0:
                                                min_val=max(0,ref-(f_min*R*T_min*np.log(10)))
                                                max_val=ref+(f_min*R*T_min*np.log(10))
                                            if ref<0:
                                                min_val=ref-(f_min*R*T_min*np.log(10))
                                                max_val=min(0,ref+(f_min*R*T_min*np.log(10)))
                                rand_dir = random.randrange(0,2)
                                if ref != 0:
                                    if rand_dir == 0: # random modif of a kinetic constant
                                        change=(max_val-val)*(random.random()*(1-ratio))**shape
                                        new_val = self.population[child1].mech.react.kin[r][j][k]+change
                                        mod_factor.append(new_val/ref)
                                    elif rand_dir == 1:
                                        change=(val-min_val)*(random.random()*(1-ratio))**shape
                                        new_val = self.population[child1].mech.react.kin[r][j][k]-change
                                        mod_factor.append(new_val/ref)
                                else:
                                    mod_factor.append(0)
                            # Introduction of the new Arrhenius parameters
                            for j in range(len(self.population[0].mech.react.kin[r])):
                                for k in range(3):
                                    ref = self.population[0].mech.react.ref_kin[r][j][k]
                                    self.population[child1].mech.react.kin[r][j][k]=ref*mod_factor[k]
                                    
                        elif self.population[0].mech.react.type[r] =='chebyshev':
                            for j in range(len(self.population[0].mech.react.kin[r])):
                                for k in range(3):
                                    ref = self.population[0].mech.react.ref_kin[r][j][k]
                                    val = self.population[child1].mech.react.kin[r][j][k]                                    
                                    if ref>0:
                                        min_val = ref*(1-incert_r[k]*(damping**try_r))
                                        max_val = ref*(1+incert_r[k]*(damping**try_r))
                                    else:
                                        max_val = ref*(1-incert_r[k]*(damping**try_r))
                                        min_val = ref*(1+incert_r[k]*(damping**try_r))
                                    rand_dir = random.randrange(0,2)
                                    if rand_dir == 0: # random modif of a kinetic constant
                                        change=(max_val-val)*(random.random()*(1-ratio))**shape
                                        self.population[child1].mech.react.kin[r][j][k]+=change
                                    elif rand_dir == 1:
                                        change=(val-min_val)*(random.random()*(1-ratio))**shape
                                        self.population[child1].mech.react.kin[r][j][k]-=change                                    
                                    
                # Validation of the new reaction rate
                valid = check_k(self.population[child1].mech.react,r)
                try_r += 1

            if not valid:
                self.population[child1].mech.react.kin[r] = copy.deepcopy(self.population[parent1].mech.react.kin[r])


#        self.population[child1].fitness_eval(conditions_list,optim_param,ref_results_list)



    def mut_3_bound(self,conditions_list,optim_param,ref_results_list,created_ind,opt):

        random.seed()
        size_react = len(self.population[0].mech.react.kin)
        size_pop   = int(optim_param.n_ind)
        parent1    = random.randrange(0, size_pop)
        child1     = int(size_pop+created_ind)
        probaMut   = optim_param.mut_intensity
        damping    = 0.75
        self.population[child1] = copy.deepcopy(self.population[parent1])

        for r in range(size_react):
            try_r = 0 ; valid = False 
            rand = random.random()*100
            if rand < probaMut:      
                incert_r = [u/100 for u in self.population[0].mech.react.incert[r]]                
                while not (valid or try_r > 10):
                    if self.population[0].mech.react.modif[r]:
                        if self.population[0].mech.react.type[r] =='three_body_reaction' \
                        or self.population[0].mech.react.type[r] =='reaction':
                            for k in range(3):
                                ref = self.population[0].mech.react.ref_kin[r][k]
                                rand_dir = random.randrange(0,2)           
                                if optim_param.Arrh_var:
                                    min_val = ref*(1-incert_r[k]*(damping**try_r))
                                    max_val = ref*(1-incert_r[k]*(damping**try_r))
                                else:
                                    f_min = self.population[0].mech.react.f_min[r]*(damping**try_r)
                                    T_min = self.population[0].mech.react.f_Tit[0]
                                    T_max = self.population[0].mech.react.f_Tit[1]            
                                    R = 8.314/4.1868 #(cal/mol)
                                    if ref==0:
                                        min_val, max_val = 0,0
                                    else:
                                        if k==0: # A
                                            if ref>0:
                                                min_val=ref/(10**f_min)
                                                max_val=ref*(10**f_min)
                                            if ref<0:
                                                min_val=ref*(10**f_min)
                                                max_val=ref/(10**f_min)
                                        if k==1: # n
                                            if ref>0:
                                                min_val=max(0,ref-(f_min/np.log10(T_max)))
                                                max_val=ref+(f_min/np.log10(T_max))
                                            if ref<0:
                                                min_val=ref-(f_min/np.log10(T_max))
                                                max_val=min(0,ref+(f_min/np.log10(T_max)))
                                        if k==2: # Ea
                                            if ref>0:
                                                min_val=max(0,ref-(f_min*R*T_min*np.log(10)))
                                                max_val=ref+(f_min*R*T_min*np.log(10))
                                            if ref<0:
                                                min_val=ref-(f_min*R*T_min*np.log(10))
                                                max_val=min(0,ref+(f_min*R*T_min*np.log(10)))                                
                                if rand_dir == 0: # random modif of a kinetic constant to the boundary value        
                                    self.population[child1].mech.react.kin[r][k] = min_val
                                elif rand_dir == 1:
                                    self.population[child1].mech.react.kin[r][k] = max_val
                        elif self.population[0].mech.react.type[r] =='falloff_reaction' \
                        or   self.population[0].mech.react.type[r] =='chemically_activated_reaction'\
                        or   self.population[0].mech.react.type[r] =='pdep_arrhenius':
                            mod_factor = [] ; j=0     
                            # calculate the variation for the three first Arrhenius parameters
                            # and then apply the same factors to the other Arrhenius parameters
                            # to keep consistancy on the modifications (see Bertolino et al. Comb and Flame 229 (2021) 111366)
                            for k in range(3):
                                ref = self.population[0].mech.react.ref_kin[r][j][k]
                                rand_dir = random.randrange(0,2)
                                if optim_param.Arrh_var:
                                    min_val = ref*(1-incert_r[k]*(damping**try_r))
                                    max_val = ref*(1-incert_r[k]*(damping**try_r))
                                else:
                                    f_min = self.population[0].mech.react.f_min[r]*(damping**try_r)
                                    T_min = self.population[0].mech.react.f_Tit[0]
                                    T_max = self.population[0].mech.react.f_Tit[1]            
                                    R = 8.314/4.1868 #(cal/mol)    
                                    if ref==0:
                                        min_val, max_val = 0,0
                                    else:
                                        if k==0: # A
                                            if ref>0:
                                                min_val=ref/(10**f_min)
                                                max_val=ref*(10**f_min)
                                            if ref<0:
                                                min_val=ref*(10**f_min)
                                                max_val=ref/(10**f_min)
                                        if k==1: # n
                                            if ref>0:
                                                min_val=max(0,ref-(f_min/np.log10(T_max)))
                                                max_val=ref+(f_min/np.log10(T_max))
                                            if ref<0:
                                                min_val=ref-(f_min/np.log10(T_max))
                                                max_val=min(0,ref+(f_min/np.log10(T_max)))
                                        if k==2: # Ea
                                            if ref>0:
                                                min_val=max(0,ref-(f_min*R*T_min*np.log(10)))
                                                max_val=ref+(f_min*R*T_min*np.log(10))
                                            if ref<0:
                                                min_val=ref-(f_min*R*T_min*np.log(10))
                                                max_val=min(0,ref+(f_min*R*T_min*np.log(10)))
                                if ref != 0:
                                    if rand_dir == 0: # random modif of a kinetic constant to the boundary value
                                        mod_factor.append(min_val/ref)                                                                    
                                    elif rand_dir == 1:                            
                                        mod_factor.append(max_val/ref)            
                                else:
                                    mod_factor.append(0)
                            # Introduction of the new Arrhenius parameters
                            for j in range(len(self.population[0].mech.react.kin[r])):
                                for k in range(3):
                                    ref = self.population[0].mech.react.ref_kin[r][j][k]
                                    self.population[child1].mech.react.kin[r][j][k]=ref*mod_factor[k]
                                    
                        elif self.population[0].mech.react.type[r] =='chebyshev':
                            for j in range(len(self.population[0].mech.react.kin[r])):
                                for k in range(3):
                                    ref = self.population[0].mech.react.ref_kin[r][j][k]
                                    rand_dir = random.randrange(0,2)
                                    min_val = ref*(1-incert_r[k]*(damping**try_r))
                                    max_val = ref*(1-incert_r[k]*(damping**try_r))
                                if rand_dir == 0: # random modif of a kinetic constant to the boundary value        
                                    self.population[child1].mech.react.kin[r][j][k] = min_val
                                elif rand_dir == 1:
                                    self.population[child1].mech.react.kin[r][j][k] = max_val
                
                    # Validation of the new reaction rate
                    valid = check_k(self.population[child1].mech.react,r)
                    try_r += 1

            if not valid:
                self.population[child1].mech.react.kin[r] = copy.deepcopy(self.population[parent1].mech.react.kin[r])

#        self.population[child1].fitness_eval(conditions_list,optim_param,ref_results_list)





#%% Others

    def popBestInd(self):
        index=0
        max_fit  = 0
        for i in range(len(self.population)):
            fit = self.population[i].fitness
            if fit > max_fit:
                index = i
                max_fit = fit
        bestIndCoeff = self.population[index].value
        bestIndFit = self.population[index].fitness

        return bestIndCoeff, bestIndFit

    def sort_fitness(self):
        """ Sort population according to chromosome fitness """
        self.population =sorted(self.population, key=operator.attrgetter('fitness'), reverse=False)

    def sort_fitness_rev(self):
        """ Sort population according to chromosome fitness """
        self.population =sorted(self.population, key=operator.attrgetter('fitness'), reverse=True)


#%% Functions

def check_k(_react,_r):

    T_min  = _react.f_Tit[0]
    T_max  = _react.f_Tit[1]
    T_incr = _react.f_Tit[2]
    T_list = np.array(range(T_min,T_max+T_incr,T_incr))
    R = 8.314/4.1868 #(cal/mol)
    P =  1e5  # (Pa)
    lp = 5e4  # (Pa)
    hp = 5e6  # (Pa)

    valid_kin = True
    # https://cantera.org/science/reactions.html
    import warnings


    Y_r = []

    if _react.type[_r] == 'reaction' or _react.type[_r] == 'three_body_reaction':
        A = _react.kin[_r][0] ; n = _react.kin[_r][1] ; Ea = _react.kin[_r][2]
        for i,T in enumerate(T_list):
            Y_r.append(A*T**n*np.exp(-Ea/(R*T)))
            if not _react.k_min[_r][i] < A*T**n*np.exp(-Ea/(R*T)) < _react.k_max[_r][i]:
                valid_kin = False

    elif _react.type[_r] == 'pdep_arrhenius':
        all_P = [lp, P, hp]
        list_P = _react.param[_r][0]

        for Pi,_pi in enumerate(all_P):
            # pressure limits of the reaction law
            idx_Pinf,idx_Psup = 0,0
            for _p in range(len(list_P)):
                if P <= float(list_P[_p]):
                    idx_Psup = _p
                    if idx_Psup==0: idx_Pinf = 0
                    else:           idx_Pinf = _p-1
                    break
            P_inf = float(list_P[idx_Pinf])
            P_sup = float(list_P[idx_Psup])

            for i,T in enumerate(T_list):
                # k_pinf:
                i=0 ; k_pinf=0
                while float(list_P[idx_Pinf-i]) == float(list_P[idx_Pinf]):
                    A = _react.kin[_r][idx_Pinf-i][0]
                    n = _react.kin[_r][idx_Pinf-i][1]
                    Ea = _react.kin[_r][idx_Pinf-i][2]
                    k_pinf +=  A * (T)**n * np.exp(-Ea/(R*T))
                    i+=1
                # k_psup:
                i=0 ; k_psup=0
                while float(list_P[idx_Psup+i]) == float(list_P[idx_Psup]):
                    A = _react.kin[_r][idx_Psup+i][0]
                    n = _react.kin[_r][idx_Psup+i][1]
                    Ea = _react.kin[_r][idx_Psup+i][2]
                    k_psup +=  A * (T)**n * np.exp(-Ea/(R*T))
                    i+=1
                    if idx_Psup+i>=len(list_P):
                        break
                # k
                if P_inf == P_sup:
                    if _pi==0: k_lp = k_pinf
                    if _pi==1: k    = k_pinf
                    if _pi==2: k_hp = k_pinf
                else:
                    if k_pinf < 1e-150: k_pinf = 1e-150  # avoid log(0) calculation
                    if k_psup < 1e-150: k_psup = 1e-150  # avoid log(0) calculation
                    if _pi==0:
                        k_lp = np.exp(np.log(k_pinf)+(np.log(k_psup)-np.log(k_pinf))*((np.log(Pi)-np.log(P_inf))/(np.log(P_sup)-np.log(P_inf))))
                        if not _react.k_min_lp[_r][i] < k_lp < _react.k_max_lp[_r][i]:
                            valid_kin = False
                            break
                    if _pi==1:
                        k = np.exp(np.log(k_pinf)+(np.log(k_psup)-np.log(k_pinf))*((np.log(Pi)-np.log(P_inf))/(np.log(P_sup)-np.log(P_inf))))
                        if not _react.k_min[_r][i] < k < _react.k_max[_r][i]:
                            valid_kin = False
                            break
                    if _pi==2:
                        k_hp = np.exp(np.log(k_pinf)+(np.log(k_psup)-np.log(k_pinf))*((np.log(Pi)-np.log(P_inf))/(np.log(P_sup)-np.log(P_inf))))
                        if not _react.k_min_hp[_r][i] < k_hp < _react.k_max_hp[_r][i]:
                            valid_kin = False

    elif _react.type[_r] == 'falloff_reaction' or _react.type[_r] == 'chemically_activated_reaction':
        # get Troe param
        try:
            Troe_func = True
            At = np.float32(_react.param[_r]['A'])
            T3t = np.float32(_react.param[_r]['T3'])
            T1t = np.float32(_react.param[_r]['T1'])
            T2t = np.float32(_react.param[_r]['T2'])
        except:
            Troe_func = False
        # calculation of k
        Af = _react.kin[_r][0][0] ; nf = _react.kin[_r][0][1] ; Eaf = _react.kin[_r][0][2]
        A0 = _react.kin[_r][1][0] ; n0 = _react.kin[_r][1][1] ; Ea0 = _react.kin[_r][1][2]
        for i,T in enumerate(T_list):
            M    = P/(8.314*T)           # mol/m3
            M    = M*1e-6                # mol/cm3
            M_lp  = (lp*1e-6)/(8.314*T)  # mol/cm3
            M_hp  = (hp*1e-6)/(8.314*T)  # mol/cm3
            kinf = Af * T**nf * np.exp(-Eaf/(R*T))
            k0   = A0 * T**n0 * np.exp(-Ea0/(R*T))
            Pr    = (k0*M)/kinf
            Pr_lp = (k0*M_lp)/kinf
            Pr_hp = (k0*M_hp)/kinf
            warnings.filterwarnings("ignore") # suppress warnings with exponential calculation
            if Troe_func:
                if T2t: Fcent = (1-At)*np.exp(-T/T3t) + At*np.exp(-T/T1t) + np.exp(-T2t/T1t)
                else:   Fcent = (1-At)*np.exp(-T/T3t) + At*np.exp(-T/T1t)
                C  = -0.4 - 0.67*np.log10(Fcent)
                N  = 0.75 - 1.27*np.log10(Fcent)
                f1 = (np.log10(Pr) + C) / (N-0.14*(np.log10(Pr)+C))
                F  = 10**(np.log10(Fcent) / (1+f1**2))
                f1_lp = (np.log10(Pr_lp) + C) / (N-0.14*(np.log10(Pr_lp)+C))
                F_lp  = 10**(np.log10(Fcent) / (1+f1_lp**2))
                f1_hp = (np.log10(Pr_hp) + C) / (N-0.14*(np.log10(Pr_hp)+C))
                F_hp  = 10**(np.log10(Fcent) / (1+f1_hp**2))
            else:
                F,F_lp,F_hp = 1,1,1
            if _react.type[_r] == 'falloff_reaction':
                k    = kinf*(Pr/(1+Pr))*F
                k_lp = kinf*(Pr_lp/(1+Pr_lp))*F_lp
                k_hp = kinf*(Pr_hp/(1+Pr_hp))*F_hp
            elif _react.type[_r] == 'chemically_activated_reaction':
                k    = k0*(1/(1+Pr))*F
                k_lp = k0*(1/(1+Pr_lp))*F
                k_hp = k0*(1/(1+Pr_hp))*F
            warnings.resetwarnings()
            if not   _react.k_min[_r][i] < k    < _react.k_max[_r][i]:
                valid_kin = False
            elif not _react.k_min_lp[_r][i] < k_lp < _react.k_max_lp[_r][i]:
                valid_kin = False
            elif not _react.k_min_hp[_r][i] < k_hp < _react.k_max_hp[_r][i]:
                valid_kin = False
    elif 'chebyshev' in _react.type[_r]:
        # get param
        Tmin = np.float32(_react.param[_r]['Tmin'])
        Tmax = np.float32(_react.param[_r]['Tmax'])
        Pmin = np.float32(_react.param[_r]['Pmin'])
        Pmax = np.float32(_react.param[_r]['Pmin'])
        if _react.param[_r]['Pmin_unit'] is not False:
            if "'Pa'" in _react.param[_r]['Pmin_unit']:    # conversion to atm
                Pmin = Pmin/101325
            elif "'Bar'" in _react.param[_r]['Pmin_unit']: # conversion to atm
                Pmin = Pmin/1.01325
        if _react.param[_r]['Pmax_unit'] is not False:
            if "'Pa'" in _react.param[_r]['Pmax_unit']:    # conversion to atm
                Pmax = Pmax/101325
            elif "'Bar'" in _react.param[_r]['Pmax_unit']: # conversion to atm
                Pmax = Pmax/1.01325


        # k calculation (Carstensen H-H. Reaction Rate Representation Using Chebyshev Polynomials)
        # size of NT Np coefficients
        N_T = len(_react.kin[_r]) ; N_P = len(_react.kin[_r][0])
        for i,T in enumerate(T_list):
            T_tilde = (2/T - 1/Tmin - 1/Tmax)/(1/Tmax - 1/Tmin)
            P_tilde = (2*np.log10(P) - np.log10(Pmin) - np.log10(Pmax))/(np.log10(Pmax) - np.log10(Pmin))

            log_k = 0
            for t in range(N_T):
                phi_tT = np.cos(t*np.arccos(T_tilde))
                for p in range(N_P):
                    phi_pP = np.cos(p*np.arccos(P_tilde))
                    a_tp = _react.kin[_r][t][p]
                    log_k += a_tp*phi_tT*phi_pP
            k = 10**log_k
            if not _react.k_min[_r][i] < k < _react.k_max[_r][i]:
                valid_kin = False

    else:
        k = list(np.array([0]*len(T_list)))


#
#    if _react.equation[_r] == 'CH3O + O2 <=> CH2O + HO2':
##    or _react.equation[_r] == 'CH3O + H (+ M) <=> CH3OH (+ M)'\
##    or _react.equation[_r] == 'H2 + O <=> H + OH'\
##    or _react.equation[_r] == 'HO2 + OH <=> H2O + O2'\
##    or _react.equation[_r] == 'C2H4 + H (+M) <=> C2H5 (+M)':
#        import tool_brooks as tb
#        print(str(valid_kin))
#
#        if valid_kin:
#            tb.plot_graph___Title_X_Y_Legend_kwargs(_react.equation[_r],
#                                                1000/T_list,
#                                                [Y_r],
#                                                ['k0'],
#                                                fill_between_X = 1000/T_list,
#                                                fill_between_Y = [_react.k_min[_r],_react.k_max[_r]],
#                                                fig_save=True,
#                                                fig_show=True,
#                                                Y_log   =True)
#        print(str(valid_kin))

    return valid_kin




def compute_kref_and_klim(mech_data,p_min,p_max):

    T_min  = mech_data.react.f_Tit[0]
    T_max  = mech_data.react.f_Tit[1]
    T_incr = mech_data.react.f_Tit[2]
    T_list = np.array(range(T_min,T_max+T_incr,T_incr))

    R = 8.314/4.1868 #(cal/mol)

    P =  1e5  # (Pa)
    lp = 5e4  # (Pa)
    hp = 5e6  # (Pa)

    # https://cantera.org/science/reactions.html
    import warnings
    mech_data.react.k_ref = []
    mech_data.react.k_min = []
    mech_data.react.k_max = []
    mech_data.react.k_ref_lp = []
    mech_data.react.k_min_lp = []
    mech_data.react.k_max_lp = []
    mech_data.react.k_ref_hp = []
    mech_data.react.k_min_hp = []
    mech_data.react.k_max_hp = []
    mech_data.react.f_min    = [] 
    mech_data.react.f_lp_min = []
    mech_data.react.f_hp_min = []

    for _r in range(len(mech_data.react.type)):
        k_ref,k_ref_lp,k_ref_hp = [],[],[]
        
        # atm pressure
        f_T,k0 = [],[]
        cT  = mech_data.react.f_T_fit[_r]

        if mech_data.react.k0_fit[_r] is not False:
            ck0 = np.poly1d(mech_data.react.k0_fit[_r])
            k0  = 10**ck0(T_list)
        else:
            k0  = False

        if type(cT) is list:
            cTf   = np.poly1d(cT)
            f_T   = cTf(T_list)
        else:
            f_T   = np.array([cT]*len(T_list))
        mech_data.react.f_min.append(min(f_T))


        # low pressure
        f_T_lp,k0_lp = [],[]
        cT  = np.poly1d(mech_data.react.f_T_fit_lp[_r])
        ck0 = np.poly1d(mech_data.react.k0_fit_lp[_r])
        f_T_lp   = cT(T_list)
        k0_lp    = 10**ck0(T_list)
        mech_data.react.f_lp_min.append(min(f_T_lp))

        # high pressure
        f_T_hp,k0_hp = [],[]
        cT  = np.poly1d(mech_data.react.f_T_fit_hp[_r])
        ck0 = np.poly1d(mech_data.react.k0_fit_hp[_r])
        f_T_hp   = cT(T_list)
        k0_hp    = 10**ck0(T_list)
        mech_data.react.f_hp_min.append(min(f_T_hp))

        if k0 is not False:
            if mech_data.react.type[_r] == 'three_body_reaction':
                M = P/(8.314/T_list)
                k0 = np.array(k0)/M
            mech_data.react.k_ref.append(np.array(k0))
            # mech_data.react.k_min.append(np.array(k0)/(10**np.array(f_T)))
            # mech_data.react.k_max.append(np.array(k0)*(10**np.array(f_T)))
            mech_data.react.k_min.append(np.array(k0)*np.exp(np.array(-f_T) * np.log(10)))
            mech_data.react.k_max.append(np.array(k0)*np.exp(np.array(f_T) * np.log(10)))

            mech_data.react.k_ref_lp.append(np.array(k0_lp))
            # mech_data.react.k_min_lp.append(np.array(k0_lp)/(10**np.array(f_T_lp)))
            # mech_data.react.k_max_lp.append(np.array(k0_lp)*(10**np.array(f_T_lp)))
            mech_data.react.k_min_lp.append(np.array(k0_lp)*np.exp(np.array(-f_T_lp) * np.log(10)))
            mech_data.react.k_max_lp.append(np.array(k0_lp)*np.exp(np.array(f_T_lp) * np.log(10)))

            mech_data.react.k_ref_hp.append(np.array(k0_hp))
            # mech_data.react.k_min_hp.append(np.array(k0_hp)/(10**np.array(f_T_hp)))
            # mech_data.react.k_max_hp.append(np.array(k0_hp)*(10**np.array(f_T_hp)))
            mech_data.react.k_min_hp.append(np.array(k0_hp)*np.exp(np.array(-f_T_hp) * np.log(10)))
            mech_data.react.k_max_hp.append(np.array(k0_hp)*np.exp(np.array(f_T_hp) * np.log(10)))

        else:
            # ---  k_ref calculation  ---
            if mech_data.react.type[_r] == 'reaction' or mech_data.react.type[_r] == 'three_body_reaction':
                A = mech_data.react.ref_kin[_r][0] ; n = mech_data.react.ref_kin[_r][1] ; Ea = mech_data.react.ref_kin[_r][2]
                for i,T in enumerate(T_list):
                    k_ref.append(A*T**n*np.exp(-Ea/(R*T)))


            elif mech_data.react.type[_r] == 'pdep_arrhenius':
                all_P = [p_min/101325, P/101325, p_max/101325] #atm

                # pressure limits of the reaction law
                list_P = mech_data.react.param[_r][0]
                list_P, list_P_unit = [],[]
                for p in range(len(mech_data.react.param[_r][0])):
                    if mech_data.react.param[_r][1][p] is not False:
                        list_P_unit.append(mech_data.react.param[_r][1][p])
                        if   list_P_unit[-1] == 'Pa':
                            list_P.append(mech_data.react.param[_r][0][p]/101325)
                        elif list_P_unit[-1] == 'bar':
                            list_P.append(mech_data.react.param[_r][0][p]/1.01325)
                        else: #atm
                            list_P.append(mech_data.react.param[_r][0][p])
                    else:
                        list_P_unit.append('atm')
                        list_P.append(mech_data.react.param[_r][0][p])

                for _pi,Pi in enumerate(all_P):
                    # pressure limits of the reaction law
                    idx_Pinf,idx_Psup = 0,0
                    for _p in range(len(list_P)):
                        if Pi <= float(list_P[_p]):
                            idx_Psup = _p
                            if idx_Psup==0: idx_Pinf = 0
                            else:           idx_Pinf = _p-1
                            break
                    P_inf = float(list_P[idx_Pinf])
                    P_sup = float(list_P[idx_Psup])

                    for T in T_list:
                        # k_pinf:
                        i=0 ; k_pinf=0
                        try:
                            while float(list_P[idx_Pinf-i]) == float(list_P[idx_Pinf]):
                                A  = mech_data.react.ref_kin[_r][idx_Pinf-i][0]
                                n  = mech_data.react.ref_kin[_r][idx_Pinf-i][1]
                                Ea = mech_data.react.ref_kin[_r][idx_Pinf-i][2]
                                k_pinf +=  A * (T)**n * np.exp(-Ea/(R*T))
                                i+=1
                        except:
                            bug = True
                            k_pinf = 1

                        # k_psup:
                        i=0 ; k_psup=0
                        try:
                            while float(list_P[idx_Psup+i]) == float(list_P[idx_Psup]):
                                A  = mech_data.react.ref_kin[_r][idx_Psup+i][0]
                                n  = mech_data.react.ref_kin[_r][idx_Psup+i][1]
                                Ea = mech_data.react.ref_kin[_r][idx_Psup+i][2]
                                k_psup +=  A * (T)**n * np.exp(-Ea/(R*T))
                                i+=1
                                if idx_Psup+i>=len(list_P):
                                    break
                        except:
                            k_psup=1
                        # k
                        if P_inf == P_sup:
                            if _pi==0: k_ref_lp.append(k_pinf)
                            if _pi==1: k_ref.append(k_pinf)
                            if _pi==2: k_ref_hp.append(k_pinf)
                        else:
                            if k_pinf < 1e-150: k_pinf = 1e-150  # avoid log(0) calculation
                            if k_psup < 1e-150: k_psup = 1e-150  # avoid log(0) calculation
                            if _pi==0: k_ref_lp.append(np.exp(np.log(k_pinf)+(np.log(k_psup)-np.log(k_pinf))*((np.log(P)-np.log(P_inf))/(np.log(P_sup)-np.log(P_inf)))))
                            if _pi==1: k_ref.append(np.exp(np.log(k_pinf)+(np.log(k_psup)-np.log(k_pinf))*((np.log(P)-np.log(P_inf))/(np.log(P_sup)-np.log(P_inf)))))
                            if _pi==2: k_ref_hp.append(np.exp(np.log(k_pinf)+(np.log(k_psup)-np.log(k_pinf))*((np.log(P)-np.log(P_inf))/(np.log(P_sup)-np.log(P_inf)))))

            elif mech_data.react.type[_r] == 'falloff_reaction' or mech_data.react.type[_r] == 'chemically_activated_reaction':

                # get Troe param
                if mech_data.react.param[_r] is not False:
                    Troe_func = True
                    At = float(mech_data.react.param[_r]['A'])
                    T3t = float(mech_data.react.param[_r]['T3'])
                    T1t = float(mech_data.react.param[_r]['T1'])
                    T2t = float(mech_data.react.param[_r]['T1'])
                else:
                    Troe_func = False

                # calculation of k
                Af = mech_data.react.ref_kin[_r][0][0] ; nf = mech_data.react.ref_kin[_r][0][1] ; Eaf = mech_data.react.ref_kin[_r][0][2]
                A0 = mech_data.react.ref_kin[_r][1][0] ; n0 = mech_data.react.ref_kin[_r][1][1] ; Ea0 = mech_data.react.ref_kin[_r][1][2]
                for T in T_list:
                    M     = P/(8.314*T)           # mol/m3
                    M     = M*1e-6                # mol/cm3
                    M_lp  = p_min        / (8.314*T)*1e-6    # mol/cm3
                    M_hp  = p_max        / (8.314*T)*1e-6    # mol/cm3
                    kinf  = Af * T**nf * np.exp(-Eaf/(R*T))
                    k0    = A0 * T**n0 * np.exp(-Ea0/(R*T))
                    Pr    = (k0*M)/kinf
                    Pr_lp = (k0*M_lp)/kinf
                    Pr_hp = (k0*M_hp)/kinf
                    warnings.filterwarnings("ignore") # suppress warnings with exponential calculation
                    if Troe_func:
                        if T2t: Fcent = (1-At)*np.exp(-T/T3t) + At*np.exp(-T/T1t) + np.exp(-T2t/T1t)
                        else:   Fcent = (1-At)*np.exp(-T/T3t) + At*np.exp(-T/T1t)
                        C      = -0.4 - 0.67*np.log10(Fcent)
                        N      = 0.75 - 1.27*np.log10(Fcent)
                        f1     = (np.log10(Pr) + C) / (N-0.14*(np.log10(Pr)+C))
                        f1_lp = (np.log10(Pr_lp) + C) / (N-0.14*(np.log10(Pr_lp)+C))
                        f1_hp = (np.log10(Pr_hp) + C) / (N-0.14*(np.log10(Pr_hp)+C))
                        F      = 10**(np.log10(Fcent) / (1+f1**2))
                        F_lp   = 10**(np.log10(Fcent) / (1+f1_lp**2))
                        F_hp   = 10**(np.log10(Fcent) / (1+f1_hp**2))
                    else:
                        F,F_lp,F_hp = 1,1,1
                    if mech_data.react.type[_r] == 'falloff_reaction':
                        k_ref.append(kinf*(Pr/(1+Pr))*F)
                        k_ref_lp.append(kinf*(Pr_lp/(1+Pr_lp))*F_lp)
                        k_ref_hp.append(kinf*(Pr_hp/(1+Pr_hp))*F_hp)
                    elif mech_data.react.type[_r] == 'chemically_activated_reaction':
                        k_ref.append(k0*(1/(1+Pr))*F)
                        k_ref_lp.append(k0*(1/(1+Pr_lp))*F_lp)
                        k_ref_hp.append(k0*(1/(1+Pr_hp))*F_hp)

                    warnings.resetwarnings()

            elif 'chebyshev' in mech_data.react.type[_r]:

                # get param
                Tmin = float(mech_data.react.param[_r]['Tmin'])
                Tmax = float(mech_data.react.param[_r]['Tmax'])
                Pmin = float(mech_data.react.param[_r]['Pmin_value'])
                Pmax = float(mech_data.react.param[_r]['Pmax_value'])
                if mech_data.react.param[_r]['Pmin_unit'] is not False:
                    if mech_data.react.param[_r]['Pmin_unit'].lower() == "pa":
                        Pmin = Pmin/101325
                    if mech_data.react.param[_r]['Pmin_unit'].lower() == "bar":
                        Pmin = Pmin/1.01325
                if mech_data.react.param[_r]['Pmax_unit'] is not False:
                    if mech_data.react.param[_r]['Pmax_unit'].lower() == "pa":
                        Pmax = Pmax/101325
                    if mech_data.react.param[_r]['Pmax_unit'].lower() == "bar":
                        Pmax = Pmax/1.01325

                # k calculation (Carstensen H-H. Reaction Rate Representation Using Chebyshev Polynomials)
                # size of NT Np coefficients
                N_T = len(mech_data.react.ref_kin[_r]) ; N_P = len(mech_data.react.ref_kin[_r][0])
                for T in T_list:
                    T_tilde = (2/T - 1/Tmin - 1/Tmax)/(1/Tmax - 1/Tmin)
                    P_tilde = (2*np.log10(P) - np.log10(Pmin) - np.log10(Pmax))/(np.log10(Pmax) - np.log10(Pmin))

                    log_k = 0
                    for t in range(N_T):
                        phi_tT = np.cos(t*np.arccos(T_tilde))
                        for p in range(N_P):
                            phi_pP = np.cos(p*np.arccos(P_tilde))
                            a_tp = mech_data.react.ref_kin[_r][t][p]
                            log_k += a_tp*phi_tT*phi_pP
                    k_ref.append(10**log_k)

            else:
                k_ref = list(np.array([0]*len(T_list)))

            mech_data.react.k_ref.append(np.array(k_ref))
            # use logarithm to avoid overflow
            mech_data.react.k_min.append(np.array(k_ref)*np.exp(-np.array(f_T)*np.log(10)))
            mech_data.react.k_max.append(np.array(k_ref)*np.exp(np.array(f_T)*np.log(10)))
            if len(k_ref_lp)>0:
                mech_data.react.k_ref_lp.append(np.array(k_ref_lp))
                mech_data.react.k_min_lp.append(np.array(k_ref_lp)*np.exp(-np.array(f_T_lp)*np.log(10)))
                mech_data.react.k_max_lp.append(np.array(k_ref_lp)*np.exp(np.array(f_T_lp)*np.log(10)))

                mech_data.react.k_ref_hp.append(np.array(k_ref_hp))
                mech_data.react.k_min_hp.append(np.array(k_ref_hp)*np.exp(-np.array(f_T_hp)*np.log(10)))
                mech_data.react.k_max_hp.append(np.array(k_ref_hp)*np.exp(np.array(f_T_hp)*np.log(10)))

            else:
                mech_data.react.k_min_lp.append(False)
                mech_data.react.k_max_lp.append(False)
                mech_data.react.k_ref_hp.append(np.array(k_ref_hp))
                mech_data.react.k_min_hp.append(False)
                mech_data.react.k_max_hp.append(False)


    return mech_data





def get_diff_refopt(mech_data_ref,mech_data_opt,Tmin=300,Tmax=2000):

    # get diff
    r_opt=0 ; diff_k=[] ; diff_A=[] ; diff_n=[] ; diff_Ea=[]
    for r_ref in range(len(mech_data_ref.react.number)):
        diff_r = [] ; diff_rA=[] ; diff_rn=[] ; diff_rEa=[]
        for r_opt in range(len(mech_data_opt.react.number)):
            if mech_data_ref.react.number[r_ref]==mech_data_opt.react.number[r_opt]:
                if mech_data_ref.react.type[r_ref] == "reaction" \
                or mech_data_ref.react.type[r_ref] == "three_body_reaction":
                    for i in range(len(mech_data_ref.react.kin[r_ref])):
                        if mech_data_ref.react.kin[r_ref][i]==0:
                            diff_r=0
                            if   i==0: diff_rA.append(0)
                            elif i==1: diff_rn.append(0)
                            elif i==2: diff_rEa.append(0)
                        else:
                            diff_r = abs((mech_data_ref.react.kin[r_ref][i]-mech_data_opt.react.kin[r_opt][i])\
                                  /mech_data_ref.react.kin[r_ref][i])*100
                            if   i==0: diff_rA.append(diff_r)
                            elif i==1: diff_rn.append(diff_r)
                            elif i==2: diff_rEa.append(diff_r)
                    _rnumber = mech_data_opt.react.number[r_opt]

                    A  = mech_data_ref.react.kin[r_ref][0]
                    n  = mech_data_ref.react.kin[r_ref][1]
                    Ea = mech_data_ref.react.kin[r_ref][2]
                    A_opt  = mech_data_opt.react.kin[r_opt][0]
                    n_opt  = mech_data_opt.react.kin[r_opt][1]
                    Ea_opt = mech_data_opt.react.kin[r_opt][2]
                    incr = 50
                    T_it = np.arange(Tmin, Tmax+incr/2, incr)
                    diff_rk = 0
                    for T in T_it:
                        k     = A*T**n*np.exp(-Ea/(34.786*T))
                        k_opt = A_opt*T**n_opt*np.exp(-Ea_opt/(34.786*T))
                        diff_rkT = max(abs(k_opt/k),abs(k/k_opt))
                        if diff_rkT > diff_rk:   diff_rk=diff_rkT

                    diff_k.append((mech_data_ref.react.number[r_ref],diff_rk))

                if mech_data_ref.react.type[r_ref] == "falloff_reaction" \
                or mech_data_ref.react.type[r_ref] == "chemically_activated_reaction"\
                or mech_data_ref.react.type[r_ref] == "chebyshev"\
                or mech_data_ref.react.type[r_ref] == "pdep_arrhenius":
                    diff_m_l = 0
                    for l in range(len(mech_data_ref.react.kin[r_ref])):
                        for i in range(len(mech_data_ref.react.kin[r_ref][l])):
                            if mech_data_ref.react.kin[r_ref][l][i]==0:
                                diff_r=0
                                if   i==0: diff_rA.append(0)
                                elif i==1: diff_rn.append(0)
                                elif i==2: diff_rEa.append(0)
                            else:

                                diff_r=abs((mech_data_ref.react.kin[r_ref][l][i]-mech_data_opt.react.kin[r_opt][l][i])\
                                          /mech_data_ref.react.kin[r_ref][l][i])*100
                                if   i==0: diff_rA.append(diff_r)
                                elif i==1: diff_rn.append(diff_r)
                                elif i==2: diff_rEa.append(diff_r)

                        A  = mech_data_ref.react.kin[r_ref][l][0]
                        n  = mech_data_ref.react.kin[r_ref][l][1]
                        Ea = mech_data_ref.react.kin[r_ref][l][2]
                        A_opt  = mech_data_opt.react.kin[r_opt][l][0]
                        n_opt  = mech_data_opt.react.kin[r_opt][l][1]
                        Ea_opt = mech_data_opt.react.kin[r_opt][l][2]
                        incr = 50
                        T_it = np.arange(Tmin, Tmax+incr/2, incr)
                        diff_rk = 0
                        for T in T_it:
                            k        = A*T**n*np.exp(-Ea/(34.786*T))
                            k_opt    = A_opt*T**n_opt*np.exp(-Ea_opt/(34.786*T))
                            diff_rkT = max(abs(k_opt/k),abs(k/k_opt))
                            if diff_rkT > diff_rk:   diff_rk=diff_rkT
                        diff_m_l = np.max([diff_m_l,diff_rk])

                    diff_k.append((mech_data_ref.react.number[r_ref],diff_m_l))

                diff_A.append((mech_data_ref.react.number[r_ref],np.mean(diff_rA)))
                diff_n.append((mech_data_ref.react.number[r_ref],np.mean(diff_rn)))
                diff_Ea.append((mech_data_ref.react.number[r_ref],np.mean(diff_rEa)))
                break

        if diff_r == []:
            diff_k.append((mech_data_ref.react.number[r_ref],False))
            diff_A.append((mech_data_ref.react.number[r_ref],False))
            diff_n.append((mech_data_ref.react.number[r_ref],False))
            diff_Ea.append((mech_data_ref.react.number[r_ref],False))

    return diff_k , diff_A , diff_n , diff_Ea



def get_kin_coeffs(mech_data_ref,mech_data_opt,Tmin=300,Tmax=2000):

    # get diff
    r_opt=0 ; kT=[] ; A=[] ; n=[] ; Ea=[]
    kT_ref=[] ; A_ref=[] ; n_ref=[] ; Ea_ref=[]
    for r_ref in range(len(mech_data_ref.react.number)):
        reac_found = False
        for r_opt in range(len(mech_data_opt.react.number)):
            if mech_data_ref.react.number[r_ref]==mech_data_opt.react.number[r_opt]:
                reac_found = True

                if mech_data_ref.react.type[r_ref] == "reaction" \
                or mech_data_ref.react.type[r_ref] == "three_body_reaction":
                    #red
                    A.append((mech_data_ref.react.number[r_ref],mech_data_opt.react.kin[r_opt][0]))
                    n.append((mech_data_ref.react.number[r_ref],mech_data_opt.react.kin[r_opt][1]))
                    Ea.append((mech_data_ref.react.number[r_ref],mech_data_opt.react.kin[r_opt][2]))
                    incr = 50
                    T_it = np.arange(Tmin, Tmax+incr/2, incr)
                    kT.append([])
                    for T in T_it:
                        k = A[-1][1]*T**n[-1][1]*np.exp(-Ea[-1][1]/(34.786*T))
                        kT[-1].append((mech_data_ref.react.number[r_ref],k))
                    #ref
                    A_ref.append((mech_data_ref.react.number[r_ref],mech_data_ref.react.kin[r_ref][0]))
                    n_ref.append((mech_data_ref.react.number[r_ref],mech_data_ref.react.kin[r_ref][1]))
                    Ea_ref.append((mech_data_ref.react.number[r_ref],mech_data_ref.react.kin[r_ref][2]))
                    incr = 50
                    T_it = np.arange(Tmin, Tmax+incr/2, incr)
                    kT_ref.append([])
                    for T in T_it:
                        k_ref = A_ref[-1][1]*T**n_ref[-1][1]*np.exp(-Ea_ref[-1][1]/(34.786*T))
                        kT_ref[-1].append((mech_data_ref.react.number[r_ref],k_ref))

                if mech_data_ref.react.type[r_ref] == "falloff_reaction" \
                or mech_data_ref.react.type[r_ref] == "chemically_activated_reaction"\
                or mech_data_ref.react.type[r_ref] == "chebyshev"\
                or mech_data_ref.react.type[r_ref] == "pdep_arrhenius":
                    A_l     = [] ; n_l     = [] ; Ea_l     = [] ; kT_l     = []
                    A_l_ref = [] ; n_l_ref = [] ; Ea_l_ref = [] ; kT_l_ref = []
                    for l in range(len(mech_data_ref.react.kin[r_ref])):
                        A_l.append(mech_data_opt.react.kin[r_opt][l][0])
                        n_l.append(mech_data_opt.react.kin[r_opt][l][1])
                        Ea_l.append(mech_data_opt.react.kin[r_opt][l][2])
                        incr = 50
                        T_it = np.arange(Tmin, Tmax+incr/2, incr)
                        kT_l.append([])
                        for T in T_it:
                            kT_l[-1].append(A_l[-1]*T**n_l[-1]*np.exp(-Ea_l[-1]/(34.786*T)))
                        A_l_ref.append(mech_data_ref.react.kin[r_ref][l][0])
                        n_l_ref.append(mech_data_ref.react.kin[r_ref][l][1])
                        Ea_l_ref.append(mech_data_ref.react.kin[r_ref][l][2])
                        incr = 50
                        T_it = np.arange(Tmin, Tmax+incr/2, incr)
                        kT_l_ref.append([])
                        for T in T_it:
                            kT_l_ref[-1].append(A_l_ref[-1]*T**n_l_ref[-1]*np.exp(-Ea_l_ref[-1]/(34.786*T)))

                    A.append((mech_data_ref.react.number[r_ref],np.mean(A_l)))
                    n.append((mech_data_ref.react.number[r_ref],np.mean(n_l)))
                    Ea.append((mech_data_ref.react.number[r_ref],np.mean(Ea_l)))
                    kT.append([])
                    for T in range(len(T_it)):
                        kT_lT = []
                        for l in range(len(mech_data_ref.react.kin[r_ref])):
                            kT_lT.append(kT_l[l][T])
                        kT[-1].append((mech_data_ref.react.number[r_ref],np.mean(kT_lT)))

                    A_ref.append((mech_data_ref.react.number[r_ref],np.mean(A_l_ref)))
                    n_ref.append((mech_data_ref.react.number[r_ref],np.mean(n_l_ref)))
                    Ea_ref.append((mech_data_ref.react.number[r_ref],np.mean(Ea_l_ref)))
                    kT_ref.append([])
                    for T in range(len(T_it)):
                        kT_lT_ref = []
                        for l in range(len(mech_data_ref.react.kin[r_ref])):
                            kT_lT_ref.append(kT_l_ref[l][T])
                        kT_ref[-1].append((mech_data_ref.react.number[r_ref],np.mean(kT_lT_ref)))

                break

        if not reac_found:
            A.append((mech_data_ref.react.number[r_ref],False))
            n.append((mech_data_ref.react.number[r_ref],False))
            Ea.append((mech_data_ref.react.number[r_ref],False))
            kT.append((mech_data_ref.react.number[r_ref],False))
            A_ref.append((mech_data_ref.react.number[r_ref],False))
            n_ref.append((mech_data_ref.react.number[r_ref],False))
            Ea_ref.append((mech_data_ref.react.number[r_ref],False))
            kT_ref.append((mech_data_ref.react.number[r_ref],False))


    return A, n, Ea, kT, A_ref, n_ref, Ea_ref, kT_ref


# def mina(a, b):
#     return a if abs(a) < abs(b) else b

# def maxa(a, b):
#     return a if abs(a) < abs(b) else b



def plot_diff_refopt(mech_data_ref,mech_data_opt,diff_all,reaction_label="number",\
                     reaction_2plot="all",fig_width=7,verbose=5):

    #diff_all = (diff_m , diff_A , diff_n , diff_Ea)

    import matplotlib.pyplot as plt; plt.rcdefaults()
#    from mpltools import style; style.use('jfm')
    from matplotlib.ticker import (MultipleLocator,FormatStrFormatter,AutoMinorLocator)


    d=0
    for diff in diff_all:
        diff.sort(reverse=True, key=lambda col: col[1])

        if reaction_2plot == "all": reaction_2plot = len(diff)
        diff_vect = [] ; reactions = []
        for r in range(reaction_2plot):
            if diff[r][1]!=False:
                for r_num in range(len(mech_data_ref.react.number)):
                    if mech_data_ref.react.number[r_num]==diff[r][0]:
                        break
                diff_vect.append(diff[r][1])
                if reaction_label == "number":
                    reactions.append(str(mech_data_ref.react.number[r_num]))
                elif reaction_label == "equation":
                    reactions.append(str(mech_data_ref.react.equation[r_num]))
                elif reaction_label == "reaction+equation":
                    reactions.append("("+str(mech_data_ref.react.number[r_num])\
                                     +") "+mech_data_ref.react.equation[r_num])

        y_pos = np.arange(len(reactions))

        fig_title = 'Comparison_'+mech_data_ref.name[0:-4]+'_'+\
                     mech_data_opt.name.split('/')[-1][0:-4]
        if   d==0:
            fig_title+='_k.png'  ; xlabel = r'$\left( k_\mathrm{opt}/k_\mathrm{ref} \right)_\mathrm{max}$'
        elif d==1:
            fig_title+='_A.png'  ; xlabel = 'differences (%)'
        elif d==2:
            fig_title+='_n.png'  ; xlabel = 'differences (%)'
        elif d==3:
            fig_title+='_Ea.png' ; xlabel = 'differences (%)'

        fig, ax = plt.subplots()
        ax.barh(y_pos, diff_vect, align='center', alpha=1)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(reactions)
        ax.set_xlabel(xlabel)
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        fig.set_size_inches(7, 2+len(reactions)/6)
        fig.subplots_adjust(left = 0.5, bottom = 0.15,
                       right = 0.9, top = 0.95, wspace = 0, hspace = 0.5)

        if verbose>5:  plt.show()

        fig.savefig(fig_title,dpi=300)


        # Diff between submech
        diff.sort(key=lambda col: col[0])

        max_subm_C = max(mech_data_ref.react.subm_C)
        sub_C_diff = []
        for _i in range(max_subm_C+1): sub_C_diff.append([])
        for _r in range(len(mech_data_ref.react.subm_C)):
            if type(diff[_r][1]) is not bool:
                if   mech_data_ref.react.subm_C[_r]  == 0 \
                and  mech_data_ref.react.subm_N[_r]  == 0 \
                and  mech_data_ref.react.subm_S[_r]  == 0 \
                and  mech_data_ref.react.subm_Si[_r] == 0 \
                and not mech_data_ref.react.subm_CO[_r]:
                    sub_C_diff[0].append(diff[_r][1])
                else:
                    sub_C_diff[mech_data_ref.react.subm_C[_r]].append(diff[_r][1])

        max_subm_N = max(mech_data_ref.react.subm_N)
        sub_N_diff = []
        for _i in range(max_subm_N-1): sub_N_diff.append([])
        for _r in range(len(mech_data_ref.react.subm_N)):
            if type(diff[_r][1]) is not bool and mech_data_ref.react.subm_N[_r]>0:
                sub_N_diff[mech_data_ref.react.subm_N[_r]-1].append(diff[_r][1])

        max_subm_S = max(mech_data_ref.react.subm_S)
        sub_S_diff = []
        for _i in range(max_subm_S-1): sub_S_diff.append([])
        for _r in range(len(mech_data_ref.react.subm_S)):
            if type(diff[_r][1]) is not bool and mech_data_ref.react.subm_S[_r]>0:
                sub_S_diff[mech_data_ref.react.subm_S[_r]-1].append(diff[_r][1])

        max_subm_Si = max(mech_data_ref.react.subm_Si)
        sub_Si_diff = []
        for _i in range(max_subm_Si-1): sub_Si_diff.append([])
        for _r in range(len(mech_data_ref.react.subm_Si)):
            if type(diff[_r][1]) is not bool and mech_data_ref.react.subm_Si[_r]>0:
                sub_Si_diff[mech_data_ref.react.subm_Si[_r]-1].append(diff[_r][1])

        sub_CO_diff = []
        for _r in range(len(mech_data_ref.react.subm_CO)):
            if type(diff[_r][1]) is not bool and mech_data_ref.react.subm_CO[_r]:
                sub_CO_diff.append(diff[_r][1])

        y_label = [] ; mean_v = [] ; max_v = []
        for _subC in range(len(sub_C_diff)):
            if _subC == 0:
                y_label.append('H2')
                y_label.append('CO')
                if len(sub_C_diff[_subC])>0:
                    mean_v.append(np.mean(sub_C_diff[_subC]))
                    max_v.append(np.max(sub_C_diff[_subC]))
                else:
                    mean_v.append(0)
                    max_v.append(0)
                if len(sub_CO_diff)>0:
                    mean_v.append(np.mean(sub_CO_diff))
                    max_v.append(np.max(sub_CO_diff))
                else:
                    mean_v.append(0)
                    max_v.append(0)
            else:
                y_label.append('C'+str(_subC))
                if len(sub_C_diff[_subC])>0:
                    mean_v.append(np.mean(sub_C_diff[_subC]))
                    max_v.append(np.max(sub_C_diff[_subC]))
                else:
                    mean_v.append(0)
                    max_v.append(0)


        for _subN in range(len(sub_N_diff)):
            if _subN == 0:
                y_label.append('') ; mean_v.append(0) ; max_v.append(0)
            y_label.append('N'+str(_subN+1))
            if len(sub_N_diff[_subN])>0:
                mean_v.append(np.mean(sub_N_diff[_subN]))
                max_v.append(np.max(sub_N_diff[_subN]))
            else:
                mean_v.append(0)
                max_v.append(0)


        for _subS in range(len(sub_S_diff)):
            if _subS == 0:
                y_label.append('') ; mean_v.append(0) ; max_v.append(0)
            y_label.append('N'+str(_subS+1))
            if len(sub_S_diff[_subS])>0:
                mean_v.append(np.mean(sub_S_diff[_subS]))
                max_v.append(np.max(sub_S_diff[_subS]))
            else:
                mean_v.append(0)
                max_v.append(0)



        for _subSi in range(len(sub_Si_diff)):
            if _subSi == 0:
                y_label.append('') ; mean_v.append(0) ; max_v.append(0)
            y_label.append('N'+str(_subSi+1))
            if len(sub_Si_diff[_subSi])>0:
                mean_v.append(np.mean(sub_Si_diff[_subSi]))
                max_v.append(np.max(sub_Si_diff[_subSi]))
            else:
                mean_v.append(0)
                max_v.append(0)


        y_pos =  np.arange(len(y_label))

        ## Plot mean values
        #fig_title = fig_title.replace('.png','_mean.png')
        #fig, ax = plt.subplots()
        #ax.barh(y_pos, mean_v, align='center', alpha=1)
        #ax.set_yticks(y_pos)
        #ax.set_yticklabels(y_label)
        #ax.set_xlabel(xlabel)
        #fig.set_size_inches(7, 2+len(y_label)/6)
        #fig.subplots_adjust(left = 0.5, bottom = 0.15,
                       #right = 0.9, top = 0.95, wspace = 0, hspace = 0.5)
        #if verbose>5:  plt.show()
        #fig.savefig(fig_title,dpi=300)

        ## Plot max values
        #fig_title = fig_title.replace('mean.png','max.png')
        #fig, ax = plt.subplots()
        #ax.barh(y_pos, max_v, align='center', alpha=1)
        #ax.set_yticks(y_pos)
        #ax.set_yticklabels(y_label)
        #ax.set_xlabel(xlabel)
        #fig.set_size_inches(7, 2+len(y_label)/6)
        #fig.subplots_adjust(left = 0.5, bottom = 0.15,
                       #right = 0.9, top = 0.95, wspace = 0, hspace = 0.5)
        #if verbose>5:  plt.show()
        #fig.savefig(fig_title,dpi=300)

        # Plot mean and max values
        width = 0.5
        fig_title = fig_title.replace('.png','mean_max.png')
        fig, ax = plt.subplots()
        ax.barh(y_pos-width/2, mean_v, width, align='center', alpha=1,label='mean')
        ax.barh(y_pos+width/2, max_v,  width, align='center', alpha=1,label='max')
        ax.set_yticks(y_pos)
        ax.set_yticklabels(y_label)
        ax.set_xlabel(xlabel)
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.legend()
        fig.set_size_inches(fig_width, 2+len(y_label)/6)
        fig.subplots_adjust(left = 0.5, bottom = 0.15,
                       right = 0.9, top = 0.95, wspace = 0, hspace = 0.5)
        if verbose>5:  plt.show()

        fig.savefig(fig_title,dpi=300)

        d+=1



def plot_stats_kin(mech_data_ref,mech_data_opt,kin_list,verbose=5):

    import matplotlib.pyplot as plt; plt.rcdefaults()
#    from mpltools import style; style.use('jfm')


    max_subm_C = max(mech_data_ref.react.subm_C)
    sub_C_A_std  = [] ; sub_C_A_std_ref  = []
    sub_C_n_std  = [] ; sub_C_n_std_ref  = []
    sub_C_Ea_std = [] ; sub_C_Ea_std_ref = []
    sub_C_k_std  = [] ; sub_C_k_std_ref  = []
    for _i in range(max_subm_C+1):
        sub_C_A_std.append([])  ; sub_C_A_std_ref.append([])
        sub_C_n_std.append([])  ; sub_C_n_std_ref.append([])
        sub_C_Ea_std.append([]) ; sub_C_Ea_std_ref.append([])
        sub_C_k_std.append([])  ; sub_C_k_std_ref.append([])

    max_subm_N = max(mech_data_ref.react.subm_N)
    sub_N_A_std  = [] ; sub_N_A_std_ref  = []
    sub_N_n_std  = [] ; sub_N_n_std_ref  = []
    sub_N_Ea_std = [] ; sub_N_Ea_std_ref = []
    sub_N_k_std  = [] ; sub_N_k_std_ref  = []
    for _i in range(max_subm_N-1):
        sub_N_A_std.append([])  ; sub_N_A_std_ref.append([])
        sub_N_n_std.append([])  ; sub_N_n_std_ref.append([])
        sub_N_Ea_std.append([]) ; sub_N_Ea_std_ref.append([])
        sub_N_k_std.append([])  ; sub_N_k_std_ref.append([])

    max_subm_S = max(mech_data_ref.react.subm_S)
    sub_S_A_std  = [] ; sub_S_A_std_ref  = []
    sub_S_n_std  = [] ; sub_S_n_std_ref  = []
    sub_S_Ea_std = [] ; sub_S_Ea_std_ref = []
    sub_S_k_std  = [] ; sub_S_k_std_ref  = []
    for _i in range(max_subm_S-1):
        sub_S_A_std.append([])  ; sub_S_A_std_ref.append([])
        sub_S_n_std.append([])  ; sub_S_n_std_ref.append([])
        sub_S_Ea_std.append([]) ; sub_S_Ea_std_ref.append([])
        sub_S_k_std.append([])  ; sub_S_k_std_ref.append([])

    max_subm_Si = max(mech_data_ref.react.subm_Si)
    sub_Si_A_std  = [] ; sub_Si_A_std_ref  = []
    sub_Si_n_std  = [] ; sub_Si_n_std_ref  = []
    sub_Si_Ea_std = [] ; sub_Si_Ea_std_ref = []
    sub_Si_k_std  = [] ; sub_Si_k_std_ref  = []
    for _i in range(max_subm_Si-1):
        sub_Si_A_std.append([])  ; sub_Si_A_std_ref.append([])
        sub_Si_n_std.append([])  ; sub_Si_n_std_ref.append([])
        sub_Si_Ea_std.append([]) ; sub_Si_Ea_std_ref.append([])
        sub_Si_k_std.append([])  ; sub_Si_k_std_ref.append([])

    sub_CO_A_std  = [] ; sub_CO_A_std_ref  = []
    sub_CO_n_std  = [] ; sub_CO_n_std_ref  = []
    sub_CO_Ea_std = [] ; sub_CO_Ea_std_ref = []
    sub_CO_k_std  = [] ; sub_CO_k_std_ref  = []

    A_std_mec    =[]; n_std_mec    =[]; Ea_std_mec    =[]; kT_std_mec    =[]
    A_std_ref_mec=[]; n_std_ref_mec=[]; Ea_std_ref_mec=[]; kT_std_ref_mec=[]
    for _r in range(len(mech_data_ref.react.number)):
        if kin_list[0][0][_r][1]!=False:

            # Computation of reaction stats
            A = [] ; n = [] ; Ea = [] ; k_T = []
            for T in range(len(kin_list[0][3][1])): k_T.append([])
            A_ref  = kin_list[0][4][_r][1]
            n_ref  = kin_list[0][5][_r][1]
            Ea_ref = kin_list[0][6][_r][1]
            k_ref  = kin_list[0][7][_r][1]

            for _l in range(len(kin_list)):
                if A_ref==0: A.append(0)
                else: A.append(kin_list[_l][0][_r][1]/A_ref)
                if n_ref==0: n.append(0)
                else: n.append(kin_list[_l][1][_r][1]/n_ref)
                if Ea_ref==0: Ea.append(0)
                else: Ea.append(kin_list[_l][2][_r][1]/Ea_ref)
                for T in range(len(kin_list[_l][3][_r][1])):
                    if k_ref[T]==0: k_T[T].append(0)
                    else: k_T[T].append(kin_list[_l][3][_r][1][T]/k_ref[T])  # k_T : reaction rate at T

            A_std_mec.append(np.std(A))
            n_std_mec.append(np.std(n))
            Ea_std_mec.append(np.std(Ea))
            A_std_ref_l = [] ; n_std_ref_l = [] ; Ea_std_ref_l = []
            for _i in range(len(A)):
                A_std_ref_l.append(np.abs(A[_i]-1)**2)
                n_std_ref_l.append(np.abs(n[_i]-1)**2)
                Ea_std_ref_l.append(np.abs(Ea[_i]-1)**2)

            A_std_ref_mec.append(np.mean(A_std_ref_l))
            n_std_ref_mec.append(np.mean(n_std_ref_l))
            Ea_std_ref_mec.append(np.mean(Ea_std_ref_l))

            kT_std = [] ; kT_std_ref = []
            for T in range(len(kin_list[_l][3][_r][1])):
                kT_std.append(np.std(k_T[T]))
                kT_std_ref_l = []
                for _i in range(len(k_T[T])):
                    kT_std_ref_l.append(np.abs(k_T[T][_i]-1)**2)
                kT_std_ref = np.mean(kT_std_ref_l)

            kT_std_mec.append(np.mean(kT_std))
            kT_std_ref_mec.append(np.mean(kT_std_ref))


            # sort stats in function of the submechanism
            if   mech_data_ref.react.subm_C[_r]  == 0 \
            and  mech_data_ref.react.subm_N[_r]  == 0 \
            and  mech_data_ref.react.subm_S[_r]  == 0 \
            and  mech_data_ref.react.subm_Si[_r] == 0 \
            and not mech_data_ref.react.subm_CO[_r]:
                sub_C_A_std[0].append(A_std_mec[-1])
                sub_C_n_std[0].append(n_std_mec[-1])
                sub_C_Ea_std[0].append(Ea_std_mec[-1])
                sub_C_k_std[0].append(kT_std_mec[-1])
                sub_C_A_std_ref[0].append(A_std_ref_mec[-1])
                sub_C_n_std_ref[0].append(n_std_ref_mec[-1])
                sub_C_Ea_std_ref[0].append(Ea_std_ref_mec[-1])
                sub_C_k_std_ref[0].append(kT_std_ref_mec[-1])
            elif mech_data_ref.react.subm_C[_r] != 0:
                sub_C_A_std[mech_data_ref.react.subm_C[_r]].append(A_std_mec[-1])
                sub_C_n_std[mech_data_ref.react.subm_C[_r]].append(n_std_mec[-1])
                sub_C_Ea_std[mech_data_ref.react.subm_C[_r]].append(Ea_std_mec[-1])
                sub_C_k_std[mech_data_ref.react.subm_C[_r]].append(kT_std_mec[-1])
                sub_C_A_std_ref[mech_data_ref.react.subm_C[_r]].append(A_std_ref_mec[-1])
                sub_C_n_std_ref[mech_data_ref.react.subm_C[_r]].append(n_std_ref_mec[-1])
                sub_C_Ea_std_ref[mech_data_ref.react.subm_C[_r]].append(Ea_std_ref_mec[-1])
                sub_C_k_std_ref[mech_data_ref.react.subm_C[_r]].append(kT_std_ref_mec[-1])
            elif mech_data_ref.react.subm_CO:
                sub_CO_A_std.append(A_std_mec[-1])
                sub_CO_n_std.append(n_std_mec[-1])
                sub_CO_Ea_std.append(Ea_std_mec[-1])
                sub_CO_k_std.append(kT_std_mec[-1])
                sub_CO_A_std_ref.append(A_std_ref_mec[-1])
                sub_CO_n_std_ref.append(n_std_ref_mec[-1])
                sub_CO_Ea_std_ref.append(Ea_std_ref_mec[-1])
                sub_CO_k_std_ref.append(kT_std_ref_mec[-1])
            elif mech_data_ref.react.subm_N[_r]  != 0:
                sub_N_A_std[mech_data_ref.react.subm_N[_r]-1].append(A_std_mec[-1])
                sub_N_n_std[mech_data_ref.react.subm_N[_r]-1].append(n_std_mec[-1])
                sub_N_Ea_std[mech_data_ref.react.subm_N[_r]-1].append(Ea_std_mec[-1])
                sub_N_k_std[mech_data_ref.react.subm_N[_r]-1].append(kT_std_mec[-1])
                sub_N_A_std_ref[mech_data_ref.react.subm_N[_r]-1].append(A_std_ref_mec[-1])
                sub_N_n_std_ref[mech_data_ref.react.subm_N[_r]-1].append(n_std_ref_mec[-1])
                sub_N_Ea_std_ref[mech_data_ref.react.subm_N[_r]-1].append(Ea_std_ref_mec[-1])
                sub_N_k_std_ref[mech_data_ref.react.subm_N[_r]-1].append(kT_std_ref_mec[-1])
            elif mech_data_ref.react.subm_S[_r]  != 0:
                sub_S_A_std[mech_data_ref.react.subm_S[_r]-1].append(A_std_mec[-1])
                sub_S_n_std[mech_data_ref.react.subm_S[_r]-1].append(n_std_mec[-1])
                sub_S_Ea_std[mech_data_ref.react.subm_S[_r]-1].append(Ea_std_mec[-1])
                sub_S_k_std[mech_data_ref.react.subm_S[_r]-1].append(kT_std_mec[-1])
                sub_S_A_std_ref[mech_data_ref.react.subm_S[_r]-1].append(A_std_ref_mec[-1])
                sub_S_n_std_ref[mech_data_ref.react.subm_S[_r]-1].append(n_std_ref_mec[-1])
                sub_S_Ea_std_ref[mech_data_ref.react.subm_S[_r]-1].append(Ea_std_ref_mec[-1])
                sub_S_k_std_ref[mech_data_ref.react.subm_S[_r]-1].append(kT_std_ref_mec[-1])
            elif mech_data_ref.react.subm_Si[_r] != 0:
                sub_Si_A_std[mech_data_ref.react.subm_Si[_r]-1].append(A_std_mec[-1])
                sub_Si_n_std[mech_data_ref.react.subm_Si[_r]-1].append(n_std_mec[-1])
                sub_Si_Ea_std[mech_data_ref.react.subm_Si[_r]-1].append(Ea_std_mec[-1])
                sub_Si_k_std[mech_data_ref.react.subm_Si[_r]-1].append(kT_std_mec[-1])
                sub_Si_A_std_ref[mech_data_ref.react.subm_Si[_r]-1].append(A_std_ref_mec[-1])
                sub_Si_n_std_ref[mech_data_ref.react.subm_Si[_r]-1].append(n_std_ref_mec[-1])
                sub_Si_Ea_std_ref[mech_data_ref.react.subm_Si[_r]-1].append(Ea_std_ref_mec[-1])
                sub_Si_k_std_ref[mech_data_ref.react.subm_Si[_r]-1].append(kT_std_ref_mec[-1])



    # Standart deviation between submech

    y_label = [] ; A_std = [] ; n_std = [] ; Ea_std = [] ; k_std = []
    A_std_ref = [] ; n_std_ref = [] ; Ea_std_ref = [] ; k_std_ref = []
    for _subC in range(max_subm_C):
        if _subC == 0:
            y_label.append('H2')
            y_label.append('CO')
            if len(sub_C_A_std[_subC])>0:
                A_std.append(np.mean(sub_C_A_std[_subC]))
                n_std.append(np.mean(sub_C_n_std[_subC]))
                Ea_std.append(np.mean(sub_C_Ea_std[_subC]))
                k_std.append(np.mean(sub_C_k_std[_subC]))
                A_std_ref.append(np.mean(sub_C_A_std_ref[_subC]))
                n_std_ref.append(np.mean(sub_C_n_std_ref[_subC]))
                Ea_std_ref.append(np.mean(sub_C_Ea_std_ref[_subC]))
                k_std_ref.append(np.mean(sub_C_k_std_ref[_subC]))
            else:
                A_std.append(0)  ;  A_std_ref.append(0)
                n_std.append(0)  ;  n_std_ref.append(0)
                Ea_std.append(0) ;  Ea_std_ref.append(0)
                k_std.append(0)  ;  k_std_ref.append(0)
            if len(sub_CO_A_std)>0:
                A_std.append(np.mean(sub_CO_A_std))
                n_std.append(np.mean(sub_CO_n_std))
                Ea_std.append(np.mean(sub_CO_Ea_std))
                k_std.append(np.mean(sub_CO_k_std))
                A_std_ref.append(np.mean(sub_CO_A_std_ref))
                n_std_ref.append(np.mean(sub_CO_n_std_ref))
                Ea_std_ref.append(np.mean(sub_CO_Ea_std_ref))
                k_std_ref.append(np.mean(sub_CO_k_std_ref))
            else:
                A_std.append(0)  ;  A_std_ref.append(0)
                n_std.append(0)  ;  n_std_ref.append(0)
                Ea_std.append(0) ;  Ea_std_ref.append(0)
                k_std.append(0)  ;  k_std_ref.append(0)
        else:
            y_label.append('C'+str(_subC))
            if len(sub_C_A_std[_subC])>0:
                A_std.append(np.mean(sub_C_A_std[_subC]))
                n_std.append(np.mean(sub_C_n_std[_subC]))
                Ea_std.append(np.mean(sub_C_Ea_std[_subC]))
                k_std.append(np.mean(sub_C_k_std[_subC]))
                A_std_ref.append(np.mean(sub_C_A_std_ref[_subC]))
                n_std_ref.append(np.mean(sub_C_n_std_ref[_subC]))
                Ea_std_ref.append(np.mean(sub_C_Ea_std_ref[_subC]))
                k_std_ref.append(np.mean(sub_C_k_std_ref[_subC]))
            else:
                A_std.append(0)  ;  A_std_ref.append(0)
                n_std.append(0)  ;  n_std_ref.append(0)
                Ea_std.append(0) ;  Ea_std_ref.append(0)
                k_std.append(0)  ;  k_std_ref.append(0)


    for _subN in range(max_subm_N-1):
        if _subN == 0:
            y_label.append('')
            A_std.append(0)  ;  A_std_ref.append(0)
            n_std.append(0)  ;  n_std_ref.append(0)
            Ea_std.append(0) ;  Ea_std_ref.append(0)
            k_std.append(0)  ;  k_std_ref.append(0)
        y_label.append('N'+str(_subN+1))
        if len(sub_N_A_std[_subN])>0:
                A_std.append(np.mean(sub_N_A_std[_subN]))
                n_std.append(np.mean(sub_N_n_std[_subN]))
                Ea_std.append(np.mean(sub_N_Ea_std[_subN]))
                k_std.append(np.mean(sub_N_k_std[_subN]))
                A_std_ref.append(np.mean(sub_N_A_std_ref[_subN]))
                n_std_ref.append(np.mean(sub_N_n_std_ref[_subN]))
                Ea_std_ref.append(np.mean(sub_N_Ea_std_ref[_subN]))
                k_std_ref.append(np.mean(sub_N_k_std_ref[_subN]))
        else:
            A_std.append(0)  ;  A_std_ref.append(0)
            n_std.append(0)  ;  n_std_ref.append(0)
            Ea_std.append(0) ;  Ea_std_ref.append(0)
            k_std.append(0)  ;  k_std_ref.append(0)


    for _subS in range(max_subm_S-1):
        if _subS == 0:
            y_label.append('')
            A_std.append(0)  ;  A_std_ref.append(0)
            n_std.append(0)  ;  n_std_ref.append(0)
            Ea_std.append(0) ;  Ea_std_ref.append(0)
            k_std.append(0)  ;  k_std_ref.append(0)
        y_label.append('N'+str(_subS+1))
        if len(sub_S_A_std[_subS])>0:
            A_std.append(np.mean(sub_S_A_std[_subS]))
            n_std.append(np.mean(sub_S_n_std[_subS]))
            Ea_std.append(np.mean(sub_S_Ea_std[_subS]))
            k_std.append(np.mean(sub_S_k_std[_subS]))
            A_std_ref.append(np.mean(sub_S_A_std_ref[_subS]))
            n_std_ref.append(np.mean(sub_S_n_std_ref[_subS]))
            Ea_std_ref.append(np.mean(sub_S_Ea_std_ref[_subS]))
            k_std_ref.append(np.mean(sub_S_k_std_ref[_subS]))
        else:
            A_std.append(0)  ;  A_std_ref.append(0)
            n_std.append(0)  ;  n_std_ref.append(0)
            Ea_std.append(0) ;  Ea_std_ref.append(0)
            k_std.append(0)  ;  k_std_ref.append(0)

    for _subSi in range(max_subm_Si-1):
        if _subSi == 0:
            y_label.append('')
            A_std.append(0)  ;  A_std_ref.append(0)
            n_std.append(0)  ;  n_std_ref.append(0)
            Ea_std.append(0) ;  Ea_std_ref.append(0)
            k_std.append(0)  ;  k_std_ref.append(0)
        y_label.append('N'+str(_subSi+1))
        if len(sub_Si_A_std[_subSi])>0:
            A_std.append(np.mean(sub_Si_A_std[_subSi]))
            n_std.append(np.mean(sub_Si_n_std[_subSi]))
            Ea_std.append(np.mean(sub_Si_Ea_std[_subSi]))
            k_std.append(np.mean(sub_Si_k_std[_subSi]))
            A_std_ref.append(np.mean(sub_Si_A_std_ref[_subSi]))
            n_std_ref.append(np.mean(sub_Si_n_std_ref[_subSi]))
            Ea_std_ref.append(np.mean(sub_Si_Ea_std_ref[_subSi]))
            k_std_ref.append(np.mean(sub_Si_k_std_ref[_subSi]))
        else:
            A_std.append(0)  ;  A_std_ref.append(0)
            n_std.append(0)  ;  n_std_ref.append(0)
            Ea_std.append(0) ;  Ea_std_ref.append(0)
            k_std.append(0)  ;  k_std_ref.append(0)

    y_label.append('')
    A_std.append(0)  ;  A_std_ref.append(0)
    n_std.append(0)  ;  n_std_ref.append(0)
    Ea_std.append(0) ;  Ea_std_ref.append(0)
    k_std.append(0)  ;  k_std_ref.append(0)
    y_label.append('Full mechanism')
    A_std.append(np.mean(A_std_mec))  ; A_std_ref.append(np.mean(A_std_ref_mec))
    n_std.append(np.mean(n_std_mec))  ; n_std_ref.append(np.mean(n_std_ref_mec))
    Ea_std.append(np.mean(Ea_std_mec)); Ea_std_ref.append(np.mean(Ea_std_ref_mec))
    k_std.append(np.mean(kT_std_mec)) ; k_std_ref.append(np.mean(kT_std_ref_mec))



#    def Plot_bar(y_label,var,fig_n,x_label):
#
#        y_pos =  np.arange(len(y_label))
#
#        # Plot A_std values
#        fig_title = fig_n + '.png'
#        fig, ax = plt.subplots()
#        ax.barh(y_pos, var, align='center', alpha=1)
#        ax.set_yticks(y_pos)
#        ax.set_yticklabels(y_label)
#        ax.set_xlabel(x_label)
#        fig.set_size_inches(7, 2+len(y_label)/6)
#        fig.subplots_adjust(left = 0.5, bottom = 0.15,
#                       right = 0.9, top = 0.95, wspace = 0, hspace = 0.5)
#        if verbose>5:  plt.show()
#        fig.savefig(fig_title,dpi=300)
#
#    x_label = 'Normalized standart deviation (%)'
#    Plot_bar(y_label,A_std,'A_std',x_label)
#    Plot_bar(y_label,n_std,'n_std',x_label)
#    Plot_bar(y_label,Ea_std,'Ea_std',x_label)
#    Plot_bar(y_label,k_std,'k_std',x_label)
#    x_label = 'Normalized standart deviation compared to reference mechanism (%)'
#    Plot_bar(y_label,A_std_ref,'A_std_ref',x_label)
#    Plot_bar(y_label,n_std_ref,'n_std_ref',x_label)
#    Plot_bar(y_label,Ea_std_ref,'Ea_std_ref',x_label)
#    Plot_bar(y_label,k_std_ref,'k_std_ref',x_label)

    dev_std = (A_std,n_std,Ea_std,k_std,y_label)
    dev_std_ref = (A_std_ref,n_std_ref,Ea_std_ref,k_std_ref,y_label)

    return dev_std, dev_std_ref


def plot_bar(y_label,var,fig_n,x_label,label=False):

    import matplotlib.pyplot as plt; plt.rcdefaults()
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'g']

    y_pos =  np.arange(len(y_label))

    # Plot A_std values
    fig_title = fig_n + '.png'
    fig, ax = plt.subplots()

    n_var = len(var)
    width = 1/n_var
    if not label : label_ = ['']*n_var
    else:  label_ = label
    for _i in range(len(var)):
        switch_pos = -(n_var-1)/2*width + width*_i
        ax.barh(y_pos+switch_pos, var[_i], width, align='center', alpha=0.5,label=label_[_i])

    ax.set_yticks(y_pos)
    ax.set_yticklabels(y_label)
    ax.set_xlabel(x_label)
    if label: ax.legend()
    fig.set_size_inches(7, 2+len(y_label)/6)
    fig.subplots_adjust(left = 0.5, bottom = 0.15,
                   right = 0.9, top = 0.95, wspace = 0, hspace = 0.5)
#    if verbose>5:  plt.show()
    fig.savefig(fig_title,dpi=300)




# optimisation en fonction des incertitudes :
#    * pour les réactions standard, three body, et pdep, falloff à patm:
#       -> fait par rapport au k(T) moyen (noté k0) à patm
#       les limites kmin et kmax sont calculées p/r à f(T)
#
#    * pour les réactions dépendantes de la pression (pdep, falloff, chemically activated):
#        -> fait par rapport aux k_ref calculés  borné à P=0.5 et 50bars pour éviter de dégrader le mécanisme si les plages de pressions sont plus faibles
#               P_min < 0.5 bar et Pmax > 50 bars si les conditions pré-définies impliquent de telles pressions
#
#        -> variation dans les limites des incertitudes f(T) (définies à p atm)
#
#    * chebychev : pas de modif des coeffs
#
#    * si pas de f(T) alors f = 0.7
#
