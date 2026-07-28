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

import brookesia.Class_def      as cdef
from  brookesia.Class_def import print_
import brookesia.Opt_tools      as ot

import os, sys, multiprocessing
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
    print_('    ' + str(optim_param.n_gen)+ ' generations  /  '+\
                        str(optim_param.n_ind)+' individuals\n',mp)

    os.chdir(mp)
    if not os.path.exists("GA"): os.mkdir("GA")
    os.chdir("GA")

    gen=0


    # Reference ind
    mech_data = ot.get_uncertainty(mech_data,optim_param,conditions_list,'GA')
    ref_ind = ot.Individual(conditions_list,mech_data,\
                         ref_results_list,red_data_list,'GA',False)
    ref_ind.fitness = ref_ind.fitness_eval(conditions_list,optim_param,ref_results_list,'GA',0,True)
    if verbose >= 1 :
        print_("Non optimized reduced mechanism fitness: "+"%.3f"%(ref_ind.fitness),mp)


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

    best_ind = copy.deepcopy(ref_ind)


    # Population evaluation:
    ot.fitness_eval_newchilds(pop, optim_param,conditions_list,ref_results_list, gen, 'GA')

    # Find new best ind and save the mech,
    # if not, replace worst ind of the current pop by the previous best ind
    pop,best_ind,new_best_ind = ot.compare_best_ind(pop,best_ind,optim_param,'GA',verbose)

    # save and display convergence informations
    warnings.filterwarnings("ignore", category=ResourceWarning)

    if verbose >= 1 : print_('Initial population:',mp)
    pop.convergence_information(gen,optim_param,verbose)
    pop.check_early_conv(optim_param,gen,verbose)


    time_1 = timer.time()

    for gen in range(1, optim_param.n_gen+1):
        print_("\n\nGeneration:" + str(gen),mp)
        pop.Xover(optim_param,conditions_list,ref_results_list,verbose)
        pop.mutation(optim_param,conditions_list,ref_results_list,gen,verbose)
        ot.fitness_eval_newchilds(pop, optim_param,conditions_list,ref_results_list, gen, 'GA')
        pop.selection(optim_param,verbose)
        pop.check_early_conv(optim_param,gen,verbose)
        pop,best_ind,new_best_ind = ot.compare_best_ind(pop,best_ind,optim_param,'GA',verbose)
        if new_best_ind and optim_param.exp_data:                              # optimization of the time iteration for reactor models
            best_ind.fitness = best_ind.fitness_eval(conditions_list,optim_param,ref_results_list,'GA')
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
        best_ind.export_data(conditions_list,optim_param,ref_results_list,'GA')

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





class Population:
    def __init__(self,conditions_list,mech_data,red_data_list,ref_results_list,size_pop):

        optim_param = red_data_list[0].optim_param
        mp          = optim_param.main_path


        self.individual = []

        mech_data = ot.get_uncertainty(mech_data,optim_param,conditions_list,'GA')

        # option to import modified mecanisms in the first individual
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

                    list_mech[-1].react.uncert      = copy.deepcopy(mech_data.react.uncert)

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
            # if import_mech, get the mecanism to import in the first individual
            try:    rand_kin = list_mech[ind]
            except: rand_kin = True

            self.individual.append(ot.Individual(conditions_list,mech_data,\
                                   ref_results_list,red_data_list,'GA',rand_kin))



    def __getitem__(self, i):
        return self.individual[i]





    def check_early_conv(self,optim_param,gen,verbose):

        mp      = optim_param.main_path
        MaxIt   = optim_param.n_gen
        if gen>0:
            mean_cv_1st_gen = optim_param.mean_cv_1st_gen

        # calculation of the fitness stats
        fitness_list = []
        for p in range(len(self.individual)):
            fitness_list.append(self.individual[p].fitness)
        best_fitness = np.max(fitness_list) ; std_fit = np.std(fitness_list)

        fitness_list = np.array((fitness_list))

        # Calculate the mean coefficients of variation for the kinetic parameters
        # a) construct the list of kinetic coeffs
        kinmech = []
        for p in range(len(self.individual)):
            kinmech.append([])
            # add every kin coeff in one flat list
            for _r in range(len(self.individual[p].mech.react.kin)):
                if self.individual[p].mech.react.modif[_r]:
                    kin_r = self.individual[p].mech.react.kin[_r]
                    if kin_r and isinstance(kin_r[0], list):
                        kinmech[-1].extend([x for sublist in kin_r for x in sublist])
                    else:
                        kinmech[-1].extend(kin_r.copy())
        
        
        # b) std/mean (if mean(kin=0) -> std/mean=0)
        std_d_mean = np.divide(
                np.std(kinmech, axis=0),
                np.mean(kinmech, axis=0),
                out=np.zeros_like(np.mean(kinmech, axis=0)),
                where=np.mean(kinmech, axis=0) != 0
                )
        
        # c) mean coefficient of variation
        mean_cv = np.mean(std_d_mean)
        if gen==0:
            optim_param.mean_cv_1st_gen = mean_cv
        
        if gen>0:
            if gen < .9*MaxIt and mean_cv<(0.1*mean_cv_1st_gen):#std_fit<(.05*best_fitness):
                print_('Early convergence detected, create new random individuals',mp)
                self.sort_fitness()
    
                for p_i in range(round(len(self.individual)/2)):
                     self.individual[p_i].randomize_kin(optim_param)



    def display(self, optim_param, ind = -3):
        mp = optim_param.main_path
        """ Affichage d'un ou de tout les individus de la population"""
        if ind == -3:
            print_("Fitness of the new population:\n",mp)
            fit_disp=""
            for i in range(optim_param.n_ind):
                fit_disp = fit_disp + '%5.3f' %(self.individual[i].fitness) + ' '
                if (i+1)%10==0 and i!=0:
                    print_(fit_disp,mp)
                    fit_disp=""
            print_(fit_disp,mp)
        else:
            print_("",mp)
            print_('individu'+str(ind)+': '+"%.3f" %(self.individual[ind].fitness),mp)
            print_("\n",mp)

    # def compare_best_ind(self, best_ind, optim_param, main_path, verbose=0):
    #     mp = optim_param.main_path

    #     mp = optim_param.main_path
    #     n_ind = optim_param.n_ind
    #     try:
    #         global_best_idx = self.find_best()
    #     except:
    #         print('toto_compare_best_ind')

    #     new_best_ind = False
    #     # compare the current best population ind to the previous best ind
    #     if self.individual[global_best_idx].fitness > best_ind.fitness:
    #         best_ind = copy.deepcopy(self.individual[global_best_idx])
    #         os.chdir(mp+'/GA')
    #         if '.cti' in self.individual[0].mech.name:
    #             best_ind.mech.write_new_mech("optim_mech.cti")
    #         else:
    #             best_ind.mech.write_yaml_mech("optim_mech.yaml")

    #         if verbose >= 3:
    #             print_("New best_ind: "+"%.3f" %(best_ind.fitness),mp)
    #         new_best_ind = True

    #     next_pop_best_idx = self.find_best(n_ind)
    #     if self.individual[next_pop_best_idx].fitness < best_ind.fitness:
    #         worst_idx=self.find_worst(next_pop_best_idx,n_ind)
    #         self.individual[worst_idx]=copy.deepcopy(best_ind)

    #     return best_ind, new_best_ind


    def def_pop2keep(self):
        pop2keep = []
        for ind in self.individual:
            if ind.mech.keep4opt == True:
                pop2keep.append(copy.deepcopy(ind))
        return pop2keep


    def insert_pop2keep(self, optim_param, pop2keep):
        n_ind = optim_param.n_ind

        # Fitness list
        list_fit = []
        for _p in range(n_ind):
            list_fit.append(self.individual[_p].fitness)

        # Get sorted element indices
        sort_index = sorted(range(len(list_fit)), key=lambda k: list_fit[k])

        # Create the order list from the sorted indices
        order = [0] * len(list_fit)
        for i, _index in enumerate(sort_index, start=1):
            order[_index] = i

        for _p2k, mech_2keep in enumerate(pop2keep):
            self.individual[order[_p2k]]=copy.deepcopy(mech_2keep)


    def find_best(self,n_ind=False):
        if not n_ind: n_ind=len(self.individual)
        best_fit=-1
        for p in range(n_ind):
            if self.individual[p].fitness > best_fit:
                best_idx = p ; best_fit = self.individual[p].fitness

        try:     A = best_idx+1-1
        except:  best_idx = 0

        return best_idx

    def find_worst(self,best_idx,n_ind=False):
        if not n_ind: n_ind=len(self.individual)
        best_fit = self.individual[best_idx].fitness
        worst_idx=0
        for p in range(n_ind):
            if self.individual[p].fitness < best_fit:
                worst_idx = p ; best_fit = self.individual[p].fitness
        return worst_idx

#     def fitness_eval_par(self,fit_eval_inp):
#
#         optim_param         = fit_eval_inp[0]
#         conditions_list     = fit_eval_inp[1]
#         ref_results_list    = fit_eval_inp[2]
#         bar                 = fit_eval_inp[3]
#         ind                 = fit_eval_inp[4]
#         is_child            = fit_eval_inp[5]
#         if is_child:
#             ind               = optim_param.n_ind + ind
#             title="New ind evaluation  "
#         else:
#             title="New pop evaluation  "
#
#         mech = conditions_list[0].mech
#
#         os.chdir(conditions_list[0].main_path)
#
#         # gas = cdef.get_gas_ct(mech)
#
#         # for _c in range(len(conditions_list)):
#         #     conditions_list[_c].composition.gas     = gas
#         #     conditions_list[_c].composition.gas_ref = gas
#
#         os.chdir("GA")
#
#
#         # try:
#         fitness = self.individual[ind].fitness_eval(conditions_list,optim_param,ref_results_list,'GA',ind)
#         # except:
#         #     fitness = 0
#
#         bar.update(ind,title)
#
#         return (fitness,ind)
#
#
#     def fitness_eval_newpop(self,optim_param,conditions_list,ref_results_list):
#
#         mp      = optim_param.main_path
#
#         n_ind   = optim_param.n_ind
#
#         gas     = conditions_list[0].composition.gas
#         gas_ref = conditions_list[0].composition.gas_ref
#
#         # saving and suppression of unpickable variables on fitness eval inputs
#         for cond in range(len(conditions_list)):
#             del conditions_list[cond].composition.gas
#             del conditions_list[cond].composition.gas_ref
#         f=[]
#         simul_time_limit = 0
#         for res in range(len(ref_results_list)):
#             f.append(ref_results_list[res].f)
#             del ref_results_list[res].gas
#             del ref_results_list[res].f
#             simul_time_limit += ref_results_list[res].simul_time
#         simul_time_limit = simul_time_limit*3
#
#         bar = cdef.ProgressBar(n_ind, '')
#         title="New pop evaluation  "
#         bar.update(0,title)
#
#         fit_eval_inp = []
#         for _i in range(n_ind):
#             fit_eval_inp.append([optim_param,conditions_list,ref_results_list,bar,_i,False])
#
#         # saving of the ref mechanism
#         mech = conditions_list[0].mech
#         os.chdir(conditions_list[0].main_path)
#         copyfile(mech,'GA/'+mech)
#         os.chdir("GA")
#
#
#         if sys.gettrace() is not None : is_debug_mode = True
#         else:                           is_debug_mode = False
#
#         if is_debug_mode:
#             print('Evaluation of new individuals in debug mode (i.e. no parallelization)')
#             # bypass parallelisation for debugging
#             fit_i = []
#             mech = conditions_list[0].mech
#
#             # gas = cdef.get_gas_ct(mech)
#             # for _c in range(len(conditions_list)):
#             #     conditions_list[_c].composition.gas     = gas
#             #     conditions_list[_c].composition.gas_ref = gas
#             for ind in range(n_ind):
#                 fit_i.append([self.individual[ind].fitness_eval(conditions_list,optim_param,ref_results_list,'GA',ind),ind])
#
#         else:
#             # Parallelized Fitness calculation
#             num_cores   = multiprocessing.cpu_count()
#             if os.name == 'nt': multiprocessing.get_context('spawn')
#             with multiprocessing.Pool(num_cores) as p:
#                 fit_i = p.map(self.fitness_eval_par, fit_eval_inp)
#             # sort in ascending order of individual index values
#             fit_i.sort(reverse=False, key=lambda col: col[1])
#
#         for _i in range(len(fit_i)):
#             self.individual[_i].fitness = fit_i[_i][0]
#
#         bar.update(n_ind,title)
#         print('\n')
#
#         for i in range(len(conditions_list)):
#             conditions_list[i].composition.gas     = gas
#             conditions_list[i].composition.gas_ref = gas_ref
#         for i in range(len(ref_results_list)):
#             ref_results_list[i].gas = gas
#             ref_results_list[i].f   = f[i]
#
#
#
#     def fitness_eval_newchilds(self,optim_param,conditions_list,ref_results_list):
#
#         mp          = optim_param.main_path
#         num_cores   = multiprocessing.cpu_count()
#         child_nb    = optim_param.total_Xover + optim_param.total_mut
#
#         gas         = conditions_list[0].composition.gas
#         gas_ref     = conditions_list[0].composition.gas_ref
#
#         # saving and suppression of unpickable variables on fitness eval inputs
#         for cond in range(len(conditions_list)):
#             del conditions_list[cond].composition.gas
#             del conditions_list[cond].composition.gas_ref
#         f=[]
#         simul_time_limit = 0
#         for res in range(len(ref_results_list)):
#             f.append(ref_results_list[res].f)
#             del ref_results_list[res].gas
#             del ref_results_list[res].f
#             simul_time_limit += ref_results_list[res].simul_time
# #        simul_time_limit = simul_time_limit*(6+np.random.uniform()*16)#*(child_nb/num_cores)
#         simul_time_limit = simul_time_limit*2*(child_nb/num_cores) + 30
# #        simul_time_limit = simul_time_limit*5 + 30
#
#
#         bar = cdef.ProgressBar(child_nb, '')
#         title="New ind evaluation  "
#         bar.update(0,title)
#
#         fit_eval_inp = []
#         for ch in range(child_nb):
#             fit_eval_inp.append([optim_param,conditions_list,ref_results_list,bar,ch,True])
#
#         # saving of the ref mechanism
#         mech = conditions_list[0].mech
#         os.chdir(conditions_list[0].main_path)
#         copyfile(mech,'GA/'+mech)
#         os.chdir("GA")
#
#
#         # Fitness calculation
#
#         if sys.gettrace() is not None : is_debug_mode = True
#         else:                           is_debug_mode = False
#
#         # no parallelisation --------------------------------------------------
#         if is_debug_mode:
#             fit_list = []
#             mech = conditions_list[0].mech
#
#             gas = cdef.get_gas_ct(mech)
#             for _c in range(len(conditions_list)):
#                 conditions_list[_c].composition.gas     = gas
#                 conditions_list[_c].composition.gas_ref = gas
#             for ch in range(child_nb):
#                 fit_list.append([self.individual[ch].fitness_eval(conditions_list,optim_param,ref_results_list,'GA',ch),ch])
#         else:
#
#             # Parallelisation 1 ---------------------------------------------------
#     #        with Pool(processes=num_cores) as pool:
#     #            fit_list = [ pool.apply_async(self.fitness_eval_par, args) for args in fit_eval_inp ]
#
#             # Parallelisation 2 ---------------------------------------------------
#     #        fit_list = []
#     #        def log_result(fit_i):
#     #            # This is called whenever foo_pool(i) returns a result.
#     #            # result_list is modified only by the main process, not the pool workers.
#     #            fit_list.append(fit_i)
#     #
#     #        pool=Pool(processes=num_cores)
#     #        for args in fit_eval_inp:
#     #            pool.apply_async(self.fitness_eval_par, args, callback=log_result)
#     #        pool.close()
#
#             # Parallelisation 3 ---------------------------------------------------
#     #        if os.name == 'nt': multiprocessing.get_context('spawn')
#     #        with multiprocessing.Pool(num_cores) as p:
#     #            fit_list = p.map(self.fitness_eval_par, fit_eval_inp)
#
#
#             # Parallelisation 4   (with simulation time check) --------------------
#             #  https://pythonhosted.org/Pebble/#pools
#             fit_list = []
#             with ProcessPool() as pool:
#                 sim_results = pool.map(self.fitness_eval_par, fit_eval_inp, timeout=simul_time_limit)
#                 try:
#                     for fit in sim_results.result():
#                         fit_list.append(fit)
#                 except TimeoutError:
#                     print_('\n\nWarning : simulation time > ' + '%.0f' %simul_time_limit + 's (> 1.5 x ref simulation time)',mp)
#                     print_("TimeoutError: aborting remaining computations",mp)
#                     fit_list.sort(reverse=False, key=lambda col: col[1])
#                     fitness_incomplete = copy.deepcopy(fit_list)
#                     list_ind_eval = []
#                     for ind_fit_inc in fitness_incomplete:
#                         list_ind_eval.append(ind_fit_inc[1])
#                     n_sim=len(fitness_incomplete)
#                     for _i in range(child_nb):
#                         try:
#                             if _i+optim_param.n_ind not in list_ind_eval:
#                                 fit_list.insert(_i,(0,_i))
#                         except:
#                             fit_list.append((0,_i))
#                     print_('Number of individuals evaluated: ' + str(n_sim) +'\n\n',mp)
#                     sim_results.cancel()
#                 except:
#                     print_('\n\nWarning : error in simulation',mp)
#                     print_("aborting remaining computations",mp)
#                     fit_list.sort(reverse=False, key=lambda col: col[1])
#                     fitness_incomplete = copy.deepcopy(fit_list)
#                     list_ind_eval = []
#                     for ind_fit_inc in fitness_incomplete:
#                         list_ind_eval.append(ind_fit_inc[1])
#                     n_sim=len(fitness_incomplete)
#                     for _i in range(child_nb):
#                         try:
#                             if _i+optim_param.n_ind not in list_ind_eval:
#                                 fit_list.insert(_i,(0,_i))
#                         except:
#                             fit_list.append((0,_i))
#                     print_('Number of individuals evaluated: ' + str(n_sim) +'\n\n',mp)
#                     sim_results.cancel()
#
#
#
# #        for _i in range(len(fit_list)):
# #            self.individual[_i].fitness = fit_list[_i][0]
#
# #        if os.name == 'nt': multiprocessing.get_context('spawn')
# #        with multiprocessing.Pool(num_cores) as p:
# #            fit_list = p.map(self.fitness_eval_par, fit_eval_inp, timeout=5)
# #            p.start()
# #            p.join(1)
# #
# #            # If thread is active
# #            if p.is_alive():
# #                print_("foo is running... let's kill it...")
# #                # Terminate foo
# #                p.terminate()
# #                # Cleanup
# #                p.join()
#
#
# #        if __name__ == '__main__':
# #            # Start foo as a process
# #        #    p = multiprocessing.Process(target=simul_flame, name="Foo", args=())
# #            p = multiprocessing.Process(target=simul_flame, args=())
# #
# #            p.start()
# #
# #            # Wait a maximum of 10 seconds for foo
# #            # Usage: join([timeout in seconds])
# #            p.join(10)
# #
# #            # If thread is active
# #            if p.is_alive():
# #                print("foo is running... let's kill it...")
# #
# #                # Terminate foo
# #                p.terminate()
# #                # Cleanup
# #                p.join()
#
#
#
#         fit_list.sort(reverse=False, key=lambda col: col[1])
#
#
#         for _i in range(len(fit_list)):
#             ind = optim_param.n_ind + _i
#             self.individual[ind].fitness = fit_list[_i][0]
#
#         bar.update(child_nb,title)
#         print('\n')
#
#
#
#
#         for i in range(len(conditions_list)):
#             conditions_list[i].composition.gas     = gas
#             conditions_list[i].composition.gas_ref = gas_ref
#         for i in range(len(ref_results_list)):
#             ref_results_list[i].gas = gas
#             ref_results_list[i].f   = f[i]




    def convergence_information(self,gen,optim_param,verbose=0):
        mp = optim_param.main_path

        vec = []
        for p in range(optim_param.n_ind):
            vec.append(self.individual[p].fitness)
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
        for p in range(len(pop_copy.individual)):
            # check if fitness is not nan
#            if pop_copy.individual[p].fitness==pop_copy.individual[p].fitness:
            fit.append(pop_copy.individual[p].fitness)
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
                self.individual[new_ind]=copy.deepcopy(pop_copy.individual[i_prob])
                new_ind += 1 #; i_prob = 0
            else:
                i_prob += 1
        #self.sort_fitness_rev()

    def select_3_rank(self, optim_param):

        pop_copy = copy.deepcopy(self)
        pop_copy.sort_fitness()

        random.seed() ; rank = [] ; proba = [] ; rand_sel_vect=[]

        # build rank vector
        for p in range(len(pop_copy.individual)):
            rank.append(p+1)
        # calculate selection probability vector
        for p in range(len(pop_copy.individual)):
            proba.append(rank[p]/np.sum(rank))
        proba = np.cumsum(proba)

        # build the random vector for the selection
        for p in range(len(pop_copy.individual)):
            rand_sel_vect.append(random.random())
        rand_sel_vect.sort()

        new_ind = 0 ; i_prob =0 ;  proba[-1]=1
        size_pop = optim_param.n_ind
        while new_ind<size_pop:
            if rand_sel_vect[new_ind]<=proba[i_prob]:
                self.individual[new_ind]=copy.deepcopy(pop_copy.individual[i_prob])
                new_ind +=1 #; i_prob = 0
            else:
                i_prob+=1

    def select_4_geomNorm(self,optim_param):


        q = optim_param.selection_options[0]
        pop_copy = copy.deepcopy(self)
        pop_copy.sort_fitness

        random.seed() ; proba = [] ; rand_sel_vect=[]

        # calculate selection probability vector
        for r in range(len(pop_copy.individual)):
            p = (q/(1-(1-q)**len(pop_copy.individual)) )*(1-q)**r
            proba.append(p)
        proba = np.cumsum(proba)/np.sum(proba)
        proba.sort()

        # build the random vector for the selection
        for i in range(len(pop_copy.individual)):
            rand_sel_vect.append(random.random())
        rand_sel_vect.sort()

        new_ind = 0 ; i_prob =0 ; proba[-1]=1
        size_pop = optim_param.n_ind
        while new_ind<size_pop:
            if rand_sel_vect[new_ind]<=proba[i_prob]:
                self.individual[new_ind]=copy.deepcopy(pop_copy.individual[i_prob])
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
        size_react = len(self.individual[0].mech.react.kin)
        size_pop   = optim_param.n_ind
        parent1    = random.randrange(0, size_pop)
        parent2    = random.randrange(0, size_pop)
        child1     = size_pop+created_ind
        child2     = size_pop+created_ind+1
        while parent1 == parent2:
            parent2 = random.randrange(0, size_pop)

        self.individual[child1] = copy.deepcopy(self.individual[parent1])
        self.individual[child2] = copy.deepcopy(self.individual[parent2])

        a = random.randrange(0,size_react)
        for i in range(a,size_react):
            self.individual[child1].mech.react.kin[i] = copy.deepcopy(self.individual[parent2].mech.react.kin[i])
            self.individual[child2].mech.react.kin[i] = copy.deepcopy(self.individual[parent1].mech.react.kin[i])
                

    def Xover_2_multiple(self,conditions_list,optim_param,ref_results_list,created_ind):

        random.seed()
        size_react = len(self.individual[0].mech.react.kin)
        size_pop   = optim_param.n_ind
        parent1    = random.randrange(0, size_pop)
        parent2    = random.randrange(0, size_pop)
        child1     = size_pop+created_ind
        child2     = size_pop+created_ind+1
        while parent1 == parent2:
            parent2 = random.randrange(0, size_pop)

        self.individual[child1] = copy.deepcopy(self.individual[parent1])
        self.individual[child2] = copy.deepcopy(self.individual[parent2])

        for i in range(size_react):
            a = random.randrange(0,2)
            if a ==1:
                self.individual[child1].mech.react.kin[i] = copy.deepcopy(self.individual[parent2].mech.react.kin[i])
                self.individual[child2].mech.react.kin[i] = copy.deepcopy(self.individual[parent1].mech.react.kin[i])


    def Xover_3_arith(self,conditions_list,optim_param,ref_results_list,created_ind):

        mp = optim_param.main_path
        random.seed()
        size_react = len(self.individual[0].mech.react.kin)
        size_pop   = optim_param.n_ind
        p1    = random.randrange(0, size_pop)
        p2    = random.randrange(0, size_pop)
        child1     = size_pop+created_ind
        child2     = size_pop+created_ind+1
        while p1 == p2:
            p2 = random.randrange(0, size_pop)

        self.individual[child1] = copy.deepcopy(self.individual[p1])
        self.individual[child2] = copy.deepcopy(self.individual[p2])


        for r in range(size_react):
            try_r = 0 ; valid1,valid2 = False,False
            while not ((valid1 and valid2) or try_r > 3):
                uncert_r = [u/100 for u in self.individual[0].mech.react.uncert[r]]
                if self.individual[0].mech.react.modif[r]:
                    mix = random.random()
                    if self.individual[0].mech.react.type[r] =='three_body_reaction' \
                    or self.individual[0].mech.react.type[r] =='reaction':
                        nTry = 0; maxRetry = 15
                        while nTry <= maxRetry:
                            val1=[];val2=[]
                            for k in range(3):
                                ref = self.individual[0].mech.react.ref_kin[r][k]
                                if ref!= 0.0:
                                    new_val1  = self.individual[p1].mech.react.kin[r][k]*mix \
                                               +self.individual[p2].mech.react.kin[r][k]*(1-mix)
                                    new_val2  = self.individual[p2].mech.react.kin[r][k]*mix \
                                               +self.individual[p1].mech.react.kin[r][k]*(1-mix)
                                    if optim_param.Arrh_var:
                                        a = ref*(1-uncert_r[k])
                                        b = ref*(1-uncert_r[k])
                                        min_val = min(a,b)
                                        max_val = max(a,b)
                                    else:
                                        f_min = np.minimum(self.individual[0].mech.react.f_min[r],5)    # f<5 to avoid overflow
                                        T_min = self.individual[0].mech.react.f_Tit[0]
                                        T_max = self.individual[0].mech.react.f_Tit[1]
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
                                self.individual[child1].mech.react.kin[r] = list(val1)
                                self.individual[child2].mech.react.kin[r] = list(val2)
                                break
                            else:
                                nTry +=1 ; mix = random.random()

                    if self.individual[0].mech.react.type[r] =='falloff_reaction' \
                    or self.individual[0].mech.react.type[r] =='chemically_activated_reaction'\
                    or self.individual[0].mech.react.type[r] =='chebyshev'\
                    or self.individual[0].mech.react.type[r] =='pdep_arrhenius':
                        mix = random.random()                      
                        nTry = 0; maxRetry = 15
                        while nTry <= maxRetry:
                            success = True                        
                            for j in range(len(self.individual[0].mech.react.kin[r])):
                                val1=[];val2=[]
                                for k in range(3):
                                    ref=self.individual[0].mech.react.ref_kin[r][j][k]
                                    if ref!= 0.0:
                                        new_val1  = self.individual[p1].mech.react.kin[r][j][k]*mix\
                                                   +self.individual[p2].mech.react.kin[r][j][k]*(1-mix)
                                        new_val2  = self.individual[p2].mech.react.kin[r][j][k]*mix \
                                                   +self.individual[p1].mech.react.kin[r][j][k]*(1-mix)
                                        if optim_param.Arrh_var or self.individual[0].mech.react.type[r] =='chebyshev':
                                            a = ref*(1-uncert_r[k])
                                            b = ref*(1-uncert_r[k])
                                            min_val = min(a,b)
                                            max_val = max(a,b)
                                        else:
                                            f_min = np.minimum(self.individual[0].mech.react.f_min[r],5)
                                            T_min = self.individual[0].mech.react.f_Tit[0]
                                            T_max = self.individual[0].mech.react.f_Tit[1]
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
                                for j in range(len(self.individual[0].mech.react.kin[r])):
                                    for k in range(3):
                                        if ref!= 0.0:
                                            new_val1  = self.individual[p1].mech.react.kin[r][j][k]*mix\
                                                       +self.individual[p2].mech.react.kin[r][j][k]*(1-mix)
                                            new_val2  = self.individual[p2].mech.react.kin[r][j][k]*mix \
                                                       +self.individual[p1].mech.react.kin[r][j][k]*(1-mix)
                                            self.individual[child1].mech.react.kin[r][j][k] = new_val1
                                            self.individual[child2].mech.react.kin[r][j][k] = new_val2
                                break
                            else:
                                nTry +=1 ; 
                    valid1 = ot.check_k(self.individual[child1].mech.react,r)
                    valid2 = ot.check_k(self.individual[child2].mech.react,r)
                    try_r += 1
                else:
                    valid1,valid2 = True,True

            if not valid1 or not valid2:
                if not valid1:
                    self.individual[child1].mech.react.kin[r] = copy.deepcopy(self.individual[p1].mech.react.kin[r])
                if not valid2:
                    self.individual[child2].mech.react.kin[r] = copy.deepcopy(self.individual[p2].mech.react.kin[r])



    def Xover_4_heuri(self,conditions_list,optim_param,ref_results_list,created_ind):

        random.seed()
        size_react = len(self.individual[0].mech.react.kin)
        size_pop   = optim_param.n_ind
        p1    = random.randrange(0, size_pop)
        p2    = random.randrange(0, size_pop)
        damping = 0.75
        while p1 == p2:
            p2 = random.randrange(0, size_pop)
        if self.individual[p1].fitness>self.individual[p2].fitness:
            best=p1; worse=p2
        else:
            best=p2; worse=p1
        child1     = size_pop+created_ind
        child2     = size_pop+created_ind+1

        self.individual[child1] = copy.deepcopy(self.individual[worse])
        self.individual[child2] = copy.deepcopy(self.individual[best])

        for r in range(size_react):
            try_r = 0 ; valid1,valid2 = False,False
            while not ((valid1 and valid2) or try_r > 3):
                uncert_r = [u/100 for u in self.individual[0].mech.react.uncert[r]]
                mix = random.random()
                if self.individual[0].mech.react.modif[r]:
                    if self.individual[0].mech.react.type[r] =='three_body_reaction' \
                    or self.individual[0].mech.react.type[r] =='reaction':
                        nTry = 0; maxRetry = 15
                        while nTry <= maxRetry:
                            success = True
                            val1=[]
                            for k in range(3):
                                ref = self.individual[0].mech.react.ref_kin[r][k]
                                if ref!= 0.0:
                                    bestVal  = self.individual[best].mech.react.kin[r][k]
                                    worseVal = self.individual[worse].mech.react.kin[r][k]
                                    new_val  = mix*(damping**try_r)*(bestVal-worseVal)+bestVal
                                    if optim_param.Arrh_var:
                                        a = ref*(1-uncert_r[k])
                                        b = ref*(1-uncert_r[k])
                                        min_val = min(a,b)
                                        max_val = max(a,b)
                                    else:
                                        f_min = np.minimum(self.individual[0].mech.react.f_min[r],5) # f<5 to avoid overflow
                                        T_min = self.individual[0].mech.react.f_Tit[0]
                                        T_max = self.individual[0].mech.react.f_Tit[1]
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
                                self.individual[child1].mech.react.kin[r] = list(val1)
                                break
                            else:
                                nTry +=1 ; mix = random.random()
                    if self.individual[0].mech.react.type[r] =='falloff_reaction' \
                    or self.individual[0].mech.react.type[r] =='chemically_activated_reaction'\
                    or self.individual[0].mech.react.type[r] =='chebyshev'\
                    or self.individual[0].mech.react.type[r] =='pdep_arrhenius':
                        nTry = 0; maxRetry = 15
                        while nTry <= maxRetry:
                            success = True                        
                            for j in range(len(self.individual[0].mech.react.kin[r])):
                                val1=[]
                                for k in range(3):
                                    ref=self.individual[0].mech.react.ref_kin[r][j][k]
                                    if ref!= 0.0:
                                        bestVal  = self.individual[best].mech.react.kin[r][j][k]
                                        worseVal = self.individual[worse].mech.react.kin[r][j][k]
                                        new_val  = mix*(damping**try_r)*(bestVal-worseVal)+bestVal
                                        if optim_param.Arrh_var:
                                            a = ref*(1-uncert_r[k])
                                            b = ref*(1-uncert_r[k])
                                            min_val = min(a,b)
                                            max_val = max(a,b)
                                        else:
                                            f_min = np.minimum(self.individual[0].mech.react.f_min[r],5)
                                            T_min = self.individual[0].mech.react.f_Tit[0]
                                            T_max = self.individual[0].mech.react.f_Tit[1]
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
                                for j in range(len(self.individual[0].mech.react.kin[r])):
                                    for k in range(3):
                                        bestVal  = self.individual[best].mech.react.kin[r][j][k]
                                        worseVal = self.individual[worse].mech.react.kin[r][j][k]
                                        new_val  = mix*(damping**try_r)*(bestVal-worseVal)+bestVal
                                        self.individual[child1].mech.react.kin[r][j][k]=new_val
                                break
                            else:
                                nTry +=1 ; mix = random.random()
                    valid1 = ot.check_k(self.individual[child1].mech.react,r)
                    valid2 = ot.check_k(self.individual[child2].mech.react,r)
                    try_r += 1
                else:
                    valid1,valid2 = True,True
            if not valid1 and valid2:
                if not valid1:
                    self.individual[child1].mech.react.kin[r] = copy.deepcopy(self.individual[p1].mech.react.kin[r])
                if not valid2:
                    self.individual[child2].mech.react.kin[r] = copy.deepcopy(self.individual[p2].mech.react.kin[r])


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
        size_react = len(self.individual[0].mech.react.kin)
        size_pop   = int(optim_param.n_ind)
        parent1    = random.randrange(0, size_pop)
        child1     = int(size_pop+created_ind)
        probaMut   = optim_param.mut_intensity
        damping    = 0.75

        self.individual[child1] = copy.deepcopy(self.individual[parent1])

        for r in range(size_react):
            try_r = 0 ; valid = False
            while not (valid or try_r > 10):
                uncert_r = [u/100 for u in self.individual[0].mech.react.uncert[r]]
                if self.individual[0].mech.react.modif[r]:
                    if self.individual[0].mech.react.type[r] =='three_body_reaction' \
                    or self.individual[0].mech.react.type[r] =='reaction':
                        val=[]
                        for k in range(3):                            
                            ref = self.individual[0].mech.react.ref_kin[r][k]
                            rand = random.random()*100
                            if ref==0:
                                rand = probaMut+1
                            if rand < probaMut: # random modif of a kinetic constant
                                if optim_param.Arrh_var:
                                    val.append(random.uniform(ref*(1-uncert_r[k]*(damping**try_r))\
                                                            ,ref*(1+uncert_r[k]*(damping**try_r))))
                                
                                else:
                                    f_min = np.minimum(self.individual[0].mech.react.f_min[r]*(damping**try_r),5) # f<5 to avoid overflow
                                    T_min = self.individual[0].mech.react.f_Tit[0]
                                    T_max = self.individual[0].mech.react.f_Tit[1]
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
                                val.append(self.individual[child1].mech.react.kin[r][k])
                        self.individual[child1].mech.react.kin[r]=list(val)

                    if self.individual[0].mech.react.type[r] =='falloff_reaction' \
                    or self.individual[0].mech.react.type[r] =='chemically_activated_reaction'\
                    or self.individual[0].mech.react.type[r] =='pdep_arrhenius':
                        mod_factor = [] ; j=0
                        # calculate the variation for the three Arrhenius parameters
                        # and then apply the same factors to the other Arrhenius parameters
                        # to keep consistancy on the modifications (see Bertolino et al. Comb and Flame 229 (2021) 111366)
                        for k in range(3):
                            ref = self.individual[0].mech.react.ref_kin[r][j][k]
                            rand = random.random()*100
                            if ref==0:
                                rand = probaMut+1
                            if rand < probaMut: # random modif of a kinetic constant
                                if optim_param.Arrh_var:
                                    val = random.uniform(ref*(1-uncert_r[k]*(damping**try_r))\
                                                            ,ref*(1+uncert_r[k]*(damping**try_r)))
                                    if ref!=0:  mod_factor.append(val/ref)
                                    else:       mod_factor.append(0)
                                else:
                                    f_min = self.individual[0].mech.react.f_min[r]*(damping**try_r)
                                    T_min = self.individual[0].mech.react.f_Tit[0]
                                    T_max = self.individual[0].mech.react.f_Tit[1]
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
                        for j in range(len(self.individual[0].mech.react.kin[r])):
                            for k in range(3):
                                ref = self.individual[0].mech.react.ref_kin[r][j][k]
                                self.individual[child1].mech.react.kin[r][j][k]=ref*mod_factor[k]
                                
                    if self.individual[0].mech.react.type[r] =='chebyshev':
                        for j in range(len(self.individual[0].mech.react.kin[r])):
                            val=[]
                            for k in range(3):
                                ref = self.individual[0].mech.react.ref_kin[r][j][k]
                                rand = random.random()*100
                                if ref==0:
                                    rand = probaMut+1
                                if rand < probaMut: # random modif of a kinetic constant
                                    val.append(random.uniform(ref*(1-uncert_r[k]*(damping**try_r))\
                                                            ,ref*(1+uncert_r[k]*(damping**try_r))))
                                else:
                                    val.append(self.individual[child1].mech.react.kin[r][j][k])
                            self.individual[child1].mech.react.kin[r][j]=list(val)
                                
                # Validation of the new reaction rate                                
                valid = ot.check_k(self.individual[child1].mech.react,r)
                try_r += 1

            if not valid:
                self.individual[child1].mech.react.kin[r] = copy.deepcopy(self.individual[parent1].mech.react.kin[r])



    def mut_2_nonUnif(self,conditions_list,optim_param,ref_results_list,created_ind,\
                      gen,option):

        random.seed()
        size_react = len(self.individual[0].mech.react.kin)
        size_pop   = int(optim_param.n_ind)
        parent1    = random.randrange(0, size_pop)
        child1     = int(size_pop+created_ind)
        probaMut   = optim_param.mut_intensity
        ratio      = gen/optim_param.n_gen
        shape      = option
        damping    = 0.75

        self.individual[child1] = copy.deepcopy(self.individual[parent1])

        for r in range(size_react):
            try_r = 0 ; valid = False
            while not (valid or try_r > 10):
                rand = random.random()*100
                uncert_r = [u/100 for u in self.individual[0].mech.react.uncert[r]]
                if rand < probaMut:
                    if self.individual[0].mech.react.modif[r]:
                        if self.individual[0].mech.react.type[r] =='three_body_reaction' \
                        or self.individual[0].mech.react.type[r] =='reaction':
                            for k in range(3):
                                ref = self.individual[0].mech.react.ref_kin[r][k]
                                val = self.individual[child1].mech.react.kin[r][k]
                                if optim_param.Arrh_var:
                                    if ref>0:
                                        min_val = ref*(1-uncert_r[k]*(damping**try_r))
                                        max_val = ref*(1+uncert_r[k]*(damping**try_r))
                                    else:
                                        max_val = ref*(1-uncert_r[k]*(damping**try_r))
                                        min_val = ref*(1+uncert_r[k]*(damping**try_r))
                                else:
                                    f_min = np.minimum(self.individual[0].mech.react.f_min[r],5)*(damping**try_r)
                                    T_min = self.individual[0].mech.react.f_Tit[0]
                                    T_max = self.individual[0].mech.react.f_Tit[1]
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
                                    self.individual[child1].mech.react.kin[r][k]+=change
                                elif rand_dir == 1:
                                    change=(val-min_val)*(random.random()*(1-ratio))**shape
                                    self.individual[child1].mech.react.kin[r][k]-=change
                                    
                        elif self.individual[0].mech.react.type[r] =='falloff_reaction' \
                        or self.individual[0].mech.react.type[r] =='chemically_activated_reaction'\
                        or self.individual[0].mech.react.type[r] =='pdep_arrhenius':
                            mod_factor = [] ; j=0     
                            # calculate the variation for the three Arrhenius parameters
                            # and then apply the same factors to the other Arrhenius parameters
                            # to keep consistancy on the modifications (see Bertolino et al. Comb and Flame 229 (2021) 111366)
                            for k in range(3):
                                ref = self.individual[0].mech.react.ref_kin[r][j][k]
                                val = self.individual[child1].mech.react.kin[r][j][k]
                                if optim_param.Arrh_var:
                                    if ref>0:
                                        min_val = ref*(1-uncert_r[k]*(damping**try_r))
                                        max_val = ref*(1+uncert_r[k]*(damping**try_r))
                                    else:
                                        max_val = ref*(1-uncert_r[k]*(damping**try_r))
                                        min_val = ref*(1+uncert_r[k]*(damping**try_r))
                                else:
                                    f_min = np.minimum(self.individual[0].mech.react.f_min[r],5)*(damping**try_r)
                                    T_min = self.individual[0].mech.react.f_Tit[0]
                                    T_max = self.individual[0].mech.react.f_Tit[1]
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
                                        new_val = self.individual[child1].mech.react.kin[r][j][k]+change
                                        mod_factor.append(new_val/ref)
                                    elif rand_dir == 1:
                                        change=(val-min_val)*(random.random()*(1-ratio))**shape
                                        new_val = self.individual[child1].mech.react.kin[r][j][k]-change
                                        mod_factor.append(new_val/ref)
                                else:
                                    mod_factor.append(0)
                            # Introduction of the new Arrhenius parameters
                            for j in range(len(self.individual[0].mech.react.kin[r])):
                                for k in range(3):
                                    ref = self.individual[0].mech.react.ref_kin[r][j][k]
                                    self.individual[child1].mech.react.kin[r][j][k]=ref*mod_factor[k]
                                    
                        elif self.individual[0].mech.react.type[r] =='chebyshev':
                            for j in range(len(self.individual[0].mech.react.kin[r])):
                                for k in range(3):
                                    ref = self.individual[0].mech.react.ref_kin[r][j][k]
                                    val = self.individual[child1].mech.react.kin[r][j][k]
                                    if ref>0:
                                        min_val = ref*(1-uncert_r[k]*(damping**try_r))
                                        max_val = ref*(1+uncert_r[k]*(damping**try_r))
                                    else:
                                        max_val = ref*(1-uncert_r[k]*(damping**try_r))
                                        min_val = ref*(1+uncert_r[k]*(damping**try_r))
                                    rand_dir = random.randrange(0,2)
                                    if rand_dir == 0: # random modif of a kinetic constant
                                        change=(max_val-val)*(random.random()*(1-ratio))**shape
                                        self.individual[child1].mech.react.kin[r][j][k]+=change
                                    elif rand_dir == 1:
                                        change=(val-min_val)*(random.random()*(1-ratio))**shape
                                        self.individual[child1].mech.react.kin[r][j][k]-=change
                                    
                # Validation of the new reaction rate
                valid = ot.check_k(self.individual[child1].mech.react,r)
                try_r += 1

            if not valid:
                self.individual[child1].mech.react.kin[r] = copy.deepcopy(self.individual[parent1].mech.react.kin[r])


#        self.individual[child1].fitness_eval(conditions_list,optim_param,ref_results_list)



    def mut_3_bound(self,conditions_list,optim_param,ref_results_list,created_ind,opt):

        random.seed()
        size_react = len(self.individual[0].mech.react.kin)
        size_pop   = int(optim_param.n_ind)
        parent1    = random.randrange(0, size_pop)
        child1     = int(size_pop+created_ind)
        probaMut   = optim_param.mut_intensity
        damping    = 0.75
        self.individual[child1] = copy.deepcopy(self.individual[parent1])

        for r in range(size_react):
            try_r = 0 ; valid = False 
            rand = random.random()*100
            if rand < probaMut:      
                uncert_r = [u/100 for u in self.individual[0].mech.react.uncert[r]]
                while not (valid or try_r > 10):
                    if self.individual[0].mech.react.modif[r]:
                        if self.individual[0].mech.react.type[r] =='three_body_reaction' \
                        or self.individual[0].mech.react.type[r] =='reaction':
                            for k in range(3):
                                ref = self.individual[0].mech.react.ref_kin[r][k]
                                rand_dir = random.randrange(0,2)           
                                if optim_param.Arrh_var:
                                    min_val = ref*(1-uncert_r[k]*(damping**try_r))
                                    max_val = ref*(1-uncert_r[k]*(damping**try_r))
                                else:
                                    f_min = self.individual[0].mech.react.f_min[r]*(damping**try_r)
                                    T_min = self.individual[0].mech.react.f_Tit[0]
                                    T_max = self.individual[0].mech.react.f_Tit[1]
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
                                    self.individual[child1].mech.react.kin[r][k] = min_val
                                elif rand_dir == 1:
                                    self.individual[child1].mech.react.kin[r][k] = max_val
                        elif self.individual[0].mech.react.type[r] =='falloff_reaction' \
                        or   self.individual[0].mech.react.type[r] =='chemically_activated_reaction'\
                        or   self.individual[0].mech.react.type[r] =='pdep_arrhenius':
                            mod_factor = [] ; j=0     
                            # calculate the variation for the three first Arrhenius parameters
                            # and then apply the same factors to the other Arrhenius parameters
                            # to keep consistancy on the modifications (see Bertolino et al. Comb and Flame 229 (2021) 111366)
                            for k in range(3):
                                ref = self.individual[0].mech.react.ref_kin[r][j][k]
                                rand_dir = random.randrange(0,2)
                                if optim_param.Arrh_var:
                                    min_val = ref*(1-uncert_r[k]*(damping**try_r))
                                    max_val = ref*(1-uncert_r[k]*(damping**try_r))
                                else:
                                    f_min = self.individual[0].mech.react.f_min[r]*(damping**try_r)
                                    T_min = self.individual[0].mech.react.f_Tit[0]
                                    T_max = self.individual[0].mech.react.f_Tit[1]
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
                            for j in range(len(self.individual[0].mech.react.kin[r])):
                                for k in range(3):
                                    ref = self.individual[0].mech.react.ref_kin[r][j][k]
                                    self.individual[child1].mech.react.kin[r][j][k]=ref*mod_factor[k]
                                    
                        elif self.individual[0].mech.react.type[r] =='chebyshev':
                            for j in range(len(self.individual[0].mech.react.kin[r])):
                                for k in range(3):
                                    ref = self.individual[0].mech.react.ref_kin[r][j][k]
                                    rand_dir = random.randrange(0,2)
                                    min_val = ref*(1-uncert_r[k]*(damping**try_r))
                                    max_val = ref*(1-uncert_r[k]*(damping**try_r))
                                if rand_dir == 0: # random modif of a kinetic constant to the boundary value        
                                    self.individual[child1].mech.react.kin[r][j][k] = min_val
                                elif rand_dir == 1:
                                    self.individual[child1].mech.react.kin[r][j][k] = max_val
                
                    # Validation of the new reaction rate
                    valid = ot.check_k(self.individual[child1].mech.react,r)
                    try_r += 1

            if not valid:
                self.individual[child1].mech.react.kin[r] = copy.deepcopy(self.individual[parent1].mech.react.kin[r])

#        self.individual[child1].fitness_eval(conditions_list,optim_param,ref_results_list)





#%% Others

    def popBestInd(self):
        index=0
        max_fit  = 0
        for i in range(len(self.individual)):
            fit = self.individual[i].fitness
            if fit > max_fit:
                index = i
                max_fit = fit
        bestIndCoeff = self.individual[index].value
        bestIndFit = self.individual[index].fitness

        return bestIndCoeff, bestIndFit

    def sort_fitness(self):
        """ Sort population according to chromosome fitness """
        self.individual =sorted(self.individual, key=operator.attrgetter('fitness'), reverse=False)

    def sort_fitness_rev(self):
        """ Sort population according to chromosome fitness """
        self.individual =sorted(self.individual, key=operator.attrgetter('fitness'), reverse=True)


#%% Functions

