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

import numpy as np
import os as os
import sys as sys
import time as timer
import cantera as ct
import copy
import platform


#import re
#from operator import xor


class State_var :
    def __init__(self,T=300,P=10325):
        self.T  = T
        self.P  = P
        self.T2 = T


class Clock:
    def __init__ (self, title = 'execution time'):
        self.title = title
        self.start_time = 0
        self.stop_time   = 0

    def start(self):
        self.start_time = timer.time()

    def stop(self):
        self.stop_time = timer.time()

    def display(self):
        difftime = self.stop_time - self.start_time
        difftuple = timer.gmtime(difftime)
        txt_disp = self.title +' computation time: %i h %i min %i sec' \
        %(difftuple.tm_hour, difftuple.tm_min, difftuple.tm_sec)
        length_txt = len(txt_disp)+2
        print('\n\n' + '=' *(length_txt) + '\n ' + txt_disp)
        print('=' *length_txt + '\n\n')

    def display_(self,mp):
        difftime = self.stop_time - self.start_time
        difftuple = timer.gmtime(difftime)
        txt_disp = self.title +' computation time: %i h %i min %i sec' \
        %(difftuple.tm_hour, difftuple.tm_min, difftuple.tm_sec)
        length_txt = len(txt_disp)+2
        print_('\n\n' + '=' *(length_txt) + '\n ' + txt_disp,mp)
        print_('=' *length_txt + '\n\n',mp)

    def text(self):
        difftime = self.stop_time - self.start_time
        difftuple = timer.gmtime(difftime)
        txt_disp = self.title +' computation time: %i h %i min %i sec' \
        %(difftuple.tm_hour, difftuple.tm_min, difftuple.tm_sec) + '\n'
        return txt_disp




class Composition :
    def __init__(self,gas=False,fuel='',oxidant='',diluent='',         \
               phi='',diluent_ratio=''):
        self.gas         = gas
        self.gas_ref     = gas
        self.fuel        = fuel
        self.oxidant     = oxidant
        self.diluent     = diluent
        self.phi         = phi
        self.diluent_ratio = diluent_ratio

        # 2nd burner for cflow_flame:
        self.fuel2        = fuel
        self.oxidant2     = oxidant
        self.diluent2     = diluent
        self.phi2         = phi
        self.diluent_ratio2 = diluent_ratio
        self.mixt         = False
        self.X2           = False


    def molarFraction(self, gas,composition_2=False):

        if not gas:  gas = self.gas

        if not composition_2:
            if type(self.phi) is float:    # i.e. if not diffusion flame
                # X_oxidant calculation
                try:
                    X_oxidant = gas.n_atoms(self.fuel, 'C')            \
                                + 1/4 * gas.n_atoms(self.fuel, 'H')    \
                                - 1/2 * gas.n_atoms(self.fuel,'O')
                except:  # non carbonated fuel
                    try:
                        X_oxidant = 1/4 * gas.n_atoms(self.fuel, 'H')      \
                                    - 1/2 * gas.n_atoms(self.fuel,'O')
                    except: # no oxygen in the mechanism
                        X_oxidant = 0

                # X_diluent calculation
                if type(self.diluent_ratio) is float :
                    X_diluent = (self.diluent_ratio/100) * (X_oxidant + self.phi)   \
                                 /(1-self.diluent_ratio/100)
                else:
                    X_diluent = X_oxidant * self.diluent_ratio[1]

                X = self.fuel    + ':' + str(self.phi)   + ',' + \
                    self.oxidant + ':' + str(X_oxidant)  + ',' + \
                    self.diluent + ':' + str(X_diluent)

            else:    # i.e. if diffusion flame
                # X_diluent calculation
                if type(self.diluent_ratio) is float :
                    X_diluent = self.diluent_ratio /(1-self.diluent_ratio)
                else:
                    X_diluent = X_oxidant * self.diluent_ratio[1]

                X = self.fuel    + ':' + str(1)          + ',' + \
                    self.diluent + ':' + str(X_diluent)

        else:
            if type(self.phi2) is float:    # i.e. if not diffusion flame
                try:
                    X_oxidant = gas.n_atoms(self.fuel2, 'C')            \
                                + 1/4 * gas.n_atoms(self.fuel2, 'H')    \
                                - 1/2 * gas.n_atoms(self.fuel2,'O')
                except:  # non carbonated fuel
                    try: # no oxygen in the mechanism
                        X_oxidant = 1/4 * gas.n_atoms(self.fuel2, 'H')      \
                                    - 1/2 * gas.n_atoms(self.fuel2,'O')
                    except: # no oxygen in the mechanism
                        X_oxidant = 0


                if type(self.diluent_ratio2) is float :
                    X_diluent = self.diluent_ratio2 * (X_oxidant + self.phi2) \
                                 / (1-self.diluent_ratio2)
                else:
                    X_diluent = X_oxidant * self.diluent_ratio2[1]

                X = self.fuel2    + ':' + str(self.phi2)  + ',' + \
                    self.oxidant2 + ':' + str(X_oxidant)  + ',' + \
                    self.diluent2 + ':' + str(X_diluent)

            else:    # i.e. if diffusion flame
                try:
                    X_oxidant = gas.n_atoms(self.fuel2, 'C')            \
                                + 1/4 * gas.n_atoms(self.fuel2, 'H')    \
                                - 1/2 * gas.n_atoms(self.fuel2,'O')
                except:  # non carbonated fuel
                    try:
                        X_oxidant = 1/4 * gas.n_atoms(self.fuel2, 'H')      \
                                    - 1/2 * gas.n_atoms(self.fuel2,'O')
                    except: # no oxygen in the mechanism
                        X_oxidant = 0


                if type(self.diluent_ratio2) is float :
                    X_diluent = (X_oxidant*self.diluent_ratio2/100)/(1-self.diluent_ratio2/100)
                else:
                    X_diluent = X_oxidant * self.diluent_ratio2[1]/100

                X = self.oxidant2    + ':' + str(1)       + ',' + \
                    self.diluent2 + ':' + str(X_diluent)

        return X

    # ========================================================================
    # deprecated option

    def molarFraction_fromvalues(self, gas,composition_2=False):
        X = np.zeros(gas.n_species)

        if not composition_2:
            # remove spaces
            self.fuel = self.fuel.replace(" ","")
            self.oxidant = self.oxidant.replace(" ","")
            self.diluent = self.diluent.replace(" ","")

            fuel_sp = self.fuel.split("(")[0]
            fuel_conc = self.fuel.split("(")[1].split(")")[0]
            ox_sp = self.oxidant.split("(")[0]
            ox_conc = self.oxidant.split("(")[1].split(")")[0]
            dil_sp = self.diluent.split("(")[0]
            dil_conc = self.diluent.split("(")[1].split(")")[0]
        else:
            # remove spaces
            self.fuel2 = self.fuel2.replace(" ","")
            self.oxidant2 = self.oxidant2.replace(" ","")
            self.diluent2 = self.diluent2.replace(" ","")

            fuel_sp = self.fuel2.split("(")[0]
            fuel_conc = self.fuel2.split("(")[1].split(")")[0]
            ox_sp = self.oxidant2.split("(")[0]
            ox_conc = self.oxidant2.split("(")[1].split(")")[0]
            dil_sp = self.diluent2.split("(")[0]
            dil_conc = self.diluent2.split("(")[1].split(")")[0]


        for i in range(len(fuel_sp.split("/"))):
            ind_fuel = gas.species_index(fuel_sp.split("/")[i])
            X[ind_fuel] = fuel_conc.split("/")[i]
        for i in range(len(ox_sp.split("/"))):
            ind_ox = gas.species_index(ox_sp.split("/")[i])
            X[ind_ox] = ox_conc.split("/")[i]
        for i in range(len(dil_sp.split("/"))):
            ind_dil = gas.species_index(dil_sp.split("/")[i])
            X[ind_dil] = dil_conc.split("/")[i]
        return X
    # ========================================================================


class conditions :
    def __init__(self,mech,config,fuel,oxidant,diluent,phi,diluent_ratio,T,P, \
                 tspc,error_calculation='points',error_coupling='mean',       \
                 error_type_fit='mean',                                       \
                 T_check=False,ig_check=False,Sl_check=False,K_check=False,   \
                 sp_T=['CO2'],sp_ig=['H'],sp_Sl=['H']):

        if mech:
            # --------------------------------------------------------------------------------
            # interpretation of the new mech
            try:
                # supress console output during the interpretation
                old_stdout = sys.stdout ; old_stderr = sys.stderr
                with open(os.devnull, "w") as devnull: sys.stdout = devnull ; sys.stderr = devnull

                ct.suppress_thermo_warnings()
                try:
                    gas = ct.Solution(mech)
                except:
                    gas = ct.Solution('_kinetic_mech/'+mech)

                # restore console output
                sys.stdout = old_stdout ; sys.stderr = old_stderr
                # --------------------------------------------------------------------------------
            except:
                # restore console output
                sys.stdout = old_stdout ; sys.stderr = old_stderr
                print('\n\n\n ! ! ! ! !\n\nmech: '+mech+' not found\n\n\n\n\n')

        else:    gas=False
        self.mech        = mech
        self.mech_ext    = mech
        self.mech_prev_red = False
        self.config      = config
        self.composition = Composition(gas,fuel,oxidant,diluent,phi,\
                                       diluent_ratio)
        self.state_var   = State_var(T,P)
        self.simul_param = Simul_param()
        self.error_param = Error_param(tspc)
        self.main_path   = ''
        self.main_folder = False
        self.exp_data    = False
        self.import_data = False
        self.conc_unit   = "volumic_concentration" # molar_fraction
        self.mech_prev_red = False



class Red_operator :
    def __init__(self,target_species,optimization=False,n_points=20,          \
                 max_error_sp=30, max_error_T=10, max_error_ig = 10,          \
                 max_error_Sl=10, max_error_K=20, inter_sp_inter = True,                      \
                 eps_init=0.06,delta_eps_init=0.005,eps_max=2.0,              \
                 r_withdraw_intensity = 10,\
#                 sp_interaction_coeffs=[],r_interaction_coeffs=[]\
                 tol_ts = [False, False]):

        self.r_withdraw_intensity = r_withdraw_intensity+1e-5
        self.opt                 = optimization

        if max_error_sp is list: self.max_error_sp = max_error_sp
        else:     self.max_error_sp = [max_error_sp]*len(target_species)

        self.max_error_T        = max_error_T
        self.max_error_ig       = max_error_ig
        self.max_error_Sl       = max_error_Sl
        self.max_error_K        = max_error_K

        if eps_init is list: self.eps_init  = eps_init
        else:                self.eps_init  = [eps_init]*len(target_species)

        if delta_eps_init is list: self.delta_eps_init = delta_eps_init
        else:     self.delta_eps_init = [delta_eps_init]*len(target_species)

        if eps_max is list:  self.eps_max   = eps_max[0:len(target_species)]
        else:                self.eps_max   = [eps_max]*len(target_species)

        self.n_points           = n_points
        self.gas                = 'no_gas_yet'

        self.inter_sp_inter     = inter_sp_inter
        # drg
        self.graph_search       = 'Dijkstra'  # Dijkstra DFS
        self.new_targets_4_DRG_r= copy.deepcopy(target_species)

        # sa
        self.tol_ts             = tol_ts

        # csp
        self.csp_method         = 'Valorani_2006'       # (e.g. 1-1e-5)
        self.tol_csp            = .01       # (e.g. 1-1e-5)
        self.nTime              = 10
        self.TimeResolution     = 1e-7         # (e.g. 1e0)
        self.epsilon_rel        = False
        self.epsilon_abs        = 1e-12
        self.I_fast             = []
        self.I_slow             = []
        self.radIdx             = []
        self.csp_refin_iter     = 0

        self.OIC_sp             = False





class Errors:
    def __init__(self,conditions,ref_results,red_results,red_data,\
                 red_data_meth=False):

        self.qoi_s = \
              self.er_estim_s(conditions,red_data,ref_results,red_results)
        if conditions.error_param.T_check:
            self.qoi_T   = \
                  self.er_estim_T(conditions,ref_results,red_results)
        else: self.qoi_T=0

        if 'free_flame' in conditions.config and conditions.error_param.Sl_check:
            self.qoi_Sl  = \
                  self.er_estim_Sl(conditions,ref_results,red_results)
        else : self.qoi_Sl = 0

        if 'reactor' in conditions.config and conditions.error_param.ig_check:
            try: self.qoi_ig = \
                      self.er_estim_ig(conditions,ref_results,red_results)
            except: self.qoi_ig = 0
        else: self.qoi_ig = 0

        if 'diff_flame' in conditions.config and conditions.error_param.K_check:
            self.qoi_K = \
                  self.er_estim_K(conditions,ref_results,red_results)
        else: self.qoi_K = 0

        if red_data_meth:
            self.under_tol_s,self.under_tol_T,self.under_tol_Sl,\
            self.under_tol_ig,self.under_tol_K,self.under_tol,self.above_tol \
                = self.qoi_test(self.qoi_s,self.qoi_T,self.qoi_Sl,self.qoi_ig,\
                  self.qoi_K,conditions,red_data_meth)


    def er_estim_s(self,conditions,red_data,ref_results,red_results):

        mp = conditions.main_path

        # main variables
        gas_ref           = conditions.composition.gas
#        gas_red           = red_data.gas_loop
        tspc              = red_data.tspc
        n_tspc            = red_data.n_tspc
        error_calculation = conditions.error_param.error_calculation
        error_type        = conditions.error_param.error_coupling
        verbose           = conditions.simul_param.verbose

        pts_scatter       = ref_results.pts_scatter
        if conditions.conc_unit == 'volumic_concentration':
            conc_ref          = ref_results.conc
            conc_red          = red_results.conc
        else:                                      # 'molar_fraction'
            conc_ref          = ref_results.X
            conc_red          = red_results.X
        n_points          = len(pts_scatter)
        QoI = []

        if error_calculation =="points":

            for k in range(n_tspc):

                index = gas_ref.species_index(tspc[k])

                data1 = np.zeros(n_points)
                data2 = np.zeros(n_points)
                for i in range(n_points):
                    data1[i] = conc_ref[i][index]
                    data2[i] = conc_red[i][index]
                if sum(data1)>0:      # absence of data if experimental optimization
                    sumDiff = 0
                    if error_type == "all":
                        sumDiff = []
                        for j in range(len(pts_scatter)):
                            sumDiff.append(abs(data1[j] - data2[j])/np.amax(data1))
                    elif error_type == "mean":
                        sumdata=0
                        for j in range(len(pts_scatter)):
                            if j==0:
                                sumDiff += abs(data1[j] - data2[j])*.5*(pts_scatter[j+1]-pts_scatter[j])
                                sumdata += abs(data1[j])*.5*(pts_scatter[j+1]-pts_scatter[j])
                            elif j==(len(pts_scatter)-1):
                                sumDiff += abs(data1[j] - data2[j])*.5*(pts_scatter[j]-pts_scatter[j-1])
                                sumdata += abs(data1[j])*.5*(pts_scatter[j]-pts_scatter[j-1])
                            else:
                                sumDiff += abs(data1[j] - data2[j])*.5*(pts_scatter[j+1]-pts_scatter[j-1])
                                sumdata += abs(data1[j])*.5*(pts_scatter[j+1]-pts_scatter[j-1])
                                          #/np.amax(data1)
                        #sumDiff = sumDiff/len(pts_scatter)
                        #sumDiff = sumDiff/np.amax(data1)
                        sumDiff/=sumdata
                    elif error_type == "max":
                        diff = abs(np.array([data1]) - np.array([data2]))
                        sumDiff = np.amax(diff)/np.amax(data1)
                    QoI.append(sumDiff)
                else:
                    QoI.append(False)


        elif error_calculation == "QoI":

            for k in range(len(tspc)):
                index = gas_ref.species_index(tspc[k])
                data1 = []
                data2 = []
                for i in range(n_points-1):
                    data1.append(conc_ref[i][index])
                    data2.append(conc_red[i][index])

                if sum(data1)>0:     # absence of data if experimental optimization

                    data1=np.array(data1) ; data2=np.array(data2)
                    ref_var_grad = []

                    # Analysis of the curve shape
                    for i in range(len(data1)-1):
                        grad = (data1[i+1]-data1[i])/(pts_scatter[i+1]-pts_scatter[i])
                        if abs(grad) < 0.1:
                            grad = 0
                        ref_var_grad.append(grad)
                    max_grad = max(ref_var_grad)
                    min_grad  = min(ref_var_grad)

                    if min_grad < 0 and max_grad > 0 : # bell curve
                        thd=1.2 #threshold
                        if data1[0]>data1[-1]:
                            if (np.max(data1)-data1[-1])/(data1[0]-data1[-1])>thd:
                                shape='bell'
                            elif (data1[0]-np.min(data1))/(data1[0]-data1[-1])>thd:
                                shape='bell (inv)'
                            else:
                                shape = 'downward'
                        elif data1[0]<data1[-1]:
                            if (np.max(data1)-data1[0])/(data1[-1]-data1[0])>thd:
                                shape='bell'
                            elif (data1[-1]-np.min(data1))/(data1[-1]-data1[0])>thd:
                                shape='bell (inv)'
                            else:
                                shape='upward'
                        else :
                            if np.max(data1)>data1[0]: shape = 'bell'
                            else: shape = 'bell (inv)'
                    elif max_grad <= 0: shape='downward'  # downward curve
                    else: shape='upward'                  # upward curve

                    if verbose>5:print_(tspc[k]+' curve shape: '+shape,mp)
                    QoI_S= self.qoi_computation(pts_scatter, data1, data2, shape)

                    if error_type == "all":
                        QoI.append(QoI_S)
                    elif error_type=="mean":
                        QoI.append(sum(QoI_S)/len(QoI_S))
                    elif error_type == "max":
                        QoI.append(max(QoI_S))
                else:
                    QoI.append(False)
        else:
            QoI = False

        # Display QoI
        if verbose >=6 and error_type != 'all' and False not in QoI:
            QoI_info = '      Errors spec (' + error_type + '): '
            for i in range(n_tspc):
                QoI_info = QoI_info + str(tspc[i]) + ':' + '%.1f' %(QoI[i]*100) + '%  '
            print_('\n\n'+QoI_info,mp)

        return QoI


    def er_estim_T(self,conditions, ref_results, red_results):
        mp = conditions.main_path

        # main variables
        error_calculation = conditions.error_param.error_calculation
        error_type        = conditions.error_param.error_coupling
        verbose           = conditions.simul_param.verbose


        pts_scatter       = ref_results.pts_scatter
        T_ref             = ref_results.T
        T_red             = red_results.T

        n_points = len(pts_scatter)

        T_data = True
        if 'False' in T_ref:     T_data = False
        elif np.min(T_ref)<273:  T_data = False

        if error_calculation=="points" and T_data:

            # ========   2 - Temperature error   ==========
            sumDiff = 0
            if error_type == "all" :
                sumDiff = []
                for j in range(len(pts_scatter)):
                    sumDiff.append(abs(T_ref[j] - T_red[j])/np.amax(T_ref))
            elif error_type == "mean":
                sumdata=0
                for j in range(len(pts_scatter)):
                    if j==0:
                        sumDiff += abs(T_ref[j] - T_red[j])*.5*(pts_scatter[j+1]-pts_scatter[j])
                        sumdata += abs(T_ref[j])*.5*(pts_scatter[j+1]-pts_scatter[j])
                    elif j==(len(pts_scatter)-1):
                        sumDiff += abs(T_ref[j] - T_red[j])*.5*(pts_scatter[j]-pts_scatter[j-1])
                        sumdata += abs(T_ref[j])*.5*(pts_scatter[j]-pts_scatter[j-1])
                    else:
                        sumDiff += abs(T_ref[j] - T_red[j])*.5*(pts_scatter[j+1]-pts_scatter[j-1])
                        sumdata += abs(T_ref[j])*.5*(pts_scatter[j+1]-pts_scatter[j-1])
                                  #/np.amax(data1)
                #sumDiff = sumDiff/len(pts_scatter)
                #sumDiff = sumDiff/np.amax(data1)
                sumDiff/=sumdata
#                for i in range(n_points):
#                    sumDiff += abs(T_ref[i] - T_red[i])
#                sumDiff = sumDiff/n_points
#                sumDiff = sumDiff/np.amax(T_ref)
            elif error_type == "max":
                diff = abs(np.array([T_ref]) - np.array([T_red]))
                sumDiff = np.amax(diff)/np.amax(T_ref)
            QoI = sumDiff


        elif error_calculation=="QoI" and T_data:
            # ========   2 - Temperature QoI   ==========
            T_variation = (max(T_ref)-min(T_ref))/min(T_ref)
            if T_variation > 0.05:
                if T_variation < 0.1:
                    QoI_min = abs(min(T_ref) - min(T_red))/min(T_ref)
                    QoI_max = abs(max(T_ref) - max(T_red))/max(T_ref)
                    QoI_T = [QoI_min, QoI_max]
                else :
                    QoI_T= self.qoi_computation(pts_scatter, T_ref, T_red)
                if error_type == "all":
                    QoI = QoI_T
                elif error_type == "mean":
                    QoI = sum(QoI_T)/len(QoI_T)
                elif error_type == "max":
                    QoI = max(QoI_T)
            else :
                QoI = 0
        else:
            QoI = False


        # Display QoI
        if verbose >=6 and error_type != 'all' and QoI:
            QoI_info = '      Errors T (' + error_type + '):    '\
                        +'%0.1f' %(QoI*100) + '%'
            print_(QoI_info,mp)

        return QoI


    def er_estim_Sl(self,conditions, ref_results, red_results):
        # main variables
        mp = conditions.main_path
        verbose           = conditions.simul_param.verbose
        if ref_results.Sl>1e-5 :
            QoI = abs(ref_results.Sl - red_results.Sl)/ref_results.Sl
            # Display QoI
            if verbose >=6 :
                print_('      Error Sl:' + '%0.1f' %(QoI*100) + '%',mp)
        else:
            QoI = False

        return QoI



    def er_estim_ig(self,conditions, ref_results, red_results):
        # main variables
        mp = conditions.main_path
        verbose           = conditions.simul_param.verbose

        if ref_results.ign_time_hr == False or red_results.ign_time_hr == False:
            ref_results.ign_time = ref_results.ign_time_sp
            red_results.ign_time = red_results.ign_time_sp
        else:
            ref_results.ign_time = ref_results.ign_time_hr
            red_results.ign_time = red_results.ign_time_hr

        diff=abs(ref_results.ign_time-red_results.ign_time)/ref_results.ign_time

        # Display QoI
        if verbose >=6 :
            print_('      Error ignition time: ' + '%0.1f' %(diff*100) +'% \n ',mp)
        return diff


    def er_estim_K(self,conditions, ref_results, red_results):
        # main variables
        mp = conditions.main_path
        verbose           = conditions.simul_param.verbose

        diff=abs(ref_results.K_ext-red_results.K_ext)/ref_results.K_ext

        # Display QoI
        if verbose >=6 :
            print_('      Error extinction : ' + '%0.1f' %(diff*100) +'% \n ',mp)
        return diff


    def qoi_computation(self,pts_scatter, ref_var, red_var, curve_type='upward'):
        if curve_type == "upward":

    #      +:
    #      +:             Ipr      .-:/+++///////.   <-- Ip
    #      +:           <---->`/o+
    #      +:                /s:`                                   <-- I75p (deprecated)
    #      +:              `so`
    #      +:              so
    #      +:             os
    #      +:           `so
    #      +:         `:y:  <-- I25p
    #      +:       :+o+`
    #      +o/////:-.`
    #      -/:::::::::::::::::::::::::::::::::::-



            Dmin = [np.amin( ref_var), np.amin( red_var)]
            Dmax = [np.amax( ref_var), np.amax( red_var)]


            ## QoI : time at 25% of max data
            D25 = Dmin+(np.array(Dmax)-np.array(Dmin))*25/100
            D25_ref, ind25_ref = searchNearest( ref_var, D25[0], 1)
            D25_red, ind25_red = searchNearest( red_var, D25[1], 1)

            I25p = abs( pts_scatter[ind25_red]- pts_scatter[ind25_ref])/ pts_scatter[ind25_ref]


#            ## QoI : time at 75% of max data
            D75 = Dmin+(np.array(Dmax)-np.array(Dmin))*75/100
            D75_ref, ind75_ref = searchNearest( ref_var, D75[0], ind25_ref)
            D75_red, ind75_red = searchNearest( red_var, D75[1], ind25_red)
#            I75p = abs( pts_scatter[ind75_red]- pts_scatter[ind75_ref])/ pts_scatter[ind75_ref]

            ## QoI : maximum value comparison
            Ipa = abs(Dmax[1]-Dmax[0])/Dmax[0]

            ## QoI : diff time 90% / 5%  (allow to catch the gradient of concentration)
            D90 = Dmin+(np.array(Dmax)-np.array(Dmin))*90/100
            D05 = Dmin+(np.array(Dmax)-np.array(Dmin))*5/100
            D90_ref, ind90_ref = searchNearest( ref_var, D90[0], ind75_ref)
            D05_ref, ind05_ref = searchNearest( ref_var, D05[0], 1)
            D90_red, ind90_red = searchNearest( red_var, D90[1], ind75_red)
            D05_red, ind05_red = searchNearest( red_var, D05[1], 1)
            D_ref =  pts_scatter[ind90_ref] -  pts_scatter[ind05_ref]
            D_red =  pts_scatter[ind90_red] -  pts_scatter[ind05_red]
            if D_ref != 0:
                rate_ref = (D90_ref-D05_ref) / (pts_scatter[ind90_ref]-pts_scatter[ind05_ref])
                rate_red = (D90_red-D05_red) / (pts_scatter[ind90_red]-pts_scatter[ind05_red])
                Ipr = 1-abs(rate_red/rate_ref)
                #Ipr = abs((D_ref-D_red)/D_ref)
            elif D_ref == D_red:
                Ipr = 0
            else:
                Ipr = 1.
            dt90 =  pts_scatter[ind90_ref+1]- pts_scatter[ind90_ref-1]
            dt05 =  pts_scatter[ind05_ref+1]- pts_scatter[ind05_ref-1]
            if dt90>0.05*(D_ref) and dt05>0.05*(D_ref) : # if mesh is not sufficiently resolved
                Ipr = np.mean([I25p, Ipa])

            QoI = [I25p, Ipr, Ipa]


        # ========================================================================
        elif curve_type == "downward":

    #              o.
    #              o.
    #              oo//////+o/`                        ___
    #              o.       `:**   <-- D75c             \
    #              o.         `so                       \
    #              o.           h-                      \
    #              o.           -h                      \  --> Dc
    #              o.            os`                    \
    #              o.       <---->+**                   \                        <-- R25c         \ (deprecated)
    #              o.          D   -y                   \
    #              o.                `-:++///::::::    ___
    #              o.
    #              :::::::::::::::::::::::::::::::::::::::::::::


            Dmin = [np.amin( ref_var), np.amin( red_var)]
            Dmax = [np.amax( ref_var), np.amax( red_var)]


            # R75c : time at 75% of max data
            D75 = Dmin+ (np.array(Dmax)-np.array(Dmin))*75/100
            D75_ref, ind75_ref = searchNearest( ref_var, D75[0], 1)
            D75_red, ind75_red = searchNearest( red_var, D75[1], 1)

            ## R25c : time at 25% of max data
            D25 = Dmin+(np.array(Dmax)-np.array(Dmin))*25/100
            D25_ref, ind25_ref = searchNearest( ref_var, D25[0], ind75_ref)
            D25_red, ind25_red = searchNearest( red_var, D25[1], ind75_red)
            if  pts_scatter[ind25_ref]!=0:
                D25c = abs( pts_scatter[ind25_red]- pts_scatter[ind25_ref])/ pts_scatter[ind25_ref]
            elif  pts_scatter[ind25_ref] ==  pts_scatter[ind25_red] :
                D25c = 0
            else:
                D25c = 1


            ## consumption amplitude
            Dc_ref = Dmax[0]-Dmin[0]
            Dc_red = Dmax[1]-Dmin[1]
            Dca = abs((Dc_ref-Dc_red)/Dc_ref)*2


            ## D : diff time 95% / 5%  (allow to catch the gradient of concentration)
            D95 = Dmin+(np.array(Dmax)-np.array(Dmin))*95/100
            D05 = Dmin+(np.array(Dmax)-np.array(Dmin))*5/100
            D95_ref, ind95_ref = searchNearest( ref_var, D95[0], 1)
            D05_ref, ind05_ref = searchNearest( red_var, D05[0], ind75_ref)
            D95_red, ind95_red = searchNearest( ref_var, D95[1], 1)
            D05_red, ind05_red = searchNearest( red_var, D05[1], ind75_ref)
            D_ref = ( pts_scatter[ind95_ref]   -  pts_scatter[ind05_ref])
            D_red = ( pts_scatter[ind95_red] -  pts_scatter[ind05_red])
            if D_ref != 0:
                rate_ref = (D05_ref-D95_ref) / (pts_scatter[ind05_ref]-pts_scatter[ind95_ref])
                rate_red = (D05_red-D95_red) / (pts_scatter[ind05_red]-pts_scatter[ind95_red])
                Ipr = 1-abs(rate_red/rate_ref)
            if D_ref != 0:
                Dcr = abs((D_ref-D_red)/D_ref)
            elif D_ref == D_red:
                Dcr = 0
            else:
                Dcr = 1.
            dt95 =  pts_scatter[ind95_ref+1]- pts_scatter[ind95_ref-1]
            dt05 =  pts_scatter[ind05_ref+1]- pts_scatter[ind05_ref-1]
            if dt95>0.05*(D_ref) and dt05>0.05*(D_ref): # if mesh is not sufficiently resolved
                Dcr = np.mean([D25c, Dca])


            QoI = [D25c, Dcr,  Dca]


        # ========================================================================
        elif curve_type == "bell":

    #        o
    #        y                     .+++++-   <--- Rmaxp + RM_mf   ____
    #        y                    sy.   `/h.                       \
    #        y                   y+       .d`                      \
    #        y                  :h         /o                      \
    #        y                 .d`         `m                      \  --> D50_mf
    #        y                 h-           s/   <--- B75d         \
    #        y       R25p --> /y            `yo-                   \
    #        y               -d`              .-::::::://-        ____
    #        y              :h.
    #        h:==========//+o+`
    #        ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::`


            Dmin = [np.amin(ref_var), np.amin(red_var)]
            Dinit = [ ref_var[0],  red_var[0]]
            Dend = [ ref_var[-1],  red_var[-1]]
            Dmax = [np.amax( ref_var), np.amax( red_var)]


            # search values and time of max
            Dmax_ref, indmax_ref = searchNearest(ref_var, Dmax[0])
            Dmax_red, indmax_red = searchNearest(red_var, Dmax[1])
            if indmax_ref==0: indmax_ref=1 # 0 exception
            # search value and time of 25% of increase
            D25i = np.array(Dinit)+(np.array(Dmax)-np.array(Dinit))*25/100
            D25_ref, ind25_ref = searchNearest(ref_var, D25i[0], 0, indmax_ref)
            try:
                D25_red, ind25_red = searchNearest(red_var, D25i[1], 0, indmax_red+1)
            except:
                D25_red, ind25_red = searchNearest(red_var, D25i[1], 0, indmax_red+1)
            if ind25_ref==0: ind25_ref=1 # 0 exception


            # B25p : time to reach 25% of max value (increasing part)
            B25p = abs( pts_scatter[ind25_red]- pts_scatter[ind25_ref])/ pts_scatter[ind25_ref]


            ## Bpr : rate of increase between 25% and the pic
            rate_red = (Dmax_red-D25_red)/(pts_scatter[indmax_red]-pts_scatter[ind25_red])
            rate_ref = (Dmax_ref-D25_ref)/(pts_scatter[indmax_ref]-pts_scatter[ind25_ref])
            Bpr = 1-abs(rate_red/rate_ref)
            Bm = abs(Dmax[1]-Dmax[0])/Dmax[0]


            ## R25f : time at 25% of the diff between the concentration pic and final value
            D25 =  np.array(Dend)+(np.array(Dmax)-np.array(Dend))*25/100
            D25f_ref, ind25f_ref = searchNearest( ref_var, D25[0], indmax_ref)
            D25f_red, ind25f_red = searchNearest( red_var, D25[1], indmax_red)

            rate_red = (D25f_red-Dmax_red)/(pts_scatter[ind25f_red]-pts_scatter[indmax_red])
            rate_ref = (D25f_ref-Dmax_ref) /(pts_scatter[ind25f_ref]-pts_scatter[indmax_ref])
            B75d = 1-abs(rate_red/rate_ref)


            ## D50_mf_ref : Diff between the ref maximum and final value of each case
            D50_mf_ref = Dmax[0]-Dend[0]
            D50_mf_red = Dmax[0]-Dend[1]
            Bca = abs((D50_mf_ref-D50_mf_red)/D50_mf_ref)


            QoI = [B25p, Bpr, Bm, B75d, Bca]

        # ========================================================================
        elif curve_type == "bell (inv)":


    #        h:==========//+o+`
    #        y              :h.
    #        y               -d`              .-::::::://-        ____
    #        y       R25p --> /y            `yo-                   \
    #        y                 h-           s/   <--- R25f         \
    #        y                 .d`         `m                      \  --> D50_mf
    #        y                  :h         /o                      \
    #        y                   y+       .d`                      \
    #        y                    sy.   `/h.                       \
    #        y                     .+++++-   <--- Rminp + RM_mf   ____
    #        o
    #        ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::`




            Dmin = [np.amin(ref_var), np.amin( red_var)]
            Dinit = [ref_var[0],  red_var[0]]
            Dend = [ ref_var[-1],  red_var[-1]]
            Dmax = [np.amax( ref_var), np.amax( red_var)]


            ## search time and value of the min
            Dmin_ref, indmin_ref = searchNearest(ref_var, Dmin[0])
            Dmin_red, indmin_red = searchNearest(red_var, Dmin[1])
            if indmin_ref==0: indmin_ref=1 # 0 exception
            # search time and value of 25% decrease
            D25i = Dinit-(np.array(Dinit)-np.array(Dmin))*25/100
            D25_ref, ind25_ref = searchNearest(ref_var, D25i[0], 0, indmin_ref)
            D25_red, ind25_red = searchNearest(red_var, D25i[1], 0, indmin_red)
            if ind25_ref==0: ind25_ref=1 # 0 exception

            #Rminp
            Bi_cr = abs(((pts_scatter[ind25_red]-pts_scatter[indmin_red])-(pts_scatter[ind25_ref]-pts_scatter[indmin_ref]))\
                        /(pts_scatter[ind25_ref]-pts_scatter[indmin_ref]))

            ## RM_mf : value of the molar fraction at the pic
            Bi_m = abs(Dmin[1]-Dmin[0])/Dmin[0]


            # R25p : time at25% of max value (increasing part)
            D25i = Dinit-(np.array(Dinit)-np.array(Dmin))*25/100
            D25_ref, ind25_ref = searchNearest(ref_var, D25i[0], 0, indmin_ref)
            D25_red, ind25_red = searchNearest(red_var, D25i[1], 0, indmin_red)
            if ind25_ref==0: ind25_ref=1 # 0 exception
            Bi_25c = abs(pts_scatter[ind25_red]- pts_scatter[ind25_ref])/ pts_scatter[ind25_ref]



            ## R25f : time at 25% of the diff between the concentration pic and final value
            D25 =  np.array(Dend)-(np.array(Dend)-np.array(Dmin))*25/100
            D25f_ref, ind25f_ref = searchNearest( ref_var, D25[0], indmin_ref)
            D25f_red, ind25f_red = searchNearest( red_var, D25[1], indmin_red)
            Bi_pr = abs(((pts_scatter[ind25f_red]-pts_scatter[indmin_red])-(pts_scatter[ind25f_ref]-pts_scatter[indmin_ref]))\
                        /(pts_scatter[ind25f_ref]-pts_scatter[indmin_ref]))


            ## D50_mf_ref : Diff between the ref maximum and final value of each case
            D50_mf_ref = Dend[0]-Dmin[0]
            D50_mf_red = Dend[0]-Dmin[1]
            Bi_pa = abs((D50_mf_ref-D50_mf_red)/D50_mf_ref)


            QoI = [Bi_25c, Bi_cr, Bi_m, Bi_pr, Bi_pa]

        QoI = np.abs(QoI)

        return QoI







    def qoi_test(self,qoi_s,qoi_T,qoi_Sl,qoi_ig,qoi_K,conditions,red_data_meth,\
                 fileNameExt='temp.cti'):

        # Species
        max_error_sp = red_data_meth.max_error_sp
        under_tol_s = []
        for i in range(conditions.error_param.n_tspc):
            tolerance = max_error_sp[i]
            if qoi_s:
                if np.amax(qoi_s[i])*100>tolerance:
                    under_tol_s.append(False)
                else:
                    under_tol_s.append(True)

        # Temperature
        under_tol_T = True
        if qoi_T:
            if conditions.error_param.T_check:
                if qoi_T*100>red_data_meth.max_error_T:
                    under_tol_T = False


        # Flame speed
        under_tol_Sl = True
        if qoi_Sl:
            if 'flame' in conditions.config and conditions.error_param.Sl_check:
                if qoi_Sl*100>red_data_meth.max_error_Sl:
                    under_tol_Sl = False


        # Ignition time
        under_tol_ig = True
        if qoi_ig:
            if 'reactor' in conditions.config and conditions.error_param.ig_check:
                if qoi_ig*100>red_data_meth.max_error_ig:
                    under_tol_ig = False

        # Extinction strain rate
        under_tol_K = True
        if qoi_K:
            if 'diff_flame' in conditions.config and conditions.error_param.K_check:
                if qoi_K*100>red_data_meth.max_error_K:
                    under_tol_K = False

        # all tol under_tol
        under_tol=True
        for under_tol_s_i in under_tol_s:
            if not under_tol_s_i:
                under_tol=False
                break
        if not under_tol_T or not under_tol_Sl or not under_tol_ig or not under_tol_K:
            under_tol=False

        # all tol above_tol
        above_tol=True
        for under_tol_s_i in under_tol_s:
            if under_tol_s_i:
                above_tol=False
                break
        if under_tol_T or under_tol_Sl or under_tol_ig or under_tol_K:
            above_tol=False


        return under_tol_s, under_tol_T, under_tol_Sl, under_tol_ig, under_tol_K,\
               under_tol, above_tol






class Error_param :
    def __init__(self, tspc,error_calculation='points',error_coupling='mean', \
                 error_type_fit='mean',                                       \
                 T_check=False,ig_check=False,Sl_check=False,K_check=False,   \
                 sp_T=['CO2'],sp_ig=['H'],sp_Sl=['H'],sp_K=['H']):
        self.error_calculation = error_calculation
        self.error_coupling    = error_coupling
        self.error_type_fit    = error_type_fit
        self.T_check           = T_check
        self.ig_check          = ig_check
        self.Sl_check          = Sl_check
        self.K_check           = K_check
        self.strain_accuracy   = 0.05
        self.sp_T              = sp_T
        self.sp_ig             = sp_ig
        self.sp_Sl             = sp_Sl
        self.sp_K              = sp_K
        self.calculated_error  = 0
        self.tspc              = tspc
        self.n_tspc            = 0



class Mech_data:
    def __init__(self,mech,verbose=0):

        self.name      = mech
        self.gas_prop  = []
        self.spec      = Species()
        self.react     = Reactions(mech)
        self.react_ref = Reactions(mech)
        self.get_data(mech,verbose)


    def get_data(self,mech,verbose=0):

        fs = open(mech, 'r')
        self.duplicate_list = []


# =============================================================================
# Gas properties
# =============================================================================
        txt = fs.readline()
        while "Species data" not in txt:
            self.gas_prop.append(txt)
            txt = fs.readline()

# =============================================================================
# Thermo
# =============================================================================

        while "Reaction data" not in txt and txt != "":
            if txt == "":
                print("problem during the new mechanism writing")
                break
            elif "species(" in txt:
                self.spec.name.append([])
                self.spec.atoms.append([])
                self.spec.thermo.append([])
                self.spec.trans.append([])
                self.spec.note.append([])

                if "'" in txt:
                  self.spec.name[-1] = txt.split("'")[1]
                else:
                  self.spec.name[-1] = txt.split("\"")[1]
                while txt.replace(' ','') != "\n" and txt != "":
                    txt = fs.readline()
                    if "atoms" in txt:
                        if "'" in txt:
                            self.spec.atoms[-1] = txt.split("'")[1]
                        else:
                            self.spec.atoms[-1] = txt.split("\"")[1]
                    if "thermo" in txt:
                        if not ("NASA" in txt):
                          txt = fs.readline()
                        # Low temp coeffs
                        Low_T = txt.split("[")[1].split("]")[0].split(",")
                        Low_T = [float(i) for i in Low_T]
                        if txt.count("[")<2: txt = fs.readline()
                        L_T_coeffs=txt.split("[")[-1].split(",")[0:-1]
                        txt = fs.readline()
                        L_T_coeffs.extend(txt.split(",")[0:-1])
                        txt = fs.readline()
                        L_T_coeffs.extend(txt.split("]")[0].split(","))
                        L_T_coeffs = [float(i) for i in L_T_coeffs]
                        txt = fs.readline()
                        # Hight temp coeffs
                        Hight_T = txt.split("[")[1].split("]")[0].split(",")
                        Hight_T = [float(i) for i in Hight_T]
                        if txt.count("[")<2: txt = fs.readline()
                        H_T_coeffs=txt.split("[")[-1].split(",")[0:-1]
                        txt = fs.readline()
                        H_T_coeffs.extend(txt.split(",")[0:-1])
                        txt = fs.readline()
                        H_T_coeffs.extend(txt.split("]")[0].split(","))
                        H_T_coeffs = [float(i) for i in H_T_coeffs]
                        self.spec.thermo[-1]=Low_T+L_T_coeffs+Hight_T+H_T_coeffs
                    if "transport" in txt:
                        self.spec.trans[-1].append(txt)
                        while ")" not in txt:
                            txt = fs.readline()
                            self.spec.trans[-1].append(txt)
                    if "note" in txt:
                        while ")" not in txt:
                            self.spec.note[-1].append(txt)
                            txt = fs.readline()
                        self.spec.note[-1].append(txt)
                    if txt.strip() == ')': # if no note
                        self.spec.note[-1].append(txt)

            txt = fs.readline()

# =============================================================================
# Reactions
# =============================================================================

        txt = fs.readline()

        while txt != "":
            if "# Reaction " in txt:
                self.react.number.append(int(txt.split(" ")[-1]))
                txt = fs.readline()
                while txt[0]=='#':
                    txt = fs.readline()
                if "three_body_reaction" in txt:
                    # 1st line data
                    self.react.type.append("three_body_reaction")
                    if "'" in txt:
                        self.react.formula.append(txt.split("'")[1])
                    else:
                        self.react.formula.append(txt.split("\"")[1])
                    kin = txt.split("[")[1].split("]")[0].split(",")
                    self.react.kin.append([float(kin[0]),float(kin[1]),float(kin[2])])
                    # efficiencies (text)
                    self.react.txt.append([])
                    txt = fs.readline()
                    self.react.eff.append([])
                    while txt.replace(' ','') != "\n" and txt != "":
                        if "efficiencies" in txt:
                            if "'" in txt:
                                self.react.eff[-1] = txt.split("'")[1].split(" ")
                            else:
                                self.react.eff[-1] = txt.split("\"")[1].split(" ")
                        else:
                            self.react.txt[-1].append(txt)
                            if 'duplicate' in txt:
                                self.duplicate_list.append(len(self.react.number)-1)
                        txt = fs.readline()
                elif "falloff_reaction" in txt:
                    # 1st line data
                    self.react.type.append("falloff_reaction")
                    if "'" in txt:
                        self.react.formula.append(txt.split("'")[1])
                    else:
                        self.react.formula.append(txt.split("\"")[1])
                    self.react.txt.append([])
                    self.react.eff.append([])
                    txt = fs.readline()
                    while txt.replace(' ','') != "\n" and txt != "":
                        if "kf0" in txt:
                            kf0=txt.split("[")[1].split("]")[0].split(",")
                            kf0=[float(kf0[0]),float(kf0[1]),float(kf0[2])]
                        elif "kf" in txt:
                            kf=txt.split("[")[1].split("]")[0].split(",")
                            kf=[float(kf[0]),float(kf[1]),float(kf[2])]
                        elif "efficiencies" in txt:
                            if "'" in txt:
                                self.react.eff[-1] = txt.split("'")[1].split(" ")
                            else:
                                self.react.eff[-1] = txt.split("\"")[1].split(" ")
                        else: # efficiencies & falloff method
                            self.react.txt[-1].append(txt)
                            if 'duplicate' in txt:
                                self.duplicate_list.append(len(self.react.number)-1)
                        txt = fs.readline()
                    try:
                        self.react.kin.append([kf,kf0])
                    except:
                        print("Error while writing falloff coeffs of reaction:")
                        print(self.react.number[-1])

                elif "chemically_activated_reaction" in txt:
                    # 1st line data
                    self.react.type.append("chemically_activated_reaction")
                    if "'" in txt:
                        self.react.formula.append(txt.split("'")[1])
                    else:
                        self.react.formula.append(txt.split("\"")[1])
                    self.react.txt.append([])
                    self.react.eff.append([])
                    txt = fs.readline()
                    while txt.replace(' ','') != "\n" and txt != "":
                        if "kLow" in txt:
                            kLow=txt.split("[")[1].split("]")[0].split(",")
                            kLow=[float(kLow[0]),float(kLow[1]),float(kLow[2])]
                        elif "kHigh" in txt:
                            kHigh=txt.split("[")[1].split("]")[0].split(",")
                            kHigh=[float(kHigh[0]),float(kHigh[1]),float(kHigh[2])]
                        elif "efficiencies" in txt:
                            if "'" in txt:
                                self.react.eff[-1] = txt.split("'")[1].split(" ")
                            else:
                                self.react.eff[-1] = txt.split("\"")[1].split(" ")
                        else: # efficiencies & falloff method
                            self.react.txt[-1].append(txt)
                            if 'duplicate' in txt:
                                self.duplicate_list.append(len(self.react.number)-1)
                        txt = fs.readline()
                    try:
                        self.react.kin.append([kLow,kHigh])
                    except:
                        print("Error while writing coeffs of reaction:")
                        print(self.react.number[-1])

                elif "pdep_arrhenius" in txt:
                    # 1st line data
                    self.react.type.append("pdep_arrhenius")
                    if "'" in txt:
                        self.react.formula.append(txt.split("'")[1])
                    else:
                        self.react.formula.append(txt.split("\"")[1])
                    self.react.txt.append([])
                    self.react.kin.append([])
                    txt = fs.readline()
                    pdep_duplicate = False
                    while txt.replace(' ','') != "\n" and txt != "" and txt[0]!='#':
                        if 'duplicate' in txt:
                            self.react.eff.append(txt+'\n')
                            pdep_duplicate=True
                            self.duplicate_list.append(len(self.react.number)-1)
                        else:
                            try:
                                kf=txt.split("[")[1].split("]")[0].split(",")
                                self.react.kin[-1].append\
                                ([float(kf[-3]),float(kf[-2]),float(kf[-1])])
                                self.react.txt[-1].append\
                                (txt.split("(")[1].split(",")[0])
                            except:
                                print("Error while writing pdep coeffs of reaction:")
                                print(self.react.number[-1])
                        txt = fs.readline()
                    if not pdep_duplicate:
                        self.react.eff.append("")

                elif "chebyshev" in txt:
                    # 1st line data
                    self.react.type.append("chebyshev")
                    if "'" in txt:
                        self.react.formula.append(txt.split("'")[1])
                    else:
                        self.react.formula.append(txt.split("\"")[1])
                    self.react.txt.append([])
                    self.react.kin.append([])
                    txt = fs.readline()
                    chebyshev_duplicate = False
                    while txt.replace(' ','') != "\n" and txt != "" and txt[0]!='#':
                        if 'duplicate' in txt:
                            self.react.eff.append(txt+'\n')
                            chebyshev_duplicate=True
                            self.duplicate_list.append(len(self.react.number)-1)
                        elif 'Tmin' in txt:
                            self.react.txt[-1].append(txt)
                        elif 'Pmin' in txt:
                            self.react.txt[-1].append(txt)
                        else:
                            txt = txt.replace('[','')
                            txt = txt.replace(']','')
                            txt = txt.replace(')','')
                            txt = txt.replace('\n','')
                            txt = txt.replace('coeffs','')
                            txt = txt.replace('=','')
                            self.react.kin[-1].append([])
                            for i in range(len(txt.split(","))):
                                if i==len(txt.split(","))-1:
                                    try:
                                        self.react.kin[-1][-1].append(float(txt.split(",")[i]))
                                    except:
                                        a=2
                                else:
                                    self.react.kin[-1][-1].append(float(txt.split(",")[i]))
#
#                            try:
#                                kf=txt.split(",")
#                                self.react.kin[-1].append\
#                                ([float(kf[0]),float(kf[1]),float(kf[2]),float(kf[3])])
#                            except:
#                                print("Error while writing chebyshev coeffs of reaction:")
#                                print(self.react.number[-1])
                        txt = fs.readline()
                    if not chebyshev_duplicate:
                        self.react.eff.append("")



                elif "reaction" in txt:
                    # 1st line data
                    self.react.type.append("reaction")
                    if "'" in txt:
                        self.react.formula.append(txt.split("'")[1])
                    else:
                        self.react.formula.append(txt.split("\"")[1])
                    try:
                        kin = txt.split("[")[1].split("]")[0].split(",")
                    except:
                        print('toto')
                    self.react.kin.append([float(kin[0]),float(kin[1]),float(kin[2])])
                    # other lines data (text)
                    self.react.txt.append([])
                    txt = fs.readline()
                    while txt.replace(' ','') != "\n" and txt != "":
                        self.react.txt[-1].append(txt)
                        if 'duplicate' in txt:
                            self.duplicate_list.append(len(self.react.number)-1)
                        txt = fs.readline()
                    self.react.eff.append("none")

            txt = fs.readline();



        # check threebody exception (+AR) (+HE) etc.
        for r in range(len(self.react.formula)):
            self.react.tbe.append(False)
            if self.react.type[r]=="falloff_reaction":
                for sp in self.spec.name:
                    if '(+'  + sp +')' in self.react.formula[r] \
                    or '(+ ' + sp +')' in self.react.formula[r]:
                        self.react.tbe[-1]=sp
                        break

        # check submechanism family:
        ct.suppress_thermo_warnings()
        gas = ct.Solution(mech)
        nu_f = gas.reactant_stoich_coeffs()
        nu_r = gas.product_stoich_coeffs()
        self.react.subm_C  = [0]*len(self.react.formula)
        self.react.subm_N  = [0]*len(self.react.formula)
        self.react.subm_S  = [0]*len(self.react.formula)
        self.react.subm_Si = [0]*len(self.react.formula)
        self.react.subm_CO = [False]*len(self.react.formula)
        for r in range(len(self.react.formula)):
            for sp in range(len(self.spec.name)):
                if (nu_f[sp][r]!=0 or nu_r[sp][r]!=0) and nu_f[sp][r]!=nu_r[sp][r]:
                    n_at = self.spec.atoms[sp].split(' ')
                    for at in n_at:
                        if 'C' in at or 'c' in at:
                            n_C = int(at.split(':')[-1])
                            self.react.subm_C[r] = max(self.react.subm_C[r],n_C)
                        if 'N' in at or 'n' in at:
                            n_N = int(at.split(':')[-1])
                            self.react.subm_N[r] = max(self.react.subm_N[r],n_N)
                        if 'S' in at or 's' in at:
                            n_S = int(at.split(':')[-1])
                            self.react.subm_S[r] = max(self.react.subm_S[r],n_S)
                        if 'Si' in at or 'si' in at:
                            n_Si = int(at.split(':')[-1])
                            self.react.subm_Si[r] = max(self.react.subm_Si[r],n_Si)
            if self.react.subm_C[r]==1:
                self.react.subm_CO[r] = True
                for sp in range(len(self.spec.name)):
                    if (nu_f[sp][r]!=0 or nu_r[sp][r]!=0) and nu_f[sp][r]!=nu_r[sp][r]:
                        n_at = self.spec.atoms[sp].split(' ')
                        at = str(n_at)
                        if ('C' in at or 'c' in at) and ('H' in at or 'h' in at):
                            self.react.subm_CO[r] = False


        self.spec.activ_p  = [False]*len(self.spec.name)       # act spec for cases loop
        self.spec.activ_m  = [True]*len(self.spec.name)       # act spec for methods loop
        self.react.activ_p = [False]*len(self.react.formula)   # act react for cases loop
        self.react.activ_m = [True]*len(self.react.formula)   # act react for methods loop
        self.react.modif   = [True] *len(self.react.formula)
        self.react_ref     = self.react
        self.react.ref_kin = copy.deepcopy(self.react.kin)


    def write_new_mech(self,filename="temp.cti",act_sp="no_arg",act_r="no_arg"):

        fd = open(filename, 'w')

        # active species and reactions lists
        if act_sp=="no_arg": act_sp = self.spec.activ_p
        if act_r=="no_arg" : act_r  = self.react.activ_p
        if not act_sp:
            act_sp=self.spec.activ_m ; act_r =self.react.activ_m
        elif True not in act_sp:
            act_sp=self.spec.activ_m ; act_r =self.react.activ_m

        sp_list=[];remove_sp_list=[]
        for sp in range(len(act_sp)):
            if act_sp[sp]:
                sp_list.append(sp)
            else:
                remove_sp_list.append(sp)
        r_list=[]
        for r in range(len(act_r)):
            if act_r[r]:
                r_list.append(r)


# =============================================================================
#         Units ans gas properties
# =============================================================================

#        self.gas_prop = []

        # Units writing
        l=0
        while self.gas_prop[l]!="\n":
            fd.write(self.gas_prop[l])
            l+=1


        # Gas characteristics writing
        while "elements" not in self.gas_prop[l]:
            fd.write(self.gas_prop[l])
            l+=1

        #elements
        txt_elements = self.gas_prop[l].split('"')[0]+'"'
        elements = ["C","H","O","N","AR","HE","Ar","He", "c", "h", "o", "n", "ar", "he"]
        for el in elements:
            if self.find_element(el,act_sp):
#                if el == "AR":
#                    txt_elements += "Ar" + " "
#                elif el == "HE":
#                    txt_elements += "He" + " "
#                else:
                txt_elements += el.capitalize() + " "
        txt_elements += '",\n'
        fd.write(txt_elements)

        #species
        txt_spec = '          species="""'
        i=0
        for sp in sp_list:
            if i<5 :
                txt_spec+=self.spec.name[sp]+"  "
                i+=1
            else:
                fd.write(txt_spec+"\n")
                txt_spec="                     "+self.spec.name[sp]+"  "
                i=1
        txt_spec+='""",\n'
        fd.write(txt_spec)

        while "reactions" not in self.gas_prop[l]:
            l+=1

        while self.gas_prop[l]!="\n":
            fd.write(self.gas_prop[l])
            l+=1

        fd.write("\n")


#-------------------------------------------------------------------------------
# Species data
#-------------------------------------------------------------------------------


# =============================================================================
#         Species
# =============================================================================

        # ---------------------------------------------------------------------
        #         CSP radicals
        # ---------------------------------------------------------------------
        if len(self.spec.CSP_radicals)>0:
            fd.write("#-------------------------------------------------------------------------------\n")
            fd.write("# CSP radicals\n")
            fd.write("#-------------------------------------------------------------------------------\n")
            fd.write("#\n#  ")
            # check if species have been identified as radicals at every conditions
            CSP_radicals_always = [True]*len(self.spec.name)
            for _l in range(len(self.spec.CSP_radicals)):
                for sp in range(len(self.spec.CSP_radicals[0])):
                    if not self.spec.CSP_radicals[_l][sp]:
                        CSP_radicals_always[sp] = False
            # write radicals
            _n_rad = 0
            for sp in range(len(CSP_radicals_always)):
                if _n_rad == 10:
                    fd.write('\n#  ') ; _n_rad = 0
                if CSP_radicals_always[sp] and self.spec.activ_m[sp]:
                    fd.write(self.spec.name[sp] + '  ') ; _n_rad +=1
            fd.write("\n#\n#\n")

        fd.write("#-------------------------------------------------------------------------------\n")
        fd.write("# Species data\n")
        fd.write("#-------------------------------------------------------------------------------\n")
        fd.write("\n")

        for sp in sp_list: # in activated species list
            l_name  = "species(name='"+self.spec.name[sp]+"',\n"
            l_atoms ="        atoms='"+self.spec.atoms[sp]+"',\n"
            l_thermo="        thermo=(NASA(["\
                + str('%0.2e' %self.spec.thermo[sp][0])+" ,"     \
                + str('%0.2e' %self.spec.thermo[sp][1])          \
                + "],\n                     [ "                  \
                + str('%0.8e' %self.spec.thermo[sp][2])+", "     \
                + str('%0.8e' %self.spec.thermo[sp][3])+", "     \
                + str('%0.8e' %self.spec.thermo[sp][4])+",\n"    \
                + "                       "                      \
                + str('%0.8e' %self.spec.thermo[sp][5])+", "     \
                + str('%0.8e' %self.spec.thermo[sp][6])+", "     \
                + str('%0.8e' %self.spec.thermo[sp][7])+",\n"    \
                + "                       "                      \
                + str('%0.8e' %self.spec.thermo[sp][8])+"]),\n"  \
                + "                NASA(["                       \
                + str('%0.2e' %self.spec.thermo[sp][9])+" ,"     \
                + str('%0.2e' %self.spec.thermo[sp][10])         \
                + "],\n                     [ "                  \
                + str('%0.8e' %self.spec.thermo[sp][11])+", "    \
                + str('%0.8e' %self.spec.thermo[sp][12])+", "    \
                + str('%0.8e' %self.spec.thermo[sp][13])+",\n"   \
                + "                       "                      \
                + str('%0.8e' %self.spec.thermo[sp][14])+", "    \
                + str('%0.8e' %self.spec.thermo[sp][15])+", "    \
                + str('%0.8e' %self.spec.thermo[sp][16])+",\n"   \
                + "                       "

            if len(self.spec.trans[sp])!=0 or len(self.spec.note[sp])!=0:
                l_thermo+= str('%0.8e' %self.spec.thermo[sp][17])+"])),\n"
            else:
                l_thermo+= str('%0.8e' %self.spec.thermo[sp][17])+"])))\n"

            l_trans=""
            for line in self.spec.trans[sp]:
                l_trans+=line
            l_note=""
            for line in self.spec.note[sp]:
                l_note+=line

            fd.write(l_name)
            fd.write(l_atoms)
            fd.write(l_thermo)
            fd.write(l_trans)
            fd.write(l_note)
            fd.write("\n")


# =============================================================================
#         Reactions
# =============================================================================

        fd.write("#-------------------------------------------------------------------------------\n")
        fd.write("# Reaction data\n")
        fd.write("#-------------------------------------------------------------------------------\n")
        fd.write("\n")

        r_list_ = copy.deepcopy(r_list)
        for r in r_list_:
            if r in self.duplicate_list:
                for d in self.duplicate_list:
                        if self.react.formula[d]==self.react.formula[r]:
                            if d not in r_list:
                                r_list.append(d)
        r_list.sort()

        for r in r_list: # in activated reactions list

            fd.write("\n# Reaction "+str(self.react.number[r])+"\n")

            if self.react.type[r]=="reaction":
                #    # ReactionXX
                #    reaction('H2O2 + OH <=> H2O + HO2', [1.740000e+12, 0.0, 318.0],
                #             options='duplicate')
                r_line = "reaction('"+self.react.formula[r]+"', ["            \
                         + str('%0.6e' %self.react.kin[r][0]) + ", "          \
                         + str('%0.4f' %self.react.kin[r][1]) + ", "          \
                         + str('%0.3f' %self.react.kin[r][2]) + "]"
                if self.react.txt[r]==[] \
                or (len(self.react.txt[r])==0 and self.react.txt[r][0].replace(' ',''))=='':
                    r_txt=")\n"
                else:
                    r_txt=",\n"
                    for line in self.react.txt[r]:
                        r_txt += line
                fd.write(r_line)
                fd.write(r_txt)

            if self.react.type[r]=="three_body_reaction":
                #    # Reaction XX
                #    three_body_reaction('O + O + M <=> O2 + M', [6.165000e+15, -0.5, 0.0],
                #    efficiencies='CO2:3.8 CO:1.9 H2:2.5 H2O:12.0 AR:0.83 CH4:2.0 C2H6:3.0 HE:0.83')
                r_line = "three_body_reaction('"+self.react.formula[r]+"', [" \
                         + str('%0.6e' %self.react.kin[r][0]) + ", "          \
                         + str('%0.4f' %self.react.kin[r][1]) + ", "          \
                         + str('%0.3f' %self.react.kin[r][2]) + "]"
                if self.react.eff[r]!=[] or self.react.txt[r]!=[]:
                    if self.react.eff[r]!=[]:
                        r_efficiencies=",\n"
                        r_efficiencies += "                    efficiencies='"
                        for sp in sp_list:
                            for coll_sp in self.react.eff[r]:
                                if self.spec.name[sp] == coll_sp.split(":")[0]:
                                    r_efficiencies+=coll_sp+" "
                    else: r_efficiencies=""
                    if self.react.txt[r]!=[]:
                        r_efficiencies+="',\n"
                        r_txt=''
                        for line in self.react.txt[r]:
                            r_txt += line
                    else: r_txt="')\n"
                else: r_efficiencies="";r_txt=")\n"
                fd.write(r_line)
                fd.write(r_efficiencies)
                fd.write(r_txt)

            if self.react.type[r]=="falloff_reaction":
                #    # Reaction 71
                #    falloff_reaction('CO + H2 (+ M) <=> CH2O (+ M)',
                #                     kf=[4.300000e+07, 1.5, 79600.0],
                #                     kf0=[5.070000e+27, -3.42, 84348.0],
                #                     efficiencies='CO2:2.0 CO:1.5 H2:2.0 H2O:6.0 AR:0.7 CH4:2.0 C2H6:3.0 HE:0.7',
                #                     falloff=Troe(A=0.932, T3=197.0, T1=1540.0, T2=10300.0))
                r_line = "falloff_reaction('"+self.react.formula[r]+"',\n"    \
                         + "                 kf=["                            \
                         + str('%0.6e' %self.react.kin[r][0][0]) + ", "       \
                         + str('%0.4f' %self.react.kin[r][0][1]) + ", "       \
                         + str('%0.3f' %self.react.kin[r][0][2]) + "],\n"     \
                         + "                 kf0=["                           \
                         + str('%0.6e' %self.react.kin[r][1][0]) + ", "       \
                         + str('%0.4f' %self.react.kin[r][1][1]) + ", "       \
                         + str('%0.3f' %self.react.kin[r][1][2]) + "]"
                r_efficiencies=""
                if self.react.eff[r]!=[] or self.react.txt[r]!=[]:
                    if self.react.eff[r]!=[]:
                        r_efficiencies += ",\n                 efficiencies='"
                        for sp in sp_list:
                            for coll_sp in self.react.eff[r]:
                                if self.spec.name[sp] == coll_sp.split(":")[0]:
                                    r_efficiencies+=coll_sp+" "
                        if self.react.txt[r]!=[]:
                            r_efficiencies+="',\n"
                    elif self.react.txt[r]!=[]:
                        r_efficiencies+=",\n"
                    if self.react.txt[r]!=[]:
                        r_txt=""
                        for line in self.react.txt[r]:
                            r_txt += line
                    else:
                        r_txt="')\n"
                else:
                    r_txt=")\n"
                fd.write(r_line)
                fd.write(r_efficiencies)
                fd.write(r_txt)

            if self.react.type[r]=="chemically_activated_reaction":
                ## Reaction 134
                #chemically_activated_reaction('ch3 + oh (+ M) <=> ch2(s) + h2o (+ M)',
                #                              kLow=[1.128000e+15, -0.63327, -493.15597],
                #                              kHigh=[2.394000e-03, 4.096, -1241.875],
                #                              falloff=Troe(A=2.122, T3=837.667, T1=2326.05, T2=4432.0))
                r_line = "chemically_activated_reaction('"+self.react.formula[r]+"',\n"    \
                         + "                 kLow=["                          \
                         + str('%0.6e' %self.react.kin[r][0][0]) + ", "       \
                         + str('%0.4f' %self.react.kin[r][0][1]) + ", "       \
                         + str('%0.3f' %self.react.kin[r][0][2]) + "],\n"     \
                         + "                 kHigh=["                         \
                         + str('%0.6e' %self.react.kin[r][1][0]) + ", "       \
                         + str('%0.4f' %self.react.kin[r][1][1]) + ", "       \
                         + str('%0.3f' %self.react.kin[r][1][2]) + "]"
                r_efficiencies=""
                if self.react.eff[r]!=[] or self.react.txt[r]!=[]:
                    if self.react.eff[r]!=[]:
                        r_efficiencies += ",\n                 efficiencies='"
                        for sp in sp_list:
                            for coll_sp in self.react.eff[r]:
                                if self.spec.name[sp] == coll_sp.split(":")[0]:
                                    r_efficiencies+=coll_sp+" "
                        if self.react.txt[r]!=[]:
                            r_efficiencies+="',\n"
                    elif self.react.txt[r]!=[]:
                        r_efficiencies+=",\n"
                    if self.react.txt[r]!=[]:
                        r_txt=""
                        for line in self.react.txt[r]:
                            r_txt += line
                    else:
                        r_txt="')\n"
                else:
                    r_txt=")\n"
                fd.write(r_line)
                fd.write(r_efficiencies)
                fd.write(r_txt)

            if self.react.type[r]=="pdep_arrhenius":
                #    # Reaction 391
                #    pdep_arrhenius('CH3COCH3 <=> CH3CO + CH3',
                #                   [(0.01, 'atm'), 2.050000e+58, -12.796, 100030.1],
                #                   [(0.1, 'atm'), 3.300000e+51, -10.574, 98221.2],
                #                   [(1.0, 'atm'), 1.310000e+42, -7.657, 94660.6],
                #                   [(10.0, 'atm'), 2.160000e+33, -4.989, 90916.5],
                #                   [(100.0, 'atm'), 9.400000e+28, -3.669, 89022.8])
                r_line = "pdep_arrhenius('"+self.react.formula[r]+"'"
                k_line = ""
                for l in range(len(self.react.txt[r])):
                    try:
                        k_line += ",\n               [("                           \
                                 + self.react.txt[r][l] + ", 'atm'), "            \
                                 + '%0.6e' %self.react.kin[r][l][0] + ", "  \
                                 + '%0.4f' %self.react.kin[r][l][1] + ", "  \
                                 + '%0.3f' %self.react.kin[r][l][2] + "]"
                    except:
                        print('toto')
                if self.react.eff[r]=="":  k_line+=")\n"
                else :                     k_line+=",\n"
                fd.write(r_line)
                fd.write(k_line)
                fd.write(self.react.eff[r])

            if self.react.type[r]=="chebyshev":
                #    # Reaction 391
                #    pdep_arrhenius('CH3COCH3 <=> CH3CO + CH3',
                #                   [(0.01, 'atm'), 2.050000e+58, -12.796, 100030.1],
                #                   [(0.1, 'atm'), 3.300000e+51, -10.574, 98221.2],
                #                   [(1.0, 'atm'), 1.310000e+42, -7.657, 94660.6],
                #                   [(10.0, 'atm'), 2.160000e+33, -4.989, 90916.5],
                #                   [(100.0, 'atm'), 9.400000e+28, -3.669, 89022.8])
                r_line = "chebyshev_reaction('"+self.react.formula[r]+"',\n"
                k_line = ""
                for l in range(len(self.react.txt[r])):
                    k_line += self.react.txt[r][l]
                for l in range(len(self.react.kin[r])):
                    if l==0:
                        k_line+='                   coeffs=[['
                    else:
                        k_line+='                           ['
                    for c in range(len(self.react.kin[r][l])):
                        k_line+='%0.5e' %self.react.kin[r][l][c]
                        if c<len(self.react.kin[r][l])-1:
                            k_line+=", "
                        else:
                            k_line+="]"
#                    k_line+= '%0.5e' %self.react.kin[r][l][0] + ", "  \
#                           + '%0.5f' %self.react.kin[r][l][1] + ", "  \
#                           + '%0.5f' %self.react.kin[r][l][2] + ","   \
#                           + '%0.5f' %self.react.kin[r][l][3] + "]"
                    if l==len(self.react.kin[r])-1:
                        k_line+=']'
                    else:
                        k_line+=',\n'

                if self.react.eff[r]=="":  k_line+=")\n"
                else :                     k_line+=",\n"
                fd.write(r_line)
                fd.write(k_line)
                fd.write(self.react.eff[r])




        timer.sleep(1.5)

        fd.close()





    def write_chemkin_mech(self,filename="temp.ck",act_sp="no_arg",act_r="no_arg"):

        if len(filename)>3:
            if filename[-4:]=='.cti':
                filename = filename[:-4] + '.ck'

        fd = open(filename, 'w')

        # active species and reactions lists
        if act_sp=="no_arg": act_sp = self.spec.activ_p
        if act_r=="no_arg" : act_r  = self.react.activ_p
        if True not in act_sp:
            act_sp=self.spec.activ_m ; act_r =self.react.activ_m

        sp_list=[];remove_sp_list=[]
        for sp in range(len(act_sp)):
            if act_sp[sp]:
                sp_list.append(sp)
            else:
                remove_sp_list.append(sp)
        r_list=[]
        for r in range(len(act_r)):
            if act_r[r]:
                r_list.append(r)


# =============================================================================
#         Units ans gas properties
# =============================================================================

#        self.gas_prop = []

#        # Units writing
#        l=0
#        while self.gas_prop[l]!="\n":
#            fd.write(self.gas_prop[l])
#            l+=1
#
#
#        # Gas characteristics writing
#        while "elements" not in self.gas_prop[l]:
#            fd.write(self.gas_prop[l])
#            l+=1

        #elements
        txt_elements = 'ELEMENTS\n'#self.gas_prop[l].split('"')[0]+'"'
        elements = ["C","H","O","N","AR","HE","Ar","He"]
        for el in elements:
            if self.find_element(el,act_sp):
                #if el == "AR":
                    #txt_elements += "Ar" + " "
                #elif el == "HE":
                    #txt_elements += "He" + " "
                #else:
                txt_elements += str.upper(el) + " "
        txt_elements += '\nEND\n'
        fd.write(txt_elements)


        #species
        txt_spec = 'SPECIES\n'
        i=0
        for sp in sp_list:
            if i<5 :
                txt_spec+=self.spec.name[sp]+"  "
                i+=1
            else:
                fd.write(txt_spec+"\n")
                txt_spec=self.spec.name[sp]+"  "
                i=1
        txt_spec+='\nEND\n'

        fd.write(txt_spec)

#        while "reactions" not in self.gas_prop[l]:
#            l+=1
#
#        while self.gas_prop[l]!="\n":
#            fd.write(self.gas_prop[l])
#            l+=1

#        fd.write("\n")
        # ---------------------------------------------------------------------
        #         CSP radicals
        # ---------------------------------------------------------------------
        if len(self.spec.CSP_radicals)>0:
            fd.write("!-------------------------------------------------------------------------------\n")
            fd.write("! CSP radicals\n")
            fd.write("!-------------------------------------------------------------------------------\n")
            fd.write("!\n!  ")
            # check if species have been identified as radicals at every conditions
            CSP_radicals_always = [True]*len(self.spec.name)
            for _l in range(len(self.spec.CSP_radicals)):
                for sp in range(len(self.spec.CSP_radicals[0])):
                    if not self.spec.CSP_radicals[_l][sp]:
                        CSP_radicals_always[sp] = False
            # write radicals
            _n_rad = 0
            for sp in range(len(CSP_radicals_always)):
                if _n_rad == 10:
                    fd.write('\n!  ') ; _n_rad = 0
                if CSP_radicals_always[sp] and self.spec.activ_m[sp]:
                    fd.write(self.spec.name[sp] + '  ') ; _n_rad +=1
            fd.write("\n!\n!\n")



# =============================================================================
#         Reactions
# =============================================================================

        fd.write("!-------------------------------------------------------------------------------\n")
        fd.write("! Reaction data\n")
        fd.write("!-------------------------------------------------------------------------------\n")
        fd.write("REACTIONS\n")

        r_list_ = copy.deepcopy(r_list)
        for r in r_list_:
            if r in self.duplicate_list:
                for d in self.duplicate_list:
                        if self.react.formula[d]==self.react.formula[r]:
                            if d not in r_list:
                                r_list.append(d)
        r_list.sort()

        for r in r_list: # in activated reactions list
            dup_txt = ''
            txt_space = ''
            for ts in range(50-len(self.react.formula[r])):
                txt_space += ' '

            # comments
            comments = ''
            for line in self.react.txt[r]:
                if "options='duplicate'" in line:
                    dup_txt = ' DUPLICATE\n'
                elif "falloff=Troe(" in line or self.react.type[r]=="pdep_arrhenius" or self.react.type[r]=="chebyshev":
                    if "options='duplicate'" in self.react.eff[r]:
                        dup_txt = ' DUPLICATE\n'
                    comments += ''
                else:
                    comments +=  '!' + line

            #fd.write("\n# Reaction "+str(self.react.number[r])+"\n")
            if self.react.type[r]=="reaction":
                # -------- Cantera --------
                #    # ReactionXX
                #    reaction('H2O2 + OH <=> H2O + HO2', [1.740000e+12, 0.0, 318.0],
                #             options='duplicate')
#                r_line = "reaction('"+self.react.formula[r]+"', ["            \
#                         + str('%0.6e' %self.react.kin[r][0]) + ", "          \
#                         + str('%0.4f' %self.react.kin[r][1]) + ", "          \
#                         + str('%0.3f' %self.react.kin[r][2]) + "]"
                r_line = self.react.formula[r] + txt_space              \
                         + str('%0.6e' %self.react.kin[r][0]) + "  "    \
                         + str('%0.4f' %self.react.kin[r][1]) + "  "    \
                         + str('%0.3f' %self.react.kin[r][2]) + "\n"
#                r_txt = ''

                fd.write(r_line)
                fd.write(comments)

            if self.react.type[r]=="three_body_reaction":
                # -------- Cantera --------
                #    # Reaction XX
                #    three_body_reaction('O + O + M <=> O2 + M', [6.165000e+15, -0.5, 0.0],
                #    efficiencies='CO2:3.8 CO:1.9 H2:2.5 H2O:12.0 AR:0.83 CH4:2.0 C2H6:3.0 HE:0.83')
                r_line = self.react.formula[r]+txt_space \
                         + str('%0.6e' %self.react.kin[r][0]) + ", "          \
                         + str('%0.4f' %self.react.kin[r][1]) + ", "          \
                         + str('%0.3f' %self.react.kin[r][2]) + "\n"
                if self.react.eff[r]!=[] or self.react.txt[r]!=[]:
                    r_efficiencies=""
                    if self.react.eff[r]!=[]:
                        r_efficiencies=""
                        #r_efficiencies += "                    efficiencies='"
                        for sp in sp_list:
                            for coll_sp in self.react.eff[r]:
                                if self.spec.name[sp] == coll_sp.split(":")[0]:
                                    r_efficiencies+=coll_sp.split(":")[0]+'/'\
                                                    +coll_sp.split(":")[1]+'/ '
                        r_efficiencies+="\n" ;
                    if self.react.txt[r]!=[]:
                        r_txt = ''
                        for line in self.react.txt[r]:
                            r_txt += '!' + line
#                            if "options='duplicate'" in line:
#                                dup_txt = ' DUPLICATE\n'
                        r_txt += '\n'
                    else: r_txt="\n"
                else: r_efficiencies="";r_txt="\n"

                fd.write(r_line)
                fd.write(r_efficiencies)
                fd.write(comments)

            if self.react.type[r]=="falloff_reaction":
                # -------- Cantera --------
                #    # Reaction 71
                #    falloff_reaction('CO + H2 (+ M) <=> CH2O (+ M)',
                #                     kf=[4.300000e+07, 1.5, 79600.0],
                #                     kf0=[5.070000e+27, -3.42, 84348.0],
                #                     efficiencies='CO2:2.0 CO:1.5 H2:2.0 H2O:6.0 AR:0.7 CH4:2.0 C2H6:3.0 HE:0.7',
                #                     falloff=Troe(A=0.932, T3=197.0, T1=1540.0, T2=10300.0))
                # -------- Chemkin --------
                #H+CH2(+M)<=>CH3(+M)                      6.000E+14     .000        .00
                #     LOW  /  1.040E+26   -2.760   1600.00/
                #     TROE/   .5620  91.00  5836.00  8552.00/
                #H2/2.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/ AR/ .70/
                r_line = self.react.formula[r]+txt_space        \
                + str('%0.6e' %self.react.kin[r][0][0]) + '  '  \
                + str('%0.4e' %self.react.kin[r][0][1]) + '  '  \
                + str('%0.3e' %self.react.kin[r][0][2]) + '\n'  \
                + '     LOW /  '                                \
                + str('%0.6e' %self.react.kin[r][1][0]) + '  '  \
                + str('%0.4e' %self.react.kin[r][1][1]) + '  '  \
                + str('%0.3e' %self.react.kin[r][1][2]) + ' /\n'

                if self.react.eff[r]!=[] or self.react.txt[r]!=[]:
#                    else: r_efficiencies=""
                    if self.react.txt[r]!=[]:
                        r_txt = ''
                        for line in self.react.txt[r]:
                            if 'falloff=Troe' in line:
                                r_txt += '     TROE/   '                     \
                                + line.split('A=')[1].split(',')[0]  + '  ' \
                                + line.split('T3=')[1].split(',')[0] + '  ' \
                                + line.split('T1=')[1].split(',')[0].split(')')[0]
                            if 'T2=' in line:
                                r_txt += '  ' + line.split('T2=')[1].split(')')[0] + '/\n'
                            else:
                                r_txt += '/\n'
                    else: r_txt=""
                    if self.react.eff[r]!=[]:
                        r_efficiencies=""
                        #r_efficiencies += "                    efficiencies='"
                        for sp in sp_list:
                            for coll_sp in self.react.eff[r]:
                                if self.spec.name[sp] == coll_sp.split(":")[0]:
                                    r_efficiencies+=coll_sp.split(":")[0]+'/'\
                                                    +coll_sp.split(":")[1]+'/ '
                        r_efficiencies+="\n"
                    else: r_efficiencies=""
                else: r_txt="" ; r_efficiencies=""
#                for line in self.react.txt[r]:
#                    r_txt += '!' + line

                fd.write(r_line)
                fd.write(r_txt)
                fd.write(r_efficiencies)
                fd.write(comments)
#                    if "options='duplicate'" in line:
#                        dup_txt = ' DUPLICATE\n'
#                fd.write(dup_txt)

            if self.react.type[r]=="pdep_arrhenius":
                # -------- Cantera --------
                #    # Reaction 391
                #    pdep_arrhenius('CH3COCH3 <=> CH3CO + CH3',
                #                   [(0.01, 'atm'), 2.050000e+58, -12.796, 100030.1],
                #                   [(0.1, 'atm'), 3.300000e+51, -10.574, 98221.2],
                #                   [(1.0, 'atm'), 1.310000e+42, -7.657, 94660.6],
                #                   [(10.0, 'atm'), 2.160000e+33, -4.989, 90916.5],
                #                   [(100.0, 'atm'), 9.400000e+28, -3.669, 89022.8])
                # -------- Chemkin --------
                #CH3+OH<=>CH2(S)+H2O                               4.936E+014    -0.669    -445.8
                #PLOG/      0.0100     4.936E+014    -0.669      -445.8/
                #PLOG/      0.1000     1.207E+015    -0.778      -175.6/
                #PLOG/      1.0000     5.282E+017    -1.518      1772.0/
                #PLOG/     10.0000     4.788E+023    -3.155      7003.0/
                #PLOG/    100.0000     8.433E+019    -1.962      8244.0/

                r_line = self.react.formula[r]+txt_space \
                             + '%0.3e' %self.react.kin[r][0][0] + "  "  \
                             + '%0.3f' %self.react.kin[r][0][1] + "  "  \
                             + '%0.1f' %self.react.kin[r][0][2] + "\n"
                k_line = ""
                for l in range(len(self.react.txt[r])):
                    if float(self.react.txt[r][l])<9:   txt_space = '      ';
                    elif float(self.react.txt[r][l])<99: txt_space = '     ';
                    elif float(self.react.txt[r][l])<999: txt_space = '    ';
                    elif float(self.react.txt[r][l])<9999: txt_space = '   ';
                    elif float(self.react.txt[r][l])<99999: txt_space = '  ';

                    k_line += "PLOG/" + txt_space                             \
                             + '%0.4f' %float(self.react.txt[r][l]) + "     " \
                             + '%0.3e' %self.react.kin[r][l][0]     + "     " \
                             + '%0.3f' %self.react.kin[r][l][1]     + "     " \
                             + '%0.1f' %self.react.kin[r][l][2]     + "/\n"
#                if self.react.eff[r]=="":  k_line+=")\n"
#                else :                     k_line+=",\n"
                fd.write(r_line)
                fd.write(k_line)
                fd.write(comments)
#                fd.write(self.react.eff[r])

            if self.react.type[r]=="chemically_activated_reaction":
                # -------- Cantera --------
                ## Reaction 134
                #chemically_activated_reaction('ch3 + oh (+ M) <=> ch2(s) + h2o (+ M)',
                #                              kLow=[1.128000e+15, -0.63327, -493.15597],
                #                              kHigh=[2.394000e-03, 4.096, -1241.875],
                #                              falloff=Troe(A=2.122, T3=837.667, T1=2326.05, T2=4432.0))
                # -------- Chemkin --------
                #ch3+oh(+m)<=>ch2(s)+h2o(+m) 1.128E+15 -0.63327 -493.15597
                #high/ 2.394E-03 4.096 -1241.875/
                #troe/ 2.122E+00 837.667 2326.05 4.432E+03/

                r_line = self.react.formula[r]+txt_space        \
                + str('%0.6e' %self.react.kin[r][0][0]) + '  '  \
                + str('%0.4e' %self.react.kin[r][0][1]) + '  '  \
                + str('%0.3e' %self.react.kin[r][0][2]) + '\n'  \
                + '     HIGH/  '                                \
                + str('%0.6e' %self.react.kin[r][1][0]) + '  '  \
                + str('%0.4e' %self.react.kin[r][1][1]) + '  '  \
                + str('%0.3e' %self.react.kin[r][1][2]) + ' /\n'

                if self.react.eff[r]!=[] or self.react.txt[r]!=[]:
#                    else: r_efficiencies=""
                    if self.react.txt[r]!=[]:
                        r_txt = ''
                        for line in self.react.txt[r]:
                            if 'falloff=Troe' in line:
                                r_txt += '     TROE/   '                     \
                                + line.split('A=')[1].split(',')[0]  + '  ' \
                                + line.split('T3=')[1].split(',')[0] + '  ' \
                                + line.split('T1=')[1].split(',')[0].split(')')[0]
                            if 'T2=' in line:
                                r_txt += '  ' + line.split('T2=')[1].split(')')[0] + '/\n'
                            else:
                                r_txt += '/\n'
                    else: r_txt=""
                    if self.react.eff[r]!=[]:
                        r_efficiencies=""
                        #r_efficiencies += "                    efficiencies='"
                        for sp in sp_list:
                            for coll_sp in self.react.eff[r]:
                                if self.spec.name[sp] == coll_sp.split(":")[0]:
                                    r_efficiencies+=coll_sp.split(":")[0]+'/'\
                                                    +coll_sp.split(":")[1]+'/ '
                        r_efficiencies+="\n"
                    else: r_efficiencies=""
                else: r_txt="" ; r_efficiencies=""
#                for line in self.react.txt[r]:
#                    r_txt += '!' + line

                fd.write(r_line)
                fd.write(r_txt)
                fd.write(r_efficiencies)
                fd.write(comments)
#                    if "options='duplicate'" in line:
#                        dup_txt = ' DUPLICATE\n'
#                fd.write(dup_txt)

            if self.react.type[r]=="chebyshev":
                # -------- Cantera --------
                ## Reaction 71
                #chebyshev_reaction('NC3H7O (+ M) <=> C3H6OH-3 (+ M)',
                #                   Tmin=290.0, Tmax=3000.0,
                #                   Pmin=(0.00986923267, 'atm'), Pmax=(98.6923267, 'atm'),
                #                   coeffs=[[-5.58710e+00,  2.82470e+00, -3.80420e-01],
                #                           [ 1.64120e+00,  5.05000e-01,  2.10230e-01],
                #                           [-4.36070e+00, -1.78570e-01,  5.11070e-02],
                #                           [-2.03440e+00, -3.36440e-01, -2.58330e-02],
                #                           [-7.42790e-01, -2.11990e-01, -4.01300e-02]])
                # -------- Chemkin --------
                # NC3H7O(+m)=C3H6OH-3(+m) 1.0E0 0.0 0.0                                                          !C4H9O(+m)=CC4H8OH(+m) from n-butanol  Cai et al. E&F 26(2012)5550-5568
                # TCHEB / 290.0 3000.0 / PCHEB / 0.009869232667160128 98.69232667160128 /
                # CHEB / 5 3 /
                # CHEB / -5.5871000e+00  2.8247000e+00 -3.8042000e-01 /
                # CHEB /  1.6412000e+00  5.0500000e-01  2.1023000e-01 /
                # CHEB / -4.3607000e+00 -1.7857000e-01  5.1107000e-02 /
                # CHEB / -2.0344000e+00 -3.3644000e-01 -2.5833000e-02 /
                # CHEB / -7.4279000e-01 -2.1199000e-01 -4.0130000e-02 /

                r_line = self.react.formula[r]+txt_space + '1.0E0 0.0 0.0\n'
                k_line = ''
                for _l in range(len(self.react.txt[r])):
                    if 'Tmin' in self.react.txt[r][_l]:
                        k_line +=   'TCHEB / ' \
                                  + self.react.txt[r][_l].split('=')[1].split(',')[0] + ' '\
                                  + self.react.txt[r][_l].split('=')[2].split(',')[0] + ' / '
                    if 'Pmin' in self.react.txt[r][_l]:
                        k_line +=   'PCHEB / ' \
                                  + self.react.txt[r][_l].split('=(')[1].split(',')[0] + ' '\
                                  + self.react.txt[r][_l].split('=(')[2].split(',')[0] + ' /'

                k_line += "\nCHEB / "
                cheby_l = len(self.react.kin[r])
                cheby_c = len(self.react.kin[r][0])
                k_line += str(int(cheby_l)) + ' ' + str(int(cheby_c)) + ' / \n'

                for l in range(cheby_l):
                    k_line += 'CHEB / '
                    for c in range(cheby_c):
                        k_line += '%0.6e' %float(self.react.kin[r][l][c]) + ' '
                    k_line += '/\n'

                fd.write(r_line)
                fd.write(k_line)
                fd.write(comments)
#                fd.write(self.react.eff[r])



            fd.write(dup_txt)

        fd.write('\n END')
        timer.sleep(1.5)

        fd.close()



    def find_element(self,element,act_sp):
        for sp in range(len(self.spec.name)):
            if act_sp[sp]:
                if isinstance(self.spec.name[sp], str) \
                and element in self.spec.name[sp]:
                    return True
        return False

    def compare_red(self,mech_data_red,mp,verbose=0):

        self.spec.activ_m  = [False]*len(self.spec.name)       # act spec for methods loop
        self.react.activ_m = [False]*len(self.react.formula)   # act react for methods loop

        for sp in range(len(self.spec.name)):
            for sp_red in mech_data_red.spec.name:
                if self.spec.name[sp]==sp_red:
                    self.spec.activ_m[sp]=True
                    break

        r_red = 0
        for r in range(len(self.react.number)):
            for r_red in range(len(mech_data_red.react.number)):
                if self.react.number[r]==mech_data_red.react.number[r_red]:
                    self.react.activ_m[r]=True

                    self.react.kin[r] = copy.deepcopy(mech_data_red.react.kin[r_red])
                    break

        print_('\n\n Previous species reduction :  '+    \
               '%4.0f' %len(self.spec.activ_m)+ '  ->'+  \
               '%5.0f' %self.spec.activ_m.count(True),mp)
        print_('\n\n Previous reaction reduction : '+    \
               '%4.0f' %len(self.react.activ_m)+ '  ->'+ \
               '%5.0f' %self.react.activ_m.count(True),mp)



class Optim_param :
    def __init__(self, tsp, n_tspc, n_gen=20,n_ind=50, \
         selection_operator=1,      selection_options=[],   \
         Xover_operator=[1,2,3,4], Xover_pct = [15,30,30,30],   \
         mut_operator=[2,3,4], mut_pct=[15,30,30], mut_intensity=[30], mut_option=[3],\
         optim_on_meth='False', Arrh_max_variation=10,nb_r2opt=30):


        self.n_gen                = n_gen
        self.n_ind                = n_ind
        # GA options
        self.selection_operator   = selection_operator
        self.selection_options    = selection_options
        self.Xover_operator       = Xover_operator
        self.Xover_pct            = Xover_pct
        self.mut_operator         = mut_operator
        self.mut_pct              = mut_pct
        self.mut_intensity        = mut_intensity
        self.mut_option           = mut_option
        # PSO options
        self.inertia_score        = True
        self.inertia_min          = 0.2
        self.inertia_i            = 1.5
        self.inertia_end          = 0.4
        self.cognitive_accel_i    = 1.7
        self.cognitive_accel_end  = 0.9
        self.social_accel_i       = 0
        self.social_accel_end     = 2


        self.optim_on_meth        = optim_on_meth
        self.target_r             = []
        if type(Arrh_max_variation) is list:
            self.Arrh_max_variation      = Arrh_max_variation
        else:
            self.Arrh_max_variation      = [Arrh_max_variation]*3
        self.diff_r_meth          = False
        self.r2opt                = False
        self.tspc                 = tsp
        self.n_tspc               = n_tspc
        self.coeff_s              = [1]*len(tsp)
        self.coeff_T              = 1
        self.coeff_ig             = 1
        self.coeff_Sl             = 1
        self.coeff_K              = 1
        self.genVec               = []
        self.best_fitness         = []
        self.mean_fitness         = []
        self.worst_fitness         = []
        self.main_path            = ''
        self.exp_data             = False
        self.nb_r2opt             = nb_r2opt
        self.opt_subm_C           = [True]*30
        self.opt_subm_CO          = True
        self.opt_subm_N           = [True]*30
        self.opt_subm_S           = [True]*30
        self.opt_subm_Si          = [True]*30



    def count_Xover(self):
        Xover_num = []
        for Xover_pct in self.Xover_pct:
            Xover_num.append(int(np.ceil(Xover_pct*self.n_ind/100)))
        #check if the number of Xover is an even number
        for i in range(len(Xover_num)):
            if (Xover_num[i]%2)!=0:
                Xover_num[i] +=1
        self.Xover_num   = Xover_num
        self.total_Xover = int(np.sum(Xover_num))

    def count_mut(self):
        mut_num = []
        for mut_pct in self.mut_pct:
            mut_num.append(int(np.ceil(mut_pct*self.n_ind/100)))
        self.mut_num   = mut_num
        self.total_mut = int(np.sum(mut_num))




class Reactions:
    def __init__(self,mech,verbose=0):

        self.number    = []
        self.type      = []
        self.formula   = []
        self.kin       = []
        self.ref_kin   = []
        self.txt       = []
        self.eff       = []
        self.tbe       = [] # three body exception in falloff reactions (+ AR) / (+HE) ...
        self.activ_p   = []
        self.activ_m   = []
        self.modif     = []
        self.incert    = []



class Red_data :
    def __init__(self,gas_ref,mech_name,tspc,n_tspc,reduction_operator='DRG_sp',   \
                 optim = False,verbose=6,targetSpeciesIdx=[],gas_loop='gas_ref',gas_act='gas_ref'):

        self.gas_ref                        = gas_ref
        self.gas_loop                       = gas_ref
        self.gas_red                        = gas_ref
        self.verbose                        = verbose


#        self.verbose                        = verbose
        self.targetSpeciesIdx=[]
        while '' in tspc: tspc.remove('')
        for i in range(len(tspc)):
            	self.targetSpeciesIdx.append(gas_ref.species_index(tspc[i]))

        self.tspc                           = tspc
        self.n_tspc                         = n_tspc
        self.reduction_operator             = reduction_operator


        if gas_act=='gas_ref': self.gas_act = gas_ref
        else :                 self.gas_act = gas_act

        self.red_op = Red_operator(self.targetSpeciesIdx)

        self.optim = optim
        self.write_results = False
        self.sens_method = 'adjoint'


class Species:
    def __init__(self,verbose=0):

        self.name           = []
        self.atoms          = []
        self.thermo         = []
        self.trans          = []
        self.note           = []
        self.activ_p        = []
        self.activ_m        = []
        self.CSP_radicals   = []


class Simul_param :
    def __init__(self, pts_scatter = [0.0, 0.03, 0.3, 0.5, 0.7, 1.0],        \
                 end_sim = 0.02,                                             \
                 tol_ss = [1.0e-5, 1.0e-8],                                  \
                 tol_ts = [1.0e-4, 1.0e-8],                                  \
                 verbose = 0,                                                \
                 tign_nPoints = 450, tign_dt = 1e-09,                    \
                 n_pts = 250, delta_npts = 20, t_max_coeff = 5,              \
                 Scal_ref = 'H2O', grad_curv_ratio = 0.5):

        self.pts_scatter      = pts_scatter      # time stepping  or  grid
        self.end_sim          = end_sim          # tmax           or  xmax
        self.tol_ss           = tol_ss
        self.tol_ts           = tol_ts
        self.tign_nPoints     = tign_nPoints
        self.tign_dt          = tign_dt
        self.n_pts            = n_pts
        self.delta_npts       = delta_npts
        self.t_max_coeff      = t_max_coeff
        self.t_max_react      = 1
        self.Scal_ref         = Scal_ref
        self.grad_curv_ratio  = grad_curv_ratio
        self.mdot             = 1   # kg/m^2/s
        self.mdot2            = 3   # kg/m^2/s
        self.verbose          = verbose
        self.show_plots       = False
        self.write_ck         = True
        self.transport_model  = "Mix"
        self.slope_ff         = 0.05
        self.curve_ff         = 0.05
        self.ratio_ff         = 2.0
        self.prune_ff         = 0.01
        self.T_lim            = 800
        self.u_0              = .006   # PFR opt - inflow velocity [m/s]
        self.area             = 1.e-4  # PFR opt - cross-sectional area [m**2]
        self.PFR_auto_time    = True  # PFR opt - automatic time discretization
        self.par_ind          = False
        self.restore_flame    = False
        self.flame_res_folder = 'False'


class Sim_Results :
    def __init__(self,conditions,gas=[], pts_scatter=[],         \
                 T=[],P=[],conc=[],kf=[],kr=[],ign_time="",Sl=[]):
#        error_param=['QoI','max','max',True,True,True]

        self.conditions  = conditions
        self.gas         = gas
        self.pts_scatter = pts_scatter
        self.T           = T
        self.P           = P
        self.conc        = conc
        self.kf          = kf
        self.kr          = kr
        self.ign_time_hr = ign_time
        self.ign_time_sp = ign_time
        self.ign_time    = ign_time
        self.Sl          = Sl
        self.K_ext       = 0
        self.f           = False
#        self.error_param = Error_param(error_param)
        self.X           = 0
        self.Y           = 0
        self.conc2X()
        self.res_fname   = 'reduction_results.csv'


    def conc2X(self):
        self.X=[]
        for step in range(len(self.conc)):
            self.X.append([])
            n_tot = sum(self.conc[step])
            for sp in range(len(self.conc[step])):
                if n_tot!=0:
                    self.X[step].append(self.conc[step][sp]/n_tot)
                else:
                    self.X[step].append(0)

    def conc2Y(self):
        gas = self.conditions.composition.gas_ref
        self.Y=[]
        for step in range(len(self.conc)):
            self.Y.append([])
            n_tot = sum(self.conc[step])
            mean_molecular_weight = 0
            for sp in range(len(self.conc[step])):
                M_sp = gas.molecular_weights[sp]
                mean_molecular_weight += self.X[step][sp]*M_sp
            for sp in range(len(self.conc[step])):
                M_sp = gas.molecular_weights[sp]
                if n_tot!=0:
                    self.Y[step].append((self.conc[step][sp]*M_sp)/(n_tot*mean_molecular_weight))
                else:
                    self.Y[step].append(0)

    def X2conc(self):
        conc=[]
#        if sum(self.X[0])>0.001:
        for step in range(len(self.X)):
            if type(self.P) is float: P = self.P
            else: P = self.P[step]
            try:
                if type(self.T[step]) is str: T = 300
                else: T =  self.T[step]
            except:
                T = 300
            ntot_V = P/(8.314*T)/1000  # kmol/m3
            conc.append([])
            for idx in range(len(self.X[step])):
                conc[step].append(self.X[step][idx]*ntot_V)
        self.conc=conc


    def plotData_opt(self,spec2plot):

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
             (0, (1, 10))]                  #'loosely dotted'                   12
        scatter_styles = ["o", "^", "*",  "+", "s", "d", "v", "x", "p", "."]
        colors_styles = ["blue","red","green","grey","purple","black", "darkturkoise",\
                         "sandybrown","saddlebrown","magenta"]

        fig, ax = plt.subplots()
        if "reactor" in self.conditions.config:
            ax.set_xlabel('t (s)')
        elif "flame" in self.conditions.config:
            ax.set_xlabel('Z (mm)')

        sp2plt_idx=[]
        for sp in spec2plot:
            sp2plt_idx.append(self.conditions.composition.gas.species_index(sp))

        for i in range(len(sp2plt_idx)):
            Xi=[];idx=sp2plt_idx[i]
            for step_X in self.X:
                Xi.append(step_X[idx])
            ax.plot(self.pts_scatter, Xi, \
                           linestyle=linestyles[i], color=colors_styles[i], \
                           linewidth=2, label=spec2plot[i])
        lines, labels = ax.get_legend_handles_labels()
        ax.legend(lines, labels, loc=0)

        plt.show()



    def write_mech_info(self,num_meth=""):

        fichier_data = open(str(num_meth)+'_'+self.res_fname, 'w')

        if self.conditions.simul_param.verbose>3: print("Writting csv file...")
        fichier_data.write("Reference mechanism;"+self.conditions.mech+\
        "\nNumber of species:;"+str(self.conditions.composition.gas.n_species)\
        +"\nNumber of reactions:;"+str(self.conditions.composition.gas.n_reactions)+"\n")

    def write_method(self,clock,method,num_meth="",\
                     fitness=False,clock_opt=False,optim_param=False):
        fichier_data = open(str(num_meth)+'_'+self.res_fname, 'a')
        fichier_data.write(method+'\n')
        fichier_data.write(clock.text())
        if fitness:
            fichier_data.write(clock_opt.text())
            fichier_data.write('Fitness : '+'%0.2f' %(fitness)+';')
            fichier_data.write(self.conditions.error_param.error_calculation+\
                ' '+self.conditions.error_param.error_coupling+\
                ' '+self.conditions.error_param.error_type_fit+';;')
            if optim_param:
                txt = 'generation:'
                for gen in range(len(optim_param.genVec)):
                    txt+=';'+str(optim_param.genVec[gen])
                txt+='\n;;;best fitness:'
                for gen in range(len(optim_param.genVec)):
                    txt+=';'+str(optim_param.best_fitness[gen])
                txt+='\n;;;mean fitness:'
                for gen in range(len(optim_param.genVec)):
                    txt+=';'+str(optim_param.mean_fitness[gen])
                txt+='\n'
                fichier_data.write(txt)


    def write_case_data(self,step,num_meth="",nb_case=False,errors=False):

        fichier_data = open(str(num_meth)+'_'+self.res_fname, 'a')

        gas = self.conditions.composition.gas_ref

        if nb_case:
            fichier_data.write("Case nb: "+str(nb_case)+"\n")
            l1_conditions = "config;fuel;oxidant;diluent;phi;diluent_ratio;mixt;"\
                +"P(Pa)"
            if "free_flame" in self.conditions.config:
                l1_conditions += ";Ti(K);rtol_ts;atol_ts;rtol_ss;atol_ss;transport_model;Sl0(m/s)"
            if "burner_flame" in self.conditions.config:
                l1_conditions += ";Ti(K);rtol_ts;atol_ts;rtol_ss;atol_ss;transport_model;mdot(kg/m2/s)"
            if "tp_flame" in self.conditions.config:
                l1_conditions += ";Ti(K);rtol_ts;atol_ts;rtol_ss;atol_ss;transport_model"
            if 'diff_' in self.conditions.config\
            or 'pp_'   in self.conditions.config:
                l1_conditions += ";mdot;fuel2;oxidant2;diluent2;phi2;diluent_ratio2;mixt2;mdot2;Ti(K);rtol_ts;atol_ts;rtol_ss;atol_ss;transport_model;K_max(1/s)"
            elif "reactor" in self.conditions.config:
                l1_conditions += ";Ti(K);rtol_ts;atol_ts;ig_time(s)"
            elif "JSR" in self.conditions.config:
                l1_conditions += ";time(s);rtol_ts;atol_ts"
            elif "PFR" in self.conditions.config:
                l1_conditions += ";Ti(K);rtol_ts;atol_ts;length(m);u_0(m/s);area(m**2)"
            l2_conditions = "\n"+self.conditions.config+";"\
                    +self.conditions.composition.fuel+";"\
                    +str(self.conditions.composition.oxidant)+";"\
                    +self.conditions.composition.diluent+";"\
                    +str(self.conditions.composition.phi)+";"\
                    +str(self.conditions.composition.diluent_ratio)+";"\
                    +str(self.conditions.composition.X)+";"\
                    +str(self.conditions.state_var.P)+";"
            if 'diff_' in self.conditions.config\
            or 'pp_'   in self.conditions.config:
                l2_conditions += str(self.conditions.simul_param.mdot)+";"\
                                +str(self.conditions.composition.fuel2)+";"\
                                +self.conditions.composition.oxidant2+";"\
                                +self.conditions.composition.diluent2+";"\
                                +str(self.conditions.composition.phi2)+";"\
                                +str(self.conditions.composition.diluent_ratio2)+";"\
                                +str(self.conditions.composition.X2)+";"\
                                +str(self.conditions.simul_param.mdot2)+";"
            if "free_flame" in self.conditions.config\
            or 'diff_' in self.conditions.config\
            or 'tp_'   in self.conditions.config\
            or 'pp_'   in self.conditions.config:
                l2_conditions+=str(self.conditions.state_var.T)+";"\
                        +str(self.conditions.simul_param.tol_ts[0])+";"\
                        +str(self.conditions.simul_param.tol_ts[1])+";"\
                        +str(self.conditions.simul_param.tol_ss[0])+";"\
                        +str(self.conditions.simul_param.tol_ss[1])+";"\
                        +self.conditions.simul_param.transport_model+";"
                if "free_" in self.conditions.config:
                    l2_conditions+=str(self.Sl)
                if "diff_" in self.conditions.config or "pp_" in self.conditions.config:
                    l2_conditions+=str(self.K_max)
                if "burner" in self.conditions.config:
                    l2_conditions+=str(self.conditions.simul_param.mdot)
            elif "reactor" in self.conditions.config:
                l2_conditions+=str(self.conditions.state_var.T)+";"\
                        +str(self.conditions.simul_param.tol_ts[0])+";"\
                        +str(self.conditions.simul_param.tol_ts[1])+";"
                if self.ign_time == '':
                    if self.ign_time_hr != '' and self.ign_time_hr != False:
                        l2_conditions+=str(self.ign_time_hr)
                    else :
                        l2_conditions+=str(self.ign_time_sp)
                else:
                    l2_conditions+=str(self.ign_time)
            elif "JSR" in self.conditions.config:
                l2_conditions+=str(self.conditions.simul_param.end_sim)+";"\
                        +str(self.conditions.simul_param.tol_ts[0])+";"\
                        +str(self.conditions.simul_param.tol_ts[1])
            elif "PFR" in self.conditions.config:
                l2_conditions+=str(self.conditions.state_var.T)+";"\
                        +str(self.conditions.simul_param.tol_ts[0])+";"\
                        +str(self.conditions.simul_param.tol_ts[1])+";"\
                        +str(self.conditions.simul_param.end_sim)+";"\
                        +str(self.conditions.simul_param.u_0)+";"\
                        +str(self.conditions.simul_param.area)
            fichier_data.write(l1_conditions)
            fichier_data.write(l2_conditions)

        fichier_data.write("\n* Step: "+step)

        if errors:
            txt_error='\n'
            if "JSR" not in self.conditions.config:
                if errors.qoi_T!=0:
                    txt_error+='Temperature error: '+'%0.1f'%(errors.qoi_T*100)+'%\n'
            if "reactor" in self.conditions.config:
                if errors.qoi_ig!=0:
                    txt_error+='Ignition delay time error: '+'%2.1f'%(errors.qoi_ig*100)+'%\n'
            elif "free_flame" in self.conditions.config:
                if errors.qoi_Sl!=0:
                    txt_error+='Flame speed error: '+'%0.1f' %(errors.qoi_Sl*100)+'%\n'
            elif "diff_flame" in self.conditions.config:
                if errors.qoi_K!=0:
                    txt_error+='Extinction strain rate error: '+'%0.1f' %(errors.qoi_K*100)+'%\n'
            txt_error+='Target species errors:;'
            for sp in range(self.conditions.error_param.n_tspc):
                    txt_error+=self.conditions.error_param.tspc[sp]\
                              +' error: '+'%0.1f' %(errors.qoi_s[sp]*100)+'%;'
            fichier_data.write(txt_error)

        if "reactor" in self.conditions.config:
            fichier_data.write("\nIgnition delay time(s):;")
            if self.ign_time == '':
                if self.ign_time_hr != '' and self.ign_time_hr != False:
                    fichier_data.write(str(self.ign_time_hr))
                else :
                    fichier_data.write(str(self.ign_time_sp))
            else:
                fichier_data.write(str(self.ign_time))
        elif "free_flame" in self.conditions.config:
            fichier_data.write("\nSl0(cm/s):;"+str(self.Sl*100))
        elif "diff_" in self.conditions.config or "pp_" in self.conditions.config:
            if self.K_ext!=0:
                fichier_data.write("\nK_ext(1/s):;"+str(self.K_ext))

        # headers
        if "reactor" in self.conditions.config:
            fichier_data.write("\nTime(s);T(K)")
        elif "JSR" in self.conditions.config:
            fichier_data.write("\nTi(K);Tf(K)")
        elif "flame" in self.conditions.config:
            fichier_data.write("\nZ(m);T(K)")
        elif "PFR" in self.conditions.config:
            fichier_data.write("\nt(s);Z(m);T(K)")
        for n_sp in range(gas.n_species):
##            # ----------------------------------
            #if gas.species_name(n_sp)=='nc7h16' \
            #or gas.species_name(n_sp)=='co' \
            #or gas.species_name(n_sp)=='co2' \
            #or gas.species_name(n_sp)=='h2' \
            #or gas.species_name(n_sp)=='o2' \
            #or gas.species_name(n_sp)=='n2' \
            #or gas.species_name(n_sp)=='c2h4':
                #fichier_data.write(";"+str(gas.species_name(n_sp)))
##            # ----------------------------------
            fichier_data.write(";"+str(gas.species_name(n_sp)))

        fichier_data.write("\n")

        # ----------------------------------
        self.conc2Y()
        # ----------------------------------

        # data
        for nb_line in range(len(self.pts_scatter)):
            # time / z / Ti
            fichier_data.write(str(self.pts_scatter[nb_line])+";")
            # for PFR only : specify the position at the time t
            if "PFR" in self.conditions.config:
                fichier_data.write(str(self.z1[nb_line])+";")
            # temperature
            fichier_data.write(str(self.T[nb_line]))
            for n_sp in range(gas.n_species):
                ## ----------------------------------
                ## concentration
                #if gas.species_name(n_sp)=='nc7h16' \
                #or gas.species_name(n_sp)=='co' \
                #or gas.species_name(n_sp)=='co2' \
                #or gas.species_name(n_sp)=='h2' \
                #or gas.species_name(n_sp)=='o2' \
                #or gas.species_name(n_sp)=='n2' \
                #or gas.species_name(n_sp)=='c2h4':
                    #fichier_data.write(";"+str(self.X[nb_line][n_sp]))
            # ----------------------------------
                fichier_data.write(";"+str(self.X[nb_line][n_sp]))
            fichier_data.write("\n")





def searchNearest(data, search_value, start =0, end_ind=-1) :
    vect = list(abs(np.array(data) - search_value))
    index = vect.index(np.amin(vect[start:end_ind]))
    value = data[index]

    return value, index


def suppress_stdout():
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
#        try:
#            yield
#        finally:
#            sys.stdout = old_stdout
    return old_stdout

def restore_stdout(old_stdout):
    sys.stdout = old_stdout

def print_(text,main_path,file_n='red_info.txt'):
    print(text)
    cur_path = os.getcwd()
    if main_path!='': os.chdir(main_path)
    fichier_data = open(file_n, 'a')
    fichier_data.write(text+'\n')
    os.chdir(cur_path)


class ProgressBar:
    """
    This class allows you to make easily a progress bar.
    """
    def __init__ (self, valmax, title='', maxbar=50):
        if valmax == 0:  valmax = 1
        if maxbar > 200: maxbar = 200
        self.valmax = valmax
        self.maxbar = maxbar
        self.title  = title

    def update(self, val,val2=''):

        # format
        if val > self.valmax: val = self.valmax

        # process
        perc  = np.round((float(val) / float(self.valmax)) * 100)
        scale = 100.0 / float(self.maxbar)
        bar   = int(perc / scale)

        # render
        out = '\r%2s %s [%s%s] %3d %% %s' % (self.title, str(val2),'=' * bar, ' ' * (self.maxbar - bar), perc,'      ')

        try:
            get_ipython().profile
            sys.stdout.write(out)
            #Extinction du curseur
            # not on Mac
            if platform.system()!='Darwin':
                os.system('setterm -cursor off')
            #Rafraichissement de la barre
            sys.stdout.flush()
        except Exception:
            pass



def get_screen_size():
    from PyQt5 import QtWidgets

    _height = QtWidgets.QDesktopWidget().screenGeometry(-1).height()

    return _height,_height
