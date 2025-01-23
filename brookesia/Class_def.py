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
import warnings
import re

warnings.filterwarnings("ignore") # suppress warnings with exponential calculation
warnings.filterwarnings("ignore", category=DeprecationWarning)


global version


#from operator import xor


class State_var :
    def __init__(self,T=300,P=101325):
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

                X = self.oxidant2    + ':' + str(X_oxidant)  + ',' + \
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
        
        
        path_lm = '/dev/null' ; path_w = 'nul'
        # # --------------------------------------------------------------------------------
        # # interpretation of the new mech
        try:
            # with open(path_lm, 'w') as fnull:  
            #     with redirect_stdout(fnull):
            try:
                gas = get_gas_ct('_kinetic_mech/'+mech)
            except:
                gas = get_gas_ct(mech)
        except:
            print('\n\n\n ! ! ! ! !\n\nmech: '+mech+' not found\n\n\n\n\n')
#
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
        self.fit_coeff     = 1
        self.write_conc  = True



class Red_operator :
    def __init__(self,target_species,optimization=False,n_points=20,          \
                 max_error_sp=30, max_error_T=10, max_error_ig = 10,          \
                 max_error_Sl=10, max_error_K=20, inter_sp_inter = True,                      \
                 eps_init=0.06,delta_eps_init=0.005,eps_max=2.0,              \
                 r_withdraw_intensity = 10,\
#                 sp_interaction_coeffs=[],r_interaction_coeffs=[]\
                 atol_ts = [1e-4, 1e-6]):

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
        self.write_results      = False


        # drg
        self.graph_search       = 'Dijkstra'  # Dijkstra DFS
        self.new_targets_4_DRG_r= copy.deepcopy(target_species)

        # sa
        self.rtol_ts            = 1e-6
        self.atol_ts            = 1e-9
        self.sens_method        = 'adjoint' # adjoint / brute_force
        self.sensi_Sl           = False
        self.sensi_T            = False
        self.sensi_igt          = False

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
            except: self.qoi_ig = 1
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
        gas_ref           = conditions.composition.gas_ref
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
            #conc_ori          = ref_results.opt_refind_conc
        else:                                      # 'molar_fraction'
            conc_ref          = ref_results.X
            conc_red          = red_results.X
            #conc_ori          = ref_results.opt_refind_X
        n_points          = len(pts_scatter)
        QoI = []

        if error_calculation =="points":

            for k in range(n_tspc):

                index = gas_ref.species_index(tspc[k])

                data1 = np.zeros(n_points)
                data2 = np.zeros(n_points)
                #data_ori = np.zeros(n_points)
                for i in range(n_points):
                    data1[i] = conc_ref[i][index]
                    data2[i] = conc_red[i][index]
                    #data_ori[i] = conc_ori[i][index]
                if sum(data1)>0:      # absence of data if experimental optimization
                    sumDiff, sumDiff_ori = 0, 0
                    if error_type == "all":
                        sumDiff = []
                        for j in range(len(pts_scatter)):
                            sumDiff.append(abs(data1[j] - data2[j])/np.amax(data1))
                    elif error_type == "mean":
                        sumdata_ref,sumdata_red=0,0
                        for j in range(len(pts_scatter)):
                            if j==0:
                                sumDiff     += abs(data1[j] - data2[j])*.5*(pts_scatter[j+1]-pts_scatter[j])
                                sumdata_ref += abs(data1[j])*.5*(pts_scatter[j+1]-pts_scatter[j])
                                sumdata_red += abs(data2[j])*.5*(pts_scatter[j+1]-pts_scatter[j])
                                # sumDiff_ori += abs(data1[j] - data_ori[j])*.5*(pts_scatter[j+1]-pts_scatter[j])
                            elif j==(len(pts_scatter)-1):
                                sumDiff     += abs(data1[j] - data2[j])*.5*(pts_scatter[j]-pts_scatter[j-1])
                                sumdata_ref += abs(data1[j])*.5*(pts_scatter[j]-pts_scatter[j-1])
                                sumdata_red += abs(data2[j])*.5*(pts_scatter[j]-pts_scatter[j-1])
                                # sumDiff_ori += abs(data1[j] - data_ori[j])*.5*(pts_scatter[j]-pts_scatter[j-1])
                            else:
                                sumDiff     += abs(data1[j] - data2[j])*.5*(pts_scatter[j+1]-pts_scatter[j-1])
                                sumdata_ref += abs(data1[j])*.5*(pts_scatter[j+1]-pts_scatter[j-1])
                                sumdata_red += abs(data2[j])*.5*(pts_scatter[j+1]-pts_scatter[j-1])
                                # sumDiff_ori += abs(data1[j] - data_ori[j])*.5*(pts_scatter[j+1]-pts_scatter[j-1])
                        sumDiff=(sumDiff)/max(sumdata_ref,sumdata_red)
                    elif error_type == "max":
                        diff = abs(np.array([data1]) - np.array([data2]))/max(max(np.array([data1])),max(np.array([data2])))
                        sumDiff = np.amax(diff)
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
        # T_ori             = ref_results.opt_refind_T

        n_points = len(pts_scatter)

        T_data = True
        if 'False' in T_ref:     T_data = False
        elif np.min(T_ref)<273:  T_data = False

        if error_calculation=="points" and T_data:

            # ========   2 - Temperature error   ==========
            sumDiff, sumDiff_ori = 0, 0
            if error_type == "all" :
                sumDiff = []
                for j in range(len(pts_scatter)):
                    sumDiff.append(abs(T_ref[j] - T_red[j])/max(T_ref[j]))
            elif error_type == "mean":
                sumdata_ref,sumdata_red=0,0
                for j in range(len(pts_scatter)):
                    if j==0:
                        sumDiff     += abs(T_ref[j] - T_red[j])*.5*(pts_scatter[j+1]-pts_scatter[j])
                        sumdata_ref += abs(T_ref[j])*.5*(pts_scatter[j+1]-pts_scatter[j])
                        sumdata_red += abs(T_red[j])*.5*(pts_scatter[j+1]-pts_scatter[j])
                        # sumDiff_ori += abs(T_ref[j] - T_ori[j])*.5*(pts_scatter[j+1]-pts_scatter[j])
                    elif j==(len(pts_scatter)-1):
                        sumDiff     += abs(T_ref[j] - T_red[j])*.5*(pts_scatter[j]-pts_scatter[j-1])
                        sumdata_ref += abs(T_ref[j])*.5*(pts_scatter[j]-pts_scatter[j-1])
                        sumdata_red += abs(T_red[j])*.5*(pts_scatter[j]-pts_scatter[j-1])
                        # sumDiff_ori += abs(T_ref[j] - T_ori[j])*.5*(pts_scatter[j]-pts_scatter[j-1])
                    else:
                        sumDiff     += abs(T_ref[j] - T_red[j])*.5*(pts_scatter[j+1]-pts_scatter[j-1])
                        sumdata_ref += abs(T_ref[j])*.5*(pts_scatter[j+1]-pts_scatter[j-1])
                        sumdata_red += abs(T_red[j])*.5*(pts_scatter[j+1]-pts_scatter[j-1])
                        # sumDiff_ori += abs(T_ref[j] - T_ori[j])*.5*(pts_scatter[j+1]-pts_scatter[j-1])
                sumDiff=(sumDiff)/max(sumdata_ref,sumdata_red)
            elif error_type == "max":
                diff = abs(np.array([T_ref]) - np.array([T_red]))/max(max(np.array([T_ref])),max(np.array([T_red])))                
                sumDiff = np.amax(diff)
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
            QoI = abs((ref_results.Sl - red_results.Sl)/max(ref_results.Sl,red_results.Sl))
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

        if ref_results.ign_time_hr == False and ref_results.ign_time_sp ==False:
            diff = False
            if verbose >=6 :
                print_('      No reference ignition time',mp)            
        else:
            if ref_results.ign_time_hr == False or red_results.ign_time_hr == False:
                ref_results.ign_time = ref_results.ign_time_sp
                red_results.ign_time = red_results.ign_time_sp
            else:
                ref_results.ign_time = ref_results.ign_time_hr
                red_results.ign_time = red_results.ign_time_hr
            diff=abs((np.log10(ref_results.ign_time)-np.log10(red_results.ign_time))\
                     /max(abs(np.log10(ref_results.ign_time)),abs(np.log10(red_results.ign_time))))
            # Display QoI
            if verbose >=6 :
                print_('      Error ignition time (in log10 scale): ' + '%0.1f' %(diff*100) +'% \n ',mp)

        return diff



    def er_estim_K(self,conditions, ref_results, red_results):
        # main variables
        mp = conditions.main_path
        verbose           = conditions.simul_param.verbose

        diff=abs((ref_results.K_ext-red_results.K_ext)/max(ref_results.K_ext,red_results.K_ext))

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

            I25p = abs( pts_scatter[ind25_red]- pts_scatter[ind25_ref]) / max(pts_scatter[ind25_ref],pts_scatter[ind25_red])


#            ## QoI : time at 75% of max data
            D75 = Dmin+(np.array(Dmax)-np.array(Dmin))*75/100
            D75_ref, ind75_ref = searchNearest( ref_var, D75[0], ind25_ref)
            D75_red, ind75_red = searchNearest( red_var, D75[1], ind25_red)
#            I75p = abs( pts_scatter[ind75_red]- pts_scatter[ind75_ref])/ pts_scatter[ind75_ref]

            ## QoI : maximum value comparison
            Ipa = abs(Dmax[1]-Dmax[0])/(max(Dmax[0],Dmax[1]))

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
                Ipr = abs(rate_red-rate_red)/max(rate_red,rate_red)
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
    #              o.       `:**   <-- D75c             |
    #              o.         `so                       |
    #              o.           h-                      |
    #              o.           -h                      |  --> Dc
    #              o.            os`                    |
    #              o.       <---->+**                   |                        <-- R25c         \ (deprecated)
    #              o.          D   -y                   |
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
                D25c = abs( pts_scatter[ind25_red]- pts_scatter[ind25_ref]) / max(pts_scatter[ind25_ref], pts_scatter[ind25_red])
            elif  pts_scatter[ind25_ref] ==  pts_scatter[ind25_red] :
                D25c = 0
            else:
                D25c = 1


            ## consumption amplitude
            Dc_ref = Dmax[0]-Dmin[0]
            Dc_red = Dmax[1]-Dmin[1]
            Dca = abs((Dc_ref-Dc_red))/max(Dc_ref,Dc_red)*2


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
                Ipr = abs(rate_red-rate_red)/max(abs(rate_red),abs(rate_red))
            if D_ref != 0:
                Dcr = abs((D_ref-D_red)/max(D_ref,D_red))
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
    #        y                    sy.   `/h.                       |
    #        y                   y+       .d`                      |
    #        y                  :h         /o                      |
    #        y                 .d`         `m                      |  --> D50_mf
    #        y                 h-           s/   <--- B75d         |
    #        y       R25p --> /y            `yo-                   |
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
            B25p = abs( pts_scatter[ind25_red]- pts_scatter[ind25_ref]) / max(pts_scatter[ind25_ref],pts_scatter[ind25_red])


            ## Bpr : rate of increase between 25% and the pic
            rate_red = (Dmax_red-D25_red)/(pts_scatter[indmax_red]-pts_scatter[ind25_red])
            rate_ref = (Dmax_ref-D25_ref)/(pts_scatter[indmax_ref]-pts_scatter[ind25_ref])
            Bpr = abs(rate_red-rate_ref)/max(rate_ref,rate_red)
            Bm = abs(Dmax[1]-Dmax[0])/max(Dmax[0],Dmax[1])


            ## R25f : time at 25% of the diff between the concentration pic and final value
            D25 =  np.array(Dend)+(np.array(Dmax)-np.array(Dend))*25/100
            D25f_ref, ind25f_ref = searchNearest( ref_var, D25[0], indmax_ref)
            D25f_red, ind25f_red = searchNearest( red_var, D25[1], indmax_red)

            rate_red = (D25f_red-Dmax_red)/(pts_scatter[ind25f_red]-pts_scatter[indmax_red])
            rate_ref = (D25f_ref-Dmax_ref) /(pts_scatter[ind25f_ref]-pts_scatter[indmax_ref])
            B75d = abs(rate_red-rate_ref)/max(abs(rate_ref),abs(rate_red))


            ## D50_mf_ref : Diff between the ref maximum and final value of each case
            D50_mf_ref = Dmax[0]-Dend[0]
            D50_mf_red = Dmax[0]-Dend[1]
            Bca = abs(D50_mf_ref-D50_mf_red)/max(D50_mf_ref,D50_mf_red)

            QoI = [B25p, Bpr, Bm, B75d, Bca]

        # ========================================================================
        elif curve_type == "bell (inv)":


    #        h:==========//+o+`
    #        y              :h.
    #        y               -d`              .-::::::://-        ____
    #        y       R25p --> /y            `yo-                   |
    #        y                 h-           s/   <--- R25f         |
    #        y                 .d`         `m                      |  --> D50_mf
    #        y                  :h         /o                      |
    #        y                   y+       .d`                      |
    #        y                    sy.   `/h.                       |
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
            Bi_cr = abs((pts_scatter[ind25_red]-pts_scatter[indmin_red])-(pts_scatter[ind25_ref]-pts_scatter[indmin_ref]))\
                   /max(abs(pts_scatter[ind25_ref]-pts_scatter[indmin_ref]),abs(pts_scatter[ind25_red]-pts_scatter[indmin_red]))

            ## RM_mf : value of the molar fraction at the pic
            Bi_m = abs(Dmin[1]-Dmin[0])/max(Dmin[0],Dmin[1])


            # R25p : time at25% of max value (increasing part)
            D25i = Dinit-(np.array(Dinit)-np.array(Dmin))*25/100
            D25_ref, ind25_ref = searchNearest(ref_var, D25i[0], 0, indmin_ref)
            D25_red, ind25_red = searchNearest(red_var, D25i[1], 0, indmin_red)
            if ind25_ref==0: ind25_ref=1 # 0 exception
            Bi_25c = abs(pts_scatter[ind25_red]- pts_scatter[ind25_ref])/ max(pts_scatter[ind25_ref],pts_scatter[ind25_red])



            ## R25f : time at 25% of the diff between the concentration pic and final value
            D25 =  np.array(Dend)-(np.array(Dend)-np.array(Dmin))*25/100
            D25f_ref, ind25f_ref = searchNearest( ref_var, D25[0], indmin_ref)
            D25f_red, ind25f_red = searchNearest( red_var, D25[1], indmin_red)
            Bi_pr = abs((pts_scatter[ind25f_red]-pts_scatter[indmin_red])-(pts_scatter[ind25f_ref]-pts_scatter[indmin_ref]))\
                        /max(abs(pts_scatter[ind25f_ref]-pts_scatter[indmin_ref]),abs(pts_scatter[ind25f_red]-pts_scatter[indmin_red]))

            ## D50_mf_ref : Diff between the ref maximum and final value of each case
            D50_mf_ref = Dend[0]-Dmin[0]
            D50_mf_red = Dend[0]-Dmin[1]
            Bi_pa = abs(D50_mf_ref-D50_mf_red)/max(abs(D50_mf_ref),abs(D50_mf_red))

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
                elif np.amax(qoi_s[i])!=np.amax(qoi_s[i]): # check nan
                    under_tol_s.append(False)
                else:
                    under_tol_s.append(True)

        # Temperature
        under_tol_T = True
        if qoi_T:
            if conditions.error_param.T_check:
                if qoi_T*100>red_data_meth.max_error_T:
                    under_tol_T = False
                elif qoi_T != qoi_T: # check nan
                    under_tol_T = False


        # Flame speed
        under_tol_Sl = True
        if qoi_Sl:
            if 'flame' in conditions.config and conditions.error_param.Sl_check:
                if qoi_Sl*100>red_data_meth.max_error_Sl:
                    under_tol_Sl = False
                elif qoi_Sl != qoi_Sl: # check nan
                    under_tol_Sl = False


        # Ignition time
        under_tol_ig = True
        if qoi_ig:
            if 'reactor' in conditions.config and conditions.error_param.ig_check:
                if qoi_ig*100>red_data_meth.max_error_ig:
                    under_tol_ig = False
                elif qoi_ig != qoi_ig: # check nan
                    under_tol_ig = False

        # Extinction strain rate
        under_tol_K = True
        if qoi_K:
            if 'diff_flame' in conditions.config and conditions.error_param.K_check:
                if qoi_K*100>red_data_meth.max_error_K:
                    under_tol_K = False
                elif qoi_K != qoi_K: # check nan
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
        if ".cti" in mech:
            self.get_data_cti(mech,verbose)
        elif ".yaml" in mech:
            self.get_data_yaml(mech,verbose)
        self.keep4opt = False

    def get_data_cti(self,mech,verbose=0):

        fs = open(mech, 'r')
        self.duplicate_list = []


# =============================================================================
# Gas properties
# =============================================================================
        txt = fs.readline()

        self.description = {'description':'', 'cantera_version':''}
        self.mech_units = {'length':'cm', 'time':'s', 'quantity':'mol', 'act_energy':'cal/mol'}

        while "units(" not in txt:
            if txt == "":
                print("problem while reading the units section in the mechanism " + mech)
                break
            self.description['description'] += txt
            txt = fs.readline()
        if self.description['description'].replace(" ","").replace("\n","") == "":
            self.description['description'] = False

        if 'length' in txt:
            self.mech_units['length'] = txt.split("length")[1].split("'")[1]
        if 'time' in txt:
            self.mech_units['time'] = txt.split("time")[1].split("'")[1]
        if 'quantity' in txt:
            self.mech_units['quantity'] = txt.split("quantity")[1].split("'")[1]
        if 'activation-energy' in txt:
            self.mech_units['act_energy'] = txt.split("act_energy")[1].split("'")[1]

        self.gas_prop = {'name':'gas', 'thermo':'ideal-gas', 'elements':False, \
                         'species':[], 'kinetics':'gas', 'reactions':'all',\
                         'transport':False, 'state_T':'300', 'state_p':'101325'}

        while "Species data" not in txt:
            if "name=" in txt:
                self.gas_prop['name'] = txt.split('name')[1].split("'")[1]
            if "gas(" in txt:
                self.gas_prop['thermo'] = txt.split('(')[0]
            elif "elements=" in txt:
                if '"' in txt:
                    self.gas_prop['elements'] = txt.split('elements')[1].split('"')[1].split(",")
                if "'" in txt:
                    self.gas_prop['elements'] = txt.split('elements')[1].split("'")[1].split(",")
            elif "species=" in txt:
                txt.replace('species="""', '')
                while '=' not in txt:
                    specs = txt.split(' ')
                    for spec_i in specs:
                        if spec_i != '':
                            self.gas_prop['species'].append(spec_i)
            # elif "- kinetics:" in txt:
            #     self.gas_prop['kinetics'] = txt.split('kinetics')[1].split("'")[1]
            elif "reactions=" in txt:
                self.gas_prop['reactions'] = txt.split('reactions')[1].split("'")[1]
            elif "transport=" in txt:
                self.gas_prop['transport'] = txt.split('transport')[1].split("'")[1]
            elif "state" in txt:
                self.gas_prop['state_T'] = txt.split('temperature=')[1].split(",")[0]
                self.gas_prop['state_p'] = txt.split('pressure=')[1].split(")")[0]

            txt = fs.readline()

# =============================================================================
# Thermo / transport
# =============================================================================

        while "Reaction data" not in txt and txt != "":
            if txt == "":
                print("problem during the new mechanism writing")
                break
            elif "species(" in txt:
                self.spec.name.append([])
                self.spec.atoms.append([])
                self.spec.thermo_temp.append([])
                self.spec.thermo_coeff_lT.append([])
                self.spec.thermo_coeff_hT.append([])
                self.spec.thermo_model.append(False)
                self.spec.trans_model.append(False)
                self.spec.trans_geom.append(False)
                self.spec.trans_wd.append(False)
                self.spec.trans_diam.append(False)
                self.spec.trans_dipole.append(False)
                self.spec.trans_polar.append(False)
                self.spec.trans_rot_relax.append(False)
                self.spec.note_therm.append(False)
                self.spec.note_trans.append(False)


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
                        self.spec.atoms[-1] = self.spec.atoms[-1].replace(' ',',')
                    if "thermo" in txt:
                        if not ("NASA" in txt):
                            txt = fs.readline()
                        if 'NASA9' in txt:
                            self.spec.thermo_model[-1] = 'NASA9'
                        else:
                            self.spec.thermo_model[-1] = 'NASA7'

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
                        self.spec.thermo_temp[-1]     = Low_T + [Hight_T[1]]
                        self.spec.thermo_coeff_lT[-1] = L_T_coeffs
                        self.spec.thermo_coeff_hT[-1] = H_T_coeffs

                        # self.spec.thermo[-1]=Low_T+L_T_coeffs+Hight_T+H_T_coeffs
                    if "transport" in txt:
                        self.spec.trans_model[-1] = 'gas'
                        if 'geom' in txt:
                            self.spec.trans_geom[-1] = txt.split("geom='")[1].split("',")[0].split("')")[0]
                        while ")" not in txt:
                            txt = fs.readline()
                            #self.spec.trans[-1].append(txt)
                            if 'geom' in txt:
                                self.spec.trans_geom[-1] = txt.split("geom='")[1].split("',")[0].split("')")[0]
                            if 'diam=' in txt:
                                self.spec.trans_diam[-1] = txt.split("diam=")[1].split(",")[0].split(")")[0]
                            if 'well_depth' in txt:
                                self.spec.trans_wd[-1] = txt.split("well_depth=")[1].split(",")[0].split(")")[0]
                            if 'dipole' in txt:
                                self.spec.trans_dipole[-1] = txt.split("dipole=")[1].split(",")[0].split(")")[0]
                            if 'polar' in txt:
                                self.spec.trans_polar[-1] = txt.split("polar=")[1].split(",")[0].split(")")[0]
                            if 'rot_relax' in txt:
                                self.spec.trans_rot_relax[-1] = txt.split("rot_relax=")[1].split(",")[0].split(")")[0]


                    if "note" in txt:
                        self.spec.note_trans[-1] = ""
                        while ")" not in txt:
                            if 'note=' in txt:
                                self.spec.note_trans[-1] += txt.split('note=')[1]
                            else:
                                self.spec.note_trans[-1] += txt.split(')')[0]
                            txt = fs.readline()
                        self.spec.note_trans[-1] += txt.split('note=')[1].split(')')[0]
                    if txt.strip() == ')': # if no note
                        self.spec.note_trans[-1].append(txt)

            txt = fs.readline()

# =============================================================================
# Reactions
# =============================================================================

        gas = get_gas_ct(mech)
        self.react.equation = gas.reaction_equations()

        txt = fs.readline()

        while txt != "":
            k_units = {'A':False, 'Ea':False}
            if "# Reaction " in txt:

                # ---------------------------------------------
                self.react.number.append(int(txt.split("Reaction ")[1]))
                self.react.type.append('reaction')
                # self.react.dup.append(False)
                self.react.kin.append('')
                self.react.kin_units.append('')
                self.react.eff.append(False)
                self.react.param.append(False)
                # self.react.opt.append(False)
                self.react.opt.append({"duplicate":False, "negative-A":False, \
                                      "orders":False, "negative-orders":False})
                # ---------------------------------------------

                txt = fs.readline()
                while txt[0]=='#':
                    txt = fs.readline()
                if "three_body_reaction" in txt:
                    # 1st line data
                    self.react.type[-1] = "three_body_reaction"
                    # kinetic parameters
                    kin = txt.split("[")[1].split("]")[0].split(",")
                    self.react.kin[-1]       = [float(kin[0]),float(kin[1]),float(kin[2])]
                    self.react.kin_units[-1] = k_units
                    txt = fs.readline()
                    while txt.replace(' ','') != "\n" and txt != "":
                        # efficiencies
                        if "efficiencies" in txt:
                            if "'" in txt:
                                self.react.eff[-1] = txt.split("'")[1].split(" ")
                            else:
                                self.react.eff[-1] = txt.split("\"")[1].split(" ")
                        else:
                            # options
                            if 'duplicate' in txt:
                                self.duplicate_list.append(len(self.react.number)-1)
                                self.react.opt[-1]['duplicate'] = True
                            if 'negative_A' in txt:
                                self.react.opt[-1]['negative-A'] = True
                        txt = fs.readline()

                elif "falloff_reaction" in txt:
                    # 1st line data
                    self.react.type[-1] = "falloff_reaction"
                    txt = fs.readline()
                    while txt.replace(' ','') != "\n" and txt != "":
                        # kinetic parameters
                        if "kf0" in txt:
                            kf0=txt.split("[")[1].split("]")[0].split(",")
                            kf0=[float(kf0[0]),float(kf0[1]),float(kf0[2])]
                        elif "kf" in txt:
                            kf=txt.split("[")[1].split("]")[0].split(",")
                            kf=[float(kf[0]),float(kf[1]),float(kf[2])]
                        # efficiencies
                        elif "efficiencies" in txt:
                            if "'" in txt:
                                self.react.eff[-1] = txt.split("'")[1].split(" ")
                            else:
                                self.react.eff[-1] = txt.split("\"")[1].split(" ")
                        # Troe parameters
                        elif "Troe" in txt:
                            A  = txt.split('A=' )[1].split(',')[0].split(')')[0]
                            T1 = txt.split('T1=')[1].split(',')[0].split(')')[0]
                            if 'T2=' in txt:
                                T2 = txt.split('T2=')[1].split(',')[0].split(')')[0]
                            else:
                                T2 = False
                            T3 = txt.split('T3=')[1].split(',')[0].split(')')[0]
                            self.react.param[-1] = {'A':A, 'T1':T1, 'T2':T2, 'T3':T3}
                        else:
                            # options
                            if 'duplicate' in txt:
                                self.duplicate_list.append(len(self.react.number)-1)
                                self.react.opt[-1]['duplicate'] = True
                            if 'negative_A' in txt:
                                self.react.opt[-1]['negative-A'] = True
                        txt = fs.readline()
                    try:
                        self.react.kin[-1]       = [kf,kf0]
                        self.react.kin_units[-1] = [k_units,k_units]
                    except:
                        print("Error while reading falloff coeffs of reaction:")
                        print(self.react.number[-1])

                elif "chemically_activated_reaction" in txt:
                    # 1st line data
                    self.react.type[-1] = "chemically_activated_reaction"
                    txt = fs.readline()
                    while txt.replace(' ','') != "\n" and txt != "":
                        # kinetic parameters
                        if "kLow" in txt:
                            kLow=txt.split("[")[1].split("]")[0].split(",")
                            kLow=[float(kLow[0]),float(kLow[1]),float(kLow[2])]
                        elif "kHigh" in txt:
                            kHigh=txt.split("[")[1].split("]")[0].split(",")
                            kHigh=[float(kHigh[0]),float(kHigh[1]),float(kHigh[2])]
                        #efficiencies
                        elif "efficiencies" in txt:
                            if "'" in txt:
                                self.react.eff[-1] = txt.split("'")[1].split(" ")
                            else:
                                self.react.eff[-1] = txt.split("\"")[1].split(" ")
                        # Troe parameters
                        elif "Troe" in txt:
                            A  = txt.split('A=' )[1].split(',')[0].split(')')[0]
                            T1 = txt.split('T1=')[1].split(',')[0].split(')')[0]
                            if 'T2=' in txt:
                                T2 = txt.split('T2=')[1].split(',')[0].split(')')[0]
                            else:
                                T2 = False
                            T3 = txt.split('T3=')[1].split(',')[0].split(')')[0]
                            self.react.param[-1] = {'A':A, 'T1':T1, 'T2':T2, 'T3':T3}
                        else:
                            # options
                            if 'duplicate' in txt:
                                self.duplicate_list.append(len(self.react.number)-1)
                                self.react.opt[-1]['duplicate'] = True
                            if 'negative_A' in txt:
                                self.react.opt[-1]['negative-A'] = True
                        txt = fs.readline()
                    try:
                        self.react.kin[-1]       = [kHigh,kLow]
                        self.react.kin_units[-1] = [k_units,k_units]
                    except:
                        print("Error while reading falloff coeffs of reaction:")
                        print(self.react.number[-1])

                elif "pdep_arrhenius" in txt:
                    # 1st line data
                    self.react.type[-1] = "pdep_arrhenius"
                    txt = fs.readline()
                    self.react.kin[-1], self.react.kin_units[-1] = [],[]
                    self.react.param[-1] = [[], []]   # [pres, p_units]
                    while txt.replace(' ','') != "\n" and txt != "" and txt[0]!='#':
                        # options
                        if 'options' in txt:
                            if 'duplicate' in txt:
                                self.duplicate_list.append(len(self.react.number)-1)
                                self.react.opt[-1]['duplicate']  = True
                            if 'negative_A' in txt:
                                self.react.opt[-1]['negative-A'] = True
                        else:
                            try:
                                # reaction rate constants
                                kf=txt.split("[")[1].split("]")[0].split(",")
                                self.react.kin[-1].append([float(kf[-3]),float(kf[-2]),float(kf[-1])])
                                self.react.kin_units[-1].append(k_units)
                                # pressures
                                pres    = txt.split("(")[1].split(",")[0]
                                p_units = txt.split("'")[1]
                                self.react.param[-1][0].append(pres)
                                self.react.param[-1][1].append(p_units)
                            except:
                                print("Error while writing pdep coeffs of reaction:")
                                print(self.react.number[-1])
                        txt = fs.readline()

                elif "chebyshev" in txt:
                    # 1st line data
                    self.react.type[-1] = "chebyshev"
                    self.react.kin[-1] = []
                    txt = fs.readline()
                    while txt.replace(' ','') != "\n" and txt != "" and txt[0]!='#':
                        if 'duplicate' in txt:
                            self.duplicate_list.append(len(self.react.number)-1)
                            self.react.opt[-1]['duplicate'] = True
                        elif 'Tmin' in txt:
                            Tmin =  txt.split("Tmin=")[1].split(',')[0]
                            Tmax =  txt.split("Tmax=")[1].split(',')[0]
                        elif "Pmin=" in txt:
                            if '(' in txt:
                                Pmin =  txt.split("Pmin=")[1].split('Pmax=')[0]
                                Pmax =  txt.split("Pmax=")[1]
                            if 'bar' in Pmin.lower() or 'atm' in Pmin.lower() or 'pa' in Pmin.lower():
                                Pmin_value = float(Pmin.split('(')[1].split(' ')[0].split(',')[0])
                                Pmin_unit  = Pmin.split(' ')[1].split(')')[0]
                                Pmin_unit = Pmin_unit.replace("'","").replace('"','')
                            else:
                                Pmin_value = Pmin.split(',')[0]
                                Pmin_unit = False
                            if 'bar' in Pmax.lower() or 'atm' in Pmax.lower() or 'pa' in Pmax.lower():
                                Pmax_value = float(Pmax.split('(')[1].split(' ')[0].split(',')[0])
                                Pmax_unit  = Pmax.split(',')[1].split(')')[0]
                                Pmax_unit = Pmax_unit.replace("'","").replace('"','')
                            else:
                                Pmax_value = Pmax
                                Pmax_unit = False
                        else:
                            txt = txt.replace('[','')
                            txt = txt.replace(']','')
                            txt = txt.replace(')','')
                            txt = txt.replace('\n','')
                            txt = txt.replace('coeffs','')
                            txt = txt.replace('=','')
                            txt = txt.replace(' ','')
                            self.react.kin[-1].append([])
                            for ki in txt.split(","):
                                if ki != '':
                                    self.react.kin[-1][-1].append(float(ki))
                        txt = fs.readline()
                    self.react.param[-1] = {'Tmin':Tmin, 'Tmax':Tmax, \
                            'Pmin_value':Pmin_value, 'Pmax_value':Pmax_value,
                            'Pmin_unit' :Pmin_unit,  'Pmax_unit' :Pmax_unit}


                elif "reaction" in txt:
                    # 1st line data
                    self.react.type[-1] = "reaction"
                    kin = txt.split("[")[1].split("]")[0].split(",")
                    self.react.kin[-1] = [float(kin[0]),float(kin[1]),float(kin[2])]
                    self.react.kin_units[-1] = k_units
                    # other lines data (text)
                    txt = fs.readline()
                    while txt.replace(' ','') != "\n" and txt != "":
                        if 'options' in txt:
                            if 'duplicate' in txt:
                                self.duplicate_list.append(len(self.react.number)-1)
                                self.react.opt[-1]['duplicate']  = True
                            if 'negative_A' in txt:
                                self.react.opt[-1]['negative-A'] = True
                        txt = fs.readline()

            txt = fs.readline();

        # check threebody exception (+AR) (+HE) etc.
        for r in range(len(self.react.equation)):
            self.react.tbe.append(False)
            react_sp  = self.react.equation[r].split('=')[0]
            react_sp  = react_sp.replace('<','').replace(" ","")
            reactants = react_sp.split('+')
            prod_sp  = self.react.equation[r].split('=')[1]
            prod_sp  = prod_sp.replace('>','').replace(" ","")
            products = prod_sp.split('+')
            for sp in self.spec.name:
                if '(+'  + sp +')' in self.react.equation[r] \
                or '(+ ' + sp +')' in self.react.equation[r]:
                    self.react.tbe[-1]=sp
                    break
                elif sp in reactants and sp in products:
                    self.react.tbe[-1]=sp
                    break

        # check submechanism family:
#        ct.suppress_thermo_warnings()
#        gas = ct.Solution(mech)
        if int(ct.__version__[0])>2:
            nu_f = gas.reactant_stoich_coeffs
            nu_r = gas.product_stoich_coeffs
        else:
            nu_f = gas.reactant_stoich_coeffs()
            nu_r = gas.product_stoich_coeffs()
        self.react.subm_C  = [0]*len(self.react.equation)
        self.react.subm_N  = [0]*len(self.react.equation)
        self.react.subm_S  = [0]*len(self.react.equation)
        self.react.subm_Si = [0]*len(self.react.equation)
        self.react.subm_CO = [False]*len(self.react.equation)
        for r in range(len(self.react.equation)):
            for sp in range(len(self.spec.name)):
                if (nu_f[sp][r]!=0 or nu_r[sp][r]!=0) and nu_f[sp][r]!=nu_r[sp][r]:
                    n_at = self.spec.atoms[sp].split(',')
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
        self.react.activ_p = [False]*len(self.react.equation)   # act react for cases loop
        self.react.activ_m = [True]*len(self.react.equation)   # act react for methods loop
        self.react.modif   = [True] *len(self.react.equation)
        self.react_ref     = self.react
        self.react.ref_kin = copy.deepcopy(self.react.kin)






    def get_data_yaml(self,mech,verbose=0):

        fs = open(mech, 'r')
        self.duplicate_list = []


# =============================================================================
# Gas properties
# =============================================================================
        txt = fs.readlines()
        l=0

        self.description = {'description':False,'cantera_version':''}
        self.mech_units = {'length':'cm', 'time':'s', 'quantity':'mol', 'act_energy':'cal/mol'}
        while "phases:" not in txt[l]:
            if "description:" in txt[l]:
                if "|-" in txt[l]:
                    self.description['description'] = txt[l].split('description: ')[1]
                    while txt[l].replace(" ","") != "\n":
                        self.description['description'] += txt[l]
                        if l == len(txt):
                            print("problem while reading the description section in the new mechanism " + mech)
                            break

            if l == len(txt):
                print("problem while reading the description section in the new mechanism " + mech)
                break


            if "units:" in txt[l]:
                if 'length' in txt[l]:
                    self.mech_units['length'] = txt[l].split("length: ")[1].split(",")[0].split("}")[0]
                if 'time' in txt[l]:
                    self.mech_units['time'] = txt[l].split("time: ")[1].split(",")[0].split("}")[0]
                if 'quantity' in txt[l]:
                    self.mech_units['quantity'] = txt[l].split("quantity: ")[1].split(",")[0].split("}")[0]
                if 'activation-energy' in txt[l]:
                    self.mech_units['act_energy'] = txt[l].split("activation-energy: ")[1].split(",")[0].split("}")[0]
            l+=1
            self.gas_prop.append(txt[l])


        self.gas_prop = {'name':'gas', 'thermo':'ideal-gas', 'elements':False, \
                         'species':False, 'kinetics':'gas', 'reactions':'all',\
                         'transport':False, 'state_T':'300', 'state_p':'101325'}
        while "species:\n" not in txt[l].replace(" ",""):
            if "- name:" in txt[l]:
                self.gas_prop['name'] = txt[l].split('name: ')[1].replace('\n','')
            elif "thermo:" in txt[l]:
                self.gas_prop['thermo'] = txt[l].split('thermo: ')[1].replace('\n','')
            elif "elements:" in txt[l]:
                self.gas_prop['elements'] = txt[l].split('elements: ')[1]
            elif "species:" in txt[l]:
                self.gas_prop['species'] = txt[l].split('species: ')[1]
            elif "kinetics:" in txt[l]:
                self.gas_prop['kinetics'] = txt[l].split('kinetics: ')[1].replace('\n','')
            elif "reactions:" in txt[l]:
                self.gas_prop['reactions'] = txt[l].split('reactions: ')[1].replace('\n','')
            elif "transport:" in txt[l]:
                if 'mixture-averaged' in txt[l].split('transport: ')[1]:
                    self.gas_prop['transport'] = 'Mix'
                if 'multicomponent' in txt[l].split('transport: ')[1]:
                    self.gas_prop['transport'] = 'Multi'
            if "T: " in txt[l]:
                self.gas_prop['state_T'] = txt[l].split('T: ')[1].split(',')[0]
            if "P: " in txt[l]:
                self.gas_prop['state_p'] = txt[l].split('P: ')[1].split('}')[0]

            l+=1


# =============================================================================
# Thermo
# =============================================================================

        while "reactions:\n" not in txt[l].replace(" ",""):
            if "- name:" in txt[l]:
                self.spec.name.append([])
                self.spec.atoms.append([])
                self.spec.thermo_temp.append([])
                self.spec.thermo_coeff_lT.append([])
                self.spec.thermo_coeff_hT.append([])
                self.spec.thermo_model.append(False)
                self.spec.trans_model.append(False)
                self.spec.trans_geom.append(False)
                self.spec.trans_wd.append(False)
                self.spec.trans_diam.append(False)
                self.spec.trans_dipole.append(False)
                self.spec.trans_polar.append(False)
                self.spec.trans_rot_relax.append(False)
                self.spec.note_therm.append(False)
                self.spec.note_trans.append(False)

                self.spec.name[-1] = txt[l].split("name: ")[1].split("\n")[0]

                while txt[l].replace(' ','') != "\n" and "name:" not in txt[l+1]\
                and "reactions:\n" not in txt[l+1].replace(" ",""):
                    l+=1

                    if "composition:" in txt[l]:
                        self.spec.atoms[-1] = txt[l].split("{")[1].split("}")[0]

                    if "thermo" in txt[l]: read_data = "thermo"
                    if "transport" in txt[l]: read_data = "transport"

                    if "model" in txt[l]:
                        if read_data == "thermo":
                            self.spec.thermo_model[-1]= txt[l].split("model:")[1].replace(' ','').replace('\n','')
                        elif read_data == "transport":
                            self.spec.trans_model[-1]= txt[l].split("model:")[1].replace(' ','').replace('\n','')

                    if "temperature-ranges" in txt[l]:
                        T_range = txt[l].split("[")[1].split("]")[0].split(",")
                        T_range = [float(T) for T in T_range]
                        self.spec.thermo_temp[-1]= T_range

                    if "data" in txt[l]:
                        # low T coeff
                        l+=1
                        while "]" not in txt[l]:
                            txt_l = txt[l].replace("\n","").replace(" ","").replace(']','')
                            if "[" in txt_l:
                                lT_coeffs = txt_l.split("[")[1].split(",")
                                while '' in lT_coeffs: lT_coeffs.remove('')
                                lT_coeffs = [float(c) for c in lT_coeffs]
                                self.spec.thermo_coeff_lT[-1] = lT_coeffs
                            else:
                                lT_coeffs = txt_l.replace("\n","").split(",")
                                while '' in lT_coeffs: lT_coeffs.remove('')
                                lT_coeffs = [float(c) for c in lT_coeffs]
                                self.spec.thermo_coeff_lT[-1] += lT_coeffs
                            l+=1
                        if "[" in txt[l]:
                            txt_l = txt[l].replace("\n","").replace(" ","").replace(']','')
                            lT_coeffs = txt_l.split("[")[1].split(",")
                            while '' in lT_coeffs: lT_coeffs.remove('')
                            lT_coeffs = [float(c) for c in lT_coeffs]
                            self.spec.thermo_coeff_lT[-1] = lT_coeffs
                        else:
                            lT_coeffs = txt[l].replace("\n","").replace(']','').split(",")
                            while '' in lT_coeffs: lT_coeffs.remove('')
                            lT_coeffs = [float(c) for c in lT_coeffs]
                            self.spec.thermo_coeff_lT[-1] += lT_coeffs
                        l+=1
                        # high T coeff
                        while "]" not in txt[l]:
                            txt_l = txt[l].replace("\n","").replace(" ","").replace(']','')
                            if "[" in txt_l:
                                hT_coeffs = txt_l.split("[")[1].split(",")
                                while '' in hT_coeffs: hT_coeffs.remove('')
                                hT_coeffs = [float(c) for c in hT_coeffs]
                                self.spec.thermo_coeff_hT[-1] = hT_coeffs
                            else:
                                hT_coeffs = txt_l.replace("\n","").replace(']','').split(",")
                                while '' in hT_coeffs: hT_coeffs.remove('')
                                hT_coeffs = [float(c) for c in hT_coeffs]
                                self.spec.thermo_coeff_hT[-1] += hT_coeffs
                            l+=1
                        if "[" in txt[l]:
                            txt_l = txt[l].replace("\n","").replace(" ","").replace(']','')
                            hT_coeffs = txt_l.split("[")[1].split(",")
                            while '' in hT_coeffs: hT_coeffs.remove('')
                            hT_coeffs = [float(c) for c in hT_coeffs]
                            self.spec.thermo_coeff_hT[-1] = hT_coeffs
                        else:
                            hT_coeffs = txt[l].replace("\n","").replace(']','').split(",")
                            while '' in hT_coeffs: hT_coeffs.remove('')
                            hT_coeffs = [float(c) for c in hT_coeffs]
                            self.spec.thermo_coeff_hT[-1] += hT_coeffs

                    if "note" in txt[l]:
                        note = txt[l].split("note: ")[1]
                        while 'transport' not in txt[l]\
                        and '- name:' not in txt[l]\
                        and 'reactions:' not in txt[l]:
                            note += txt[l]
                            l+=1
                            if l == len(txt):
                                print("problem while reading the note of species" + self.spec.name[-1])
                                break

                        if read_data == "thermo":
                            self.spec.note_therm[-1] = note.split("note: ")[1].replace("\n","")
                        elif read_data == "transport":
                            self.spec.note_trans[-1] = note.split("note: ")[1].replace("\n","")
                        l-=1

                    if "geometry:"  in txt[l]:
                        self.spec.trans_geom[-1] = txt[l].split("geometry: ")[1].replace("\n","")
                    if "well-depth:"  in txt[l]:
                        self.spec.trans_wd[-1] = txt[l].split("well-depth: ")[1].replace("\n","")
                    if "diameter:"  in txt[l]:
                        self.spec.trans_diam[-1] = txt[l].split("diameter: ")[1].replace("\n","")
                    if 'dipole' in txt[l]:
                        self.spec.trans_dipole[-1] = txt[l].split("dipole: ")[1].replace("\n","")
                    if 'polar' in txt[l]:
                        self.spec.trans_polar[-1] = txt[l].split("polarizability: ")[1].replace("\n","")
                    if "rotational-relaxation:"  in txt[l]:
                        self.spec.trans_rot_relax[-1] = txt[l].split("rotational-relaxation: ")[1].replace("\n","")

            l+=1
            if l == len(txt):
                print("problem while reading the species section in the new mechanism " + mech)
                break



# =============================================================================
# Reactions
# =============================================================================

        gas = get_gas_ct(mech)
        self.react.equation = gas.reaction_equations()

        l+=1

        def get_arrhenius_cst(txt):
            A = txt.split('A:')[1].split(',')[0]
            if 'c' in A.lower() or 'm' in A.lower() or '^' in A:
                A_value = float(A.split(' ')[0])
                A_unit  = A.split(' ')[1]
            else:
                A_value = float(A)
                A_unit  = False
            b_value  = float(txt.split('b:')[1].split(',')[0])
            Ea = txt.split('Ea:')[1].split('}')[0]
            if '/' in Ea.lower() or 'm' in Ea.lower() or '^' in Ea:
                Ea_value = float(Ea.split(' ')[0])
                Ea_unit  = Ea.split(' ')[1]
            else:
                Ea_value = float(Ea)
                Ea_unit  = False

            kf       = [A_value,b_value,Ea_value]
            kf_units = {'A':A_unit,'Ea':Ea_unit}

            return kf, kf_units

        while l != (len(txt)):
            if "equation:" in txt[l]:
                self.react.number.append(int(txt[l].split('# Reaction ')[1]))
                self.react.type.append('reaction')
                self.react.opt.append({"duplicate":False, "negative-A":False, \
                                      "orders":False, "negative-orders":False})
                self.react.kin.append('')
                self.react.kin_units.append('')
                self.react.eff.append(False)
                self.react.param.append(False)

                l+=1
                while "- equation:" not in txt[l]:
                    # General options
                    if "type: " in txt[l]:
                        if  'three-body' in txt[l]:                         self.react.type[-1] = 'three_body_reaction'
                        elif  'falloff' in txt[l]:                          self.react.type[-1] = 'falloff_reaction'
                        elif  'chemically-activated' in txt[l]:             self.react.type[-1] = 'chemically_activated_reaction'
                        elif  'pressure-dependent-Arrhenius' in txt[l]:     self.react.type[-1] = 'pdep_arrhenius'
                        elif  'Chebychev' in txt[l]:                        self.react.type[-1] = 'chebychev'

                    if "duplicate: " in txt[l]:
                        if "true" in txt[l].lower():
                            self.react.opt[-1]["duplicate"] = True

                    if "negative-A: " in txt[l]:
                        if "true" in txt[l].lower():
                            self.react.opt[-1]["negative-A"] = True

                    if "orders: " in txt[l]:
                        self.react.opt[-1]["orders"] = txt[l].split('{')[1].split('}')[0]

                    if "negative-orders: " in txt[l]:
                        if "true" in txt[l].lower():
                            self.react.opt[-1]["negative-orders"] = True


                    # elementary and three body reactions
                    if " rate-constant: " in txt[l]:
                        self.react.kin[-1], self.react.kin_units[-1] = get_arrhenius_cst(txt[l])

                    # options for falloff ans chemicallly activated reactions
                    if "low-P-rate-constant: " in txt[l]:
                        kf0, kf0_units = get_arrhenius_cst(txt[l])

                    if "high-P-rate-constant: " in txt[l]:
                        kf, kf_units = get_arrhenius_cst(txt[l])

                    if "efficiencies: " in txt[l]:
                        while "}" not in txt[l]:
                            txt_i = txt[l].replace("\n","").replace("}","")
                            if "{" in txt_i:
                                eff = txt_i.split("{")[1].split(",")
                            else: eff += txt_i.split(",")
                            l+=1
                        txt_i = txt[l].replace("\n","").replace("}","")
                        if "{" in txt_i:
                            eff = txt_i.split("{")[1].split(",")
                        else: eff += txt_i.split(",")
                        eff = [effi.replace(" ","") for effi in eff]
                        if "" in eff : eff.remove("")
                        self.react.eff[-1] = eff

                    if "Troe"  in txt[l]:
                        A  = txt[l].split('A:' )[1].split(',')[0].split('}')[0]
                        T1 = txt[l].split('T1:')[1].split(',')[0].split('}')[0]
                        if 'T2:' in txt[l]:
                            T2 = txt[l].split('T2:')[1].split(',')[0].split('}')[0]
                        else:
                            T2 = False
                        T3 = txt[l].split('T3:')[1].split(',')[0].split('}')[0]
                        self.react.param[-1] = {'A':A, 'T1':T1, 'T2':T2, 'T3':T3}


                    # options for chebychev reactions
                    if "temperature-range: " in txt[l]:
                        Tmin =  txt[l].split("[")[1].split(',')[0]
                        Tmax =  txt[l].split(",")[1].split(']')[0]
                    if "pressure-range: " in txt[l]:
                        Pmin =  txt[l].split("[")[1].split(',')[0]
                        if 'bar' in Pmin.lower() or 'atm' in Pmin.lower() or 'pa' in Pmin.lower():
                            Pmin_value = float(Pmin.split(' ')[0])
                            Pmin_unit  = Pmin.split(' ')[1]
                        else:
                            Pmin_value = Pmin
                            Pmin_unit = False
                        Pmax =  txt[l].split(",")[1].split(']')[0]
                        if 'bar' in Pmax.lower() or 'atm' in Pmax.lower() or 'pa' in Pmax.lower():
                            Pmax_value = float(Pmax.split(' ')[0])
                            Pmax_unit  = Pmax.split(' ')[1]
                        else:
                            Pmax_value = Pmax
                            Pmax_unit = False
                    if "data: " in txt[l]:
                        while "]]" not in txt[l]:
                            if "[[" in txt[l]:
                                data = txt[l].split("[[")[1].split("],")[0].split(",")
                            else:
                                data += txt[l].split("[")[1].split("],")[0].split(",")
                            l+=1
                        data += txt[l].split("[")[1].split("]]")[0].split(",")


                    # options for pressure-dependent-Arrhenius reactions
                    if "rate-constants:" in txt[l]:
                        l+=1
                        kfp, kfp_units, pres, p_units = [], [], [], []
                        while '{P:' in txt[l]:
                            p_i = txt[l].split("P:")[1].split(",")[0]
                            if 'bar' in p_i.lower() or 'pa'  in p_i.lower()\
                            or 'atm' in p_i.lower() or 'mbar' in p_i.lower()\
                            or 'mpa' in p_i.lower():
                                p_units.append(p_i.split(" ")[-1])
                                pres.append(float(p_i.split(p_units[-1])[0]))
                            else:
                                pres.append(float(p_i))
                                p_units.append(False)
                            kfpi, kfpi_units = get_arrhenius_cst(txt[l])
                            kfp.append(kfpi) ; kfp_units.append(kfpi_units)
                            if l+1 == len(txt):  break
                            else:                l+=1
                        l-=1
                    l+=1
                    if l == len(txt):
                        break

                # once reactions lines are read...

                if 'falloff' in self.react.type[-1]:
                    try:
                        self.react.kin[-1]        = [kf,kf0]
                        self.react.kin_units[-1]  = [kf_units,kf0_units]
                    except:
                        print("Error while reading falloff coeffs of reaction:")
                        print(self.react.number[-1])

                if 'chebyshev' in self.react.type[-1]:
                    self.react.param[-1] = {'Tmin':Tmin, 'Tmax':Tmax, \
                            'Pmin_value':Pmin_value, 'Pmax_value':Pmax_value,
                            'Pmin_unit':Pmin_unit, 'Pmax_unit':Pmax_unit}
                    self.react.kin[-1] = data

                if 'pdep' in self.react.type[-1]:
                    try:
                        self.react.kin[-1]        = kfp
                        self.react.kin_units[-1]  = kfp_units
                        self.react.param[-1]      = [pres, p_units]
                    except:
                        print("Error while reading falloff coeffs of reaction:")
                        print(self.react.number[-1])

                if 'kf'   in locals(): del kf
                if 'kf0'  in locals(): del kf0
                if 'eff'  in locals(): del eff
                if 'Tmin' in locals(): del Tmin
                if 'Tmax' in locals(): del Tmax
                if 'Pmin' in locals(): del Pmin
                if 'Pmax' in locals(): del Pmax
                if 'data' in locals(): del data

            if l == len(txt):
                break
            elif "- equation:" not in txt[l]:
                l+=1



        # check threebody exception (+AR) (+HE) etc.
        for r in range(len(self.react.equation)):
            self.react.tbe.append(False)
            react_sp  = self.react.equation[r].split('=')[0]
            react_sp  = react_sp.replace('<','').replace(" ","")
            reactants = react_sp.split('+')
            prod_sp  = self.react.equation[r].split('=')[1]
            prod_sp  = prod_sp.replace('>','').replace(" ","")
            products = prod_sp.split('+')
            for sp in self.spec.name:
                if '(+'  + sp +')' in self.react.equation[r] \
                or '(+ ' + sp +')' in self.react.equation[r]:
                    self.react.tbe[-1]=sp
                    break
                elif sp in reactants and sp in products:
                    self.react.tbe[-1]=sp
                    break


        # check submechanism family:
        if int(ct.__version__[0])>2:
            nu_f = gas.reactant_stoich_coeffs
            nu_r = gas.product_stoich_coeffs
        else:
            nu_f = gas.reactant_stoich_coeffs()
            nu_r = gas.product_stoich_coeffs()
        self.react.subm_C  = [0]*len(self.react.equation)
        self.react.subm_N  = [0]*len(self.react.equation)
        self.react.subm_S  = [0]*len(self.react.equation)
        self.react.subm_Si = [0]*len(self.react.equation)
        self.react.subm_CO = [False]*len(self.react.equation)
        for r in range(len(self.react.equation)):
            for sp in range(len(self.spec.name)):
                if (nu_f[sp][r]!=0 or nu_r[sp][r]!=0) and nu_f[sp][r]!=nu_r[sp][r]:
                    n_at = self.spec.atoms[sp].split(',')
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
        self.react.activ_p = [False]*len(self.react.equation)   # act react for cases loop
        self.react.activ_m = [True]*len(self.react.equation)   # act react for methods loop
        self.react.modif   = [True] *len(self.react.equation)
        self.react_ref     = self.react
        self.react.ref_kin = copy.deepcopy(self.react.kin)







    def write_new_mech(self,filename="temp.cti",act_sp="no_arg",act_r="no_arg",  **kwargs):

        if 'version' in kwargs:     version = kwargs['version']
        else:                       version = False

        if len(filename)>3:
            if filename[-5:]=='.yaml':
                filename = filename[:-5] + '.cti'


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

        import datetime
        now    = datetime.datetime.now()
        date   = now.strftime("y:%Y m:%m d:%d,  %H:%M")

        fd.write("#-------------------------------------------------------------------------------\n")
        if version is not False:
            fd.write("# Kinetic mechanism converted by Brookesia "+version+"\n")
        else:
            fd.write("# Kinetic mechanism converted by Brookesia\n")
        fd.write("# " + date + "\n")
        fd.write("# https://github.com/Brookesia-py/Brookesia\n#\n")
        fd.write("#-------------------------------------------------------------------------------\n")
        fd.write("# Species and Elements\n")
        fd.write("#-------------------------------------------------------------------------------\n\n")

        l=0

        # Comments
        if self.description['description'] is not False:
            fd.write(self.description['description'])

        # Units
        fd.write("units(length='"   + self.mech_units['length']     \
               + "', time='"        + self.mech_units['time']       \
               + "', quantity='"    + self.mech_units['quantity']   \
               + "', act_energy='"  + self.mech_units['act_energy'] \
               + "')\n\n")

        # Gas characteristics
        self.gas_prop = {'name':'gas', 'thermo':'ideal-gas', 'elements':False,\
                 'elements':False, 'kinetics':'gas', 'reactions':'all',       \
                 'transport':False, 'state_T':'300', 'state_p':'101325'}

        #elements
        elements = ["C","H","O","N","AR","HE","Ar","He", "c", "h", "o", "n", "ar", "he"]
        txt_elements = ""
        for el in elements:
            if self.find_element(el,act_sp):
                txt_elements += el.capitalize() + " "
        #species
        i=0 ; txt_spec = ""
        for sp in sp_list:
            if i<5 :
                txt_spec+=self.spec.name[sp]+"  "
                i+=1
            else:
                txt_spec+="\n                     "+self.spec.name[sp]+"  "
                i=1

        # if self.gas_prop['thermo'] is not False:
        #     lname  = "name='"      + self.gas_prop['thermo']    + "',\n          "
        # else:
        #     lname = ""
        lname      = "name='gas',\n"
        l_elements = 'elements="'  + txt_elements[:-1]          + '",\n          '
        l_species  = 'species="""' + txt_spec[:-2]              + '""",\n          '

        if self.gas_prop['reactions'] is not False:
            l_reactions = "reactions='" + self.gas_prop['reactions'] + "',\n          "
        else:
            l_reactions = ""
        if self.gas_prop['transport'] is not False:
            l_transp    = "transport='" + self.gas_prop['transport'] + "',\n          "
        else:
            l_transp = ""
        l_state = "initial_state=state(temperature=" + self.gas_prop['state_T']\
                 + ", pressure=OneAtm))\n\n\n"

        fd.write(self.gas_prop['thermo'].replace('-','_') + "("                                \
           + lname + l_elements + l_species + l_reactions + l_transp + l_state)


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
            l_atoms ="        atoms='"+self.spec.atoms[sp].replace(': ',':')+"',\n"
            if self.spec.thermo_model=='NASA9':
                l_thermo="        thermo=(NASA9(["
            else:
                l_thermo="        thermo=(NASA(["
            l_thermo += str('%0.2e' %self.spec.thermo_temp[sp][0])+", "     \
                     + str('%0.2e' %self.spec.thermo_temp[sp][1])           \
                     + "],\n                     [ "                        \
                     + str('%0.8e' %self.spec.thermo_coeff_lT[sp][0])+", "  \
                     + str('%0.8e' %self.spec.thermo_coeff_lT[sp][1])+", "  \
                     + str('%0.8e' %self.spec.thermo_coeff_lT[sp][2])+",\n" \
                     + "                       "                            \
                     + str('%0.8e' %self.spec.thermo_coeff_lT[sp][3])+", "  \
                     + str('%0.8e' %self.spec.thermo_coeff_lT[sp][4])+", "  \
                     + str('%0.8e' %self.spec.thermo_coeff_lT[sp][5])+",\n" \
                     + "                       "
            if self.spec.thermo_model=='NASA9':
                l_thermo += str('%0.8e' %self.spec.thermo_coeff_lT[sp][6])+", "    \
                          + str('%0.8e' %self.spec.thermo_coeff_lT[sp][7])+", "    \
                          + str('%0.8e' %self.spec.thermo_coeff_lT[sp][8])+"]),\n" \
                          + "               NASA9(["
            else:
                l_thermo += str('%0.8e' %self.spec.thermo_coeff_lT[sp][6])+"]),\n" \
                          + "                NASA(["

            l_thermo +=  str('%0.2e' %self.spec.thermo_temp[sp][1])+", "     \
                     + str('%0.2e' %self.spec.thermo_temp[sp][2])         \
                     + "],\n                     [ "                  \
                     + str('%0.8e' %self.spec.thermo_coeff_hT[sp][0])+", "     \
                     + str('%0.8e' %self.spec.thermo_coeff_hT[sp][1])+", "     \
                     + str('%0.8e' %self.spec.thermo_coeff_hT[sp][2])+",\n"    \
                     + "                       "                      \
                     + str('%0.8e' %self.spec.thermo_coeff_hT[sp][3])+", "     \
                     + str('%0.8e' %self.spec.thermo_coeff_hT[sp][4])+", "     \
                     + str('%0.8e' %self.spec.thermo_coeff_hT[sp][5])+",\n"    \
                     + "                       "
            if self.spec.thermo_model=='NASA9':
                l_thermo += str('%0.8e' %self.spec.thermo_coeff_hT[sp][6])+", "    \
                          + str('%0.8e' %self.spec.thermo_coeff_hT[sp][7])+", "
                if len(self.spec.trans_model[sp]) is not False \
                or len(self.spec.note_trans[sp])  is not False:
                    l_thermo += str('%0.8e' %self.spec.thermo_coeff_hT[sp][8])+"])),\n"
                else:
                    l_thermo += str('%0.8e' %self.spec.thermo_coeff_hT[sp][8])+"])))\n"
            else:
                if len(self.spec.trans_model[sp]) is not False \
                or len(self.spec.note_trans[sp])  is not False:
                    l_thermo += str('%0.8e' %self.spec.thermo_coeff_hT[sp][6])+"])),\n"
                else:
                    l_thermo += str('%0.8e' %self.spec.thermo_coeff_hT[sp][6])+"])))\n"

            space = '                                '
            if self.spec.trans_model[sp] is not False:
                l_trans="        transport=gas_transport("
                if self.spec.trans_geom[sp] is not False:
                    l_trans += "geom='" + self.spec.trans_geom[sp] + "',\n"
                if self.spec.trans_diam[sp] is not False:
                    l_trans += space + "diam=" + self.spec.trans_diam[sp] + ",\n"
                if self.spec.trans_wd[sp] is not False:
                    l_trans += space + "well_depth=" + self.spec.trans_wd[sp] + ",\n"
                if self.spec.trans_dipole[sp] is not False:
                    l_trans += space + "dipole=" + self.spec.trans_dipole[sp] + ",\n"
                if self.spec.trans_polar[sp] is not False:
                    l_trans += space + "polar=" + self.spec.trans_polar[sp] + ",\n"
                if self.spec.trans_rot_relax[sp] is not False:
                    l_trans += space + "rot_relax=" + self.spec.trans_rot_relax[sp] + ",\n"
                if self.spec.note_trans[sp] is not False:
                    l_trans = l_trans[:-1] + "),\n"
                    self.spec.note_trans[sp] = self.spec.note_trans[sp].replace("'","")
                    l_trans += "        note='" + self.spec.note_trans[sp] + "')\n"
                else:
                    l_trans = l_trans[:-1] + "))\n"
            else:
                l_trans = ""



            fd.write(l_name)
            fd.write(l_atoms)
            fd.write(l_thermo)
            fd.write(l_trans)
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
            duplicate_r, tba  = False, False
            # if self.react.number[r]==7:
            #     print('toto')
            if self.react.opt[r]['duplicate'] == True:
                reactants = self.react.equation[r].split('=')[0]
                reactants = reactants.replace('<','').replace('(+','+').replace('( +','+')
                reactants = reactants.split('+')
                reactants_a = [react.replace(' ','') for react in reactants]
                products  = self.react.equation[r].split('=')[1]
                products  = products.replace('>','').replace('(+','+').replace('( +','+')
                products  = products.split('+')
                products_a  = [prod.replace(' ','') for prod in products]
                if self.react.tbe[r] is not False:
                    tbe = self.react.tbe[r]
                    if tbe     in reactants_a:
                        reactants_a.remove(tbe)
                        tba = True
                    if tbe+")" in reactants_a:
                        reactants_a.remove(tbe + ")")
                        tba = True
                    if tbe     in products_a:
                        products_a.remove(tbe)
                        tba = True
                    if tbe+")" in products_a:
                        products_a.remove(tbe + ")")
                        tba = True
                if "M"  in reactants_a:
                    reactants_a.remove("M")
                    tba = True
                if "M)" in reactants_a:
                    reactants_a.remove("M)")
                    tba = True
                if "M"  in products_a:
                    products_a.remove("M")
                    tba = True
                if "M)" in products_a:
                    products_a.remove("M)")
                    tba = True
                for d in range(len(self.react.equation)):
                    if d != r:
                        tbb = False
                        # if self.react.number[r]==7 and self.react.number[d]==75:
                        #     print('toto')
                        if self.react.opt[d]['duplicate'] == True:
                            reactants = self.react.equation[d].split('=')[0]
                            reactants = reactants.replace('<','').replace('(+','+').replace('( +','+')
                            reactants = reactants.split('+')
                            reactants_b = [react.replace(' ','') for react in reactants]
                            products  = self.react.equation[d].split('=')[1]
                            products  = products.replace('>','').replace('(+','+').replace('( +','+')
                            products  = products.split('+')
                            products_b  = [prod.replace(' ','') for prod in products]
                            if self.react.tbe[d] is not False:
                                tbe = self.react.tbe[d]
                                if tbe     in reactants_b:
                                    reactants_b.remove(tbe) ;        tbb = True
                                if tbe+")" in reactants_b:
                                    reactants_b.remove(tbe + ")") ;  tbb = True
                                if tbe     in products_b:
                                    products_b.remove(tbe) ;         tbb = True
                                if tbe+")" in products_b:
                                    products_b.remove(tbe + ")") ;   tbb = True
                            if "M"  in reactants_b:
                                reactants_b.remove("M")  ; tbb = True
                            if "M)" in reactants_b:
                                reactants_b.remove("M)") ; tbb = True
                                tbb = True
                            if "M"  in products_b:
                                products_b.remove("M")   ; tbb = True
                            if "M)" in products_b:
                                products_b.remove("M)")  ; tbb = True
                            common_reac_1 = set(reactants_a).intersection(reactants_b)
                            common_prod_1 = set(products_a).intersection(products_b)
                            common_reac_2 = set(products_a).intersection(reactants_b)
                            common_prod_2 = set(products_b).intersection(reactants_a)
                            if  len(common_reac_1) == len(reactants_a) == len(reactants_b)\
                            and len(common_prod_1) == len(products_a)  == len(products_b)\
                            and tba==tbb:
                                if self.react.tbe[d] is not False:  # remove duplicate species if third body contribution is neglectable
                                # ... otherwise it may prevent the withdrawal of the third body species
                                    if d in r_list:  duplicate_r = True
                                else:
                                    duplicate_r = True # add duplicate species without third body
                                    if d not in r_list: r_list.append(d)
                            if  len(common_reac_2) == len(products_a)  == len(reactants_b)\
                            and len(common_prod_2) == len(reactants_a) == len(products_b)\
                            and tba==tbb:
                                if self.react.tbe[d] is not False:
                                    if d in r_list:  duplicate_r = True
                                else:
                                    duplicate_r = True
                                    if d not in r_list: r_list.append(d)

            self.react.opt[r]['duplicate'] = duplicate_r

        r_list.sort()


        for r in r_list: # in activated reactions list

            fd.write("\n# Reaction "+str(self.react.number[r])+"\n")

            if self.react.type[r]=="reaction":
                #    # ReactionXX
                #    reaction('H2O2 + OH <=> H2O + HO2', [1.740000e+12, 0.0, 318.0],
                #             options='duplicate')
                r_line = "reaction('"+self.react.equation[r]+"', ["            \
                         + str('%0.6e' %self.react.kin[r][0]) + ", "          \
                         + str('%0.4f' %self.react.kin[r][1]) + ", "          \
                         + str('%0.3f' %self.react.kin[r][2]) + "]"
                # options:
                r_opt = ""
                if self.react.opt[r]['duplicate']  is not False\
                or self.react.opt[r]['negative-A'] is not False:
                    r_line += ",\n"
                    r_opt = "         options=['"
                    if self.react.opt[r]['duplicate']  is not False:
                        r_opt += "duplicate'"
                        if self.react.opt[r]['negative-A']  is not False:
                            r_opt += ", 'negative_A'])\n\n"
                        else:
                            r_opt += "])\n\n"
                    else:
                        r_opt += "negative_A'])\n\n"
                else:
                    r_line += ")\n\n"

                fd.write(r_line)
                fd.write(r_opt)

            if self.react.type[r]=="three_body_reaction":
                #    # Reaction XX
                #    three_body_reaction('O + O + M <=> O2 + M', [6.165000e+15, -0.5, 0.0],
                #    efficiencies='CO2:3.8 CO:1.9 H2:2.5 H2O:12.0 AR:0.83 CH4:2.0 C2H6:3.0 HE:0.83')
                r_line = "three_body_reaction('"+self.react.equation[r]+"', [" \
                         + str('%0.6e' %self.react.kin[r][0]) + ", "          \
                         + str('%0.4f' %self.react.kin[r][1]) + ", "          \
                         + str('%0.3f' %self.react.kin[r][2]) + "]"

                r_efficiencies, r_opt = "", ""
                if self.react.eff[r]               is not False\
                or self.react.opt[r]['duplicate']  is not False\
                or self.react.opt[r]['negative-A'] is not False:
                    # efficiencies
                    if self.react.eff[r] is not False:
                        r_efficiencies = "                    efficiencies='"
                        for sp in sp_list:
                            for coll_sp in self.react.eff[r]:
                                if self.spec.name[sp] == coll_sp.split(":")[0]:
                                    r_efficiencies+=coll_sp+" "
                        if r_efficiencies == "                    efficiencies='":    # if all collision species are removed
                            r_efficiencies = ""
                        else:
                            if self.react.opt[r]['duplicate']  is not False\
                            or self.react.opt[r]['negative-A'] is not False:
                                r_efficiencies = r_efficiencies[:-1] + "',\n"
                            else:
                                r_efficiencies = r_efficiencies[:-1] + "')\n\n"
                    # options:
                    if self.react.opt[r]['duplicate']  is not False\
                    or self.react.opt[r]['negative-A'] is not False:
                        r_opt = "                    options=['"
                        if self.react.opt[r]['duplicate']  is not False:
                            r_opt += "duplicate'"
                            if self.react.opt[r]['negative-A']  is not False:
                                r_opt += ", 'negative_A'])\n\n"
                            else:
                                r_opt += "])\n\n"
                        else:
                            r_opt += "negative_A'])\n\n"
                if r_efficiencies + r_opt != "":
                    r_line += ",\n"
                else:
                    r_opt=")\n\n"

                fd.write(r_line)
                fd.write(r_efficiencies)
                fd.write(r_opt)

            if self.react.type[r]=="falloff_reaction" \
            or self.react.type[r]=="chemically_activated_reaction":
                #    # Reaction 71
                #    falloff_reaction('CO + H2 (+ M) <=> CH2O (+ M)',
                #                     kf=[4.300000e+07, 1.5, 79600.0],
                #                     kf0=[5.070000e+27, -3.42, 84348.0],
                #                     efficiencies='CO2:2.0 CO:1.5 H2:2.0 H2O:6.0 AR:0.7 CH4:2.0 C2H6:3.0 HE:0.7',
                #                     falloff=Troe(A=0.932, T3=197.0, T1=1540.0, T2=10300.0),
                #                     options=['duplicate', 'negative_A'])
                #
                ## Reaction 134
                #    chemically_activated_reaction('ch3 + oh (+ M) <=> ch2(s) + h2o (+ M)',
                #                                  kLow=[1.128000e+15, -0.63327, -493.15597],
                #                                  kHigh=[2.394000e-03, 4.096, -1241.875],
                #                                  falloff=Troe(A=2.122, T3=837.667, T1=2326.05, T2=4432.0))
                #
                if self.react.type[r]=="falloff_reaction" :
                    r_line = "falloff_reaction('"
                    space  = "                 "
                    ka     = "kf" ; kb = "kf0"
                else:
                    r_line = "chemically_activated_reaction('"
                    space  = "                              "
                    ka     = "kLow" ; kb = "kHigh"

                r_line += self.react.equation[r]+"',\n"   \
                         + space + ka + "=["                                     \
                         + str('%0.6e' %self.react.kin[r][0][0]) + ", "       \
                         + str('%0.4f' %self.react.kin[r][0][1]) + ", "       \
                         + str('%0.3f' %self.react.kin[r][0][2]) + "],\n"     \
                         + space + kb + "=["                           \
                         + str('%0.6e' %self.react.kin[r][1][0]) + ", "       \
                         + str('%0.4f' %self.react.kin[r][1][1]) + ", "       \
                         + str('%0.3f' %self.react.kin[r][1][2]) + "]"
                r_efficiencies, r_troe, r_opt = "", "", ""
                if self.react.eff[r]               is not False\
                or self.react.param[-1]            is not False\
                or self.react.opt[r]['duplicate']  is not False\
                or self.react.opt[r]['negative-A'] is not False:
                    # efficiencies
                    if self.react.eff[r] is not False:
                        r_efficiencies += space + "efficiencies='"
                        for sp in sp_list:
                            for coll_sp in self.react.eff[r]:
                                if self.spec.name[sp] == coll_sp.split(":")[0]:
                                    r_efficiencies+=coll_sp+" "
                        if r_efficiencies[-14:]=="efficiencies='":    # if all collision species are removed
                            r_efficiencies = ""
                        else:
                            if self.react.param[r] is not False\
                            or self.react.opt[r]['duplicate']  is not False\
                            or self.react.opt[r]['negative-A'] is not False:
                                r_efficiencies = r_efficiencies[:-1] + "',\n"
                            else:
                                r_efficiencies = r_efficiencies[:-1] + "')\n\n"
                    # troe
                    if self.react.param[r] is not False:
                        A      = 'A='  + self.react.param[r]['A']  + ','
                        T1     = 'T1=' + self.react.param[r]['T1'] + ','
                        if self.react.param[r]['T2'] is not False:
                            T2 = 'T2=' + self.react.param[r]['T2'] + ','
                        else:
                            T2 = ''
                        T3     = 'T3=' + self.react.param[r]['T3'] + ')'
                        r_troe= space + "falloff=Troe(" +A+T1+T2+T3
                        if self.react.opt[r]['duplicate']  is not False\
                        or self.react.opt[r]['negative-A'] is not False:
                            r_troe+=",\n"
                        else:
                            r_troe+=")\n\n"
                    # options:
                    if self.react.opt[r]['duplicate']  is not False\
                    or self.react.opt[r]['negative-A'] is not False:
                        r_opt = space + "options=['"
                        if self.react.opt[r]['duplicate']  is not False:
                            r_opt += "duplicate'"
                            if self.react.opt[r]['negative-A']  is not False:
                                r_opt += ", 'negative_A'])\n\n"
                            else:
                                r_opt += "])\n\n"
                        else:
                            r_opt += "negative_A'])\n\n"
                if r_efficiencies + r_troe + r_opt != "":
                    r_line += ",\n"
                else:
                    r_opt=")\n\n"

                fd.write(r_line)
                fd.write(r_efficiencies)
                fd.write(r_troe)
                fd.write(r_opt)


            if self.react.type[r]=="pdep_arrhenius":
                #    # Reaction 391
                #    pdep_arrhenius('CH3COCH3 <=> CH3CO + CH3',
                #                   [(0.01, 'atm'), 2.050000e+58, -12.796, 100030.1],
                #                   [(0.1, 'atm'), 3.300000e+51, -10.574, 98221.2],
                #                   [(1.0, 'atm'), 1.310000e+42, -7.657, 94660.6],
                #                   [(10.0, 'atm'), 2.160000e+33, -4.989, 90916.5],
                #                   [(100.0, 'atm'), 9.400000e+28, -3.669, 89022.8])
                r_line = "pdep_arrhenius('"+self.react.equation[r]+"'"
                k_line = ""
                for l in range(len(self.react.kin[r])):
                    k_line += ",\n               [("                    \
                             + str(self.react.param[r][0][l]) + ", '"        \
                             + str(self.react.param[r][1][l]) + "'), "       \
                             + '%0.6e' %self.react.kin[r][l][0] + ", "  \
                             + '%0.4f' %self.react.kin[r][l][1] + ", "  \
                             + '%0.3f' %self.react.kin[r][l][2] + "]"
                r_efficiencies, r_opt = "", ""
                if self.react.eff[r]               is not False\
                or self.react.opt[r]['duplicate']  is not False\
                or self.react.opt[r]['negative-A'] is not False:
                    # efficiencies
                    if self.react.eff[r]           is not False:
                        r_efficiencies += "               efficiencies='"
                        for sp in sp_list:
                            for coll_sp in self.react.eff[r]:
                                if self.spec.name[sp] == coll_sp.split(":")[0]:
                                    r_efficiencies+=coll_sp+" "
                        if r_efficiencies[-14:]=="efficiencies='":    # if all collision species are removed
                            r_efficiencies = ""
                        else:
                            if self.react.opt[r]['duplicate']  is not False\
                            or self.react.opt[r]['negative-A'] is not False:
                                r_efficiencies = r_efficiencies[:-1] + "',\n"
                            else:
                                r_efficiencies = r_efficiencies[:-1] + "')\n\n"
                    # options:
                    if self.react.opt[r]['duplicate']  is not False\
                    or self.react.opt[r]['negative-A'] is not False:
                        r_opt = space + "options=['"
                        if self.react.opt[r]['duplicate']  is not False:
                            r_opt += "duplicate'"
                            if self.react.opt[r]['negative-A']  is not False:
                                r_opt += ", 'negative_A'])\n\n"
                            else:
                                r_opt += "])\n\n"
                        else:
                            r_opt += "negative_A'])\n\n"
                if r_efficiencies + r_opt != "":
                    k_line += ",\n"
                else:
                    k_line+=")\n\n"

                fd.write(r_line)
                fd.write(k_line)
                fd.write(r_efficiencies)
                fd.write(r_opt)

            if self.react.type[r]=="chebyshev":
                # # Reaction 274
                # chebyshev_reaction('H + C2H2 (+ M) <=> C2H3 (+ M)',
                #                    Tmin=500.0, Tmax=2000.0,
                #                    Pmin=(0.001, 'atm'), Pmax=(100.0, 'atm'),
                #                    coeffs=[[ 1.06310e+01,  2.25240e+00, -1.51140e-01, -6.30360e-02],
                #                            [-3.35090e-01,  2.48140e-01,  1.42620e-01,  4.89960e-02],
                #                            [-3.51450e-01, -5.96700e-03,  3.35700e-03,  8.66000e-03],
                #                            [-9.37400e-02, -9.53930e-03, -5.47860e-03, -9.99650e-04],
                #                            [-2.04170e-02, -2.92750e-03, -2.25310e-03, -1.27840e-03],
                #                            [-2.22070e-03, -3.08030e-04, -3.52370e-04, -3.72070e-04],
                #                            [ 8.88590e-04,  2.44050e-04,  1.29870e-04,  1.45970e-05]])
                r_line = "chebyshev_reaction('"+self.react.equation[r]+"',\n"
                space = "                   "
                # T
                k_line  = space + 'Tmin=' + self.react.param[r]['Tmin'] + ', '\
                                  'Tmax=' + self.react.param[r]['Tmax'] + ',\n'
                # P
                k_line += space + 'Pmin=(' + self.react.param[r]['Pmin_value'] + ", '"  \
                                + self.react.param[r]['Pmin_unit']             + "'), " \
                                  'Pmax=(' + self.react.param[r]['Pmax_value'] + ", '"  \
                                + self.react.param[r]['Pmin_unit']             + "'),\n"\
                # coeffs
                for l in range(len(self.react.kin[r])):
                    if l==0:
                        k_line += space + 'coeffs=[['
                    else:
                        k_line += space + '        ['
                    for c in range(len(self.react.kin[r][l])):
                        k_line+='%0.5e' %self.react.kin[r][l][c]
                        if c<len(self.react.kin[r][l])-1:
                            k_line+=", "
                        else:
                            k_line+="]"
                    if l==len(self.react.kin[r])-1:
                        k_line+=']\n\n'
                    else:
                        k_line+=',\n'

                fd.write(r_line)
                fd.write(k_line)


        timer.sleep(1.5)

        fd.close()


    def write_yaml_mech(self,filename="temp.cti",act_sp="no_arg",act_r="no_arg", **kwargs):

        if 'version' in kwargs:     version = kwargs['version']
        else:                       version = False

        if len(filename)>3:
            if filename[-4:]=='.cti':
                filename = filename[:-4] + '.yaml'

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

        import datetime
        now    = datetime.datetime.now()
        date   = now.strftime("y:%Y m:%m d:%d,  %H:%M")

        fd.write("#-------------------------------------------------------------------------------\n")
        if version is not False:
            fd.write("# Kinetic mechanism converted by Brookesia "+version+"\n")
        else:
            fd.write("# Kinetic mechanism converted by Brookesia\n")
        fd.write("# " + date + "\n")
        fd.write("# https://github.com/Brookesia-py/Brookesia\n#\n")
        fd.write("#-------------------------------------------------------------------------------\n")
        fd.write("# Units, Species and Elements\n")
        fd.write("#-------------------------------------------------------------------------------\n\n")


        # Units writing
        fd.write("units: {length: " + self.mech_units['length']     \
               + ", time: "         + self.mech_units['time']       \
               + ", quantity: "     + self.mech_units['quantity']   \
               + ", activation-energy: "   + self.mech_units['act_energy'] \
               + "}\n\n")


        # Gas characteristics
        #elements
        elements = ["C","H","O","N","AR","HE","Ar","He", "c", "h", "o", "n", "ar", "he"]
        txt_elements = ""
        for el in elements:
            if self.find_element(el,act_sp):
                txt_elements += el.capitalize() + ", "
        #species
        i=0 ; txt_spec = ""
        for sp in sp_list:
            if i<8 :
                txt_spec+=self.spec.name[sp]+", "
                i+=1
            else:
                txt_spec+="\n"
                txt_spec+="    "+self.spec.name[sp]+", "
                i=1

        #transport
        if self.gas_prop['transport'] == 'Mix':
            transport = 'mixture-averaged'
        elif self.gas_prop['transport'] == 'Multi':
            transport = 'multicomponent'
        else:
            transport = False
        fd.write("phases:\n")
        fd.write(  "- name: "      + self.gas_prop['name']      + "\n"  \
                 + "  thermo: "    + self.gas_prop['thermo'].replace('_','-')    + "\n"  \
                 + '  elements: [' + txt_elements[:-2]          + ']\n' \
                 + '  species: ['  + txt_spec[:-2]              + ']\n' \
                 + '  kinetics: '  + self.gas_prop['kinetics']  + '\n'  \
                 + "  reactions: " + self.gas_prop['reactions'] + "\n")
        if transport is not False:
            fd.write("  transport: " + transport                + "\n")
        fd.write("  state:\n    T: " + self.gas_prop['state_T']         \
                 +       "\n    P: 101325\n\n")



#         # Gas characteristics writing
#         while "elements" not in self.gas_prop[l]:
#             fd.write(self.gas_prop[l])
#             l+=1

#         #elements
#         txt_elements = self.gas_prop[l].split('[')[0]+'['
#         elements = ["C","H","O","N","AR","HE","Ar","He", "c", "h", "o", "n", "ar", "he"]
#         for el in elements:
#             if self.find_element(el,act_sp):
# #                if el == "AR":
# #                    txt_elements += "Ar" + " "
# #                elif el == "HE":
# #                    txt_elements += "He" + " "
# #                else:
#                 txt_elements += el.capitalize() + " "
#         txt_elements += ']\n'
#         fd.write(txt_elements)

#         #species
#         txt_spec = '  species: ['
#         i=0
#         for sp in sp_list:
#             if i<10 :
#                 txt_spec+=self.spec.name[sp]+", "
#                 i+=1
#             else:
#                 fd.write(txt_spec+"\n")
#                 txt_spec="    "+self.spec.name[sp]+", "
#                 i=1
#         txt_spec[:-2]+=']\n'
#         fd.write(txt_spec)

#         while "reactions" not in self.gas_prop[l]:
#             l+=1

#         while self.gas_prop[l]!="\n":
#             fd.write(self.gas_prop[l])
#             l+=1

#         fd.write("\n")


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
        fd.write("species:\n")


        for sp in sp_list: # in activated species list
            composition = self.spec.atoms[sp].replace(' ','').replace(':', ': ')
            l_name    = "- name: "+self.spec.name[sp]+"\n"
            l_atoms   = "  composition: {"+composition+"}\n"
            # thermo
            l_thermo  = "  thermo:\n"
            l_thermo += "    model: " + self.spec.thermo_model[sp] + "\n"
            l_thermo += "    temperature-ranges: ["
            for T in self.spec.thermo_temp[sp]: l_thermo += str(T) + ", "
            l_thermo = l_thermo[:-2] + "]\n"
            l_thermo += "    data:\n    - ["
            for i,T in enumerate(self.spec.thermo_coeff_lT[sp]):
                if i==4:
                    l_thermo += str('%0.8e' %T) + ",\n      "
                elif i == len(self.spec.thermo_coeff_lT[sp])-1:
                    l_thermo += str('%0.8e' %T) + "]\n    - ["
                else:
                    l_thermo += str('%0.8e' %T) + ", "
            for i,T in enumerate(self.spec.thermo_coeff_hT[sp]):
                if i==4:
                    l_thermo += str('%0.8e' %T) + ",\n      "
                elif i == len(self.spec.thermo_coeff_hT[sp])-1:
                    l_thermo += str('%0.8e' %T) + "]\n"
                else:
                    l_thermo += str('%0.8e' %T) + ", "
            if self.spec.note_therm[sp] is not False:
                l_thermo += "  note: " + self.spec.note_therm[sp]

            # transport
            if self.spec.trans_model[sp] is not False:
                l_trans = "  transport:\n"
                l_trans += "    model: " + self.spec.trans_model[sp] + "\n"
                if self.spec.trans_geom[sp]      is not False:
                    l_trans += "    geometry: "              + self.spec.trans_geom[sp]      + "\n"
                if self.spec.trans_wd[sp]        is not False:
                    l_trans += "    well-depth: "            + self.spec.trans_wd[sp]        + "\n"
                if self.spec.trans_diam[sp]      is not False:
                    l_trans += "    diameter: "              + self.spec.trans_diam[sp]      + "\n"
                if self.spec.trans_dipole[sp]      is not False:
                    l_trans += "    dipole: "                + self.spec.trans_dipole[sp]    + "\n"
                if self.spec.trans_polar[sp]      is not False:
                    l_trans += "    polarizability: "        + self.spec.trans_polar[sp]     + "\n"
                if self.spec.trans_rot_relax[sp] is not False:
                    l_trans += "    rotational-relaxation: " + self.spec.trans_rot_relax[sp] + "\n"
                if self.spec.note_trans[sp]      is not False:
                    l_trans += "  note: "                    + self.spec.note_trans[sp]
            else:
                l_trans = ""

            fd.write(l_name)
            fd.write(l_atoms)
            fd.write(l_thermo)
            fd.write(l_trans)
            fd.write("\n")


# =============================================================================
#         Reactions
# =============================================================================

        fd.write("#-------------------------------------------------------------------------------\n")
        fd.write("# Reaction data\n")
        fd.write("#-------------------------------------------------------------------------------\n")
        fd.write("\n")
        fd.write("reactions:\n")



        r_list_ = copy.deepcopy(r_list)
        for r in r_list_:
            duplicate_r, tba  = False, False
            # if self.react.number[r]==7:
            #     print('toto')
            if self.react.opt[r]['duplicate'] == True:
                reactants = self.react.equation[r].split('=')[0]
                reactants = reactants.replace('<','').replace('(+','+').replace('( +','+')
                reactants = reactants.split('+')
                reactants_a = [react.replace(' ','') for react in reactants]
                products  = self.react.equation[r].split('=')[1]
                products  = products.replace('>','').replace('(+','+').replace('( +','+')
                products  = products.split('+')
                products_a  = [prod.replace(' ','') for prod in products]
                if self.react.tbe[r] is not False:
                    tbe = self.react.tbe[r]
                    if tbe     in reactants_a:
                        reactants_a.remove(tbe)
                        tba = True
                    if tbe+")" in reactants_a:
                        reactants_a.remove(tbe + ")")
                        tba = True
                    if tbe     in products_a:
                        products_a.remove(tbe)
                        tba = True
                    if tbe+")" in products_a:
                        products_a.remove(tbe + ")")
                        tba = True
                if "M"  in reactants_a:
                    reactants_a.remove("M")
                    tba = True
                if "M)" in reactants_a:
                    reactants_a.remove("M)")
                    tba = True
                if "M"  in products_a:
                    products_a.remove("M")
                    tba = True
                if "M)" in products_a:
                    products_a.remove("M)")
                    tba = True
                for d in range(len(self.react.equation)):
                    if d != r:
                        tbb = False
                        # if self.react.number[r]==7 and self.react.number[d]==75:
                        #     print('toto')
                        if self.react.opt[d]['duplicate'] == True:
                            reactants = self.react.equation[d].split('=')[0]
                            reactants = reactants.replace('<','').replace('(+','+').replace('( +','+')
                            reactants = reactants.split('+')
                            reactants_b = [react.replace(' ','') for react in reactants]
                            products  = self.react.equation[d].split('=')[1]
                            products  = products.replace('>','').replace('(+','+').replace('( +','+')
                            products  = products.split('+')
                            products_b  = [prod.replace(' ','') for prod in products]
                            if self.react.tbe[d] is not False:
                                tbe = self.react.tbe[d]
                                if tbe     in reactants_b:
                                    reactants_b.remove(tbe) ;        tbb = True
                                if tbe+")" in reactants_b:
                                    reactants_b.remove(tbe + ")") ;  tbb = True
                                if tbe     in products_b:
                                    products_b.remove(tbe) ;         tbb = True
                                if tbe+")" in products_b:
                                    products_b.remove(tbe + ")") ;   tbb = True
                            if "M"  in reactants_b:
                                reactants_b.remove("M")  ; tbb = True
                            if "M)" in reactants_b:
                                reactants_b.remove("M)") ; tbb = True
                                tbb = True
                            if "M"  in products_b:
                                products_b.remove("M")   ; tbb = True
                            if "M)" in products_b:
                                products_b.remove("M)")  ; tbb = True
                            common_reac_1 = set(reactants_a).intersection(reactants_b)
                            common_prod_1 = set(products_a).intersection(products_b)
                            common_reac_2 = set(products_a).intersection(reactants_b)
                            common_prod_2 = set(products_b).intersection(reactants_a)
                            if  len(common_reac_1) == len(reactants_a) == len(reactants_b)\
                            and len(common_prod_1) == len(products_a)  == len(products_b)\
                            and tba==tbb:
                                if self.react.tbe[d] is not False:  # remove duplicate species if third body contribution is neglectable
                                # ... otherwise it may prevent the withdrawal of the third body species
                                    if d in r_list:  duplicate_r = True
                                else:
                                    duplicate_r = True # add duplicate species without third body
                                    if d not in r_list: r_list.append(d)
                            if  len(common_reac_2) == len(products_a)  == len(reactants_b)\
                            and len(common_prod_2) == len(reactants_a) == len(products_b)\
                            and tba==tbb:
                                if self.react.tbe[d] is not False:
                                    if d in r_list:  duplicate_r = True
                                else:
                                    duplicate_r = True
                                    if d not in r_list: r_list.append(d)

            self.react.opt[r]['duplicate'] = duplicate_r

        r_list.sort()


        # https://cantera.org/documentation/docs-2.6/sphinx/html/cython/kinetics.html#chebyshevreaction
        for r in r_list: # in activated reactions list
            fd.write('- equation: ' + self.react.equation[r] \
                     + '  # Reaction ' + str(self.react.number[r]) + '\n')

            if   self.react.type[r] == 'three_body_reaction':
                fd.write('  type: three-body\n')
            elif self.react.type[r] == 'falloff_reaction':
                fd.write('  type: falloff\n')
            elif self.react.type[r] == 'chemically_activated_reaction':
                fd.write('  type: chemically-activated\n')
            elif self.react.type[r] == 'pdep_arrhenius':
                fd.write('  type: pressure-dependent-Arrhenius\n')
            elif self.react.type[r] == 'pdep_arrhenius':
                fd.write('  type: Chebyshev\n')

            r_kin, r_efficiencies, r_troe, r_opt, r_param = '', '', '', '', ''

            if self.react.type[r] == 'reaction' \
            or self.react.type[r] == 'three_body_reaction':
                if self.react.kin_units[r]['A'] is not False:
                      A_unit = ' ' + self.react.kin_units[r]['A']
                else: A_unit = ''
                if self.react.kin_units[r]['Ea'] is not False:
                      Ea_unit = ' ' + self.react.kin_units[r]['Ea']
                else: Ea_unit = ''
                r_kin = '  rate-constant: {'\
                  + 'A: '    + str('%0.6e' %self.react.kin[r][0]) + A_unit\
                  + ', b: '  + str('%0.4f' %self.react.kin[r][1])\
                  + ', Ea: ' + str('%0.3f' %self.react.kin[r][2]) + Ea_unit + '}\n'

            if any(opt is True for opt in self.react.opt[r].values()):
                if self.react.opt[r]['duplicate'] is True:
                    r_opt += '  duplicate: true\n'
                if self.react.opt[r]['negative-A'] is True:
                    r_opt += '  negative-A: true\n'
                if self.react.opt[r]['orders'] is True:
                    r_opt += '  orders: {' + self.react.opt[r]['orders'] + '}\n'
                if self.react.opt[r]['negative-orders'] is True:
                    r_opt += '  negative-orders: true\n'

            if self.react.type[r] == 'three_body_reaction'\
            or self.react.type[r] == 'falloff_reaction'\
            or self.react.type[r] == 'chemically_activated_reaction':
                if self.react.eff[r]!=False:
                    r_efficiencies = '  efficiencies: {'
                    count_sp = 0
                    for sp in sp_list:
                        for coll_sp in self.react.eff[r]:
                            if self.spec.name[sp] == coll_sp.split(":")[0]:
                                if count_sp == 6:
                                    r_efficiencies += '\n    '
                                    count_sp = 0
                                r_efficiencies+=coll_sp.replace(':',': ') + ", "
                                count_sp += 1
                    if r_efficiencies == '  efficiencies: {':    # if all collision species are removed
                        r_efficiencies = ""
                    else:
                        r_efficiencies = r_efficiencies[:-2] + '}\n'
                else:
                    r_efficiencies = ''

            if self.react.type[r] == 'falloff_reaction'\
            or self.react.type[r] == 'chemically_activated_reaction':
                r_kin = ''
                for p in [0,1]:
                    if self.react.kin_units[r][p]['A'] is not False:
                          A_unit = ' ' + self.react.kin_units[r][p]['A']
                    else: A_unit = ''
                    if self.react.kin_units[r][p]['Ea'] is not False:
                          Ea_unit = ' ' + self.react.kin_units[r][p]['Ea']
                    else: Ea_unit = ''
                    if p==0:
                        prefix_reac = '  low-P-rate-constant: {'
                        pk = 1
                    else :
                        prefix_reac = '  high-P-rate-constant: {'
                        pk=0
                    r_kin += prefix_reac\
                      + 'A: '    + str('%0.6e' %self.react.kin[r][pk][0]) + A_unit\
                      + ', b: '  + str('%0.4f' %self.react.kin[r][pk][1])\
                      + ', Ea: ' + str('%0.3f' %self.react.kin[r][pk][2]) + Ea_unit + '}\n'

                if self.react.param[r]:
                    r_troe = '  Troe: {' \
                        + 'A: '    + str(self.react.param[r]['A'])\
                        + ', T3: ' + str(self.react.param[r]['T3'])\
                        + ', T1: ' + str(self.react.param[r]['T1'])
                    if self.react.param[r]['T2'] is not False:
                        r_troe += ', T2: ' + str(self.react.param[r]['T2']) + '}\n'
                    else:
                        r_troe += '}\n'

            if self.react.type[r] == 'pdep_arrhenius':
                r_kin = '  rate-constants:\n'

                for p in range(len(self.react.kin[r])):
                    if self.react.param[r][1][p] is not False:
                          p_unit = ' ' + self.react.param[r][1][p]
                    else: p_unit = ''
                    if self.react.kin_units[r][p]['A'] is not False:
                          A_unit = ' ' + self.react.kin_units[r][p]['A']
                    else: A_unit = ''
                    if self.react.kin_units[r][p]['Ea'] is not False:
                          Ea_unit = ' ' + self.react.kin_units[r][p]['Ea']
                    else: Ea_unit = ''
                    r_kin += '  - {P: ' + str(self.react.param[r][0][p]) + p_unit \
                      + ', A: '  + str('%0.6e' %self.react.kin[r][p][0]) + A_unit\
                      + ', b: '  + str('%0.4f' %self.react.kin[r][p][1])\
                      + ', Ea: ' + str('%0.3f' %self.react.kin[r][p][2]) + Ea_unit + '}\n'

            if self.react.type[r] == 'chebyshev':
                # T range
                r_param = '  temperature-range: [' + str(self.react.param[r]['Tmin'])\
                    + ', ' + str(self.react.param[r]['Tmax']) + ']\n'
                # Pres range
                r_param += '  pressure-range: [' + str('%0.6e' %self.react.param[r]['Pmin_value'])
                if self.react.param[r]['Pmin_unit'] is not False:
                    r_param += ' ' + self.react.param[r]['Pmin_unit']
                r_param += ', ' + str('%0.6e' %self.react.param[r]['Pmax_value'])
                if self.react.param[r]['Pmax_unit'] is not False:
                    r_param += ' ' + self.react.param[r]['Pmax_unit']
                r_param += ']\n'
                # data
                r_param += '  data:\n'
                r_param += '  - ['
                for _d in range(len(self.react.kin[r])):
                    for _di in range(len(self.react.kin[r][_d])):
                        r_param += str('%0.6e' %self.react.kin[r][_d][_di])
                        if _di != len(self.react.kin[r][_d])-1: r_param += ', '
                        elif _d != len(self.react.kin[r])-1: r_param += ']\n  - ['
                        else: r_param += ']\n'



            fd.write(r_kin)
            fd.write(r_efficiencies)
            fd.write(r_troe)
            fd.write(r_opt)
            fd.write(r_param)

        timer.sleep(1.5)

        fd.close()






    def write_chemkin_mech(self,filename,version,act_sp="no_arg",act_r="no_arg"):

        # ----- Raise error if thermo data in NASA9 format
        def stop_program():
            raise Exception("Cannot write kinetic mechanism in chemkin format with NASA9 formalism")
        try:
            for therm_formalism in self.spec.thermo_model:
                if therm_formalism == 'NASA9':
                    stop_program()
        except Exception as e:
            print(f"Error: {e}")
        # -----

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

        fd.write("!-------------------------------------------------------------------------------\n")
        fd.write("! Kinetic mechanism converted from the Cantera format (.cti) by Brookesia "+version+"\n")
        fd.write("! https://github.com/Brookesia-py/Brookesia\n!\n")

        fd.write("!-------------------------------------------------------------------------------\n")
        fd.write("! Species and Elements\n")
        fd.write("!-------------------------------------------------------------------------------\n")

        #elements
        txt_elements = 'ELEMENTS\n'#self.gas_prop[l].split('"')[0]+'"'
        elements = ["C","H","O","N","AR","HE","Ar","He"]
        for el in elements:
            if self.find_element(el,act_sp):
                txt_elements += str.capitalize(el) + " "
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
#         Thermo
# =============================================================================


        fd.write("!-------------------------------------------------------------------------------\n")
        fd.write("! Thermophysical properties\n")
        fd.write("!-------------------------------------------------------------------------------\n")
        fd.write("THERMO ALL\n   300.0  1000.0  5000.0 \n")

        for sp in sp_list: # in activated species list

            # -------------    1st line    -------------
            l_name  = self.spec.name[sp]
            for i in range(24-len(self.spec.name[sp])):
                l_name += ' '
            atoms = self.spec.atoms[sp].split(' ')
            l_atoms = ''
            for atom in atoms:
                l_atoms += atom.split(':')[0] + ' ' + atom.split(':')[1] + '  '
            for i in range(20-len(l_atoms)):
                l_atoms += ' '
            l_atoms +='G   '
            if float(self.spec.thermo_coeff_lT[sp][0])<1000: l_atoms+='  '
            else:                                   l_atoms+=' '
            l_T =  '%0.1f' %self.spec.thermo_temp[sp][0]  + '    ' + \
                   '%0.1f' %self.spec.thermo_temp[sp][1]  + '  ' + \
                   '%0.1f' %self.spec.thermo_temp[sp][2]
            for i in range(29-len(l_T)):
                l_T += ' '
            l_T += '1\n'
            first_line = l_name+l_atoms+l_T

            # -------------    2nd line    -------------
            second_line = ''
            for i in range (5):
                if self.spec.thermo_coeff_hT[sp][i]>=0:  second_line += ' '
                second_line += str('%0.8e' %self.spec.thermo_coeff_hT[sp][i]).upper()
            second_line += '    2\n'

            # -------------    3rd line - high T coeffs  -------------
            third_line = ''
            for i in range (2):
                if self.spec.thermo_coeff_hT[sp][5+i]>=0:   third_line += ' '
                third_line += str('%0.8e' %self.spec.thermo_coeff_hT[sp][5+i]).upper()
            for i in range (3):
                if self.spec.thermo_coeff_lT[sp][i]>=0:
                    third_line += ' '
                third_line += str('%0.8e' %self.spec.thermo_coeff_lT[sp][i]).upper()
            third_line += '    3\n'

            # -------------    4th line    -------------
            fourth_line = ''
            for i in range(4):
                if self.spec.thermo_coeff_lT[sp][3+i]>=0:   fourth_line += ' '
                fourth_line += str('%0.8e' %self.spec.thermo_coeff_lT[sp][3+i]).upper()
            fourth_line += '                   4\n'

            fd.write(first_line)
            fd.write(second_line)
            fd.write(third_line)
            fd.write(fourth_line)

        fd.write('END\n!\n')


# =============================================================================
#         Transport
# =============================================================================


        if self.spec.trans_model is not False:
            fd.write("!-------------------------------------------------------------------------------\n")
            fd.write("! Transport properties\n")
            fd.write("!-------------------------------------------------------------------------------\n")
            fd.write("TRANSPORT\n")

            for sp in sp_list: # in activated species list
                line_trans  = self.spec.name[sp] + ' '
                if len(line_trans)<14:
                    for i in range(14-len(line_trans)): line_trans+=' '
                else:
                    line_trans += ' '
                trans_diam, trans_wdepth, trans_polar, trans_rrelax, trans_dipole = '0.0  ', '0.0  ', '0.0  ', '0.0', '0.0  '

                # geometry
                if 'atom' in self.spec.trans_geom[sp]:
                    trans_geom = '0  '
                elif 'nonlinear' in self.spec.trans_geom[sp]:
                    trans_geom = '2  '
                else:
                    trans_geom = '1  '
                # well_depth
                if self.spec.trans_wd[sp] is not False:
                    trans_wdepth = self.spec.trans_wd[sp]
                    if len(trans_wdepth)<7:
                        for i in range(7-len(trans_wdepth)): trans_wdepth+=' '
                    else:
                        trans_wdepth += ' '
                # diameter
                if self.spec.trans_diam[sp] is not False:
                    trans_diam = self.spec.trans_diam[sp]
                    if len(trans_diam)<7:
                        for i in range(7-len(trans_diam)): trans_diam+=' '
                    else:
                        trans_diam += ' '
                # dipole
                if self.spec.trans_dipole[sp] is not False:
                    trans_dipole = self.spec.trans_dipole[sp]
                    if len(trans_dipole)<7:
                        for i in range(7-len(trans_dipole)): trans_dipole+=' '
                    else:
                        trans_dipole += ' '
                # polar
                if self.spec.trans_polar[sp] is not False:
                    trans_polar = self.spec.trans_polar[sp]
                    if len(trans_polar)<7:
                        for i in range(7-len(trans_polar)): trans_polar+=' '
                    else:
                        trans_polar += ' '
                # rot-relax
                if self.spec.trans_rot_relax[sp] is not False:
                    trans_rrelax = self.spec.trans_rot_relax[sp]

                line_trans += trans_geom   + \
                              trans_wdepth + \
                              trans_diam   + \
                              trans_dipole + \
                              trans_polar  + \
                              trans_rrelax + '\n'
                fd.write(line_trans)
            fd.write('END\n!\n')



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
                    if self.react.equation[d]==self.react.equation[r]:
                        if d not in r_list:
                            r_list.append(d)
        r_list.sort()

        for r in r_list: # in activated reactions list
            dup_txt = ''
            txt_space = ''
            for ts in range(50-len(self.react.equation[r])):
                txt_space += ' '

            if self.react.opt[r]['duplicate'] is not False:
                dup_txt = ' DUPLICATE\n'

            if self.react.type[r]=="reaction":
                r_line = self.react.equation[r] + txt_space              \
                         + str('%0.6e' %self.react.kin[r][0]) + "  "    \
                         + str('%0.4f' %self.react.kin[r][1]) + "  "    \
                         + str('%0.3f' %self.react.kin[r][2]) + "\n"

                fd.write(r_line)

            if self.react.type[r]=="three_body_reaction":
                r_line = self.react.equation[r]+txt_space \
                         + str('%0.6e' %self.react.kin[r][0]) + "  "          \
                         + str('%0.4f' %self.react.kin[r][1]) + "  "          \
                         + str('%0.3f' %self.react.kin[r][2]) + "\n"
                if self.react.eff[r] is not False :
                    r_efficiencies=""
                    if self.react.eff[r]!=[]:
                        r_efficiencies=""
                        for sp in sp_list:
                            for coll_sp in self.react.eff[r]:
                                if self.spec.name[sp] == coll_sp.split(":")[0]:
                                    r_efficiencies+=coll_sp.split(":")[0]+'/'\
                                                    +coll_sp.split(":")[1]+'/ '
                        r_efficiencies+="\n" ;
                else: r_efficiencies=""

                fd.write(r_line)
                fd.write(r_efficiencies)

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
                r_line = self.react.equation[r]+txt_space        \
                + str('%0.6e' %self.react.kin[r][0][0]) + '  '  \
                + str('%0.4e' %self.react.kin[r][0][1]) + '  '  \
                + str('%0.3e' %self.react.kin[r][0][2]) + '\n'  \
                + '     LOW /  '                                \
                + str('%0.6e' %self.react.kin[r][1][0]) + '  '  \
                + str('%0.4e' %self.react.kin[r][1][1]) + '  '  \
                + str('%0.3e' %self.react.kin[r][1][2]) + ' /\n'
                # Troe
                r_txt = ""
                if self.react.param[r] is not False:
                    r_txt += '     TROE/   ' \
                        + self.react.param[r]['A'] + '  ' \
                        + self.react.param[r]['T3'] + '  ' \
                        + self.react.param[r]['T1'] + '  '

                    if self.react.param[r]['T2'] is not False:
                        r_txt += self.react.param[r]['T2'] + '/\n'
                    else:
                        r_txt += '/\n'
                # efficiencies
                r_efficiencies=""
                if self.react.eff[r] is not False:
                    for sp in sp_list:
                        for coll_sp in self.react.eff[r]:
                            if self.spec.name[sp] == coll_sp.split(":")[0]:
                                r_efficiencies+=coll_sp.split(":")[0]+'/'\
                                                +coll_sp.split(":")[1]+'/ '
                    r_efficiencies+="\n"

                fd.write(r_line)
                fd.write(r_txt)
                fd.write(r_efficiencies)

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

                r_line = self.react.equation[r]+txt_space \
                             + '%0.3e' %self.react.kin[r][0][0] + "  "  \
                             + '%0.3f' %self.react.kin[r][0][1] + "  "  \
                             + '%0.1f' %self.react.kin[r][0][2] + "\n"
                k_line = ""
                for p in range(len(self.react.param[r][0])):
                    if float(self.react.param[r][0][p])<9:   txt_space = '      ';
                    elif float(self.react.param[r][0][p])<99: txt_space = '     ';
                    elif float(self.react.param[r][0][p])<999: txt_space = '    ';
                    elif float(self.react.param[r][0][p])<9999: txt_space = '   ';
                    elif float(self.react.param[r][0][p])<99999: txt_space = '  ';
                    k_line += "PLOG/" + txt_space                             \
                             + '%0.4f' %float(self.react.param[r][0][p]) + "     " \
                             + '%0.3e' %self.react.kin[r][p][0]     + "     " \
                             + '%0.3f' %self.react.kin[r][p][1]     + "     " \
                             + '%0.1f' %self.react.kin[r][p][2]     + "/\n"

                fd.write(r_line)
                fd.write(k_line)

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

                r_line = self.react.equation[r]+txt_space        \
                + str('%0.6e' %self.react.kin[r][1][0]) + '  '  \
                + str('%0.4e' %self.react.kin[r][1][1]) + '  '  \
                + str('%0.3e' %self.react.kin[r][1][2]) + '\n'  \
                + '     HIGH/  '                                \
                + str('%0.6e' %self.react.kin[r][0][0]) + '  '  \
                + str('%0.4e' %self.react.kin[r][0][1]) + '  '  \
                + str('%0.3e' %self.react.kin[r][0][2]) + ' /\n'
                # Troe
                r_txt = ""
                if self.react.param[r] is not False:
                    r_txt += '     TROE/   ' \
                        + self.react.param[r]['A'] + '  ' \
                        + self.react.param[r]['T3'] + '  ' \
                        + self.react.param[r]['T1'] + '  '

                    if self.react.param[r]['T2'] is not False:
                        r_txt += self.react.param[r]['T2'] + '/\n'
                    else:
                        r_txt += '/\n'
                # efficiencies
                r_efficiencies=""
                if self.react.eff[r] is not False:
                    for sp in sp_list:
                        for coll_sp in self.react.eff[r]:
                            if self.spec.name[sp] == coll_sp.split(":")[0]:
                                r_efficiencies+=coll_sp.split(":")[0]+'/'\
                                                +coll_sp.split(":")[1]+'/ '
                    r_efficiencies+="\n"

                fd.write(r_line)
                fd.write(r_txt)
                fd.write(r_efficiencies)


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

                r_line = self.react.equation[r] + ' 1.0E0 0.0 0.0\n'
                # T / P range
                k_line =  'TCHEB / '                                                 \
                        + self.react.param[r]['Tmin']                       + ' '    \
                        + self.react.param[r]['Tmax']                       + ' / '  \
                        + 'PCHEB / '                                                 \
                        + str('%0.6e' %self.react.param[r]['Pmin_value'])   + ' '    \
                        + str('%0.6e' %self.react.param[r]['Pmax_value'])   + ' /\n'
                # data
                k_line += "CHEB / "
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

            fd.write(dup_txt)

        fd.write('END')
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
        self.react.activ_m = [False]*len(self.react.equation)   # act react for methods loop

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
         selection_operator='Geometric_norm',      selection_options=[0.2],   \
         Xover_operator=['simple_Xover', 'multiple_Xover', 'arith_Xover', 'heuristic_Xover'], \
         Xover_pct = [10,10,20,20],   \
         mut_operator=['uniform_mutation', 'non_uniform_mutation', 'boundary_mutation'], \
         mut_pct=[10,15,15], mut_intensity=30, mut_option=[0,3,0],\
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
        self.inertia_score        = False
        self.inertia_min          = 0.2
        self.inertia_i            = 0.8
        self.inertia_end          = 0.3
        self.cognitive_accel_i    = 1.6
        self.cognitive_accel_end  = 2.0
        self.social_accel_i       = 0.5
        self.social_accel_end     = 2

        self.optim_on_meth        = optim_on_meth
        self.display_react2opt    = True
        self.target_r             = []
        if type(Arrh_max_variation) is list:
            self.Arrh_max_variation      = Arrh_max_variation
        else:
            self.Arrh_max_variation      = [Arrh_max_variation]*3
        self.Arrh_var             = False
        self.f_default            = 0.7
        self.diff_r_meth          = False
        self.r2opt                = False
        self.tspc                 = tsp
        self.n_tspc               = n_tspc
        self.coeff_s              = [1]*len(tsp)
        self.coeff_T              = 1
        self.coeff_ig             = 1
        self.coeff_Sl             = 1
        self.coeff_K              = 1
        self.coeff_cond           = False
        self.genVec               = []
        self.best_fitness         = []
        self.mean_fitness         = []
        self.worst_fitness         = []
        self.main_path            = ''
        self.exp_data             = False
        self.nb_r2opt             = nb_r2opt # number of reaction to optimized when optimization based on DRG/SA is selected
        self.reactions2opt        = False
        self.opt_subm_C           = [True]*30
        self.opt_subm_CO          = True
        self.opt_subm_N           = [True]*30
        self.opt_subm_S           = [True]*30
        self.opt_subm_Si          = [True]*30
        self.import_mech          = False
        self.keep_until_gen       = int(n_gen * 2/3)



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
        self.equation  = []
        self.opt       = []
        self.kin       = []
        self.ref_kin   = []
        self.kin_units = []
        self.param     = []
        self.eff       = []
        self.tbe       = [] # three body exception in falloff reactions (+ AR) / (+HE) ...
        self.activ_p   = []
        self.activ_m   = []
        self.modif     = []
        self.incert    = []
        self.f_max     = []
        self.f_min     = []
        self.f_T_fit   = []
        self.f_Tit     = [300,2500,200]
        self.k0_fit    = []
        self.f_max_lp  = []
        self.f_T_fit_lp= []
        self.k0_fit_lp = []
        self.f_max_hp  = []
        self.f_T_fit_hp= []
        self.k0_fit_hp = []

        self.k_ref     = []
        self.k_min     = []
        self.k_max     = []
        self.k_ref_lp  = []
        self.k_min_lp  = []
        self.k_max_lp  = []
        self.k_ref_hp  = []
        self.k_min_hp  = []
        self.k_max_hp  = []





class Red_data :
    def __init__(self,gas_ref,mech_name,tspc,n_tspc,reduction_operator='DRG_sp',   \
                 optim = 'False',verbose=6,targetSpeciesIdx=[],gas_loop='gas_ref',gas_act='gas_ref'):

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



class Species:
    def __init__(self,verbose=0):

        self.name            = []
        self.atoms           = []
        self.thermo_temp     = []
        self.thermo_coeff_lT = []
        self.thermo_coeff_hT = []
        self.thermo_model    = []
        self.trans_model     = []
        self.trans_geom      = []
        self.trans_wd        = []
        self.trans_diam      = []
        self.trans_dipole    = []
        self.trans_polar     = []
        self.trans_rot_relax = []
        self.note_therm      = []
        self.note_trans      = []
        self.activ_p         = []
        self.activ_m         = []
        self.CSP_radicals    = []


class Simul_param :
    def __init__(self, pts_scatter = [0.0, 0.03, 0.3, 0.5, 0.7, 1.0],        \
                 end_sim = 0.02,                                             \
                 tol_ss = [1.0e-5, 1.0e-8],                                  \
                 verbose = 0,                                                \
                 tign_nPoints = 450, tign_dt = 1e-09,                        \
                 n_pts = 250, delta_npts = 20, t_max_coeff = 5,              \
                 Scal_ref = 'H2O', grad_curv_ratio = 0.5):

        self.pts_scatter      = pts_scatter      # time stepping  or  grid
        self.end_sim          = end_sim          # tmax           or  xmax
        self.rtol_ts           = 1e-8
        self.atol_ts           = 1e-13
        self.tol_ss           = tol_ss
        self.tol_ts_f           = [1.0e-4, 1.0e-8]
        self.tign_nPoints     = tign_nPoints
        self.tign_dt          = tign_dt
        self.n_pts            = n_pts
        self.delta_npts       = delta_npts
        self.t_max_coeff      = t_max_coeff
        self.t_max_s          = False
        self.t_max_react      = 1
        self.Scal_ref         = Scal_ref
        self.grad_curv_ratio  = grad_curv_ratio
        self.mdot             = 1   # kg/m^2/s
        self.mdot2            = 3   # kg/m^2/s
        self.verbose          = verbose
        self.show_plots       = False
        self.write_ck         = False
        self.transport_model  = 'mixture-averaged'
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
            or 'burner' in self.conditions.config\
            or 'diff_' in self.conditions.config\
            or 'tp_'   in self.conditions.config\
            or 'pp_'   in self.conditions.config:
                l2_conditions+=str(self.conditions.state_var.T)+";"\
                        +str(self.conditions.simul_param.rtol_ts)+";"\
                        +str(self.conditions.simul_param.atol_ts)+";"\
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
                        +str(self.conditions.simul_param.rtol_ts)+";"\
                        +str(self.conditions.simul_param.atol_ts)+";"
                if self.ign_time == '':
                    if self.ign_time_hr != '' and self.ign_time_hr != False:
                        l2_conditions+=str(self.ign_time_hr)
                    else :
                        l2_conditions+=str(self.ign_time_sp)
                else:
                    l2_conditions+=str(self.ign_time)
            elif "JSR" in self.conditions.config:
                l2_conditions+=str(self.conditions.simul_param.end_sim)+";"\
                        +str(self.conditions.simul_param.rtol_ts)+";"\
                        +str(self.conditions.simul_param.atol_ts)
            elif "PFR" in self.conditions.config:
                l2_conditions+=str(self.conditions.state_var.T)+";"\
                        +str(self.conditions.simul_param.rtol_ts)+";"\
                        +str(self.conditions.simul_param.atol_ts)+";"\
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
        if self.conditions.write_conc is True:
            nline_conc = len(self.pts_scatter)
        else:
            nline_conc = 1
        for nb_line in range(nline_conc):
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


def generate_plot_style(plt, fig,  plot_style,delta_polices=[0,0,0,0]) :

    from matplotlib import rcParams

    delta_police_legend = delta_polices[0]
    delta_police_axes   = delta_polices[1]
    delta_police_ticks  = delta_polices[2]
    delta_police_title  = delta_polices[3]


    if plot_style=='presentation':
        plt.style.use('seaborn-v0_8-talk') # 'seaborn-talk'
        fig.subplots_adjust(left = 0.1, bottom = 0.15,
                right = 0.95, top = 0.93, wspace = 0, hspace = 0.5)

        #rcParams['font.family'] = 'sans-serif'
        #rcParams['font.sans-serif'] = ['Verdana']   #['Tahoma', 'DejaVu Sans','Lucida Grande', 'Verdana']
        plt.rc('font', size=17)
        plt.rc('text', usetex=False)

        title_size   = 17 + delta_police_title
        axes_size    = 17 + delta_police_axes
        tick_size    = 14 + delta_police_ticks
        legend_size  = 14 + delta_police_legend

    elif 'paper' in plot_style:
        plt.style.use('seaborn-v0_8-paper')
        fig.subplots_adjust(left = 0.25, bottom = 0.2,
                right = 0.95, top = 0.93, wspace = 0, hspace = 0.5)

        plt.rc('font', family='serif', size=15)
        if 'latex' in plot_style:
            plt.rc('text', usetex=True)
        else:
            plt.rc('text', usetex=False)

        title_size   = 15 + delta_police_title
        axes_size    = 15 + delta_police_axes
        tick_size    = 12 + delta_police_ticks
        legend_size  = 12 + delta_police_legend


    elif plot_style=='poster':
        plt.style.use('seaborn-v0_8-paper')
        fig.subplots_adjust(left = 0.13, bottom = 0.15,
                right = 0.95, top = 0.93, wspace = 0, hspace = 0.5)
        plt.rc('font', size=19)
        plt.rc('text', usetex=False)

        title_size   = 19 + delta_police_title
        axes_size    = 19 + delta_police_axes
        tick_size    = 16 + delta_police_ticks
        legend_size  = 10 + delta_police_legend

    parameters = {'axes.labelsize': axes_size,
                  'axes.titlesize': title_size,
                  'xtick.labelsize': tick_size,
                  'ytick.labelsize': tick_size,
                  'legend.fontsize': legend_size}
    plt.rcParams.update(parameters)

    #axes.legend(fontsize=legend_size)

    #axes.tick_params(axis='x', labelsize=tick_size)
    #axes.tick_params(axis='y', labelsize=tick_size)

    plt.rc('xtick', direction='in')
    plt.rc('ytick', direction='in')

    return plt,fig


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
    warnings.filterwarnings("ignore", category=ResourceWarning)

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

        # progress bar deactivated on Windows
        if os.name != 'nt':
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


def get_react_species(str_react):
    # remove + M or (+ M)
    str_react_wt = re.sub(r'\(\s*\+\s*M\s*\)|(?<=\s)\+\s*M(?=\s|$)', '', str_react)
    # remove spaces
    # str_react_wt = ' '.join(str_react_wt.split())
    str_react = str_react_wt.replace(' ','')
    # split reactants
    react_list = str_react_wt.split('+')
    return react_list


def get_screen_size():
    from PyQt5 import QtWidgets

    _height = QtWidgets.QDesktopWidget().screenGeometry(-1).height()

    return _height,_height

def get_gas_ct(mech):

    ct.suppress_thermo_warnings()
    warnings.filterwarnings("ignore", category=DeprecationWarning)

    error           = True
    reaction_error  = False
    while error:
        try:
            if float(ct.__version__[0:3])>2.4:
                gas = ct.Solution(mech, transport_model = 'mixture-averaged')
            else:
                gas = ct.Solution(mech, transport_model = 'Mix')
            error = False
        except ct.CanteraError as e:

            # Catch the exception and capture the error message
            error_message = str(e)
            print(error_message)

            if 'duplicate' in error_message:

                # Extract reaction numbers from the error message
                # if '.cti' in mech:
                #     reaction_numbers = re.findall(r'Reaction (\d+):', error_message)
                #     reaction_numbers = [int(num) for num in reaction_numbers]  # Convert to integers
                # elif '.yaml' in mech:
                reaction_numbers = re.findall(r'reaction number (\d+) ', error_message)
                reaction_numbers = [int(num) for num in reaction_numbers]  # Convert to integers
                if reaction_numbers==[]:
                    reaction_numbers = re.findall(r'Reaction (\d+):', error_message)
                    reaction_numbers = [int(num) for num in reaction_numbers]  # Convert to integers

                if reaction_numbers==[]:
                    print(error_message)
                    break
                else:
                    reaction_error = True

                    # Read the file content
                    with open(mech, 'r') as f:
                        lines = f.readlines()

                    if '.cti' in mech:
                        # Regex pattern to match reactions and capture the reaction number
                        reaction_pattern = re.compile(r"# Reaction (\d+)")

                        # Go through the file and add 'options="duplicate"' to the identified reactions
                        for i, line in enumerate(lines):
                            match = reaction_pattern.match(line)
                            if match:
                                reaction_number = int(match.group(1))
                                # Check if the reaction number is in the list from the error message
                                if reaction_number in reaction_numbers:
                                    dup_declared = False
                                    # Find the next line where the actual reaction is defined (i.e., contains 'reaction' or 'three_body_reaction')
                                    for j in range(i + 1, len(lines)):
                                        if 'duplicate' in lines[j]:
                                            dup_declared = True
                                        if ')\n' in lines[j]:
                                            closure_line = np.int32(j)
                                        if '# Reaction' in lines[j] and not dup_declared \
                                            or j==len(lines) and not dup_declared :
                                            lines[closure_line] = lines[closure_line].replace(')\n',",\n         options='duplicate')\n")
                                            break

                    elif '.yaml' in mech:
                        # Go through the file and add 'options="duplicate"' to the identified reactions
                        r_nb_m = 0
                        for i, line in enumerate(lines):
                            if '- equation: ' in line:
                                r_nb_m += 1
                                if r_nb_m in reaction_numbers:
                                    # if reaction_number in reaction_numbers:
                                    dup_declared = False
                                    for j in range(i + 1, len(lines)):
                                        if 'duplicate' in lines[j]:
                                            dup_declared = True
                                        if '- equation: ' in lines[j] and not dup_declared :
                                            lines[j-1] = lines[j-1] + "  duplicate: true\n"
                                            break

                    print('(mech: ' + mech + ')  Undeclared duplicate reactions: ' + str(reaction_numbers) + ' -> corrected')
                    # Write the modified content to a new file
                    with open(mech, 'w') as f:
                        f.writelines(lines)

            elif 'negative pre-exponential' in error_message:
                if '.yaml' in mech:
                    # Extract reaction equation from the error message
                    reaction_eq = error_message.split('> - equation: ')[1].split('#')[0].split('\n')[0]

                    # Read the file content
                    with open(mech, 'r') as f:
                        lines = f.readlines()

                    # Go through the file and add 'options="duplicate"' to the identified reactions
                    r_nb_m = 0
                    for i, line in enumerate(lines):
                        if '- equation: ' + reaction_eq in line:
                            r_nb_m += 1
                            NegA_declared = False
                            for j in range(i + 1, len(lines)):
                                if 'negative-A: true' in lines[j]:
                                    NegA_declared = True
                                if 'rate-constant: {A:' in lines[j]:
                                    if float(lines[j].split('A:')[1].split(',')[0])>=0:
                                        NegA_declared = True
                                if '- equation: ' in lines[j] and not NegA_declared :
                                    lines[j-1] = lines[j-1] + '  negative-A: true\n'
                                elif '- equation: ' in lines[j]:
                                    break

                print('(mech: ' + mech + ')  Undeclared negative pre-exponential factor: ' \
                       + str(reaction_eq) + str(' -> corrected'))
                # Write the modified content to a new file
                with open(mech, 'w') as f:
                    f.writelines(lines)

            else:
                error = False

            # Pause for 1 second
            timer.sleep(1)

    return gas
