
import cantera as ct
#import pyjacob
import numpy as np
<<<<<<< HEAD
import brookesia.Class_def   as cdef
import brookesia.Computation as comp
import brookesia.SA          as sa
from  brookesia.Class_def import print_
=======
import brookesia_devLOI.Class_def   as cdef
import brookesia_devLOI.Computation as comp
import brookesia_devLOI.SA as sa
from  brookesia_devLOI.Class_def import print_
>>>>>>> origin/LOI
import pandas as pd
import copy
import sys
import os
import time as timer

#from brookesia_devLOI.SA import sensitivities_computation_SA




def LOI_computation(red_data, mech_data, red_results): 
    """
    Level of importance computation
    """


    # gas_ref = red_data.gas_ref   ## utiliser pour comparer ou normaliser
    # gas = ct.Solution('gri30.yaml')  ## utiliser pour calculer 
    # gas = red_results.gas## Objet Cantera modifié au fil du temps dans Brookesia
    conditions = red_results.conditions 
    mp = conditions.main_path
    gas_ref = conditions.composition.gas_ref
    
    # gas_red = gas_ref ; print_('gas_red = gas_ref',mp)
    gas_red = red_results.gas #; print_('gas_red_ok',mp)
    
    
    # gas = red_data.gas_ref
    n_sp = gas_ref.n_species  # n_sp_ref = red_data.gas_ref.n_species
    # J_diag = np.zeros(n_sp)
    

    T = conditions.state_var.T
    P = conditions.state_var.P
    X_red = conditions.composition.X
    gas_red.TPX = conditions.state_var.T, conditions.state_var.P, X_red    
    
    pts_scatter = red_results.pts_scatter
    n_points = len(pts_scatter)
    n_points_LOI = int(red_data.red_op.n_points)
    
    
    
    # =============================================================================
    # 1- timescales computation
    # =============================================================================
        
    timescales_all = np.ones((n_points_LOI, n_sp))
       
    
    if "reactor" in conditions.config:
        
        # parametres de tolerance
        abs_tol = conditions.simul_param.atol_ts
        rel_tol = conditions.simul_param.rtol_ts
        
        
        if conditions.config == "reactor_UV":
            reactor = ct.IdealGasReactor(gas_red)
        elif conditions.config == "reactor_HP":
            reactor = ct.IdealGasConstPressureReactor(gas_red)
        sim = ct.ReactorNet([reactor])

        # tolerances
        sim.rtol = rel_tol
        sim.atol = abs_tol

        t_i = -1
        for t in range(n_points-int(n_points/red_data.red_op.n_points+1)):
            sim.advance(pts_scatter[t])
            if t != 0 and t % int(max(n_points/red_data.red_op.n_points, 1)) == 0:
                t_i += 1

    
            #on initialise le vecteur timescales avec des 1 (pour pouvoir diviser par Jii ensuite)
            timescales = np.ones(n_sp)
            omega_ref  = np.zeros(n_sp)
            conc_ref = gas_red.concentrations
            
            for idx_mech1 in range(n_sp) : 
                try:
                    # ----- if spec i in gas_red:
                    
                    idx_red = gas_red.species_index(gas_ref.species_name(idx_mech1))
                    # conc_ref[i]   = gas_red.concentrations[idx_sp_in_red]
                    omega_ref[idx_mech1]  = gas_red.net_production_rates[idx_red]
                    conc_i        = conc_ref[idx_red]
                    
                    #si la concentration est trop faible tau=1 par defaut
                    if conc_i > abs_tol:
                        #calcul du pas de perturbation
                        dc = abs_tol + conc_i * rel_tol
                        
                        #etat perturbé
                        conc_pert = conc_ref.copy()
                        conc_pert[idx_red] += dc
                        
                        #mise a jour de l'etat cantera
                        gas_red.concentrations = conc_pert
                        omega_pert_sp = gas_red.net_production_rates[idx_red]
                        
                        
                        #approximation de la dérivée (Jacobien diagonal)
                        J_diag_sp = (omega_pert_sp - omega_ref[idx_mech1]) / dc
                        
                        #calcul implicite du timescale avec borne max a 1.0
                        if abs(J_diag_sp) > 1 : 
                            timescales[idx_mech1] = 1.0 / abs(J_diag_sp)
                        else : 
                            timescales[idx_mech1] = 1.0                        
                        
                        
                        # retrait de la perturbation
                        gas_red.concentrations = conc_ref
                            
                    
                except ValueError: #if spec not in red mech
                    pass
                                        
                    
                timescales_all[t_i] = timescales.copy()
            
    
    elif "flame" in conditions.config:
        
        if conditions.mech_prev_red:
            # in this configuration, the flame object red_results.f come from 
            # the reference simulation (as defined in the first lines of 
            # the function reduction() in GeneralFunction.py
            # => calculation of the flame object from the reduced mech
            red_results = comp.red_computation(
                conditions,
                gas_red,
                mech_data.spec.activ_p,
                mech_data.react.activ_p)

        f = red_results.f
        T_prof = f.T

        
        # parametres de tolerance        
        abs_tol = conditions.simul_param.tol_ss[1]
        rel_tol = conditions.simul_param.tol_ss[0]
        
        
        z_i = 0
        for z in range(n_points):
            if z % int(max(n_points/n_points_LOI, 1)) == 0 \
               and 0.01*(max(red_results.T)-min(red_results.T)) < T_prof[z]-T_prof[0] < 0.99*(max(red_results.T)-min(red_results.T)):
                #sensi_scatter.append(z)
                f.set_gas_state(z)
                
                # conc_ref = f.gas.concentrations.copy()
                # omega_ref = f.gas.net_production_rates.copy()
                timescales = np.ones(n_sp)
                conc_ref = f.gas.concentrations
                omega_ref = np.zeros(n_sp)
                
                for idx_mech1 in range(n_sp):
                    try:
                        idx_red = f.gas.species_index(gas_ref.species_name(idx_mech1))
                        omega_ref[idx_mech1]  = f.gas.net_production_rates[idx_red]                    
                        conc_i = conc_ref[idx_red]
                        
                        if conc_i > abs_tol:                        
                            dc = abs_tol + conc_i * rel_tol
                            conc_pert = conc_ref.copy()
                            conc_pert[idx_red] += dc
                            
                            #mise a jour de l'etat cantera
                            f.gas.concentrations = conc_pert
                            omega_pert_sp = f.gas.net_production_rates[idx_red]
                            J_diag = (omega_pert_sp - omega_ref[idx_mech1]) / dc
                            
                            #calcul implicite du timescale 
                            if abs(J_diag) > 0 : 
                                timescales[idx_mech1] = 1.0 / abs(J_diag)
                            else : 
                                timescales[idx_mech1] = 1.0                        
                            # if J_ii != 0 and abs(J_ii) < 1e30:
                            #     timescales[idx_mech1] = 1.0 / abs(J_ii)
                            # if timescales[idx_mech1] > 1.0:
                            #     timescales[idx_mech1] = 1.0
                            
                    except ValueError: #if spec not in red mech
                        pass
                            
                timescales_all[z_i] = timescales.copy()
                z_i += 1
   
    
    elif 'JSR' in conditions.config:
        
        abs_tol = conditions.simul_param.atol_ts
        rel_tol = conditions.simul_param.rtol_ts
        
        # =============================================================================
        # Reactor model
        # =============================================================================
        residenceTime            = conditions.simul_param.end_sim #s
        reactorVolume            = 30.5*(1e-2)**3 #m3
        pressureValveCoefficient = .05
        maxSimulationTime        = 50 # seconds
        
        fuelAirMixtureTank       = ct.Reservoir(gas_red)
        exhaust                  = ct.Reservoir(gas_red)
        stirredReactor           = ct.IdealGasReactor(gas_red, energy='off', volume=reactorVolume)
        massFlowController       = ct.MassFlowController(upstream=fuelAirMixtureTank,
                                                   downstream=stirredReactor,
                                                   mdot=stirredReactor.mass/residenceTime)
        pressureRegulator        = ct.Valve(upstream=stirredReactor,
                                     downstream=exhaust,
                                     K=pressureValveCoefficient)
        reactorNetwork           = ct.ReactorNet([stirredReactor])
        # now compile a list of all variables for which we will store data
        columnNames = [stirredReactor.component_name(item) for item in range(stirredReactor.n_vars)]
        columnNames = ['pressure'] + columnNames
        # use the above list to create a DataFrame
        timeHistory = pd.DataFrame(columns=columnNames)

        # =============================================================================
        # Modeling data
        # =============================================================================
        # Create a data frame to store values for the above points
        tempDependence            = pd.DataFrame(columns=timeHistory.columns)
        tempDependence.index.name = 'Temperature'
        inletConcentrations       = X_red
        X_outlet                   = inletConcentrations
        abs_tol = conditions.simul_param.atol_ts
        rel_tol = conditions.simul_param.rtol_ts
        
        T_list  = list(red_results.pts_scatter)

        t=-1 ; t_i=-1
        for T in T_list:
            t+=1
            #Re-initialize the gas
            reactorTemperature = T #Kelvin
            reactorVolume = 30.5*(1e-2)**3 #m3
            gas_red.TPX = reactorTemperature,conditions.state_var.P,inletConcentrations
            # Re-initialize the dataframe used to hold values
            timeHistory = pd.DataFrame(columns=columnNames)
            # Re-initialize all the reactors, reservoirs, etc
            fuelAirMixtureTank = ct.Reservoir(gas_red)
            stirredReactor = ct.IdealGasReactor(gas_red, energy='off', volume=reactorVolume)

            # use concentrations from the previous iteration to speed up convergence
            # gas_red.TPX = reactorTemperature,conditions.state_var.P,concentrations
            gas_red.X = X_outlet
            exhaust = ct.Reservoir(gas_red)
            
            massFlowController = ct.MassFlowController(upstream=fuelAirMixtureTank,
                                                       downstream=stirredReactor,
                                                       mdot=stirredReactor.mass/residenceTime)
            pressureRegulator  = ct.Valve(upstream=stirredReactor,
                                         downstream=exhaust,
                                         K=pressureValveCoefficient)
            
            reactorNetwork     = ct.ReactorNet([stirredReactor])
            reactorNetwork.rtol = conditions.simul_param.rtol_ts
            reactorNetwork.atol = conditions.simul_param.atol_ts
            # if red_data.red_op.rtol_ts:
            #     reactorNetwork.rtol_sensitivity = red_data.red_op.rtol_ts
            #     reactorNetwork.atol_sensitivity = red_data.red_op.atol_ts
            # else:
            #     reactorNetwork.rtol_sensitivity = conditions.simul_param.rtol_ts
            #     reactorNetwork.atol_sensitivity = conditions.simul_param.atol_ts

            # for r in range(gas_red.n_reactions):
            #     stirredReactor.add_sensitivity_reaction(r)

            # Re-run the isothermal simulations
            time_s = 0

            while time_s < maxSimulationTime:
                time_s = reactorNetwork.step()

            # record concentrations to speed up convergence of the next iteration 
            X_outlet = stirredReactor.thermo.X
            
            # conc_ref, omega_ref = np.zeros(n_sp), np.zeros(n_sp)
            # timescales = np.ones(n_sp)
            #on initialise le vecteur timescales avec des 1 (pour pouvoir diviser par Jii ensuite)
            timescales = np.ones(n_sp)
            omega_ref  = np.zeros(n_sp)
            conc_ref = gas_red.concentrations

            for idx_mech1 in range(n_sp):
                try:
                    idx_red = gas_red.species_index(gas_ref.species_name(idx_mech1))
                    omega_ref[idx_mech1]  = gas_red.net_production_rates[idx_red]
                    conc_i                = conc_ref[idx_red]


                    #on initialise le vecteur timescales avec des 1 (pour pouvoir diviser par Jii ensuite)
                                        
                    #si la concentration est trop faible tau=1 par defaut
                    if conc_i > abs_tol:
                        #calcul du pas de perturbation
                        dc = abs_tol + conc_i * rel_tol
                        
                        #etat perturbé
                        conc_pert = conc_ref.copy()
                        conc_pert[idx_red] += dc
                        
                        #mise a jour de l'etat cantera
                        gas_red.concentrations = conc_pert
                        omega_pert_sp = gas_red.net_production_rates[idx_red]

                        #approximation de la dérivée (Jacobien diagonal)
                        J_diag = (omega_pert_sp - omega_ref[idx_mech1]) / dc
                        
                        #calcul implicite du timescale avec borne max a 1.0
                        if abs(J_diag) > 1 : 
                            timescales[idx_mech1] = 1.0 / abs(J_diag)
                        else : 
                            timescales[idx_mech1] = 1.0                        
                    
                except ValueError: #if spec not in red mech
                    pass

                    
                
            timescales_all[t_i] = timescales.copy()
   
    
    timescales_all = np.array(timescales_all)   # taille (n_points_LOI, n_sp)
    
    
    fn = conditions.num+'tp_'+conditions.composition.fuel.replace('/','').split('(')[0]\
            +'_'+'%.2f' %conditions.composition.phi\
            +'_'+'%.0f'%conditions.state_var.T+'_'+'%.0f'%(conditions.state_var.P)
    fn = fn.replace('.','p') + '.png'

    plot_timescales_bars(timescales_all, gas_ref.species_names, list(map(str, range(len(timescales_all)))), fn)
    
    
    # =============================================================================
    # 2- sensitivity analysis
    # =============================================================================
    
    red_data = sa.sensitivities_computation_SA(red_data, mech_data,red_results,LOI_calc=True)
    S = red_data.red_op.S_AB_z
    
       
    
    # =============================================================================
    # 3- LOI computation
    # =============================================================================
    
    tsp_idx = red_data.targetSpeciesIdx
    n_tsp = len(tsp_idx)
    tsp_name = red_data.tspc
    
    LOI = np.zeros((n_points_LOI, n_tsp, n_sp))
    
    min_tsc = np.min(timescales_all[timescales_all > 0]) # min > 0
    max_tsc = np.max(timescales_all)
    
    min_S = np.min(S[S>0]) # min > 0
    max_S = np.max(S)
    
    for tsp in range(n_tsp):
        
       for i in range(n_sp):
           
           t_i = -1
           for t in range(n_points-int(n_points/red_data.red_op.n_points+1)):
               # sim.advance(pts_scatter[t])
               if t != 0 and t % int(max(n_points/red_data.red_op.n_points, 1)) == 0:
                   t_i += 1
                   
                   # ---------- Calculation of LOI : method 1
                   # definition of Lovas (2007)
                   # print_('LOI: method 1',mp)
                   LOI[t_i, tsp, i] = S[ t_i, tsp, i] * timescales_all[t_i, i]
                   
                   
                   # ---------- Calculation of LOI : method 2
                   # calculation with log-average between S(AB) and timescale : 
                   # print_('LOI: method 2',mp)
                   # # log-space normalization of timescale
                   # tsc_log_norm = (np.log(np.max([timescales_all[t_i, i],min_tsc])) - np.log(min_tsc)) /  \
                   #                (np.log(max_tsc) - np.log(min_tsc))
                   # # log-space normalization of S
                   # S_log_norm   = (np.log(np.max([S[ t_i, tsp, i],min_S])) - np.log(min_S)) / \
                   #                (np.log(max_S) - np.log(min_S))
                   # # LOI = log-average
                   # LOI[t_i, tsp, i] = 0.5*(tsc_log_norm + S_log_norm)

                   # ---------- Calculation of LOI : method 3
                   #  S(AB) x log normalized timescale : 
                   # print_('LOI: method 3',mp)
                   # log-space normalization of timescale
                   # tsc_log_norm = (np.log(np.max([timescales_all[t_i, i],min_tsc])) - np.log(min_tsc)) /  \
                   #                (np.log(max_tsc) - np.log(min_tsc))
                   # # LOI = log-average
                   # LOI[t_i, tsp, i] =  S[ t_i, tsp, i] * tsc_log_norm


    # Calcul du maximum de LOI pour chaque (tsp, i)
    LOI_max = np.max(np.abs(LOI), axis=0)   # taille = (n_tsp, n_sp)
    for _t in range(len(LOI_max)):
        loi_max_targ = max(LOI_max[_t])
        for _sp in range(len(LOI_max[_t])):
            LOI_max[_t][_sp] = LOI_max[_t][_sp]/loi_max_targ

    red_data.red_op.LOI_max = LOI_max 
        
    red_data.red_op.LOI = LOI
    
    return red_data
    



def speciesWithdrawal(conditions, red_data, red_method, mech_data, eps):
    # gas, fuel, diluant, sensi_species, eps, target_species, active_species_main, verbose = 3) :

    mp = conditions.main_path

    # main variables
    gas = red_data.gas_ref

    # sensi_AB = red_data.red_op.sensi_sp
    LOI = copy.deepcopy(red_data.red_op.LOI_max)

    tsp_idx = red_data.targetSpeciesIdx
    verbose = conditions.simul_param.verbose

    ns = gas.n_species
#    nr = gas.n_reactions
    active_species = list(mech_data.spec.activ_p)

    # Add fuel / ox /diluent to the conserved species
    init_spec = conditions.composition.X
    init_spec = init_spec.split(',')
    for spec in init_spec:
        ind_spec = gas.species_index(spec.split(":")[0].replace(' ', ''))
        if not active_species[ind_spec]:
            active_species[ind_spec] = True
    if conditions.composition.X2:
        init_spec = conditions.composition.X2
        init_spec = init_spec.split(',')
        for spec in init_spec:
            ind_spec = gas.species_index(spec.split(":")[0].replace(' ', ''))
            if not active_species[ind_spec]:
                active_species[ind_spec] = True


    # if '_sp' in red_method:
    # Species sensitivities based withdrawal
    # for t in range(len(tsp_idx)):
    #     spA = tsp_idx[t]
    #     if not active_species[spA]:
    #         active_species[spA] = True
    #     for spB in range(ns):
    #         # removed species from previous reductions
    #         if not mech_data.spec.activ_m[spB]:
    #             active_species[spB] = False
    #         else:
    #             if LOI[t][spB] > eps[t] and not active_species[spB]:
    #                 active_species[spB] = True
    sp_rank = np.zeros(ns)
    if '_sp' in red_method:
        # Species sensitivities based withdrawal
        for t in range(len(tsp_idx)):
            spA=tsp_idx[t]
            if not active_species[spA]: active_species[spA]=True
            for spB in range(ns):
                if not mech_data.spec.activ_m[spB]:
                    active_species[spB] = False
                else:
                    sp_rank[spB] = max(sp_rank[spB],LOI[t][spB])
                    if LOI[t][spB]>eps[t] and not active_species[spB]:
                        active_species[spB]=True


    red_data.red_op.sp_rank = sp_rank


    return active_species, red_data


def reactionWithdrawal(conditions, mech_data, active_species, red_data, red_method, eps_r):
    # gas, sensi_reactions, active_species , eps, target_species, active_reactions_main, verbose=3) :

    # main variables
    gas_ref = red_data.gas_ref

    ns = gas_ref.n_species
    nr = gas_ref.n_reactions

    tsp_idx = red_data.targetSpeciesIdx

    if int(ct.__version__[0]) > 2:
        nu_f = gas_ref.reactant_stoich_coeffs
        nu_r = gas_ref.product_stoich_coeffs
    else:
        nu_f = gas_ref.reactant_stoich_coeffs()
        nu_r = gas_ref.product_stoich_coeffs()

    active_reactions = list(mech_data.react.activ_p)


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
                        if nu_f[sp][r] != 0 or nu_r[sp][r] != 0:
                            active_reactions[r] = False
                            break


    # check threebody exception (+AR) (+HE) etc.
    for r in range(len(mech_data.react.equation)):
        if type(mech_data.react.tbe[r]) is str:
            for sp in range(len(mech_data.spec.name)):
                if mech_data.react.tbe[r] == mech_data.spec.name[sp]\
                        and not active_species[sp]:
                    active_reactions[r] = False

    return active_reactions, active_species




def plot_timescales_bars(timescales_all, group_labels=None, bar_labels=None,fn='fig.png'):
    import matplotlib.pyplot as plt    
    """
    Args:
        timescales_all (np.ndarray): timescales_all 2D of shape (Species names, Position)
        Species names (list):  X (optionnel)
        Position (list)
    """
    
    # Replace values equal to 1 with NaN to remove them from the graph
    timescales_all = timescales_all.astype(float)  # pour pouvoir mettre des NaN
    timescales_all[timescales_all == 1] = np.nan

    n_series, n_groups = timescales_all.shape
    x = np.arange(n_groups)  # position of groups
    width = 0.8 / n_series   # bar width

    fig, ax = plt.subplots(figsize=(n_groups/1.5, 5))  # largeur proportionnelle au nb de groupes²
    # fig, ax = plt.subplots(figsize=(10, 5))

    for i in range(n_series):
        ax.bar(x + i*width, timescales_all[i], width, label=bar_labels[i] if bar_labels else f"Série {i+1}")

    # labels
    ax.set_xlabel("Species")
    ax.set_ylabel("Timescale (s)")
    if group_labels:
        ax.set_xticks(x + width*(n_series-1)/2)
        ax.set_xticklabels(group_labels)
    else:
        ax.set_xticks(x + width*(n_series-1)/2)
        ax.set_xticklabels([f"Col {j+1}" for j in range(n_groups)])

    if bar_labels:
        ax.legend()
        
    ax.set_yscale("log")

    plt.tight_layout()
    plt.savefig(str(fn))
    
    plt.show()

