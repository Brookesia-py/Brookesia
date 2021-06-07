"""
    Brookesia
    Reduction and optimization of kinetic mechanisms

    Copyright (C) 2020  Chakravarty, Matynia
    contact : hkchakravarty@gmail.com

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


# importing all the required libraries........................................

import warnings
warnings.filterwarnings('ignore')
import cantera as ct
import numpy as np
import scipy.integrate
import sympy as sym
import matplotlib.pyplot as plt
import copy
import brookesia.Class_def as cdef

from sympy.physics.vector import dynamicsymbols


import multiprocessing
import os
#import sys

from  brookesia.Class_def import print_

# import definitions
# get_ipython().run_line_magic('matplotlib', 'inline')

#start = time.time()

#config = configparser.ConfigParser()



def CSP_analysis(red_data,mech_data,results):

# =============================================================================
#%%    1- Evaluation of global reaction rate vector
# =============================================================================
    clock = cdef.Clock('1- Evaluation of global reaction rate vector.') ; clock.start()


#    filename = 'grireduced.cti'
#
#    # custom.py (taken from cantera site)................................................................................
#    gas = ct.Solution(filename)

    conditions = results.conditions
    mp = conditions.main_path
    print_('\n' + conditions.config + " CSP analysis", mp)

    os.chdir(conditions.main_path + '/__CSP')
    mech_data.write_new_mech('csp_to.cti', False, False)

#    gas_ref = red_data.gas_ref
#    gas = red_data.gas_red

    pts_scatter = results.pts_scatter
    verbose     = conditions.simul_param.verbose
    show_plots  = conditions.simul_param.show_plots

    n_points = len(pts_scatter)



    # -------------------------

    # saving and suppression of unpickable variables for parallelization
    rd_gas     = red_data.gas_ref
    rd_gas_ref = red_data.gas_loop
    rd_gas_red = red_data.gas_red
    rd_gas_act = red_data.gas_act
    del red_data.gas_ref ; del red_data.gas_loop ; del red_data.gas_red ; del red_data.gas_act

    par_arg_comm = []
    par_arg_comm.append(red_data)
    par_arg_comm.append(mech_data)
    par_arg_comm.append(n_points)

    num_cores = multiprocessing.cpu_count()



    # Initial condition
    X = conditions.composition.X
    T = conditions.state_var.T
    P = conditions.state_var.P


    # CSP options
#    csp_method = red_data.red_op.csp_method

    # an option that might be added in future
    tr_plot_interaction = False
#    if not tr_plot_interaction:
#        TimeResolution = red_data.red_op.TimeResolution

    epsilon_rel = red_data.red_op.epsilon_rel
#    epsilon_abs_op = red_data.red_op.epsilon_abs


# =============================================================================
#     Test Bench
#     Lam SH, Goussis DA. The CSP method for simplifying kinetics. International Journal of Chemical Kinetics 1994;26:461–486.
# =============================================================================

    # print_('test bench',mp)

#    list_reaction = ['A + A = B + H1', 'A = B + H2', 'B + B = A + H3']
#
#    H1_0 = 1.1e+4  ;   H2_0 = 1.0e+5  ;   H3_0 = -2.9e+5
#    A_0  = 1.5e-4  ;   B_0  = 0.1e-6  ;   C_0  = 300
#    A_0_ = 1.46e-4 ;   B_0_ = 1.95e-6 ;   C_0_ = 300.02
#    k1_0 = 1.0e+4  ;   k2_0 = 0.1     ;   k3_0 = 1.0e+4
#    K1_0 = 1.1e-2  ;   K2_0 = 1.1e+2  ;   K3_0 = 0.8e-8
#    initial_concentration_1 = [A_0, B_0,C_0]
#    initial_concentration_2 = [A_0_, B_0_,C_0_]
#    constants = [H1_0, H2_0, H3_0, K1_0, K2_0, K3_0, k1_0, k2_0, k3_0]
#    s1_ = [-2,1,const_list[0]]
#    s2_ = [-1,1,const_list[1]]
#    s3_ = [1,-2,const_list[2]]
#    s1 = sym.Matrix(len(species_list),1,s1_)
#    s2 = sym.Matrix(len(species_list),1,s2_)
#    s3 = sym.Matrix(len(species_list),1,s3_)




    # =============================================================================
    #     # To implement :
    #     # nTime = red_data.red_op.n_points
    # =============================================================================







#    # User must enter the name of input file here (such as param_general_1.ini)
#    config.read('param_general_1.ini') # enter path of the parameter file
#
#    filename = config.get('parameters','filename')
#    t_end = config.getfloat('parameters','end_time')
#    T = config.getfloat('parameters','T')
#    P = config.getfloat('parameters','P')*ct.one_atm
#    X = config.get('parameters','X')
#
#    equiv_ratio = config.getfloat('parameters','equiv_ratio')
#    fuel = config.get('parameters','fuel')
#    oxidizer = config.get('parameters','oxidizer')
#    inert_gas = config.get('parameters','inert_gas')
#    tsp_name = config.get('parameters','tsp_name').split()


    # Combustion and Flame 146 (2006) 29–51


    # custom.py from cantera has been used only to get concentration of species as a function of time
    # and solving a ignitions problem at constant pressure where the governing equations are implemented in Python


#    gas = ct.Solution(filename)
    gas = results.gas
    gas.TPX = T, P, X
    # other conditions
    y0 = np.hstack((gas.T, gas.Y))


#    gas
#    gas.set_equivalence_ratio(equiv_ratio, fuel, oxidizer+','+inert_gas)

#    inlet = ct.Reservoir(gas)

    class ReactorOde:
        def __init__(self, gas):
            # Parameters of the ODE system and auxiliary data are stored in the
            # ReactorOde object.
            self.gas = gas
            self.P = gas.P

        def __call__(self, t, y):
            """the ODE function, y' = f(t,y) """

            # State vector is [T, Y_1, Y_2, ... Y_K]
            self.gas.set_unnormalized_mass_fractions(y[1:])
            self.gas.TP = y[0], self.P
            rho = self.gas.density

            wdot = self.gas.net_production_rates
            dTdt = - (np.dot(self.gas.partial_molar_enthalpies, wdot) /
                      (rho * self.gas.cp))
            dYdt = wdot * self.gas.molecular_weights / rho

            return np.hstack((dTdt, dYdt))


    # Set up objects representing the ODE and the solver
    ode = ReactorOde(gas)
    solver = scipy.integrate.ode(ode)
    solver.set_integrator('vode', method='bdf', with_jacobian=True)
    solver.set_initial_value(y0, 0.0)

    # Integrate the equations, keeping T(t) and Y(k,t)
    # t_end = 1e-3


#    t_end = conditions.simul_param.pts_scatter[-1]

#    states = ct.SolutionArray(gas, 1, extra={'t': [0.0]})
#    dt =  t_end/50# 1e-5
#    while solver.successful() and solver.t < t_end:
#        solver.integrate(solver.t + dt)
#        gas.TPY = solver.y[0], P, solver.y[1:]
#        states.append(gas.state, t=solver.t)

#    states = ct.SolutionArray(gas, 1, extra={'t': [0.0]})
    k_f, k_r = [], []
#    dt =  t_end/50# 1e-5
#    while solver.successful() and solver.t < t_end:
    for t in range(n_points-int(n_points/red_data.red_op.n_points+1)):
        if t!=0 and t%int(max(n_points/red_data.red_op.n_points,1))==0 :
            solver.integrate(pts_scatter[t])
            gas.TPY = solver.y[0], P, solver.y[1:]
            if len(k_f) == 0:
                states = ct.SolutionArray(gas, 1, extra={'t': [0.0]})
            else:
                states.append(gas.state, t=pts_scatter[t])
            k_f.append(gas.forward_rate_constants)
            k_r.append(gas.reverse_rate_constants)

    conc = states.Y
    T    = states.T
    P    = states.P


#    if test_bench = true:
#        conc =
#        T =
#        P =
#        k =
#        K =

    #.........................end of custom.py..............................................................................................

    # Evaluation of the global reaction rate vector using the stoichiometric vector and reaction rates.

    species_names = gas.species_names
    for s_i in range(len(species_names)):
        species_names[s_i] = species_names[s_i].replace(',','\,')
        species_names[s_i] = species_names[s_i].replace('(','ob')
        species_names[s_i] = species_names[s_i].replace(')','cb')
    nSpecies, nReaction = len(species_names), len(gas.reactions())
    if verbose>5: print('No. of species: {}'.format(nSpecies))
    if verbose>5: print('Total no. of chemical reactions: {}'.format(nReaction))

    species_list = []
    for species in species_names:
        species_list.append(dynamicsymbols(species))
#    dicts_conc_initial = dict(zip(species_list,states.Y[0,:]))

    F = [] ; g = [] ; info4kin = [] ; _t = 0
    for t in range(n_points-int(n_points/red_data.red_op.n_points+1)):
        if t!=0 and t%int(max(n_points/red_data.red_op.n_points,1))==0 :
            info4kin.append([conc[_t],T[_t],P[_t],k_f[_t],k_r[_t]])
            F_t, Sr = F_S_vector(red_data,mech_data,gas,info4kin[-1],verbose)
            F.append(F_t)
            g.append(Sr*F_t)
            _t += 1

#    g = []
#    for _t in range(len(pts_scatter)):
#        g=Sr*F
#    if verbose >7: print('physical global reaction rate vector: {}'.format(g))

#    dicts_g = {}
#    for i in range(len(species_list)):
#        dicts_g[i] = g[i].subs(dicts_conc_initial)
#        dicts_g[species_list[i].diff()]=dicts_g.pop(i)

    clock.stop() ; clock.display_(mp)


# =============================================================================
#%%    2- Evaluation of jacobian matrix by differentiating g with respect to y
# =============================================================================
    clock = cdef.Clock('2- Evaluation of jacobian matrix.') ; clock.start()

#    nTime = len(states.Y)


    # Jacobian
    J_fun = [] ; _t = 0
    for t in range(n_points-int(n_points/red_data.red_op.n_points+1)):
        if t!=0 and t%int(max(n_points/red_data.red_op.n_points,1))==0 :
            J_fun.append(sym.Matrix(nSpecies,nSpecies,[sym.diff(g[_t][i],species_list[j]) for i in range(nSpecies) for j in range(nSpecies)]))
            _t += 1
#            if verbose>8: print ('J_fun: \n{}\n'.format(J_fun))

    J_val = []
#    for iTime in range(nTime):
    _t = 0
    for t in range(n_points-int(n_points/red_data.red_op.n_points+1)):
        if t!=0 and t%int(max(n_points/red_data.red_op.n_points,1))==0 :
            dicts_conc = dict(zip(species_list,conc[_t,:]))
            J_val.append(np.array(J_fun[_t].subs(dicts_conc),dtype=float))
            _t += 1

    if verbose >7:
        print_('Jacobian: \n{}\n'.format(J_val))

    Tau = []# list containing the timeSceles at all times
    _t = 0
    for t in range(n_points-int(n_points/red_data.red_op.n_points+1)):
        if t!=0 and t%int(max(n_points/red_data.red_op.n_points,1))==0 :
            timeScales = np.reciprocal(np.sort(np.absolute(np.linalg.eigvals(J_val[_t])))[::-1])
            Tau.append(timeScales)
            _t += 1

    clock.stop() ; clock.display_(mp)


# =============================================================================
#%%    3a- Représentation alternative du problème
# =============================================================================
    a_i=[] ; _t = 0
    for t in range(n_points-int(n_points/red_data.red_op.n_points+1)):
        if t!=0 and t%int(max(n_points/red_data.red_op.n_points,1))==0 :
        # Computing the eigen values(w) and eigen vectors (v) of the Jacobian
            w, v = np.linalg.eig(J_val[_t])
            a_i.append(v)
            _t += 1

# =============================================================================
#%%    3b- Start of Decompostion of phase space
# =============================================================================
    clock = cdef.Clock('3- Decomposition of phase space.') ; clock.start()

    if verbose > 3 and show_plots:
        _t = 0
        for t in range(n_points-int(n_points/red_data.red_op.n_points+1)):
            if t!=0 and t%int(max(n_points/red_data.red_op.n_points,1))==0 :
                plt.semilogy(range(nSpecies),Tau[_t],'.b')
                _t += 1

    # add interaction procedure  -> tr_plot_interaction
        plt.grid()
        plt.ylabel('time (s)')
        plt.xlabel('# modes ')
        plt.title('Decompostion of phase space')
        plt.show()

    if tr_plot_interaction: # option that might be added in the future
        red_data.red_op.TimeResolution = float(input('Enter a time resolution; reactions with time scales less than this threshold would be considered as initial guess for the total number of fast modes \n(e.g. 1e0)'))

    arg_psd = []
    arg_psd.append(Tau)
    arg_psd.append(species_list)
    arg_psd.append(conc)
    arg_psd.append(nSpecies)
    arg_psd.append(n_points)
    arg_psd.append(show_plots)
    arg_psd.append(verbose)

    red_data.red_op.nFastMode = phase_space_decomposition(red_data,arg_psd,True)



    if verbose > 3 and show_plots:
        plt.plot(red_data.red_op.nFastMode,'s')
        plt.grid()

    red_data.red_op.nFastMode = red_data.red_op.nFastMode

    clock.stop() ; clock.display_(mp)



# =============================================================================
#%%    4- performing the refinement of intial trial basis vector (two-step recursive refinement procedure)
# =============================================================================

    clock = cdef.Clock('4- Refinement of initial trial basis vector.') ; clock.start()

    a_, b_, t_a, a_refined, b_refined = [], [], [], [], []
    _t = 0
    for t in range(n_points-int(n_points/red_data.red_op.n_points+1)):
        if t!=0 and t%int(max(n_points/red_data.red_op.n_points,1))==0 :
        # Computing the eigen values(w) and eigen vectors (v) of the Jacobian
#            w, v = np.linalg.eig(J_val[_t])
            a_.append(a_i[_t])
            b_.append(np.linalg.inv(a_i[_t]))
            t_a.append(pts_scatter[t])
            a_refined.append(sym.re(sym.Matrix(a_[-1])))
            b_refined.append(sym.re(sym.Matrix(np.linalg.inv(a_[-1]))))
            _t += 1

    _ref_i=0 #; red_data.red_op.csp_refin_iter=2
    while _ref_i < red_data.red_op.csp_refin_iter:
        dadt, dbdt = [], []
        _t = 0
        for t in range(n_points-int(n_points/red_data.red_op.n_points+1)):
            if t!=0 and t%int(max(n_points/red_data.red_op.n_points,1))==0:
                if _t < len(a_)-1:
                    dadt.append((a_[_t+1]-a_[_t])/(t_a[_t+1]-t_a[_t]))
                    dbdt.append((b_[_t+1]-b_[_t])/(t_a[_t+1]-t_a[_t]))
                else:
                    dadt.append((a_[_t]-a_[_t-1])/(t_a[_t]-t_a[_t-1]))
                    dbdt.append((b_[_t]-b_[_t-1])/(t_a[_t]-t_a[_t-1]))
                _t += 1

        tau_refined = [], [], []
        _t = 0

        Rf_V = []
        Rf_V = []
        _t = 0
        for t in range(n_points-int(n_points/red_data.red_op.n_points+1)):
            if t!=0 and t%int(max(n_points/red_data.red_op.n_points,1))==0 :
                Rf_V.append([a_[_t], b_[_t], red_data.red_op.nFastMode[_t],\
                             info4kin[_t], dadt[_t], dbdt[_t], verbose, par_arg_comm, _t])
                _t += 1
        _ref_i +=1


#    a_initial    = IIV_arg[0]
#    b_initial    = IIV_arg[1]
#    nFastMode    = IIV_arg[2]
#    conc         = IIV_arg[3]
#    dadt         = IIV_arg[4]
#    dbdt         = IIV_arg[5]
#    J_fun        = IIV_arg[6]
#    verbose      = IIV_arg[7]

        a_b_tau_ref=[]
        if os.name == 'nt': multiprocessing.get_context('spawn')
        with multiprocessing.Pool(num_cores) as p:
            a_b_tau_ref=p.map(refinement, Rf_V)

        a_b_tau_ref.sort(reverse=False, key=lambda col: col[3])

        a_refined, b_refined, tau_refined = [], [], []
        for abt in a_b_tau_ref:
                a_refined.append(sym.re(sym.Matrix(abt[0])))
                b_refined.append(sym.re(sym.Matrix(abt[1])))
                tau_refined.append(abt[2])

        print('Refinement iteration ' + str(_ref_i) + ' done')


#for iTime in range(nTime):
#    # Computing the eigen values(w) and eigen vectors (v) of the Jacobian
#    w, v = np.linalg.eig(J_val[iTime])
#    a_refined.append(sym.re(sym.Matrix(v)))
#    b_refined.append(sym.re(sym.Matrix(np.linalg.inv(v))))


    clock.stop() ; clock.display_(mp)



# =============================================================================
#%%    4bis- Evaluation of number of exhausted modes
# =============================================================================



    if epsilon_rel:

        clock = cdef.Clock('4bis- re-evaluation of number of exhausted modes with error analysis') ; clock.start()

        arg_psd.append(a_refined)
        arg_psd.append(b_refined)
        arg_psd.append(g)

        red_data.red_op.nFastMode = phase_space_decomposition(red_data,arg_psd,False)


        clock.stop() ; clock.display_(mp)



# =============================================================================
#%%    5- Evaluation of importance index from M. valorani, Combustion and Flame 146 (2006) 29–51
# =============================================================================
    clock = cdef.Clock('5- Determination of projection matrix') ; clock.start()



    IIV_arg = []
    _t = 0
    for t in range(n_points-int(n_points/red_data.red_op.n_points+1)):
        if t!=0 and t%int(max(n_points/red_data.red_op.n_points,1))==0 :
            IIV_arg.append([a_refined[_t], b_refined[_t], red_data.red_op.nFastMode[_t], info4kin[_t],par_arg_comm,_t])
            _t += 1



    # ---------------------------

    proj_mat=[]
    if os.name == 'nt': multiprocessing.get_context('spawn')
    with multiprocessing.Pool(num_cores) as p:
        proj_mat=p.map(projectionMatrix, IIV_arg)

    proj_mat.sort(reverse=False, key=lambda col: col[4])

    radicalSpecies_list = []
    for pm_i in range(len(proj_mat)):
        radicalSpecies_list.append(proj_mat[pm_i][1])
        red_data.red_op.radIdx.append(proj_mat[pm_i][3])

        #Q, radicalSpecies_list, nonRadicalSpecies_list, radIdx



#    I_fast, I_slow = [], []
#    for iTime in range(nTime):
#        I_f, I_s = importanceIndex_Valorani(a_refined[iTime], b_refined[iTime], nFastModeGuess[iTime], conc[iTime,:])
#        red_data.red_op.I_fast.append(I_f)
#        red_data.red_op.I_slow.append(I_s)
#    del I_f, I_s

#    print('Finished in {} seconds.'.format(end-start))

    clock.stop() ; clock.display_(mp)


# =============================================================================
#%%    6- Determination of projection matrix
# =============================================================================
    clock = cdef.Clock('6- Evaluation of importance index') ; clock.start()


    # -------------------------

    IIV=[]
    if os.name == 'nt': multiprocessing.get_context('spawn')
    with multiprocessing.Pool(num_cores) as p:
        IIV=p.map(importanceIndex_Valorani, IIV_arg)

    IIV.sort(reverse=False, key=lambda col: col[2])

    for IIV_i in range(len(IIV)):
        red_data.red_op.I_fast.append(IIV[IIV_i][0])
        red_data.red_op.I_slow.append(IIV[IIV_i][1])
    del IIV



    del proj_mat


    # ---------------------------

#    Q, radicalSpecies_list, nonRadicalSpecies_list = [], [], []
#    for iTime in range(nTime):
#        Q_, rad, nonRad, radIdx_ = projectionMatrix(a_refined[iTime], b_refined[iTime], nFastModeGuess[iTime], conc[iTime,:])
#        Q.append(Q_)
#        radicalSpecies_list.append(rad)
#        nonRadicalSpecies_list.append(nonRad)
#        red_data.red_op.radIdx.append(radIdx_)
#    if verbose>7:
#        print('Radical species: {}\nNon radical species: {}\n'.format(radicalSpecies_list[iTime], nonRadicalSpecies_list[iTime]))


    CSP_radicals = copy.deepcopy([True]*len(mech_data.spec.name))
    for _t in range(len(red_data.red_op.radIdx)):
        sp_red = -1
        for sp in range(len(mech_data.spec.activ_m)):
            if mech_data.spec.activ_m[sp]:
                sp_red += 1
                if sp_red not in red_data.red_op.radIdx[_t]:
                    CSP_radicals[sp] = False
            else:
                CSP_radicals[sp] = False

    radicals_name_list = []
    for sp in range(len(CSP_radicals)):
        if CSP_radicals[sp]:
            if mech_data.spec.name[sp] not in radicals_name_list:
                radicals_name_list.append(mech_data.spec.name[sp])

#                for sp in range(len(mech_data.spec.name)):
#                    if mech_data.spec.name[sp] not in radicalSpecies_list[t]:
#                        CSP_radicals[sp] = False
#    print(red_data.red_op.radIdx)
#    print(str(CSP_radicals))
    if verbose>2:
        print_('\nCSP radicals: \n', mp)
        print_(str(radicals_name_list), mp)

    clock.stop() ; clock.display_(mp)

    # recovering unpickable variables
    red_data.gas_ref    = rd_gas
    red_data.gas_loop   = rd_gas_ref
    red_data.gas_red    = rd_gas_red
    red_data.gas_act    = rd_gas_act
    del red_data.gas_ref ; del red_data.gas_loop ; del red_data.gas_red


    return red_data, CSP_radicals




def F_S_vector(red_data,mech_data,gas,info4kin,verbose=0):


#    gas         = results.gas
    nReaction   = gas.n_reactions
    nSpecies    = gas.n_species

    sym.var('t');


#    list_ = []
#    for i in range(nSpecies):
#        list_.append([0]*nReaction)
#    S_reactant = sym.Matrix(list_)
#    S_product = sym.Matrix(list_)
#    S_product.shape

    list_ = []
    for i in range(nSpecies):
        list_.append([0]*nReaction)
    S_reactant = sym.Matrix(list_)
    S_product = sym.Matrix(list_)

    for i in range(nSpecies):
        for j in range(nReaction):
            S_reactant[i,j] = gas.reactant_stoich_coeffs()[i,j]
            S_product[i,j] = gas.product_stoich_coeffs()[i,j]
    Sr = S_product - S_reactant
    red_data.red_op.csp_Sr = Sr
    if verbose>7: print('Stoichoiometric vectors: {}'.format(Sr.shape))


#    conc = info4kin[0]
    T    = info4kin[1]
    P    = info4kin[2]
    k_f  = info4kin[3]
    k_r  = info4kin[4]


    if verbose>7: print('rateconstants: {}'.format(k_f.shape))


    F = sym.Matrix([0]*nReaction)
    F_forward = sym.Matrix([0]*nReaction)
    F_backward = sym.Matrix([0]*nReaction)
    for iReaction in range(nReaction):
        reactant_species_list = gas.reactants(iReaction).split()
        #remove ( / ) / +
        for rs_i in range(len(reactant_species_list)):
            if '+' in reactant_species_list[rs_i]:
                reactant_species_list[rs_i] = reactant_species_list[rs_i].replace('+','')
                reactant_species_list[rs_i] = reactant_species_list[rs_i].replace('(','')
                reactant_species_list[rs_i] = reactant_species_list[rs_i].replace(')','')
            else:
                reactant_species_list[rs_i] = reactant_species_list[rs_i].replace('(','ob')
                reactant_species_list[rs_i] = reactant_species_list[rs_i].replace(')','cb')
                reactant_species_list[rs_i] = reactant_species_list[rs_i].replace(',','\,')
        reactant_species_list[:] = (value for value in reactant_species_list if value != '')

        for species in reactant_species_list:
            if species.isdigit():
                idx = reactant_species_list.index(species)
                reactant_species_list[idx+1] = species + reactant_species_list[idx+1]
                reactant_species_list.remove(species)

        product_species_list = gas.products(iReaction).split()
        for ps_i in range(len(product_species_list)):
            if '+' in product_species_list[ps_i]:
                product_species_list[ps_i] = product_species_list[ps_i].replace('+','')
                product_species_list[ps_i] = product_species_list[ps_i].replace('(','')
                product_species_list[ps_i] = product_species_list[ps_i].replace(')','')
            else:
                product_species_list[ps_i] = product_species_list[ps_i].replace(',','\,')
                product_species_list[ps_i] = product_species_list[ps_i].replace('(','ob')
                product_species_list[ps_i] = product_species_list[ps_i].replace(')','cb')

        product_species_list[:] = (value for value in product_species_list if value != '')

        for species in product_species_list:
            if species.isdigit():
                idx = product_species_list.index(species)
                product_species_list[idx+1] = species + product_species_list[idx+1]
                product_species_list.remove(species)

        factor = 1
        for species in reactant_species_list:
            if species[0].isdigit():
                factor*= int(species[0])
                idx = reactant_species_list.index(species)
                reactant_species_list[idx] = reactant_species_list[idx].replace(species[0],'',1)
        product_of_reactants = factor*dynamicsymbols(reactant_species_list[0])
        for i in range(1,len(reactant_species_list)):
            product_of_reactants*=dynamicsymbols(reactant_species_list[i])

        factor = 1
        for species in product_species_list:
            if species[0].isdigit():
                factor*= int(species[0])
                idx = product_species_list.index(species)
                product_species_list[idx] = product_species_list[idx].replace(species[0],'',1)

        product_of_products = factor*dynamicsymbols(product_species_list[0])
        for i in range(1,len(product_species_list)):
            product_of_products*=dynamicsymbols(product_species_list[i])
        F_forward[iReaction] = k_f[iReaction]*product_of_reactants
        F_backward[iReaction] = k_r[iReaction]*product_of_products
        F[iReaction]=(F_forward[iReaction]-F_backward[iReaction])



    # 1- Only for pressure-dependant reaction with collision partner
        #    => Calculation of M           -> omega = [M]*(...)
    third_body_species_name = ['M', '(+M)']
    third_body_species_list = []
    for tbsn in third_body_species_name:
        third_body_species_list.append(dynamicsymbols(tbsn))
#    third_body_species_list = dynamicsymbols(third_body_species_name)

#    mech_data = cdef.Mech_data(filename)
    third_body_list = np.ones(len(mech_data.react.type)) #value of M.....
    ns = gas.n_species

    for r in range(len(mech_data.react.type)):  # loop on all reactions


        if mech_data.react.type[r] == 'three_body_reaction'    or mech_data.react.type[r] == 'falloff_reaction'    or mech_data.react.type[r] == 'pdep_arrhenius':

            # if there is no efficiency specified for the reaction
            if mech_data.react.eff[r]=='':
                third_body_list[r] = P/(8.314*T)

            # if there is efficiency specified for the reaction
            else:
                third_body_list[r]=0
                        # identification of collisional species
                col_sp = [] ; col_coeff = []
                for col_ in mech_data.react.eff[r]:
                    for sp in range(ns):
                        if col_.split(':')[0] == gas.species_name(sp):
                            col_sp.append(sp)
                            col_coeff.append(float(col_.split(':')[1]))
                        # calculation of M as M = Sum(col_coeff*C_sp)
                for sp in range(ns):
                    if sp in col_sp:
                        third_body_list[r] += col_coeff[col_sp.index(sp)]*gas.concentrations[sp]
                    else:
                        third_body_list[r] += 1*gas.concentrations[sp]

    M = third_body_list
    for i in range(nReaction):
        F[i] = F[i].subs({third_body_species_list[0]:M[i], third_body_species_list[1]:M[i]})
        F_forward[i] = F_forward[i].subs({third_body_species_list[0]:M[i], third_body_species_list[1]:M[i]})
        F_backward[i] = F_backward[i].subs({third_body_species_list[0]:M[i], third_body_species_list[1]:M[i]})


    return F, Sr


def importanceIndex_Valorani(IIV_arg):

    a_refined    = IIV_arg[0]
    b_refined    = IIV_arg[1]
    nFastMode    = IIV_arg[2]
    info4kin     = IIV_arg[3]
    IIV_arg_comm = IIV_arg[4]
    red_data     = IIV_arg_comm[0]
    mech_data    = IIV_arg_comm[1]
#    n_points     = IIV_arg_comm[2]
    _t           = IIV_arg[5]

    conc         = info4kin[0]

    ct.suppress_thermo_warnings()
    gas             = ct.Solution('csp_to.cti')
    nReaction       = gas.n_reactions
    nSpecies        = gas.n_species
    species_names   = gas.species_names
    for s_i in range(len(species_names)):
        species_names[s_i] = species_names[s_i].replace(',','\,')
        species_names[s_i] = species_names[s_i].replace('(','ob')
        species_names[s_i] = species_names[s_i].replace(')','cb')

    species_list = []
    for species in species_names:
        species_list.append(dynamicsymbols(species))

#    F = [] #; g = []
#    for t in range(n_points-int(n_points/red_data.red_op.n_points+1)):
#        if t!=0 and t%int(max(n_points/red_data.red_op.n_points,1))==0 :
    F_t, Sr = F_S_vector(red_data,mech_data,gas,info4kin)
#            F.append(F_t)
#            g.append(Sr*F_t)

    dict_ = dict(zip(species_list,conc))

    R = F_t.subs(dict_)

    B = np.zeros((nSpecies,nReaction),dtype=float)
    for i in range(nSpecies):
        for j in range(nReaction):
            B[i,j] = np.array((b_refined[i,:]*Sr[:,j]).subs(dict_),dtype=float)[0,0]
    I_fast = np.zeros((nSpecies,nReaction),dtype=float)
    I_slow = np.zeros((nSpecies,nReaction),dtype=float)
    den = np.zeros((nSpecies,nReaction),dtype=float)

    for i in range(nSpecies):
        for j in range(nSpecies):
            for r in range(nFastMode):
                den[i,j] += float(a_refined[i,r])*B[r,j]*R[j]
    den_ = np.sum(np.absolute(den),axis=1)


    for i in range(nSpecies):
        for k in range(nReaction):
            num  = 0
            for r in range(nFastMode):
#                print('a_refined[i,r]: ' + str(a_refined[i,r]))
#                print('B[r,k]:         ' + str(B[r,k]))
#                print('R[k]:           ' + str(R[k]))
                num += float(a_refined[i,r])*B[r,k]*R[k]
#            print(str(num))
#            print(str(den_[i]))
            I_fast[i,k] = num/den_[i]


    den = np.zeros((nSpecies,nSpecies),dtype=float)

    for i in range(nSpecies):
        for j in range(nSpecies):
            for r in range(nFastMode,nSpecies):
                den[i,j] += float(a_refined[i,r])*B[r,j]*R[j]
    den_ = np.sum(np.absolute(den),axis=1)

    for i in range(nSpecies):
        for k in range(nReaction):
            num  = 0
            for r in range(nFastMode,nSpecies):
                num += float(a_refined[i,r])*B[r,k]*R[k]
            I_slow[i,k] = num/np.real(den_[i])



    return I_fast, I_slow, _t



def refinement(IIV_arg):

# Refinement functions: The basis vectors generated from the Jacobian matrix provide
# leading order approximation and has been used as intial trial basis vectors in the refinement procedure.
# Since basis vectors changes as a function of time, it is important to perform two step
# CSP refinement procedure to generate the correct refined set of vectors valid at all the time (time-dependent basis vector).
# In general, CSP method provides programmable two-step recursive procedure for refinement of the linearly  independent basis vectors
# obtained from the Jacobian matrix of chemical kinetic mechanism to generate the correct refined
# set of orthogonal time-dependent basis vector to develop a reduce model with higher accuracy
# Combustion Science and Technology, 89, 5-6, p. 375, 1993.

# Defining the CSP refinement function


    a_initial    = IIV_arg[0]
    b_initial    = IIV_arg[1]
    nFastMode    = IIV_arg[2]
    info4kin     = IIV_arg[3]
    dadt         = IIV_arg[4]
    dbdt         = IIV_arg[5]
    verbose      = IIV_arg[6]
    par_arg_comm = IIV_arg[7]
    red_data     = par_arg_comm[0]
    mech_data    = par_arg_comm[1]
    _t           = IIV_arg[8]
#
    ct.suppress_thermo_warnings()
    gas             = ct.Solution('csp_to.cti')
#    nReaction       = gas.n_reactions
    nSpecies        = gas.n_species
    species_names   = gas.species_names
    for s_i in range(len(species_names)):
        species_names[s_i] = species_names[s_i].replace(',','\,')
        species_names[s_i] = species_names[s_i].replace('(','ob')
        species_names[s_i] = species_names[s_i].replace(')','cb')

    species_list = []
    for species in species_names:
        species_list.append(dynamicsymbols(species))
    conc            = info4kin[0]

#    F = [] #; g = []
    F_t, Sr = F_S_vector(red_data,mech_data,gas,info4kin)
#    F.append(F_t)
    g_t = Sr*F_t
#        J_fun = [] ; _t = 0
#    for t in range(n_points-int(n_points/red_data.red_op.n_points+1)):
#        if t!=0 and t%int(max(n_points/red_data.red_op.n_points,1))==0 :
    J_fun = sym.Matrix(nSpecies,nSpecies,[sym.diff(g_t[i],species_list[j]) for i in range(nSpecies) for j in range(nSpecies)])
#            _t += 1
#    print(nFastMode)

    dicts_conc = dict(zip(species_list,conc))
    J_val = np.array(J_fun.subs(dicts_conc),dtype=float)

    Lambda_initial = np.zeros((nFastMode,nFastMode),dtype=float)  # equation 6.15, without the derivative term.....
    for i in range(nFastMode):
        for j in range(nFastMode):
            Lambda_initial[i,j] = np.matmul(np.matmul(b_initial[i,:],J_val),a_initial[:,j])
#    print(Lambda_initial)
    if verbose>8: print('Initial Lambda matrix: {}\n'.format(Lambda_initial))
    tau_initial = np.linalg.inv(Lambda_initial)
    if verbose>8: print('Initial Tau matrix: {}\n'.format(tau_initial))

    b_refined_1 = b_initial
    for i in range(nFastMode):
        sum_ = np.zeros((1,nSpecies),dtype=float)
        for j in range(nFastMode):
            sum_ += np.real(tau_initial[i,j]*np.matmul(np.reshape(b_initial[j,:],(1,nSpecies)),J_val))
            b_refined_1[i,:] = sum_
    if verbose>8: print('Rows of refined vector b (1), corresponding to fast modes (only): \n{}\n'.format(b_refined_1))

    Lambda_1 = np.zeros((nFastMode,nFastMode),dtype=float)  # equation 6.15, without the derivative term.....
    for i in range(nFastMode):
        for j in range(nFastMode):
            Lambda_1[i,j] = np.matmul(np.matmul(b_refined_1[i,:],J_val),a_initial[:,j])
    if verbose>8: print('Initial Lambda matrix: {}\n'.format(Lambda_1))
    tau_1 = np.linalg.pinv(Lambda_1)
    if verbose>8: print('Initial Tau matrix: {}\n'.format(tau_1))

    a_refined_1 = a_initial
    for i in range(nFastMode):
        sum_ = np.zeros((nSpecies,1),dtype=float)
        for j in range(nFastMode):
            sum_ += np.real(tau_1[j,i]*np.matmul(J_val,np.reshape(a_initial[:,j],(nSpecies,1))))
            a_refined_1[:,i] = np.reshape(sum_,(nSpecies,))
    if verbose>8: print('Rows of refined vector a (1), corresponding to fast modes (only): \n{}\n'.format(a_refined_1))

# step2: including the derivative term
    b_refined_2 = b_refined_1
    for i in range(nFastMode):
        sum_ = np.zeros((1,nSpecies),dtype=float)
        for j in range(nFastMode):
            sum_ += np.real(tau_1[i,j]*(np.reshape(dbdt[:,j],(1,nSpecies))+np.matmul(np.reshape(b_refined_1[j,:],(1,nSpecies)),J_val)))
            b_refined_2[i,:] = sum_
    if verbose>8: print('Rows of refined vector b (1), corresponding to fast modes (only): \n{}\n'.format(b_refined_1))

    Lambda_2 = np.zeros((nFastMode,nFastMode),dtype=float)  # equation 6.15, without the derivative term.....
    for i in range(nFastMode):
        for j in range(nFastMode):
            Lambda_2[i,j] = np.matmul(np.matmul(b_refined_2[i,:],J_val),a_refined_1[:,j])
    if verbose>8: print('Initial Lambda matrix: {}\n'.format(Lambda_1))
    tau_2 = np.linalg.pinv(Lambda_2)
    if verbose>8: print('Initial Tau matrix: {}\n'.format(tau_1))

    a_refined_2 = a_refined_1
    for i in range(nFastMode):
        sum_ = np.zeros((nSpecies,1),dtype=float)
        for j in range(nFastMode):
            sum_ += np.real(tau_1[j,i]*(-np.reshape(dadt[:,j],(nSpecies,1))+np.matmul(J_val,np.reshape(a_refined_1[:,j],(nSpecies,1)))))
            a_refined_2[:,i] = np.reshape(sum_,(nSpecies,))


    if verbose>8: print('Rows of refined vector a (1), corresponding to fast modes (only): \n{}\n'.format(a_refined_1))


#    return sym.re(sym.Matrix(a_refined_2)), sym.re(sym.Matrix(b_refined_2)), tau_2, _t
    return a_refined_2, b_refined_2, tau_2, _t




def projectionMatrix(IIV_arg):


    a_refined    = IIV_arg[0]
    b_refined    = IIV_arg[1]
    nFastMode    = IIV_arg[2]
#    conc         = IIV_arg[3]
#    par_arg_comm = IIV_arg[4]
#    red_data     = par_arg_comm[0]
#    mech_data    = par_arg_comm[1]
#    n_points     = par_arg_comm[2]
    _t           = IIV_arg[5]

    ct.suppress_thermo_warnings()
    gas             = ct.Solution('csp_to.cti')
#    nReaction       = gas.n_reactions
    nSpecies        = gas.n_species
    species_names   = gas.species_names
    for s_i in range(len(species_names)):
        species_names[s_i] = species_names[s_i].replace(',','\,')
        species_names[s_i] = species_names[s_i].replace('(','ob')
        species_names[s_i] = species_names[s_i].replace(')','cb')

#    species_list    = dynamicsymbols(species_names)

#    dicts_conc = dict(zip(species_list,conc))
    Q = np.zeros((nSpecies,nSpecies,nFastMode),dtype=float)
#    radicals = []
    radIdx = []

    for i in range(nFastMode):
#        if i == nFastMode-1 and _t >7:
#            print('\n\n\n\n\n last \n\n')
#            print(str(a_refined[:,i]*b_refined[i,:]))
#        Q[:,:,i] = np.array((a_refined[:,i]*b_refined[i,:]).subs(dicts_conc),dtype=float)
        Q[:,:,i] = np.array((a_refined[:,i]*b_refined[i,:]),dtype=float)
#
#        if i == nFastMode-1 and _t >7:
#            print('\n\n\n\n\n list  \n\n' + str(_t) + '\n\n\n')
#            print(Q)
#            print('\n\n\n\n\n Q diag :   \n\n' + str(_t) + '\n\n\n')
#            print(abs(np.diagonal(Q[:,:,i])-1))
#            print('\n\n\n\n\n min abs Q :   \n\n' + str(_t) + '\n\n\n')
#            print(min(abs(np.diagonal(Q[:,:,i])-1)))

        idx = np.where(abs(np.diagonal(Q[:,:,i])-1)==min(abs(np.diagonal(Q[:,:,i])-1)))[0][0]

#        radicals.append((i+1,idx))
        radIdx.append(idx)

#    g_fast = sym.Matrix(np.sum(Q,axis=2))*g
#    g_slow = sym.Matrix(np.eye(nSpecies)-np.sum(Q,axis=2))*g
#    radicalSpecies_list = [species_list[i] for i in radIdx]
#    nonRadicalSpecies_list = [species for species in species_list if species not in radicalSpecies_list]

    radicalSpecies_list = []
    nonRadicalSpecies_list = []
    return Q, radicalSpecies_list, nonRadicalSpecies_list, radIdx, _t




def phase_space_decomposition(red_data, arg_psd,select_tr):

    epsilon_rel    = red_data.red_op.epsilon_rel
    epsilon_abs_op = red_data.red_op.epsilon_abs
    TimeResolution = red_data.red_op.TimeResolution
    Tau            = arg_psd[0]
    species_list   = arg_psd[1]
    conc           = arg_psd[2]
    nSpecies       = arg_psd[3]
    n_points       = arg_psd[4]
    show_plots     = arg_psd[5]
    verbose        = arg_psd[6]

    nFastMode = []


    if select_tr:
        _t = 0
        for t in range(n_points-int(n_points/red_data.red_op.n_points+1)):
            if t!=0 and t%int(max(n_points/red_data.red_op.n_points,1))==0 :
                nFast = 0
                for i in range(nSpecies):
                    if Tau[_t][i] <=TimeResolution:
                        nFast+=1
                nFastMode.append(nFast)
                _t += 1

    else:

        a_refined      = arg_psd[7]
        b_refined      = arg_psd[8]
        g              = arg_psd[9]

        # The criterion used here to determine the number M of exhausted modes is based on
        # an error vector Yerror. Equation 6 from Valorani paper: M. valorani, Combustion and Flame 146 (2006) 29–51:
        f = np.zeros((len(conc),nSpecies),dtype=float)
        _t = 0
        for t in range(n_points-int(n_points/red_data.red_op.n_points+1)):
            if t!=0 and t%int(max(n_points/red_data.red_op.n_points,1))==0 :
                for i in range(nSpecies):
                    dict_ = dict(zip(species_list,conc[_t,:]))
                    f[_t, i] = np.array((b_refined[_t][i,:]*g[_t]).subs(dict_),dtype=float)
                _t += 1

        epsilon_abs = epsilon_abs_op*np.ones((nSpecies,),dtype=float)
    #    epsilon_rel = float(input('Enter a value in the range [1e-1-1e-3] (epsilon relative) ;\n(e.g. 1e-2 (typical), )'))
    #    _ = float(input('Enter a value in the range [1e-10-1e-14] (epsilon absolute) ;\n(e.g. 1e-14 (typical), )'))
    #    epsilon_abs = _*np.ones((nSpecies,),dtype=float)

        y_error = []
        _t = 0
        for t in range(n_points-int(n_points/red_data.red_op.n_points+1)):
            if t!=0 and t%int(max(n_points/red_data.red_op.n_points,1))==0 :
                dict_ = dict(zip(species_list,conc[_t,:]))
                y_vec = np.array([species.subs(dict_) for species in species_list],dtype=float)
                y_error.append(epsilon_rel*np.absolute(y_vec)+epsilon_abs)
                _t += 1
        # implementing equation 8 of Valorani paper M. valorani, Combustion and Flame 146 (2006) 29–51

        _t = 0
        for t in range(n_points-int(n_points/red_data.red_op.n_points+1)):
            if t!=0 and t%int(max(n_points/red_data.red_op.n_points,1))==0 :
        #     print('Time stamp: {}'.format(iTime+1))
                for i in range(sum(~np.isinf(Tau[_t]))-2,-1,-1):
                    sum_ = 0
                    for r in range(i):
                        sum_+= (np.array(a_refined[_t][:,r],dtype=float))*f[_t,r]
                    if(np.absolute(Tau[_t][i+1]*sum_) < y_error[_t]).all():
                        break
                    n = i+1
                nFastMode.append(n)
                _t += 1


        if verbose > 3 and show_plots:
            plt.plot(nFastMode,'s')
            plt.grid()

    return nFastMode


def reactions_withdrawal(conditions,red_data,mech_data,red_results,eps):

    # The simplifications achieved by the CSP-derived reduced model depends
    # not only on the user-specified error tolerance thresholds, but also on the
    # user’s selection of the set of primary species of interest

#    mp = conditions.main_path

    if 'Valorani' in red_data.red_op.csp_method:

# =============================================================================
#%%    7a- Simplification Algorithm Valorani, Combustion and Flame 146 (2006) 29–51
# =============================================================================
#        clock = cdef.Clock('7- Simplification Algorithm (Valorani 2006).') ; clock.start()

        active_species, active_reactions = simplificationAlgorithm(\
                                conditions,red_data,mech_data,red_results,eps)
                                                                          #nTime, tsp_name, tsp_idx, tol, mech_data)

##        for i, tol in enumerate(tolerances):
#        plt.semilogx(tolerances[0],active_species.count(True),'sb')
#        plt.grid()
#        plt.xlabel('tolerance')
#        plt.ylabel('# active species')
#        plt.show()
##        for i, tol in enumerate(tolerances):
#        plt.semilogx(tolerances[0],len(active_reactions(True)),'sg')
#        plt.grid()
#        plt.xlabel('tolerance')
#        plt.ylabel('# active reactions')
#        plt.show()

    elif 'Lam' in red_data.red_op.csp_method:

# =============================================================================
#%%    7b- Simplification Algorithm Lam and Goussis (1993)
# =============================================================================
#        clock = cdef.Clock('7- Simplification Algorithm (Lam 1993).') ; clock.start()

        a_refined       = red_data.red_op.a_refined
        b_refined       = red_data.red_op.b_refined
        nFastModeGuess  = red_data.red_op.nFastModeGuess
        conc            = red_data.red_op.csp_conc
        nTime           = len(conc)

#        start = time.time()
        P = []
        for iTime in range(nTime):
        #     print(iTime)
            _P_ = participationIndex_Goussis(b_refined[iTime], nFastModeGuess[iTime], conc[iTime,:])
            P.append(_P_)
        del _P_

        print('Summation of absolute values........\n')
        for iTime in range(nTime):
            print('at time: {}'.format(iTime))
            for iSpecies in range(nFastModeGuess[iTime]):
                print('For species {}: {}'.format(mech_data.spec.name[iSpecies], np.sum(np.absolute(P[iTime][iSpecies,:]))))

        I = []
        for iTime in range(nTime):
        #     print(iTime)
            _I_ = importanceIndex_Goussis(a_refined[iTime], b_refined[iTime], nFastModeGuess[iTime], conc[iTime,:])
            I.append(_I_)
        del _I_


        print('Summation of absolute values........\n')
        for iTime in range(nTime):
            print('at time: {}'.format(iTime))
            for iSpecies in range(nFastModeGuess[iTime]):
                print('For species {}: {}'.format(mech_data.spec.name[iSpecies], np.sum(np.absolute(I[iTime][iSpecies,:]))))
#        end = time.time()


        print('Simplification Algorithm.....................\n\n')
        _ = list((-1)*np.arange(0,30,dtype=float)[::-1])
#        tolerances = [pow(10,i) for i in _]

        tolerances = [0.01]

        S_global, R_global = [], []

        for tol in tolerances:
            S, R = simplificationAlgorithmGoussis(conditions, red_data, I, P, eps)
            S_global.append(S)
            R_global.append(R)

        for i, tol in enumerate(tolerances):
            plt.semilogx(tol,len(S_global[i]),'sb')
        plt.grid()
        plt.xlabel('tolerance')
        plt.ylabel('# active species')
        plt.show()
        for i, tol in enumerate(tolerances):
            plt.semilogx(tol,len(R_global[i]),'sg')
        plt.grid()
        plt.xlabel('tolerance')
        plt.ylabel('# active reactions')
        plt.show()

        #.............................End of Simplification algorithm based on Lam and Goussis paper 1993...................................

#    clock.stop() ; clock.display_(mp)

    return active_reactions, active_species #, CSP_radicals



def simplificationAlgorithm(conditions,red_data,mech_data,red_results,eps):

    gas_ref  = conditions.composition.gas_ref
#    gas_red  = red_results.gas
    nTime    = int(red_data.red_op.n_points)
    tsp_idx  = red_data.targetSpeciesIdx
#    tol      = eps
    radIdx   = red_data.red_op.radIdx
    I_fast   = red_data.red_op.I_fast
    I_slow   = red_data.red_op.I_slow
    verbose         = conditions.simul_param.verbose


#        list_reaction = [str(reaction) for reaction in gas.reactions()]
    r_formulas = mech_data.react.formula
    sp_name    = mech_data.spec.name
#        R_set = set()
#        S_set = set(tsp_name)
#    if verbose>6: print('Tolerance: {} and kernel species: {}\n'.format(tol,red_data.tspc))
#        S = set()
#[       S_global = set()
#        R_global = set()
    active_reactions = copy.deepcopy(mech_data.react.activ_p)
    active_species = copy.deepcopy(mech_data.spec.activ_p)
    active_species_f = copy.deepcopy(mech_data.spec.activ_p)

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

    for t_i in tsp_idx:
        active_species[t_i] = True

    for iTime in range(1,nTime-1):

        while (active_species!=active_species_f):
            active_species_f = copy.deepcopy(active_species)
            sp_red = -1
            for sp_ref in range(len(sp_name)):
                if active_species[sp_ref]:
                    sp_red +=1
                    if sp_ref in tsp_idx:
                        idx_ = tsp_idx.index(sp_ref)
                        tol = eps[idx_]
                    else:
                        tol = min(eps)
                    r_red = 0
                    for r_ref in range(len(r_formulas)):
                        if mech_data.react.activ_m[r_ref]:
                            if ((sp_name[sp_ref] in r_formulas[r_ref].split()) and (I_slow[iTime][sp_red,r_red]>tol) \
                                or ((I_fast[iTime][radIdx[iTime],r_red]>tol).all())):
        #                        if  ((species in reaction.split()) and (I_slow[iTime][species_names.index(species),list_reaction.index(reaction)]>tol) \
        #                        and ((I_fast[iTime][radIdx[iTime],list_reaction.index(reaction)]>tol).all())):
                                if verbose>8 and not active_reactions:
                                    print_('\nIncluding reaction {}, as it contains species {}, \nwith I_slow: {} and I_fast: {}'.format(r_formulas[r_ref],sp_name[sp_ref],I_slow[iTime][sp_red,r_red],I_fast[iTime][radIdx[iTime],r_red]))
        #                            active_reactions[r] = True
                                active_reactions[r_ref] = True
                            r_red +=1

            for r_ref in range(len(r_formulas)):
                if active_reactions[r_ref]:
#                    for reaction in R_set:
                    for sp_ref in range(len(sp_name)):
                        if not active_species[sp_ref]:
                            if (sp_name[sp_ref] in r_formulas[r_ref].split()):
                                if verbose>8 and not active_species[sp_ref]:
                                    print_('\nIncluding species {}, as it appears in reaction {}.'.format(sp_name[sp_ref], r_formulas[r_ref]))
    #                            S_set.add(species)
                                active_species[sp_ref] = True

#                S_global = S_global.union(S_set)

#                for sp in range(len(sp_name)):
#                    if active_species[sp]:
#                        for r in range(len(r_formulas)):
#                            if sp_name[sp] in r_formulas[r].split():
#        #                        if  ((species in reaction.split()) and (I_slow[iTime][species_names.index(species),list_reaction.index(reaction)]>tol) \
#        #                        and ((I_fast[iTime][radIdx[iTime],list_reaction.index(reaction)]>tol).all())):
#                                if verbose>8 and not active_reactions[r]:
#                                    print_('\nIncluding reaction '+ r_formulas[r] + ' as it contains species ' + sp_name[sp])
#        #                            active_reactions[r] = True
#                                active_reactions[r] = True
#
#                            reaction_list = list_reaction[r].split()
#                            reaction_list[:] = (value for value in reaction_list if value != '+')
#                            if (species in reaction_list):
#                                R_set.add(list_reaction[iReaction])
#                    R_global = R_global.union(R_set)
#        if verbose>7 and not active_reactions[r]:
#            print('\n\n Active species set: {}'.format(S_global))
#            print('\n\nActive reaction set (#{}): \n{}'.format(len(R_global),R_global))
#            print('----------------------')

    return active_species, active_reactions



#    def simplificationAlgorithm(nTime, tsp_name, tsp_idx, tol, mech_data):
##
##
#        list_reaction = [str(reaction) for reaction in gas.reactions()]
#
#        R_set = set()
#        S_set = set(tsp_name)
#        print('Tolerance: {} and kernel species: {}\n'.format(tol,tsp_name))
#        S = set()
#        S_global = set()
#        R_global = set()
#        for iTime in range(1,nTime-1):
#
#            while (not(S==S_set)):
#                S = S_set
#                for species in S_set:
#                    for reaction in list_reaction:
#                        if ((species in reaction.split()) and (I_slow[iTime][species_names.index(species),list_reaction.index(reaction)]>tol) and ((I_fast[iTime][radIdx[iTime],list_reaction.index(reaction)]>tol).all())):
#                            print('\nIncluding reaction {}, as it contains species {}, \nwith I_slow: {} and I_fast: {}'.format(reaction,species,I_slow[iTime][species_names.index(species),list_reaction.index(reaction)],I_fast[iTime][radIdx[iTime],list_reaction.index(reaction)]))
#                            R_set.add(str(reaction))
#
#                for reaction in R_set:
#                    for species in species_names:
#                        if (species in reaction.split() and not(float(Sr[species_names.index(species),list_reaction.index(reaction)])==0)):
#                            print('\nIncluding species {}, as it appears in reaction {}.'.format(species, str(reaction)))
#                            S_set.add(species)
#
#                S_global = S_global.union(S_set)
#
#                for species in S_set:
#                    for iReaction in range(nReaction):
#                        reaction_list = list_reaction[iReaction].split()
#                        reaction_list[:] = (value for value in reaction_list if value != '+')
#                        if (species in reaction_list):
#                            R_set.add(list_reaction[iReaction])
#                R_global = R_global.union(R_set)
#        print('\n\n Active species set: {}'.format(S_global))
#        print('\n\nActive reaction set (#{}): \n{}'.format(len(R_global),R_global))
#        print('----------------------')
#
#        return S_global, R_global



def simplificationAlgorithmGoussis(conditions, red_data, I, P, eps):

    gas             = conditions.composition.gas_ref
    species_names   = gas.species_names
#    from sympy.physics.vector import dynamicsymbols
#    species_list    = dynamicsymbols(species_names)
#    ns              = gas.n_species()
    tol             = eps
    tsp_name        = red_data.tspc
    nr              = gas.n_reactions()
    Sr              = red_data.red_op.csp_Sr
    nTime           = int(red_data.red_op.n_points)

    verbose         = conditions.simul_param.verbose


#        nTime, tsp_name, tol):
#    gas, species_names, I



    list_reaction = [str(reaction) for reaction in gas.reactions()]

    R_set = set()
    S_set = set(tsp_name)
    if verbose>6: print('Tolerance: {} and kernel species: {}\n'.format(tol,tsp_name))
    S = set()
    S_global = set()
    R_global = set()
    for iTime in range(1,nTime-1):

        while (not(S==S_set)):
            S = S_set
            for species in S_set:
                for reaction in list_reaction:
                    if (species in reaction and np.absolute(I[iTime][species_names.index(species),list_reaction.index(reaction)])>tol and (np.absolute(I[iTime][species_names.index(species),list_reaction.index(reaction)])>tol).all):
                        print('\nIncluding reaction {}, as it contains species {}, \nwith |I|: {} and |I_radical|: {}'.format(reaction,species,np.absolute(I[iTime][species_names.index(species),list_reaction.index(reaction)]),np.absolute(I[iTime][species_names.index(species),list_reaction.index(reaction)])))
                        R_set.add(str(reaction))

            for reaction in R_set:
                for species in species_names:
                    if (species in reaction.split() and not(float(Sr[species_names.index(species),list_reaction.index(reaction)])==0)):
                        print('\nIncluding species {}, as it appears in reaction {}.'.format(species, str(reaction)))
                        S_set.add(species)

            S_global = S_global.union(S_set)

            for species in S_set:
                for iReaction in range(nr):
                    reaction_list = list_reaction[iReaction].split()
                    reaction_list[:] = (value for value in reaction_list if value != '+')
                    if (species in reaction_list):
                        R_set.add(list_reaction[iReaction])
            R_global = R_global.union(R_set)
    if verbose>8: print('\n\n Active species set: {}'.format(S_global))
    if verbose>8: print('\n\nActive reaction set (#{}): \n{}'.format(len(R_global),R_global))
    if verbose>8: print('----------------------')

    return S_global, R_global

    #..........................End of functins........................................................

def participationIndex_Goussis(conditions, red_data):

    gas             = conditions.composition.gas_ref
    species_names   = gas.species_names
    from sympy.physics.vector import dynamicsymbols
    species_list = []
    for species in species_names:
        species_list.append(dynamicsymbols(species))
    ns              = gas.n_species()
    nr              = gas.n_reactions()
    Sr              = red_data.red_op.csp_Sr
    F               = red_data.red_op.csp_F
    b_refined       = red_data.red_op.b_refined
    nFastMode       = red_data.red_op.nFastMode
    conc            = red_data.red_op.conc

    b = np.array(b_refined,dtype=float)
    P = np.zeros((nFastMode,nr),dtype=float)
    dict_ = dict(zip(species_list,conc))

    sum_ = np.zeros((nFastMode,1))
    for m in range(nFastMode):
        for k in range(nr):
            sum_[m,0] +=np.absolute(float(np.matmul(np.reshape(b[m,:],(1,ns)),(np.array(Sr[:,k],dtype=float)*float(F[k].subs(dict_))))))

    for m in range(nFastMode):
        for j in range(nr):
            P[m,j] = float(np.matmul(np.reshape(b[m,:],(1,ns)),(np.array(Sr[:,j],dtype=float)*float(F[j].subs(dict_)))))/sum_[m,0]

    return P


def importanceIndex_Goussis(a_refined, b_refined, nFastMode, conc,  conditions, red_data):

    gas             = conditions.composition.gas_ref
    species_names   = gas.species_names
    from sympy.physics.vector import dynamicsymbols
    species_list = []
    for species in species_names:
        species_list.append(dynamicsymbols(species))
    ns              = gas.n_species()
    nr              = gas.n_reactions()
    Sr              = red_data.red_op.csp_Sr
    F               = red_data.red_op.csp_F

    I = np.zeros((ns,nr),dtype=float)
    dict_ = dict(zip(species_list,conc))

    s = nFastMode
    a = np.zeros((ns,ns))
    a[s:,s:] = np.array(a_refined[s:,s:],dtype=float)
    b = np.zeros((ns,ns))
    b[s:,s:] = np.array(b_refined[s:,s:],dtype=float)
#     A = np.matmul(a,b)
    A = np.matmul(np.array(a_refined,dtype=float), np.array(b_refined,dtype=float))

    sum_ = np.zeros((nFastMode,1))
    for n in range(nFastMode):
        for j in range(nr):
            sum_[n,0] +=np.absolute(float(np.matmul(np.reshape(A[n,:],(1,ns)),(np.array(Sr[:,j],dtype=float)*float(F[j].subs(dict_))))))

    for n in range(nFastMode):
        for j in range(nr):
            I[n,j] = float(np.matmul(np.reshape(A[n,:],(1,ns)),(np.array(Sr[:,j],dtype=float)*float(F[j].subs(dict_)))))/sum_[n,0]

    return I

