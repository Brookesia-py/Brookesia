import random, copy, os, sys, warnings, multiprocessing, csv

import numpy                    as np
import brookesia.Computation    as comp
import brookesia.Class_def      as cdef
import cantera                  as ct

from   brookesia.Class_def      import print_
from   shutil                   import copyfile
from   pebble                   import ProcessPool
from   concurrent.futures       import TimeoutError



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
        # Activer la conversion des overflows en exceptions
        for i, T in enumerate(T_list):
            Y_r.append(A * T ** n * np.exp(-Ea / (R * T)))
            if not _react.k_min[_r][i] < A * T ** n * np.exp(-Ea / (R * T)) < _react.k_max[_r][i]:
                # print(str(_r))
                # print(_react.equation[_r])
                # print(_react.number[_r])
                # print(_react.type[_r])
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

            if not   _react.k_min[_r][i]    < k    < _react.k_max[_r][i]:
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
    nr = len(mech_data.react.type)
    mech_data.react.k_ref    = [False for _ in range(nr)]
    mech_data.react.k_min    = [False for _ in range(nr)]
    mech_data.react.k_max    = [False for _ in range(nr)]
    mech_data.react.k_ref_lp = [False for _ in range(nr)]
    mech_data.react.k_min_lp = [False for _ in range(nr)]
    mech_data.react.k_max_lp = [False for _ in range(nr)]
    mech_data.react.k_ref_hp = [False for _ in range(nr)]
    mech_data.react.k_min_hp = [False for _ in range(nr)]
    mech_data.react.k_max_hp = [False for _ in range(nr)]
    mech_data.react.f_min    = [False for _ in range(nr)]
    mech_data.react.f_lp_min = [False for _ in range(nr)]
    mech_data.react.f_hp_min = [False for _ in range(nr)]

    for _r in range(nr):
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
        mech_data.react.f_min[_r] = min(min(f_T),5)


        # low pressure
        f_T_lp,k0_lp = [],[]
        cT  = np.poly1d(mech_data.react.f_T_fit_lp[_r])
        ck0 = np.poly1d(mech_data.react.k0_fit_lp[_r])
        f_T_lp   = cT(T_list)
        k0_lp    = 10**ck0(T_list)
        mech_data.react.f_lp_min[_r] = min(min(f_T_lp),5)

        # high pressure
        f_T_hp,k0_hp = [],[]
        cT  = np.poly1d(mech_data.react.f_T_fit_hp[_r])
        ck0 = np.poly1d(mech_data.react.k0_fit_hp[_r])
        f_T_hp   = cT(T_list)
        k0_hp    = 10**ck0(T_list)
        mech_data.react.f_hp_min[_r] = min(min(f_T_hp),5)

        if k0 is not False: # k0 = k(T) computed from the data of the file uncertainties_f.csv
            if mech_data.react.type[_r] == 'three_body_reaction':
                M = P/(8.314/T_list)
                k0 = np.array(k0)/M
            # f_T < 5 to avoid overflow
            mech_data.react.k_ref[_r] = np.array(k0)
            mech_data.react.k_min[_r] = np.array(k0)/(10**np.minimum(f_T,5))  - 10**-3
            mech_data.react.k_max[_r] = np.array(k0)*(10**np.minimum(f_T,5))  + 10**-3        # 10**-3 offset to prevent optimization issues caused by extremely small values
            # mech_data.react.k_min.append(np.array(k0)*np.exp(np.array(-f_T) * np.log(10)))
            # mech_data.react.k_max.append(np.array(k0)*np.exp(np.array(f_T) * np.log(10)))

            mech_data.react.k_ref_lp[_r] = np.array(k0_lp)
            mech_data.react.k_min_lp[_r] = np.array(k0_lp)/(10**np.minimum(f_T_lp,5))  - 10**-3
            mech_data.react.k_max_lp[_r] = np.array(k0_lp)*(10**np.minimum(f_T_lp,5))  + 10**-3
            # mech_data.react.k_min_lp.append(np.array(k0_lp)*np.exp(np.array(-f_T_lp) * np.log(10)))
            # mech_data.react.k_max_lp.append(np.array(k0_lp)*np.exp(np.array(f_T_lp) * np.log(10)))

            mech_data.react.k_ref_hp[_r] = np.array(k0_hp)
            mech_data.react.k_min_hp[_r] = np.array(k0_hp)/(10**np.minimum(f_T_hp,5))  - 10**-3
            mech_data.react.k_max_hp[_r] = np.array(k0_hp)*(10**np.minimum(f_T_hp,5))  + 10**-3
            # mech_data.react.k_min_hp.append(np.array(k0_hp)*np.exp(np.array(-f_T_hp) * np.log(10)))
            # mech_data.react.k_max_hp.append(np.array(k0_hp)*np.exp(np.array(f_T_hp) * np.log(10)))

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

            mech_data.react.k_ref[_r] = np.array(k_ref)
            # f_T < 5 to avoid overflow
            mech_data.react.k_min[_r] = np.array(k_ref)/(10**np.minimum(f_T,5))   - 10**-3
            mech_data.react.k_max[_r] = np.array(k_ref)*(10**np.minimum(f_T,5))   + 10**-3  # 10**-3 offset to prevent optimization issues caused by extremely small values
            if len(k_ref_lp)>0:
                mech_data.react.k_ref_lp[_r] = np.array(k_ref_lp)
                mech_data.react.k_min_lp[_r] = np.array(k_ref_lp)/(10**np.minimum(f_T,5))  - 10**-3
                mech_data.react.k_max_lp[_r] = np.array(k_ref_lp)*(10**np.minimum(f_T,5))  + 10**-3

                mech_data.react.k_ref_hp[_r] = np.array(k_ref_hp)
                mech_data.react.k_min_hp[_r] = np.array(k_ref_hp)/(10**np.minimum(f_T,5))  - 10**-3
                mech_data.react.k_max_hp[_r] = np.array(k_ref_hp)*(10**np.minimum(f_T,5))  + 10**-3

            else:
                mech_data.react.k_min_lp[_r] = False
                mech_data.react.k_max_lp[_r] = False
                mech_data.react.k_ref_hp[_r] = np.array(k_ref_hp)
                mech_data.react.k_min_hp[_r] = False
                mech_data.react.k_max_hp[_r] = False

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


def get_uncertainty(mech_data,optim_param,conditions_list,opt):

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

    mech_data.react.uncert = []
    for r in range(len(mech_data.react.equation)):
        react_found=False
        for r_inc in range(len(uncertainty)):
            if int(num_react[r_inc])==mech_data.react.number[r]\
            and uncertainty[r_inc]!='':
                mech_data.react.uncert.append([float(u) for u in uncertainty[r_inc].split(',')])
                react_found = True
                break
        if not react_found:
            mech_data.react.uncert.append(optim_param.Arrh_max_variation)


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
                    react_f['f_T_fit'][-1] = optim_param.f_default
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
                    react_f['f_T_fit_lp'][-1] = optim_param.f_default
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
                    react_f['f_T_fit_hp'][-1] = optim_param.f_default
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

    os.chdir(opt)

    return mech_data

def compare_best_ind(pop, best_ind, optim_param, opt, verbose=0):

    mp = optim_param.main_path
    n_ind = optim_param.n_ind
    global_best_idx = pop.find_best()

    new_best_ind = False
    # compare the current best population ind to the previous best ind
    if pop.individual[global_best_idx].fitness > best_ind.fitness:
        best_ind = copy.deepcopy(pop.individual[global_best_idx])
        os.chdir(mp+'/'+opt)
        if '.cti' in pop.individual[0].mech.name:
            best_ind.mech.write_new_mech("optim_mech.cti")
        else:
            best_ind.mech.write_yaml_mech("optim_mech.yaml")

        if verbose >= 3:
            print_("New best_ind: "+"%.3f" %(best_ind.fitness),mp)
        new_best_ind = True

    next_pop_best_idx = pop.find_best(n_ind)
    if pop.individual[next_pop_best_idx].fitness < best_ind.fitness:
        # best_ind of all generations replace current worst ind
            worst_idx=pop.find_worst(next_pop_best_idx,n_ind)
            pop.individual[worst_idx]=copy.deepcopy(best_ind)

    return pop, best_ind, new_best_ind


def fitness_eval_par(fit_eval_inp):

    pop                 = fit_eval_inp[0]
    optim_param         = fit_eval_inp[1]
    conditions_list     = fit_eval_inp[2]
    ref_results_list    = fit_eval_inp[3]
    bar                 = fit_eval_inp[4]
    ind                 = fit_eval_inp[5]
    is_child            = fit_eval_inp[6]
    opt                 = fit_eval_inp[7]
    if is_child:
        if opt == 'GA': ind = optim_param.n_ind + ind
        title = "New ind evaluation  "
    else:
        title = "First pop evaluation  "

    mech = conditions_list[0].mech

    os.chdir(conditions_list[0].main_path)

    # gas = cdef.get_gas_ct(mech)

    # for _c in range(len(conditions_list)):
    #     conditions_list[_c].composition.gas     = gas
    #     conditions_list[_c].composition.gas_ref = gas

    os.chdir(opt)

    try:
        fitness = pop.individual[ind].fitness_eval(conditions_list, optim_param, ref_results_list, opt, ind)
    except:
        fitness = 0

    bar.update(ind, title)

    return (fitness, ind)

# def fitness_eval_newpop(pop, optim_param, conditions_list, ref_results_list, opt):

#     mp = optim_param.main_path

#     # n_ind = optim_param.n_ind
#     n_ind = len(pop.individual)

#     gas = conditions_list[0].composition.gas
#     gas_ref = conditions_list[0].composition.gas_ref

#     # saving and suppression of unpickable variables on fitness eval inputs
#     for cond in range(len(conditions_list)):
#         del conditions_list[cond].composition.gas
#         del conditions_list[cond].composition.gas_ref
#     f = []
#     simul_time_limit = 0
#     for res in range(len(ref_results_list)):
#         f.append(ref_results_list[res].f)
#         del ref_results_list[res].gas
#         del ref_results_list[res].f
#         simul_time_limit += ref_results_list[res].simul_time
#     simul_time_limit = simul_time_limit * 3

#     bar = cdef.ProgressBar(n_ind, '')
#     title = "First pop evaluation  "
#     bar.update(0, title)

#     fit_eval_inp = []
#     for _i in range(n_ind):
#         fit_eval_inp.append([pop, optim_param, conditions_list, ref_results_list, bar, _i, False, opt])


#     # saving of the ref mechanism
#     mech = conditions_list[0].mech
#     os.chdir(conditions_list[0].main_path)
#     copyfile(mech, opt + '/' + mech)
#     os.chdir(opt)

#     if sys.gettrace() is not None:
#         is_debug_mode = True
#     else:
#         is_debug_mode = False

#     # is_debug_mode = True

#     if is_debug_mode:
#         print('Evaluation of first generation in debug mode (i.e. no parallelization)')
#         # bypass parallelisation for debugging
#         fit_i = []
#         mech = conditions_list[0].mech

#         # gas = cdef.get_gas_ct(mech)
#         # for _c in range(len(conditions_list)):
#         #     conditions_list[_c].composition.gas     = gas
#         #     conditions_list[_c].composition.gas_ref = gas
#         for ind in range(n_ind):
#             fit_i.append(
#                 [pop.individual[ind].fitness_eval(conditions_list, optim_param, ref_results_list, opt, ind), ind])

#     else:
#         # Parallelized Fitness calculation
#         num_cores = multiprocessing.cpu_count()
#         if os.name == 'nt': multiprocessing.get_context('spawn')
        
#         with ProcessPool() as pool:
#         # with multiprocessing.Pool(num_cores) as p:
#             fit_i = pool.map(fitness_eval_par, fit_eval_inp)            
#         # sort in ascending order of individual index values
#         fit_i.sort(reverse=False, key=lambda col: col[1])

#     for _i in range(len(fit_i)):
#         pop.individual[_i].fitness = fit_i[_i][0]

#     bar.update(n_ind, title)
#     print('\n')

#     for i in range(len(conditions_list)):
#         conditions_list[i].composition.gas = gas
#         conditions_list[i].composition.gas_ref = gas_ref
#     for i in range(len(ref_results_list)):
#         ref_results_list[i].gas = gas
#         ref_results_list[i].f = f[i]
#     return pop


def fitness_eval_newchilds(pop, optim_param, conditions_list, ref_results_list, gen, opt):

    mp = optim_param.main_path
    num_cores = multiprocessing.cpu_count()

    if gen == 0:
        ind_nb   = len(pop.individual)
        is_child = False
    else: 
        is_child = True
        if opt == 'GA':
            ind_nb = optim_param.total_Xover + optim_param.total_mut
        elif opt == 'PSO':
            ind_nb = optim_param.n_ind

    gas = conditions_list[0].composition.gas
    gas_ref = conditions_list[0].composition.gas_ref

    # saving and suppression of unpickable variables on fitness eval inputs
    for cond in range(len(conditions_list)):
        del conditions_list[cond].composition.gas
        del conditions_list[cond].composition.gas_ref
    f = []
    simul_time_limit = 0
    for res in range(len(ref_results_list)):
        f.append(ref_results_list[res].f)
        del ref_results_list[res].gas
        del ref_results_list[res].f
        simul_time_limit += ref_results_list[res].simul_time
    #        simul_time_limit = simul_time_limit*(6+np.random.uniform()*16)#*(child_nb/num_cores)
    simul_time_limit = simul_time_limit * 2 * (ind_nb / num_cores) + 30
    #        simul_time_limit = simul_time_limit*5 + 30

    bar = cdef.ProgressBar(ind_nb, '')
    if gen == 0:    title = "First pop evaluation  "
    else:           title = "New ind evaluation  "
    bar.update(0, title)

    fit_eval_inp = []
    for ind in range(ind_nb):
        fit_eval_inp.append([pop, optim_param, conditions_list, ref_results_list, bar, ind, is_child, opt])

    # saving of the ref mechanism
    mech = conditions_list[0].mech
    os.chdir(conditions_list[0].main_path)
    copyfile(mech, opt + '/' + mech)
    os.chdir(opt)

    # Fitness calculation
    if sys.gettrace() is not None:
        is_debug_mode = True
    else:
        is_debug_mode = False
    # is_debug_mode = True

    # no parallelisation --------------------------------------------------
    if is_debug_mode:
        fit_list = []
        mech = conditions_list[0].mech
        gas = cdef.get_gas_ct(mech)
        for _c in range(len(conditions_list)):
            conditions_list[_c].composition.gas = gas
            conditions_list[_c].composition.gas_ref = gas
        for ch in range(ind_nb):
            fit_list.append(
                [pop.individual[ch].fitness_eval(conditions_list, optim_param, ref_results_list, opt, ch), ch])
    else:
        # Parallelisation 4   (with simulation time check) --------------------
        #  https://pythonhosted.org/Pebble/#pools
        fit_list = []
        with ProcessPool() as pool:
            sim_results = pool.map(fitness_eval_par, fit_eval_inp, timeout=simul_time_limit)
            try:
                for fit in sim_results.result():
                    fit_list.append(fit)
            except TimeoutError:
                print_(
                    '\n\nWarning : simulation time > ' + '%.0f' % simul_time_limit + 's (> 2 x ref simulation time)',
                    mp)
                print_("TimeoutError: aborting remaining computations", mp)
                fit_list.sort(reverse=False, key=lambda col: col[1])
                fitness_incomplete = copy.deepcopy(fit_list)
                list_ind_eval = []
                for ind_fit_inc in fitness_incomplete:
                    list_ind_eval.append(ind_fit_inc[1])
                n_sim = len(fitness_incomplete)
                for _i in range(ind_nb):
                    if gen == 0 or opt != 'GA':
                        idx_ind = _i
                    else: 
                        idx_ind = _i + optim_param.n_ind #if GA and gen>0, evaluate new childs only
                    try:
                        if idx_ind not in list_ind_eval:
                            fit_list.insert(_i, (0, _i))
                    except:
                        fit_list.append((0, _i))
                print_('Number of individuals evaluated: ' + str(n_sim) + '\n\n', mp)
                sim_results.cancel()
            except:
                print_('\n\nWarning : error in simulation', mp)
                print_("aborting remaining computations", mp)
                fit_list.sort(reverse=False, key=lambda col: col[1])
                fitness_incomplete = copy.deepcopy(fit_list)
                list_ind_eval = []
                for ind_fit_inc in fitness_incomplete:
                    list_ind_eval.append(ind_fit_inc[1])
                n_sim = len(fitness_incomplete)
                for _i in range(ind_nb):
                    if gen == 0 or opt != 'GA':
                        idx_ind = _i
                    else: 
                        idx_ind = _i + optim_param.n_ind #if GA and gen>0, evaluate new childs only                    
                    try:
                        if idx_ind not in list_ind_eval:
                            fit_list.insert(_i, (0, _i))
                    except:
                        fit_list.append((0, _i))
                print_('Number of individuals evaluated: ' + str(n_sim) + '\n\n', mp)
                sim_results.cancel()

    fit_list.sort(reverse=False, key=lambda col: col[1])

    for _i in range(len(fit_list)):
        if gen == 0 or opt != 'GA':
            ind = _i
        else: 
            ind = optim_param.n_ind + _i
        pop.individual[ind].fitness = fit_list[_i][0]

    bar.update(ind_nb, title)
    print('\n')

    for i in range(len(conditions_list)):
        conditions_list[i].composition.gas = gas
        conditions_list[i].composition.gas_ref = gas_ref
    for i in range(len(ref_results_list)):
        ref_results_list[i].gas = gas
        ref_results_list[i].f = f[i]

    return pop







class Individual:
     def __init__(self, conditions_list, mech_data, ref_results_list, red_data_list,opt,rand_kin=True):
         optim_param = red_data_list[0].optim_param
         if type(rand_kin) is bool:
             self.mech = copy.deepcopy(mech_data)
         else:  # if the option to import kinetic mechanism have been selected
             rand_kin.react.uncert =  copy.deepcopy(mech_data.react.uncert)
             self.mech = copy.deepcopy(rand_kin)
         self.r2opt = []
         self.find_r2opt(red_data_list, optim_param, conditions_list)

         if opt=='PSO':
             self.particle_speed_init(optim_param.n_gen)
             self.best_fitness  = 0
             self.fitness       = 0
             self.best_mech = copy.deepcopy(mech_data)

         if rand_kin:
             self.randomize_kin(optim_param)

     def particle_speed_init(self, Max_it):
         size_react = len(self.mech.react.kin)
         self.particle_speed = []
         uncert_r = 0.01
         for r in range(size_react):
             # uncert_r = [u / 100 for u in self.mech.react.uncert[r]]
             r1 = random.random()  # qui permet de donner une valeure aléatoire entre [0,1]
             r2 = random.random()  # qui permet de donner une valeure aléatoire entre [0,1]
             if self.mech.react.type[r] == 'three_body_reaction' \
                     or self.mech.react.type[r] == 'reaction':
                 self.particle_speed.append([0, 0, 0])
                 for k in range(3):
                     ref = self.mech.react.ref_kin[r][k]
                     self.particle_speed[r][k] = ref * uncert_r * random.uniform(-1, 1) * (10 / Max_it)

             elif self.mech.react.type[r] == 'falloff_reaction' \
                     or self.mech.react.type[r] == 'pdep_arrhenius' \
                     or self.mech.react.type[r] == 'chemically_activated_reaction' \
                     or self.mech.react.type[r] == 'chebyshev':
                 self.particle_speed.append([])
                 for j in range(len(self.mech.react.kin[r])):
                     self.particle_speed[-1].append([0, 0, 0])
                     for k in range(3):
                         ref = self.mech.react.ref_kin[r][j][k]
                         self.particle_speed[-1][j][k] = ref * uncert_r * random.uniform(-1, 1) * (10 / Max_it)
             else:
                 print(self.mech.react.type[r])

     def find_r2opt(self, red_data_list, optim_param, conditions_list):

         if optim_param.optim_on_meth != 'False':
             n_tspc = len(red_data_list[0].tspc)
             n_r2opt = optim_param.nb_r2opt  # total number of react to opt
             self.mech.react.modif = [False] * len(self.mech.react.modif)
             mp = optim_param.main_path

             # Get list of reactions to optimize
             if 'SA' in red_data_list[0].reduction_operator \
                     or optim_param.optim_on_meth == 'SA':
                 # keep maximal sensitivities
                 if len(red_data_list[0].red_op.sensi_r) != 0:
                     add_col = 0
                     for l in range(len(red_data_list)):
                         if red_data_list[l].red_op.sensi_Sl is not False \
                                 and conditions_list[l].error_param.Sl_check:
                             if red_data_list[l].red_op.sensi_T is not False \
                                     and conditions_list[l].error_param.T_check:
                                 if red_data_list[l].red_op.sensi_igt is not False \
                                         and conditions_list[l].error_param.ig_check:
                                     # sensitivity analysis on: Sl, T, igt
                                     add_col = 3
                                     idx_Sl = -3;
                                     idx_T = -2;
                                     idx_igt = -1
                                 else:
                                     # sensitivity analysis on: Sl, T
                                     add_col = max(add_col, 2)
                                     idx_Sl = -2;
                                     idx_T = -1
                             elif red_data_list[l].red_op.sensi_igt is not False \
                                     and conditions_list[l].error_param.ig_check:
                                 # sensitivity analysis on: Sl, igt
                                 add_col = max(add_col, 2)
                                 idx_Sl = -2;
                                 idx_igt = -1
                             else:
                                 # sensitivity analysis on: Sl
                                 add_col = max(add_col, 1)
                                 idx_Sl = -1
                         elif red_data_list[l].red_op.sensi_T is not False \
                                 and conditions_list[l].error_param.T_check:
                             if red_data_list[l].red_op.sensi_igt is not False \
                                     and conditions_list[l].error_param.ig_check:
                                 # sensitivity analysis on: T, igt
                                 add_col = max(add_col, 2)
                                 idx_T = -2;
                                 idx_igt = -1
                             else:
                                 # sensitivity analysis on: T
                                 add_col = max(add_col, 1)
                                 idx_T = -1
                         elif red_data_list[l].red_op.sensi_igt is not False \
                                 and conditions_list[l].error_param.ig_check:
                             # sensitivity analysis on: igt
                             add_col = max(add_col, 1)
                             idx_igt = -1

                     # create max_sens_list
                     max_sens_list = np.zeros(
                         (len(red_data_list[0].red_op.sensi_r) + add_col, len(red_data_list[0].red_op.sensi_r[0])))

                     # fill max_sens_list   - species sensitivities
                     for l in range(len(red_data_list)):
                         for idx in range(len(red_data_list[l].red_op.sensi_r)):
                             for r in range(len(red_data_list[l].red_op.sensi_r[idx])):
                                 if abs(red_data_list[l].red_op.sensi_r[idx][r]) > max_sens_list[idx][r]:
                                     max_sens_list[idx][r] = abs(red_data_list[l].red_op.sensi_r[idx][r])

                     # fill max_sens_list...
                     for l in range(len(red_data_list)):
                         #   - Sl sensitivities
                         if red_data_list[l].red_op.sensi_Sl is not False \
                                 and conditions_list[l].error_param.Sl_check:
                             if len(max_sens_list[idx_Sl]) == 1:
                                 max_sens_list[idx_Sl] = red_data_list[l].red_op.sensi_Sl
                             else:
                                 for r in range(len(red_data_list[l].red_op.sensi_Sl)):
                                     if abs(red_data_list[l].red_op.sensi_Sl[r]) > max_sens_list[idx_Sl][r] \
                                             and conditions_list[l].error_param.Sl_check:
                                         max_sens_list[idx_Sl][r] = abs(red_data_list[l].red_op.sensi_Sl[r])
                         #   - T sensitivities
                         if red_data_list[l].red_op.sensi_T is not False \
                                 and conditions_list[l].error_param.T_check:
                             if len(max_sens_list[idx_T]) == 1:
                                 max_sens_list[idx_T] = red_data_list[l].red_op.sensi_T
                             else:
                                 for r in range(len(red_data_list[l].red_op.sensi_T)):
                                     if abs(red_data_list[l].red_op.sensi_T[r]) > max_sens_list[idx_T][r] \
                                             and conditions_list[l].error_param.T_check:
                                         max_sens_list[idx_T][r] = abs(red_data_list[l].red_op.sensi_T[r])
                         #   - igt sensitivities
                         if red_data_list[l].red_op.sensi_igt is not False \
                                 and conditions_list[l].error_param.ig_check:
                             if len(max_sens_list[idx_igt]) == 1:
                                 max_sens_list[idx_igt] = red_data_list[l].red_op.sensi_igt
                             else:
                                 for r in range(len(red_data_list[l].red_op.sensi_igt)):
                                     if abs(red_data_list[l].red_op.sensi_igt[r]) > max_sens_list[idx_igt][r] \
                                             and conditions_list[l].error_param.ig_check:
                                         max_sens_list[idx_igt][r] = abs(red_data_list[l].red_op.sensi_igt[r])
                     if len(max_sens_list) == 1:
                         print_('Warning, no sensitivity data', mp)

                 n_r2opt_sp_max = round(n_r2opt / len(
                     red_data_list[0].red_op.sensi_r))  # number of react to opt per target data (spec / Sl / ...)
                 react2mod = []
                 for idx in range(len(max_sens_list)):
                     react2mod.append([])
                     sensi_r = []
                     for r in range(len(max_sens_list[idx])):
                         sensi_r.append((max_sens_list[idx][r], r + 1))
                     sensi_r.sort()

                     n_r2opt_sp = 0
                     for r1 in range(len(sensi_r)):
                         x = -(r1 + 1)
                         for r2 in range(len(max_sens_list[idx])):
                             if max_sens_list[idx][r2] == sensi_r[x][0]:  # find most sensitive reactions
                                 if not self.mech.react.modif[r2]:
                                     # check if sub mechs can be modified:
                                     n_C = self.mech.react.subm_C[r2]
                                     sub_H = self.mech.react.subm_C[r2] == 0 and not self.mech.react.subm_CO[r2]
                                     n_N = self.mech.react.subm_N[r2]
                                     n_S = self.mech.react.subm_S[r2]
                                     n_Si = self.mech.react.subm_Si[r2]

                                     if n_C > 0 or sub_H:
                                         sub_C = optim_param.opt_subm_C[n_C]
                                     else:
                                         sub_C = False
                                     if self.mech.react.subm_CO[r2]:
                                         sub_CO = optim_param.opt_subm_CO
                                     else:
                                         sub_CO = False
                                     if n_N > 0:
                                         sub_N = optim_param.opt_subm_N[n_N]
                                     else:
                                         sub_N = False
                                     if n_S > 0:
                                         sub_S = optim_param.opt_subm_N[n_S]
                                     else:
                                         sub_S = False
                                     if n_Si > 0:
                                         sub_Si = optim_param.opt_subm_Si[n_Si]
                                     else:
                                         sub_Si = False
                                     # if sub mechs can be modified:
                                     if sub_C or sub_CO or sub_N or sub_S or sub_Si:
                                         self.mech.react.modif[r2] = True;
                                         n_r2opt_sp += 1
                                         react2mod[-1].append(r2)
                             if n_r2opt_sp == n_r2opt_sp_max: break
                         if n_r2opt_sp == n_r2opt_sp_max: break

             elif 'DRG' in red_data_list[0].reduction_operator \
                     or optim_param.optim_on_meth != 'False':
                 if n_tspc != 0: n_r2opt_sp_max = round(
                     n_r2opt / n_tspc)  # number of react to opt per target data (spec / Sl / ...)
                 max_coeffs_list = np.zeros((n_tspc, len(red_data_list[0].red_op.r_interaction_coeffs[0])))
                 for l in range(len(red_data_list)):
                     for idx in range(n_tspc):
                         for r in range(len(red_data_list[l].red_op.r_interaction_coeffs[idx])):
                             if abs(red_data_list[l].red_op.r_interaction_coeffs[idx][r]) > max_coeffs_list[idx][r]:
                                 max_coeffs_list[idx][r] = abs(red_data_list[l].red_op.r_interaction_coeffs[idx][r])

                 react2mod = []
                 for idx in range(n_tspc):
                     react2mod.append([])
                     drg_coeffs = []
                     for j in range(len(max_coeffs_list[idx])):
                         drg_coeffs.append((max_coeffs_list[idx][j], j + 1))
                     drg_coeffs.sort()

                     # DRG_sorted.sort()
                     n_r2opt_sp = 0
                     for r1 in range(len(drg_coeffs)):
                         x = -(r1 + 1)
                         for r2 in range(len(max_coeffs_list[idx])):
                             if max_coeffs_list[idx][r2] == drg_coeffs[x][0]:
                                 if not self.mech.react.modif[r2]:
                                     # check if sub mechs can be modified:
                                     n_C = self.mech.react.subm_C[r2]
                                     sub_H = self.mech.react.subm_C[r2] == 0 and not self.mech.react.subm_CO[r2]
                                     n_N = self.mech.react.subm_N[r2]
                                     n_S = self.mech.react.subm_S[r2]
                                     n_Si = self.mech.react.subm_Si[r2]

                                     if n_C > 0 or sub_H:
                                         sub_C = optim_param.opt_subm_C[n_C]
                                     else:
                                         sub_C = False
                                     if self.mech.react.subm_CO[r2]:
                                         sub_CO = optim_param.opt_subm_CO
                                     else:
                                         sub_CO = False
                                     if n_N > 0:
                                         sub_N = optim_param.opt_subm_N[n_N]
                                     else:
                                         sub_N = False
                                     if n_S > 0:
                                         sub_S = optim_param.opt_subm_N[n_S]
                                     else:
                                         sub_S = False
                                     if n_Si > 0:
                                         sub_Si = optim_param.opt_subm_Si[n_Si]
                                     else:
                                         sub_Si = False
                                     # if sub mechs can be modified:
                                     if sub_C or sub_CO or sub_N or sub_S or sub_Si:
                                         self.mech.react.modif[r2] = True;
                                         n_r2opt_sp += 1
                                         react2mod[-1].append(r2)

                             if n_r2opt_sp == n_r2opt_sp_max: break
                         if n_r2opt_sp == n_r2opt_sp_max: break

             if optim_param.display_react2opt:
                 n_r_opt = 0
                 for _r in range(len(self.mech.react.modif)):
                     if self.mech.react.modif[_r] == True: n_r_opt += 1

                 print_('-------------------------------- \n' + str(n_r_opt) + ' reactions to optimize:', mp)
                 for _tn in range(len(red_data_list[0].tspc)):
                     print_(' * for target: ' + red_data_list[0].tspc[_tn] + ':', mp)
                     for _r in react2mod[_tn]:
                         if int(self.mech.react.number[_r]) < 10:
                             spaces = '    -   '
                         elif int(self.mech.react.number[_r]) < 100:
                             spaces = '   -   '
                         elif int(self.mech.react.number[_r]) < 1000:
                             spaces = '  -   '
                         elif int(self.mech.react.number[_r]) < 10000:
                             spaces = ' -   '
                         print_(str(self.mech.react.number[_r]) + spaces + self.mech.react.equation[_r], mp)

                 # ------------ 10/07/2023
                 if 'idx_Sl' in locals():
                     print_(' * for target Flame speed' + ':', mp)
                     for _r in react2mod[idx_Sl]:
                         if int(self.mech.react.number[_r]) < 10:
                             spaces = '    -   '
                         elif int(self.mech.react.number[_r]) < 100:
                             spaces = '   -   '
                         elif int(self.mech.react.number[_r]) < 1000:
                             spaces = '  -   '
                         elif int(self.mech.react.number[_r]) < 10000:
                             spaces = ' -   '
                         print_(str(self.mech.react.number[_r]) + spaces + self.mech.react.equation[_r], mp)

                 if 'idx_T' in locals():
                     print_(' * for target Temperature' + ':', mp)
                     for _r in react2mod[idx_T]:
                         if int(self.mech.react.number[_r]) < 10:
                             spaces = '    -   '
                         elif int(self.mech.react.number[_r]) < 100:
                             spaces = '   -   '
                         elif int(self.mech.react.number[_r]) < 1000:
                             spaces = '  -   '
                         elif int(self.mech.react.number[_r]) < 10000:
                             spaces = ' -   '
                         print_(str(self.mech.react.number[_r]) + spaces + self.mech.react.equation[_r], mp)

                 print_('--------------------------------\n\n', mp)
                 optim_param.display_react2opt = False


         else:  # selection of submech in gui

             if False in optim_param.opt_subm_C \
                     or False in optim_param.opt_subm_N \
                     or False in optim_param.opt_subm_S \
                     or False in optim_param.opt_subm_Si \
                     or optim_param.opt_subm_CO == False:  # opt_subm_C : selection of submech in gui
                 self.mech.react.modif = [True] * len(self.mech.react.modif)
                 # CxHyOz sub-mechanisms (H2 submech is considered as a C0 submech)
                 for n_C in range(len(optim_param.opt_subm_C)):
                     if optim_param.opt_subm_C[n_C] == False:
                         for r in range(len(self.mech.react.equation)):
                             if self.mech.react.subm_C[r] == n_C \
                                     and not self.mech.react.subm_CO[r] \
                                     and self.mech.react.subm_N[r] == 0 \
                                     and self.mech.react.subm_S[r] == 0 \
                                     and self.mech.react.subm_Si[r] == 0:
                                 self.mech.react.modif[r] = False
                 # CO sub-mechanism
                 if optim_param.opt_subm_CO == False:
                     for r in range(len(self.mech.react.equation)):
                         if self.mech.react.subm_CO[r] == True:
                             self.mech.react.modif[r] = False
                 # N sub-mechanisms
                 for n_N_ in range(len(optim_param.opt_subm_N) - 1):
                     n_N = n_N_ + 1
                     if optim_param.opt_subm_N[n_N] == False:
                         for r in range(len(self.mech.react.equation)):
                             if self.mech.react.subm_N[r] == n_N:
                                 self.mech.react.modif[r] = False
                 # S sub-mechanism
                 if optim_param.opt_subm_S == False:
                     for r in range(len(self.mech.react.equation)):
                         if self.mech.react.subm_S[r] == True:
                             self.mech.react.modif[r] = False
                 # Si sub-mechanism
                 if optim_param.opt_subm_Si == False:
                     for r in range(len(self.mech.react.equation)):
                         if self.mech.react.subm_Si[r] == True:
                             self.mech.react.modif[r] = False
             else:
                 self.mech.react.modif = [True] * len(self.mech.react.modif)

         if type(optim_param.reactions2opt) is not bool:
             if False not in self.mech.react.modif:  # if no restriction is defined, prevent the modification of not specified reactions
                 self.mech.react.modif = [False] * len(self.mech.react.modif)
             # allow the modification of specified reactions
             for r2mod in optim_param.reactions2opt:
                 for r in range(len(self.mech.react.modif)):
                     if r2mod == r + 1:
                         self.mech.react.modif[r] = True

     def shift_flame_data(self, conditions_list, optim_param, ref_results_list, opt):
         verbose = conditions_list[0].simul_param.verbose

         os.chdir(conditions_list[0].main_path + '/' + opt)

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

         qoi_tot = [];
         qoi_tot_pond = [];
         pond = 0

         for i in range(len(conditions_list)):
             conditions = conditions_list[i]
             if conditions.config == 'free_flame':

                 T_check = conditions.error_param.T_check

                 # Simulation conditions
                 ref_results = ref_results_list[i]

                 Opt_results = comp.red_computation(ref_results.conditions, \
                                                    gas, self.mech.spec.activ_m, self.mech.react.activ_m)

                 end_sim = ref_results.conditions.simul_param.end_sim
                 shifting = -end_sim
                 original_pts_scatter = np.array(ref_results.pts_scatter)
                 pts_scatter = np.array(ref_results.pts_scatter)
                 fitness = 0
                 for shif_it in range(100):
                     shifting += end_sim / 100
                     ref_results.pts_scatter = pts_scatter + shifting
                     conditions.simul_param.pts_scatter = pts_scatter + shifting
                     errors = cdef.Errors(conditions, ref_results, Opt_results, \
                                          optim_param)
                     for sp in range(optim_param.n_tspc):
                         qoi_tot.append(errors.qoi_s[sp])
                         qoi_tot_pond.append(errors.qoi_s[sp] * optim_param.coeff_s[sp])
                         pond += optim_param.coeff_s[sp]
                     if T_check:
                         qoi_tot.append(errors.qoi_T)
                         qoi_tot_pond.append(errors.qoi_T * optim_param.coeff_T)
                         pond += optim_param.coeff_T
                     if conditions.error_param.error_type_fit == 'mean':
                         fitness_i = 1 / (np.sum(qoi_tot) / pond)
                     elif conditions.error_param.error_type_fit == 'max':
                         fitness_i = 1 / np.max(qoi_tot)
                     if fitness_i > fitness:
                         best_shift = shifting
                 ref_results.pts_scatter = original_pts_scatter + best_shift
                 conditions.simul_param.pts_scatter = original_pts_scatter + best_shift
                 ref_results.conditions.simul_param.shift = best_shift

         return conditions_list, ref_results_list

     def time_step_optim(self, conditions_list, ref_results_list, opt):
         mp = conditions_list[0].main_path
         verbose = conditions_list[0].simul_param.verbose

         print_('time step optimization', mp)
         os.chdir(conditions_list[0].main_path + '/' + opt)
         if '.cti' in conditions_list[0].mech:
             filename = 'temp.cti'
             self.mech.write_new_mech(filename)
         else:
             filename = 'temp.yaml'
             self.mech.write_yaml_mech(filename)

         # --------------------------------------------------------------------------------
         # interpretation of the new mech

         for i in range(len(conditions_list)):
             conditions = conditions_list[i]
             conditions.composition.gas = cdef.get_gas_ct(filename)
             if 'reactor' in conditions.config:
                 opt_results, conditions = comp.ref_computation(conditions)
                 ref_results = comp.red_computation(conditions,
                                                    conditions.composition.gas_ref,
                                                    self.mech.spec.activ_m, self.mech.react.activ_m)
                 ref_results_list[i] = ref_results
                 conditions_list[i] = conditions
         # --------------------------------------------------------------------------------

         return conditions_list, ref_results_list

     def randomize_kin(self, optim_param):

         for r in range(len(self.mech.react.type)):
             try_r = 0;
             valid = False;
             damping = 0.75;
             max_try = 10
             while not valid and try_r < max_try:
                 var_range = 1 - (try_r / max_try)
                 if self.mech.react.modif[r]:
                     uncert_r = [u / 100 for u in self.mech.react.uncert[r]]
                     if self.mech.react.type[r] == "three_body_reaction" \
                             or self.mech.react.type[r] == "reaction":
                         for k in range(len(self.mech.react.kin[r])):
                             self.mech.react.kin[r][k] = self.mech.react.ref_kin[r][k] \
                                                         + self.mech.react.ref_kin[r][k] * random.uniform(-var_range,
                                                                                                          var_range) * \
                                                         uncert_r[k] * (damping ** try_r)

                     elif self.mech.react.type[r] == "falloff_reaction" \
                             or self.mech.react.type[r] == "pdep_arrhenius" \
                             or self.mech.react.type[r] == "chemically_activated_reaction" \
                             or self.mech.react.type[r] == "chebyshev":
                         for k1 in range(len(self.mech.react.kin[r])):
                             for k2 in range(len(self.mech.react.kin[r][k1])):
                                 self.mech.react.kin[r][k1][k2] = self.mech.react.ref_kin[r][k1][k2] \
                                                                  + self.mech.react.ref_kin[r][k1][k2] * random.uniform(
                                     -var_range, var_range) * uncert_r[k2] * (damping ** try_r)

                     valid = check_k(self.mech.react, r)
                 else:
                     valid = True

                 try_r += 1

             if not valid:
                 self.mech.react.kin[r] = copy.deepcopy(self.mech.react.ref_kin[r])

     def fitness_eval(self, conditions_list, optim_param, ref_results_list, opt, n_par=0, ref_ind=False):

         verbose    = conditions_list[0].simul_param.verbose
         mp         = conditions_list[0].main_path
         os.chdir(conditions_list[0].main_path + '/' + opt)

         filename = 'temp_' + str(n_par)
         # if self.mech.keep4opt == True:
         #     filename = 'keep_' + filename
         if '.cti' in conditions_list[0].mech:
             filename += '.cti'
             self.mech.write_new_mech(filename)
         else:
             filename += '.yaml'
             self.mech.write_yaml_mech(filename)

         gas = cdef.get_gas_ct(filename)

         txt_f          = '\n'
         qoi_tot        = []
         pond           = 0
         qoi_tot_mean   = 0
         for i in range(len(conditions_list)):

             # -------------------------------
             # Reduction loop

             T_check = conditions_list[i].error_param.T_check
             Sl_check = conditions_list[i].error_param.Sl_check
             ig_check = conditions_list[i].error_param.ig_check
             K_check = conditions_list[i].error_param.K_check

             # Simulation conditions
             conditions = conditions_list[i]
             conditions.simul_param.par_ind = str(n_par)
             ref_results = ref_results_list[i]

             cur_path = os.getcwd()
             Opt_results = comp.red_computation(conditions, gas, \
                                                self.mech.spec.activ_m, self.mech.react.activ_m)
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

             errors = cdef.Errors(conditions, ref_results, Opt_results, optim_param)

             # Condition fitness weighting
             if optim_param.coeff_cond:
                 coeff_cond = optim_param.coeff_cond[i]
             else:
                 coeff_cond = 1

             qoi_case, pond_case = 0, 0
             cc = conditions.config
             for sp in range(optim_param.n_tspc):
                 if errors.qoi_s[sp]:
                     pond_i = optim_param.coeff_s[sp] * coeff_cond
                     txt_f += 'Fit ' + cc + '  sp - ' + optim_param.tspc[sp] + ' : ' + "%.3f" % (
                                 1 - errors.qoi_s[sp]) + '  pond = ' + "%.1f" % pond_i + '\n'
                     qoi_tot_mean += (1 - errors.qoi_s[sp]) * pond_i
                     qoi_case += (1 - errors.qoi_s[sp]) * pond_i
                     qoi_tot.append((1 - errors.qoi_s[sp]) * min(1, np.ceil(pond_i)))
                     pond += pond_i;
                     pond_case += pond_i
             if 'JSR' not in conditions.config and T_check and errors.qoi_T:
                 pond_i = optim_param.coeff_T * coeff_cond
                 txt_f += 'Fit ' + cc + ' - T: ' + "%.3f" % (1 - errors.qoi_T) + '  pond = ' + "%.1f" % pond_i + '\n'
                 qoi_tot_mean += (1 - errors.qoi_T) * pond_i
                 qoi_case += (1 - errors.qoi_T) * pond_i
                 qoi_tot.append((1 - errors.qoi_T) * min(1, np.ceil(pond_i)))
                 pond += pond_i;
                 pond_case += pond_i
             if 'reactor' in conditions.config and ig_check and errors.qoi_ig:
                 pond_i = optim_param.coeff_ig * coeff_cond
                 txt_f += 'Fit reactor - igt: ' + "%.3f" % (1 - errors.qoi_ig) + '  pond = ' + "%.1f" % pond_i + '\n'
                 qoi_tot_mean += (1 - errors.qoi_ig) * pond_i
                 qoi_case += (1 - errors.qoi_ig) * pond_i
                 qoi_tot.append((1 - errors.qoi_ig) * min(1, np.ceil(pond_i)))
                 pond += pond_i;
                 pond_case += pond_i
             if 'free_flame' in conditions.config and Sl_check and errors.qoi_Sl:
                 pond_i = optim_param.coeff_Sl * coeff_cond
                 txt_f += 'Fit free_flame - Sl: ' + "%.3f" % (1 - errors.qoi_Sl) + '  pond = ' + "%.1f" % pond_i + '\n'
                 qoi_tot_mean += (1 - errors.qoi_Sl) * pond_i
                 qoi_case += (1 - errors.qoi_Sl) * pond_i
                 qoi_tot.append((1 - errors.qoi_Sl) * min(1, np.ceil(pond_i)))
                 pond += pond_i;
                 pond_case += pond_i
             if ('diff_flame' in conditions.config or 'pp_flame' in conditions.config) and K_check and errors.qoi_K:
                 pond_i = optim_param.coeff_K * coeff_cond
                 txt_f += 'Fit diff_flame - K: ' + "%.3f" % (1 - errors.qoi_K) + '  pond = ' + "%.1f" % pond_i + '\n'
                 qoi_tot_mean += (1 - errors.qoi_K) * pond_i
                 qoi_case += (1 - errors.qoi_K) * pond_i
                 qoi_tot.append((1 - errors.qoi_K) * min(1, np.ceil(pond_i)))
                 pond += pond_i;
                 pond_case += pond_i

             # detailed informations on initial fitness (for each case)
         #            print_detailed_fit = True
         #            if print_detailed_fit:
         #                fit_i =  1/(qoi_case/pond_case)
         #                txt_fit = 'Fitness condition ' + str(i+1) + ': ' + str(fit_i) \
         #                          + '   pond: ' + str(pond_case)
         #                print_(txt_fit, mp)
         if verbose >= 7:
             print_(txt_f, mp)

         if 'no data' in qoi_tot:  qoi_tot.remove('no data')
         if conditions.error_param.error_type_fit == 'mean':
             #             self.fitness = 1/(np.sum(qoi_tot)/pond)
             # fitness = 1/(qoi_tot_mean/pond)
             fitness = qoi_tot_mean / pond
             if verbose >= 5:
                 print_('Fitness of the individual= ' + '%.3f' % fitness, mp)
         elif conditions.error_param.error_type_fit == 'max':
             fitness = 1 - np.max(qoi_tot)

         # check nan
         if fitness != fitness:
             fitness = 0

         return max(fitness, 0)

     def export_data(self, conditions_list, optim_param, ref_results_list, opt, filename='temp.cti'):
         verbose = conditions_list[0].simul_param.verbose

         errors_list = [];
         Opt_results_list = []

         qoi_tot = [];
         pond = 0

         os.chdir(conditions_list[0].main_path + '/' + opt)
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
             conditions = conditions_list[i]
             ref_results = ref_results_list[i]
             os.chdir(conditions_list[0].main_path + '/' + opt)
             if '.cti' in filename:
                 self.mech.write_new_mech(filename)
             else:
                 self.mech.write_yaml_mech(filename)

             T_check = conditions_list[i].error_param.T_check
             Sl_check = conditions_list[i].error_param.Sl_check
             ig_check = conditions_list[i].error_param.ig_check
             K_check = conditions_list[i].error_param.K_check

             Opt_results_list.append(comp.red_computation(conditions,
                    gas, self.mech.spec.activ_m, self.mech.react.activ_m))
             os.chdir(conditions_list[0].main_path + '/' + opt)
             errors_list.append(cdef.Errors(conditions, ref_results,
                    Opt_results_list[-1], optim_param))
             for sp in range(optim_param.n_tspc):
                 if errors_list[-1].qoi_s[sp]:
                     qoi_tot.append(errors_list[-1].qoi_s[sp])
                     pond += optim_param.coeff_s[sp]
             if 'JSR' not in conditions.config and T_check and errors_list[-1].qoi_T:
                 qoi_tot.append(errors_list[-1].qoi_T)
                 pond += optim_param.coeff_T
             if 'reactor' in conditions.config and ig_check and errors_list[-1].qoi_ig:
                 qoi_tot.append(errors_list[-1].qoi_ig)
                 pond += optim_param.coeff_ig
             if 'free_flame' in conditions.config and Sl_check and errors_list[-1].qoi_Sl:
                 qoi_tot.append(errors_list[-1].qoi_Sl)
                 pond += optim_param.coeff_Sl
             if ('diff_flame' in conditions.config or 'pp_flame' in conditions.config) and K_check and errors_list[
                 -1].qoi_K:
                 qoi_tot.append(errors_list[-1].qoi_K)
                 pond += optim_param.coeff_K

         if 'no data' in qoi_tot:  qoi_tot.remove('no data')
         if conditions.error_param.error_type_fit == 'mean':
             self.fitness = 1 / np.mean(qoi_tot)
         elif conditions.error_param.error_type_fit == 'max':
             self.fitness = 1 / np.max(qoi_tot)

         return Opt_results_list, errors_list, self.fitness



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

