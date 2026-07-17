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
import os, copy




def find_max_concentration(conditions,red_data,mech_data,red_results):
    """r_AB =sum(nu_iA rate_i delta_Bi) /max(prod A, conso A) """

    mp = red_results.conditions.main_path

    gas_ref = conditions.composition.gas_ref
    gas_red = red_results.gas #; print_('gas_red_ok',mp)
    ns      = gas_ref.n_species
    
    # pts_scatter = red_results.pts_scatter
    # n_points = len(pts_scatter)

    conc_max_red = np.max(red_results.conc, axis=0)
    conc_max     = np.zeros(ns)
    
    # fill conc_max vector with the vector dimension of the ref. mechanism
    for idx_red in range(gas_red.n_species):
        idx_ref = gas_ref.species_index(gas_red.species_name(idx_red)) 
        conc_max[idx_ref] = conc_max_red[idx_red]
        

    return conc_max




def speciesWithdrawal(conditions, red_data, red_method, mech_data, conc_max, step3 = False, ranking=True):

    mp = conditions.main_path

    # main variables
    gas_ref  = red_data.gas_ref

    tsp_idx = red_data.targetSpeciesIdx
    verbose=conditions.simul_param.verbose

    ns = gas_ref.n_species
    # active_sp_p = list(mech_data.spec.activ_p)
    active_species = copy.deepcopy(mech_data.spec.activ_pm)

    # Add fuel / ox /diluent to the conserved species
    init_spec = conditions.composition.X
    init_spec = init_spec.split(',')
    idx_mixt  = []
    for spec in init_spec:
        idx_mixt.append(gas_ref.species_index(spec.split(":")[0].replace(' ','')))
    if conditions.composition.X2:
        init_spec = conditions.composition.X2
        init_spec = init_spec.split(',')
        for spec in init_spec:
            idx_mixt.append(gas_ref.species_index(spec.split(":")[0].replace(' ','')))

    if step3:
        if 'DRG' in red_method:
            red_coeff    = red_data.red_op.interaction_coeffs # [point][target, sp]
            # get, for each sp, the max of interaction across all targets and inspected points.
            red_coeff = np.array(red_coeff)
            red_coeff_sp = red_coeff.max(axis=(0, 1))
        if 'SA' in red_method:
            red_coeff    = red_data.red_op.sensi_sp # [target, sp]
            # get, for each sp, the max of interaction across all targets
            red_coeff_sp = red_coeff.max(axis=(0)) 
    else:
        red_coeff_sp = np.ones(ns)

    # consider the rank:
    if ranking: 
        conc_max     = np.argsort(np.argsort(conc_max))     + 1
        red_coeff_sp = np.argsort(np.argsort(red_coeff_sp)) + 1
    
    prod_c_rc = conc_max*red_coeff_sp
    
    
    for _sp in range(ns):
        idx = list(prod_c_rc).index(np.min(prod_c_rc))
        if idx not in idx_mixt:
            if  active_species[idx] == True               \
            and mech_data.spec.activ_p[idx] == False \
            and idx not in tsp_idx:
                active_species[idx] = False
                break
        prod_c_rc[idx] = np.max(prod_c_rc)


    return active_species










def species_reactions_withdrawal(conditions,red_data,mech_data,conc_max,eps_TRO):

    # main variables
    gas = red_data.gas_ref
    mp  = conditions.main_path

    tsp_idx = red_data.targetSpeciesIdx
    verbose=conditions.simul_param.verbose


    ns = gas.n_species
    nr = gas.n_reactions

    nu_f = gas.reactant_stoich_coeffs()
    nu_r = gas.product_stoich_coeffs()


#    active_species = list([True]*ns)
    active_species = list(mech_data.spec.activ_p)


    conc_max_tab = np.array(conc_max)
    # max vector of each species for every conditions
    conc_max_4all_conditions = np.max(conc_max_tab,0)


    # find threshold concentration for withdrawal
    conc_max_sort = np.sort(conc_max_4all_conditions)
    while 0 in conc_max_sort:
        conc_max_sort = np.delete(conc_max_sort,0)

    for i in range(eps_TRO-1):
        conc_max_sort = np.delete(conc_max_sort,0)
    threshold_conc = conc_max_sort[0]

    sp_removed = []
    for sp in range(ns):
#        if not mech_data.spec.activ_m[sp]:
#            active_species[sp] = False
#        else:
        if conc_max_4all_conditions[sp] >= threshold_conc:
            active_species[sp] = True
        elif mech_data.spec.activ_m[sp]:
            sp_removed.append(gas.species_name(sp))



#    if verbose >= 4 and len(sp_removed)>1:
    if verbose >= 4:
        nb_sp_rem   = len(sp_removed)
        sp_removed  = str(sp_removed)
        sp_removed  = sp_removed.replace('[','')
        sp_removed  = sp_removed.replace(']','')
        print_(str(nb_sp_rem)+' species removed: '+sp_removed,mp)


    # ----- Make sure to keep species in the initial mixture and target species -----

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
    # Species withdrawal
    for t in range(len(tsp_idx)):
        spA=tsp_idx[t]
        if not active_species[spA]: active_species[spA]=True


   # Reactions withdrawal
    active_reactions = [True]*nr
    for r in range(nr):
        if not mech_data.react.activ_m[r]:
            active_reactions[r] = False
        else:
            for sp in range(ns):
                # remove reactions involving non active species
                if not active_species[sp]:
                    if nu_f[sp][r]!=0 or nu_r[sp][r]!=0:
                        active_reactions[r] = False
                        break

    # check threebody exception (+AR) (+HE) etc.
    for r in range(nr):
        if mech_data.react.type[r]=="falloff_reaction" \
        and type(mech_data.react.tbe[r]) is str:
            for sp in range(len(mech_data.spec.name)):
                if mech_data.react.tbe[r]==mech_data.spec.name[sp]\
                and not active_species[sp]:
                    active_reactions[r]=False



    return active_reactions, active_species

