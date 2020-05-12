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
#    Packages
#==============================================================================

import os
import sys
from shutil import copyfile
import cantera as ct
import traceback

import __packages.GeneralFunctions as genf
import __packages.Class_def as cdef
from  __packages.Class_def import print_

# set language parameters (for xml cantera files)
#import locale as lc
#lc.setlocale(lc.LC_ALL, 'en_GB.utf8')
#==============================================================================

print('\n'*100)


def run_reduction(filename):

    #==============================================================================
    #%%      Read reduction conditions
    #==============================================================================

    conditions_list, red_data_list, ref_results_list = genf.get_reduction_parameters(filename)


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

    print_('Computed with :\n * Cantera  '+ct.__version__+'\n * Brookesia 1.3\n\n',mp)

    try:
        #==============================================================================
        #%%      Mechanism informations
        #==============================================================================

        verbose = conditions_list[0].simul_param.verbose
        mech_data = cdef.Mech_data(mech,verbose)

        if mech_prev_red:
            mech_data_prev = cdef.Mech_data(mech_prev_red,verbose)
            mech_data.compare_red(mech_data_prev,mp)

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
                results,conditions = genf.computation_reference(conditions_list[i],verbose)
                conditions_list[i] = conditions
                ref_results_list.append(results)

    #    ref_results_list[0].plotData_opt(['CH4'])

        #==============================================================================
        #%%    Reduction procedure
        ##=============================================================================
        os.chdir(mp)

        ns_red,nr_red,red_results_list,red_data_list\
            =genf.reduction(conditions_list,ref_results_list,red_data_list,mech_data)


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


if len(sys.argv)>1:
    if sys.argv[1] == 'convert':
        mech = sys.argv[2]
        mech_data = cdef.Mech_data(mech)
        mech_data.write_chemkin_mech(mech)
    else:
        filename = sys.argv[1]
        run_reduction(filename)

