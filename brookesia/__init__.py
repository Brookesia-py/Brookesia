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
#    Version
#==============================================================================

__version__ = '1.5.4'



#==============================================================================
#    Packages
#==============================================================================

import os
#import sys
#from shutil import copyfile
#import cantera as ct

import brookesia.GeneralFunctions as genf
import brookesia.Class_def as cdef
#from   brookesia.Class_def import print_
import brookesia.GUI_red_mech as bkgui
import brookesia.write_conditions_input as wc
import brookesia.write_kinetic_mech     as wk
import brookesia.write_results_input    as wr

import pandas as pd

root_path   = __file__.split('__init__.py')[0]
python_path = root_path.split('/lib/')[0] + '/bin/python'
#import shutil
#path2python = shutil.which("python")

# set language parameters (for xml cantera files)
#import locale as lc
#lc.setlocale(lc.LC_ALL, 'en_GB.utf8')
#==============================================================================

def run_reduction(filename):
#    global WD_path
    genf.main_redopt_algo(filename,WD_path)
    os.chdir(WD_path)

def convert(mech):
    mech_data = cdef.Mech_data(mech)
    mech_data.write_chemkin_mech(mech)

def gui():
#    global WD_path
#    global WD_name
    bkgui.display_gui(WD_path, WD_name, root_path)
    os.chdir(WD_path)

def create_wd():
    WD_path, WD_name = bkgui.get_wd_folder()
    os.chdir(root_path)
    fd = open('wd_activ.txt', 'w')
    fd.write(WD_name + ';' + WD_path)
    os.chdir(WD_path)
    wc.write_conditions()
    os.chdir(WD_path)
    wk.write_kinetic()
    os.chdir(WD_path)
    wr.write_results()
    return WD_path, WD_name

def check_wd():
    os.chdir(root_path)

    ready4select = False
    if 'working_dir.txt' in os.listdir():
        workdir_table = pd.read_table("working_dir.txt",sep =';',header = 0, index_col=[0])
        if len(workdir_table['directory'])>=1:
            ready4select = True

    if ready4select: # is some working dir already exist
        if len(workdir_table['directory'])==1:
            WD_name = workdir_table['working_dir_name'].iloc[0]
            WD_path = workdir_table['directory'].iloc[0]
        else:
            select_wd()
            os.chdir(root_path)
            fs = open('wd_activ.txt', 'r')
            txt = fs.readline()
            WD_name = txt.split(';')[0]
            WD_path = txt.split(';')[1]
    else:
        WD_path, WD_name = create_wd()
        data = {'working_dir_name':[WD_name],
                'directory': [WD_path]}
        workdir_table = pd.DataFrame(data=data)
        os.chdir(root_path)
        workdir_table.to_csv('working_dir.txt', sep=';')

    os.chdir(WD_path)

    return WD_path, WD_name


def select_wd():
    os.chdir(root_path)

    workdir_table = pd.read_table("working_dir.txt",sep =';',header = 0, index_col=[0])

    ready4select = False
    if 'working_dir.txt' in os.listdir():
        workdir_table = pd.read_table("working_dir.txt",sep =';',header = 0, index_col=[0])
        if len(workdir_table['directory'])>=1:
            ready4select = True

    if ready4select:
        os.system(python_path + ' ' + root_path + "wd_select_gui.py " + root_path)
        fs = open('wd_activ.txt', 'r')
        txt = fs.readline()
        WD_name = txt.split(';')[0]
        WD_path = txt.split(';')[1]

    else:
        WD_path, WD_name = create_wd()
        data = {'working_dir_name':[WD_name],
                'directory': [WD_path]}
        workdir_table = pd.DataFrame(data=data)
        os.chdir(root_path)
        workdir_table.to_csv('working_dir.txt', sep=';')

    os.chdir(WD_path)



WD_path, WD_name = check_wd()

#gui()
#run_reduction('1_reactor.inp')
#run_reduction('2_JSR.inp')
#run_reduction('3_Freeflame.inp')
#run_reduction('4_diff.inp')
#run_reduction('5_pp.inp')
#run_reduction('6_import_Sl.inp')
#run_reduction('0_deg_GRI.inp')

