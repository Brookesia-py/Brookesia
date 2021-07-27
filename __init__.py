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
global version

version = '1.7.1'
__version__ = version

#==============================================================================
#    Packages
#==============================================================================

import os
import sys
#from shutil import copyfile
#import cantera as ct

import brookesia.GeneralFunctions as genf
import brookesia.Class_def as cdef
#from   brookesia.Class_def import print_
import brookesia.GUI_red_mech as bkgui
import brookesia.write_conditions_input as wc
import brookesia.write_kinetic_mech     as wk
import brookesia.write_results_input    as wr

import gc
import pandas as pd


# path
root_path   = __file__.split('__init__.py')[0]
python_path = root_path.split('/lib/')[0] + '/bin/python'
if os.name != 'nt': # for Linux or Mac
    pers_config_path = '~/.Brookesia'
    try:
        os.chdir(os.path.expanduser(pers_config_path))
    except:
        os.chdir(os.path.expanduser('~'))
        try:
            os.mkdir('.Brookesia')
        except:
            a=False
        os.chdir('.Brookesia')
        pers_config_path = os.getcwd()
else: # for windows
    pers_config_path = root_path



#import shutil
#path2python = shutil.which("python")

# set language parameters (for xml cantera files)
#import locale as lc
#lc.setlocale(lc.LC_ALL, 'en_GB.utf8')
#==============================================================================

def run_reduction(filename):
    global version
#    global WD_path
    gc.collect()
    genf.main_redopt_algo(filename,WD_path,version)
    os.chdir(WD_path)

def convert4chemkin(mech):
    global version
    mech_data = cdef.Mech_data(mech)
    mech_data.write_chemkin_mech(mech,version)

def gui():
#    global WD_path
#    global WD_name
    os.chdir(WD_path)
    bkgui.display_gui(WD_path, WD_name, root_path)
    os.chdir(WD_path)

def create_wd():
    wd_path_ok = False
    while not wd_path_ok:
        WD_path, WD_name = bkgui.get_wd_folder()
        try:    os.chdir(pers_config_path)
        except: os.chdir(os.path.expanduser(pers_config_path))
        fd = open('wd_activ.txt', 'w')
        fd.write(WD_name + ';' + WD_path)
        try:
            os.chdir(WD_path)
            wd_path_ok = True
        except:
            print('\n\n ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ')
            print('Brookesia cannot access to the working directory:' + WD_path)
            print('Please, select a new working directory')
            print(' ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! \n\n')


    os.chdir(WD_path)
    wc.write_conditions()
    os.chdir(WD_path)
    wk.write_kinetic()
    os.chdir(WD_path)
    wr.write_results()

    return WD_path, WD_name

def check_wd():

    try:    os.chdir(pers_config_path)
    except: os.chdir(os.path.expanduser(pers_config_path))

    ready4select = False
    create_wd_file = False

    if 'working_dir.txt' in os.listdir():
        workdir_table = pd.read_table("working_dir.txt",sep =';',header = 0, index_col=[0])
        if len(workdir_table['directory'])>=1:
            ready4select = True
    else:
        create_wd_file = True


    if ready4select: # is some working dir already exist
        if len(workdir_table['directory'])==1:
            WD_name = workdir_table['working_dir_name'].iloc[0]
            WD_path = workdir_table['directory'].iloc[0]
            try:
                os.chdir(WD_path)
            except:
                print('\n\n ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ')
                print('Brookesia cannot access to the working directory:' + WD_path)
                print('Please, select a new working directory')
                print(' ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! \n\n')

                create_wd_file = True
        else:
            wd_path_ok = False
            while not wd_path_ok:
                select_wd()
                try:    os.chdir(pers_config_path)
                except: os.chdir(os.path.expanduser(pers_config_path))
                fs = open('wd_activ.txt', 'r')
                txt = fs.readline()
                WD_name = txt.split(';')[0]
                WD_path = txt.split(';')[1]
                try:
                    os.chdir(WD_path)
                    wd_path_ok = True
                except:
                    print('\n\n ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ')
                    print('Brookesia cannot access to the working directory:' + WD_path)
                    print('Please, change the path of this working directory')
                    print('working directory name: ' + WD_name + ')')
                    print(' ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! \n\n')

    if create_wd_file:
        WD_path, WD_name = create_wd()
        data = {'working_dir_name':[WD_name],
                'directory': [WD_path]}
        workdir_table = pd.DataFrame(data=data)
        try:    os.chdir(pers_config_path)
        except: os.chdir(os.path.expanduser(pers_config_path))
        workdir_table.to_csv('working_dir.txt', sep=';')

    os.chdir(WD_path)

    return WD_path, WD_name


def select_wd():
    try:    os.chdir(pers_config_path)
    except: os.chdir(os.path.expanduser(pers_config_path))

    workdir_table = pd.read_table("working_dir.txt",sep =';',header = 0, index_col=[0])

    ready4select = False
    if 'working_dir.txt' in os.listdir():
        workdir_table = pd.read_table("working_dir.txt",sep =';',header = 0, index_col=[0])
        if len(workdir_table['directory'])>=1:
            ready4select = True

    if ready4select:
        wd_path_ok = False
        while not wd_path_ok:
#            select_wd()
            os.system(python_path + ' ' + root_path + "wd_select_gui.py " + root_path)
            fs = open('wd_activ.txt', 'r')
            txt = fs.readline()
            WD_name = txt.split(';')[0]
            WD_path = txt.split(';')[1]
            try:
                os.chdir(WD_path)
                wd_path_ok = True
            except:
                print('\n\n ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ')
                print('Brookesia cannot access to the working directory:' + WD_path)
                print('Please, change the path of this working directory')
                print('working directory name: ' + WD_name + ')')
                print(' ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! \n\n')

    else:
        WD_path, WD_name = create_wd()
        data = {'working_dir_name':[WD_name],
                'directory': [WD_path]}
        workdir_table = pd.DataFrame(data=data)
#        os.chdir(root_path)
        workdir_table.to_csv('working_dir.txt', sep=';')

    os.chdir(WD_path)

def plot_results_____sp_fileN_displaylegend(species = 'no_arg', fileName = ['0_reduction_results.csv'], display_legend = True):
    c_path = os.getcwd()

    if species == 'no_arg':
        os.chdir(root_path)
        print('\n\n -------- Choix des paramètres -------- \n\n')
        os.system('kwrite _toolbox_Plot_species.py')

        if 'brookesia._toolbox_Plot_species' in sys.modules:
            del sys.modules["brookesia._toolbox_Plot_species"]
        import brookesia._toolbox_Plot_species as plt
        plt.plot_data_v7(c_path)
    else:
        import brookesia._toolbox_Plot_species as plt
        plt.plot_data_v7(c_path, species, fileName, display_legend)
    os.chdir(c_path)


def plot_fitness(see_file = False):
    c_path = os.getcwd()

    if see_file:
        os.chdir(root_path)
        print('\n\n -------- Choix des paramètres -------- \n\n')
        os.system('kwrite _toolbox_Plot_fitness.py')

        if 'brookesia._toolbox_Plot_fitness' in sys.modules:
            del sys.modules["brookesia._toolbox_Plot_fitness"]
        import brookesia._toolbox_Plot_fitness as pltf
        pltf.plot_fitness(c_path)
    else:
        import brookesia._toolbox_Plot_fitness as pltf
        pltf.plot_fitness(c_path)
    os.chdir(c_path)


os.chdir(root_path)
WD_path, WD_name = check_wd()

