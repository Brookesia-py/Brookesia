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

# =============================================================================
#  default options
# =============================================================================

d_verbose               = 4
d_show_plots            = True
d_chemkin_f             = True

d_T_check               = True
d_Sl_check              = True
d_ig_check              = True
d_K_check               = False
d_sp_T                  = 'CO2'
d_sp_Sl                 = 'H'
d_sp_K                  = 'H'
d_sp_ig                 = 'CH3'
d_error_calculation     = 'points'
d_error_coupling        = 'mean'

# Main options
d_fuel                  = 'CH4'
d_oxidant               = 'O2'
d_diluent               = 'N2'
d_diluent_ratio         = ['N2/O2', 3.76]
d_P_min                 = '1e5'
d_P_max                 = '1e5'
d_P_incr                = '1e5'
d_phi_min               = '0.5'
d_phi_max               = '1.5'
d_phi_incr              = '0.5'

# options for reactors
d_T_min_r               = '1400'
d_T_max_r               = '1800'
d_T_incr_r              = '200'
d_tol_ts_r              = '1e-06,1e-12'         # rel/abs
d_n_pts_r               = 250
d_delta_npts            = 20
d_t_max_coeff           = 5
d_Scal_ref              = 'H2O'
d_grad_curv_ratio       = 0.5
d_tign_nPoints          = '130'
d_tign_dt               = '1e-9'

# options for JSR
d_T_min_jsr             = '600'
d_T_max_jsr             = '1200'
d_T_incr_jsr            = '10'
d_tol_ts_jsr            = '1e-10,1e-15'         # rel/abs
d_t_max                 = 0.2
d_diluent_ratio_jsr     = '0.99'

# options for PFR
d_T_min_pfr             = '1600'
d_T_max_pfr             = '1600'
d_T_incr_pfr            = '200'
d_tol_ts_pfr            = '1e-10,1e-15'         # rel/abs
d_length_pfr            = 1.5e-5               # m
d_u_0_pfr               = 0.006                # m/s
d_n_pts_pfr             = 2000
d_area_pfr              = 1.e-4                # m**2
d_pfr_discr             = True

# options for flames
d_T_min_ff              = '300'
d_T_max_ff              = '300'
d_T_incr_ff             = '100'
d_ratio_ff              = 2.0
d_slope_ff              = 0.05
d_curve_ff              = 0.05
d_prune_ff              = 0.01
d_n_pts_ff              = 250
d_tol_ts_ff             = '1e-05,1e-08'      # rel/abs
d_tol_ss_ff             = '1e-06,1e-08'      # rel/abs
d_transport_model       = 'Mix'
d_xmax                  = '.02'
d_pts_scatter           = '0.0, 0.03, 0.3, 0.5, 0.7, 1.0'


# options for diff_flames
d_diluent_ratio_diff1   = 0
d_diluent_ratio_diff2   = 79

# options for partially premixed flames
d_diluent_ratio_pp1   = ['N2/O2', 3.76]
d_diluent_ratio_pp2   = ['N2/O2', 3.76]
d_pts_scatter_cf      = '0.0, 0.2, 0.4, 0.6, 0.8, 1.0'

# SA
d_SA_eps                = 0.01
d_SA_deps               = 0.01
d_SA_pt_num             = 10.0
d_SA_tgt_error          = 30.0
d_SA_tol_1              = '1e-5'
d_SA_tol_2              = '1e-8'
d_SA_tol                = True
d_SA_ISI                = True

# DRG
d_DRG_eps               = 0.01
d_DRG_deps              = 0.01
d_DRG_pt_num            = 10.0
d_DRG_tgt_error         = 30.0
d_DRG_ISI               = True

# GA
d_GA_gen                = 20
d_GA_ind                = 20
d_GA_A                  = 5
d_GA_n                  = 5
d_GA_Ea                 = 5
d_GA_meth               = True
d_GA_meth_fract         = 30
d_GA_meth_pts           = 20



d_GA_selection_operator = 'Roulette' # Roulette Rank Geometric_norm Elitism
d_GA_sel_opt            = '0.2'
d_GA_Xover_op_1         = True       # Simple Xover
d_GA_Xover_op_2         = True       # Multiple Xover
d_GA_Xover_op_3         = True       # Arithmetic Xover
d_GA_Xover_op_4         = True       # Heuristic Xover
d_GA_Xover_int_1        = '10'
d_GA_Xover_int_2        = '20'
d_GA_Xover_int_3        = '20'
d_GA_Xover_int_4        = '20'
d_GA_Xover_opt_1        = ''
d_GA_Xover_opt_2        = ''
d_GA_Xover_opt_3        = ''
d_GA_Xover_opt_4        = ''
d_GA_mut_op_1           = True       # Uniform mutation
d_GA_mut_op_2           = True       # Non-uniform mutation
d_GA_mut_op_3           = True       # Boundary mutation
d_GA_mut_int_1          = '30'
d_GA_mut_int_2          = '30'
d_GA_mut_int_3          = '10'
d_GA_mut_opt_1          = ''
d_GA_mut_opt_2          = '3'
d_GA_mut_opt_3          = ''
d_GA_mut_prob           = 30
d_GA_fit                = 'mean'     #  mean / max

# PSO
d_PSO_inert_score       = True
d_PSO_inert_min         = 0.2   # default value with score integration
d_PSO_inert_i           = 1.5   # default value with score integration
d_PSO_inert_end         = 0.4   # default value with score integration
d_PSO_inert_i_ws        = 0.6   # default value without score integration
d_PSO_inert_end_ws      = 0.2   # default value without score integration
d_PSO_cogn_i            = 1.3
d_PSO_cogn_end          = .7
d_PSO_social_i          = 1
d_PSO_social_end        = 2




# =============================================================================
#%%         package import
# =============================================================================

from PyQt5 import QtCore, QtGui, QtWidgets
import cantera as ct
ct.suppress_thermo_warnings()
import copy
import numpy as np
import brookesia.GeneralFunctions as genf
import brookesia.Class_def as cdef
#import brookesia.__init__ as mred
import os
import sys
import time as timer

# set language parameters (for xml cantera files)
import locale as lc
try:
    lc.setlocale(lc.LC_ALL, 'en_US.utf8')
except:
    try:
        lc.setlocale(lc.LC_ALL, 'en_US.UTF-8')
    except:
        a=False


# tab_flag
global condition_tab_removed
condition_tab_removed = False
# manage load cond
global load_prev_cond
load_prev_cond = False

# flag opt_num
global opt_num
opt_num = 0

#
## Screen dimensions
#app = QtWidgets.QApplication(sys.argv)
#
#screen = app.primaryScreen()
#_width = QtWidgets.QDesktopWidget().screenGeometry(-1).width()
##screen_resolution = app.desktop().screenGeometry()
#_width, _height = cdef.get_screen_size()
##width, height = screen_resolution.width(), screen_resolution.height()
##size = screen.size()
#sz_w = _width/1366
#sz_h = _height/768
#sz_w=1 ; sz_h=1
#del _width ; del _height



class Ui_MainWindow(object):
    def setupUi(self, MainWindow, WD_path, WD_name, root_path):
        global sz_w
        global sz_h
        _height = QtWidgets.QDesktopWidget().screenGeometry(-1).height()
        _width = QtWidgets.QDesktopWidget().screenGeometry(-1).width()
        sz_w = _width/1366
        sz_h = _height/768

        self.root_path = root_path
        self.WD_path   = WD_path
        self.WD_name   = WD_name

        MainWindow.setObjectName("Brookesia")
        MainWindow.resize(790*sz_w, 671*sz_h)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.tablet = QtWidgets.QTabWidget(self.centralwidget)
        self.tablet.setGeometry(QtCore.QRect(10*sz_w, 20*sz_h, 771*sz_w, 610*sz_h))
        self.tablet.setUsesScrollButtons(False)
        self.tablet.setTabsClosable(False)
        self.tablet.setMovable(False)
        self.tablet.setTabBarAutoHide(False)
        self.tablet.setObjectName("tablet")



# =============================================================================
#%%         Main parameters
# =============================================================================

        self.Main_param = QtWidgets.QWidget()
        self.Main_param.setObjectName("Main_param")

        self.label_folder_name = QtWidgets.QLabel(self.Main_param)
        self.label_folder_name.setGeometry(QtCore.QRect(400*sz_w, 110*sz_h, 81*sz_w, 18*sz_h))
        self.label_folder_name.setObjectName("label_folder_name")
        self.Text_folder_name = QtWidgets.QPlainTextEdit(self.Main_param)
        self.Text_folder_name.setGeometry(QtCore.QRect(580*sz_w, 100*sz_h, 150*sz_w, 31*sz_h))
        self.Text_folder_name.setObjectName("Text_folder_name")

        self.Text_file_name = QtWidgets.QPlainTextEdit(self.Main_param)
        self.Text_file_name.setGeometry(QtCore.QRect(580*sz_w, 55*sz_h, 150*sz_w, 31*sz_h))
        self.Text_file_name.setObjectName("Text_file_name")

        self.pB_load = QtWidgets.QPushButton(self.Main_param)
        self.pB_load.setGeometry(QtCore.QRect(400*sz_w, 10*sz_h, 171*sz_w, 34*sz_h))
        self.pB_load.setObjectName("pB_load")

        self.pB_save = QtWidgets.QPushButton(self.Main_param)
        self.pB_save.setGeometry(QtCore.QRect(400*sz_w, 50*sz_h, 171*sz_w, 34*sz_h))
        self.pB_save.setObjectName("pB_save")

        self.pB_Reference_mechanism = QtWidgets.QPushButton(self.Main_param)
        self.pB_Reference_mechanism.setGeometry(QtCore.QRect(10*sz_w, 10*sz_h, 151*sz_w, 34*sz_h))
        self.pB_Reference_mechanism.setObjectName("pB_Reference_mechanism")
        self.pB_Reduced_mechanism = QtWidgets.QPushButton(self.Main_param)
        self.pB_Reduced_mechanism.setGeometry(QtCore.QRect(10*sz_w, 50*sz_h, 151*sz_w, 34*sz_h))
        self.pB_Reduced_mechanism.setObjectName("pB_Reduced_mechanism")
        self.pB_External_data = QtWidgets.QPushButton(self.Main_param)
        self.pB_External_data.setGeometry(QtCore.QRect(-1000*sz_w, 90*sz_h, 151*sz_w, 34*sz_h))
        self.pB_External_data.setObjectName("pB_External_data")

        self.label_reference_mechanism = QtWidgets.QLabel(self.Main_param)
        self.label_reference_mechanism.setGeometry(QtCore.QRect(180*sz_w, 20*sz_h, 350*sz_w, 18*sz_h))
        font = QtGui.QFont()
        font.setItalic(True)
        self.label_reference_mechanism.setFont(font)
        self.label_reference_mechanism.setObjectName("label_reference_mechanism")
        self.label_reduced_mechanism = QtWidgets.QLabel(self.Main_param)
        self.label_reduced_mechanism.setGeometry(QtCore.QRect(180*sz_w, 60*sz_h, 350*sz_w, 18*sz_h))
        font = QtGui.QFont()
        font.setItalic(True)
        self.label_reduced_mechanism.setFont(font)
        self.label_reduced_mechanism.setObjectName("label_reduced_mechanism")
        self.reduced_mechanism = False

        self.label_External_data = QtWidgets.QLabel(self.Main_param)
        self.label_External_data.setGeometry(QtCore.QRect(-1800*sz_w, 100*sz_h, 350*sz_w, 18*sz_h))
        font = QtGui.QFont()
        font.setItalic(True)
        self.label_External_data.setFont(font)
        self.label_External_data.setObjectName("label_External_data")
        self.external_results = False

        self.Gb_MP_reduction = QtWidgets.QGroupBox(self.Main_param)
        self.Gb_MP_reduction.setGeometry(QtCore.QRect(10*sz_w, 160*sz_h, 520*sz_w, 131*sz_h))
        self.Gb_MP_reduction.setAlignment(QtCore.Qt.AlignCenter)
        self.Gb_MP_reduction.setObjectName("Gb_MP_reduction")
        # label
        self.label_Ref_method_r = QtWidgets.QLabel(self.Gb_MP_reduction)
        self.label_Ref_method_r.setGeometry(QtCore.QRect(10*sz_w, 90*sz_h, 131*sz_w, 18*sz_h))
        self.label_Ref_method_r.setObjectName("label_Ref_method_r")
        # DRG_sp
        self.pB_DRG_sp = QtWidgets.QPushButton(self.Gb_MP_reduction)
        self.pB_DRG_sp.setGeometry(QtCore.QRect(140*sz_w, 40*sz_h, 88*sz_w, 34*sz_h))
        self.pB_DRG_sp.setObjectName("pB_DRG_sp")
        # DRGEP_sp
        self.pB_DRGEP_sp = QtWidgets.QPushButton(self.Gb_MP_reduction)
        self.pB_DRGEP_sp.setGeometry(QtCore.QRect(230*sz_w, 40*sz_h, 88*sz_w, 34*sz_h))
        self.pB_DRGEP_sp.setObjectName("pB_DRGEP_sp")
        # SAR_sp
        self.pB_SAR_sp = QtWidgets.QPushButton(self.Gb_MP_reduction)
        self.pB_SAR_sp.setGeometry(QtCore.QRect(330*sz_w, 40*sz_h, 88*sz_w, 34*sz_h))
        self.pB_SAR_sp.setObjectName("pB_SAR_sp")
        # SARGEP_sp
        self.pB_SARGEP_sp = QtWidgets.QPushButton(self.Gb_MP_reduction)
        self.pB_SARGEP_sp.setGeometry(QtCore.QRect(420*sz_w, 40*sz_h, 88*sz_w, 34*sz_h))
        self.pB_SARGEP_sp.setObjectName("pB_SARGEP_sp")
        # label
        self.label_Red_method_sp = QtWidgets.QLabel(self.Gb_MP_reduction)
        self.label_Red_method_sp.setGeometry(QtCore.QRect(10*sz_w, 50*sz_h, 131*sz_w, 18*sz_h))
        self.label_Red_method_sp.setObjectName("label_Red_method_sp")
        self.tablet.addTab(self.Main_param, "")
        # DRG_r
        self.pB_DRG_r = QtWidgets.QPushButton(self.Gb_MP_reduction)
        self.pB_DRG_r.setGeometry(QtCore.QRect(140*sz_w, 80*sz_h, 88*sz_w, 34*sz_h))
        self.pB_DRG_r.setObjectName("pB_DRG_r")
        # SAR_r
        self.pB_SAR_r = QtWidgets.QPushButton(self.Gb_MP_reduction)
        self.pB_SAR_r.setGeometry(QtCore.QRect(230*sz_w, 80*sz_h, 88*sz_w, 34*sz_h))
        self.pB_SAR_r.setObjectName("pB_SAR_r")


        self.Gb_MP_error = QtWidgets.QGroupBox(self.Main_param)
        self.Gb_MP_error.setGeometry(QtCore.QRect(540*sz_w, 160*sz_h, 181*sz_w, 131*sz_h))
        self.Gb_MP_error.setAlignment(QtCore.Qt.AlignCenter)
        self.Gb_MP_error.setObjectName("Gb_MP_error")

        self.Gb_MP_qoi = QtWidgets.QFrame(self.Gb_MP_error)
        self.Gb_MP_qoi.setGeometry(QtCore.QRect(10*sz_w, 30*sz_h, 161*sz_w, 40*sz_h))
        self.Gb_MP_qoi.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.Gb_MP_qoi.setFrameShadow(QtWidgets.QFrame.Raised)
        self.Gb_MP_qoi.setObjectName("Gb_MP_qoi")

        self.rB_MP_qoi_1 = QtWidgets.QRadioButton(self.Gb_MP_qoi)
        self.rB_MP_qoi_1.setGeometry(QtCore.QRect(15*sz_w, 10*sz_h, 60*sz_w, 18*sz_h))
        self.rB_MP_qoi_1.setObjectName("rB_MP_qoi_1")
        self.rB_MP_qoi_2 = QtWidgets.QRadioButton(self.Gb_MP_qoi)
        self.rB_MP_qoi_2.setGeometry(QtCore.QRect(85*sz_w, 10*sz_h, 60*sz_w, 18*sz_h))
        self.rB_MP_qoi_2.setObjectName("rB_MP_qoi_2")

        self.Gb_MP_errmean = QtWidgets.QFrame(self.Gb_MP_error)
        self.Gb_MP_errmean.setGeometry(QtCore.QRect(10*sz_w, 80*sz_h, 161*sz_w, 40*sz_h))
        self.Gb_MP_errmean.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.Gb_MP_errmean.setFrameShadow(QtWidgets.QFrame.Raised)
        self.Gb_MP_errmean.setObjectName("Gb_MP_errmean")

        self.rB_MP_errmean_1 = QtWidgets.QRadioButton(self.Gb_MP_errmean)
        self.rB_MP_errmean_1.setGeometry(QtCore.QRect(15*sz_w, 10*sz_h, 60*sz_w, 18*sz_h))
        self.rB_MP_errmean_1.setObjectName("rB_MP_errmean_1")
        self.rB_MP_errmean_2 = QtWidgets.QRadioButton(self.Gb_MP_errmean)
        self.rB_MP_errmean_2.setGeometry(QtCore.QRect(85*sz_w, 10*sz_h, 60*sz_w, 18*sz_h))
        self.rB_MP_errmean_2.setObjectName("rB_MP_errmean_2")

        self.Gb_MP_targets = QtWidgets.QGroupBox(self.Main_param)
        self.Gb_MP_targets.setGeometry(QtCore.QRect(10*sz_w, 300*sz_h, 711*sz_w, 212*sz_h))
        self.Gb_MP_targets.setAlignment(QtCore.Qt.AlignCenter)
        self.Gb_MP_targets.setObjectName("Gb_MP_targets")

        self.label_MP_target = QtWidgets.QLabel(self.Gb_MP_targets)
        self.label_MP_target.setGeometry(QtCore.QRect(165*sz_w, 30*sz_h, 81*sz_w, 30*sz_h))
        self.label_MP_target.setObjectName("label_MP_target")
        self.text_MP_T = QtWidgets.QPlainTextEdit(self.Gb_MP_targets)
        self.text_MP_T.setGeometry(QtCore.QRect(170*sz_w, 65*sz_h, 50*sz_w, 31*sz_h))
        self.text_MP_T.setObjectName("text_MP_T")
        self.text_MP_ig = QtWidgets.QPlainTextEdit(self.Gb_MP_targets)
        self.text_MP_ig.setGeometry(QtCore.QRect(170*sz_w, 100*sz_h, 50*sz_w, 31*sz_h))
        self.text_MP_ig.setObjectName("text_MP_ig")
        self.text_MP_Sl = QtWidgets.QPlainTextEdit(self.Gb_MP_targets)
        self.text_MP_Sl.setGeometry(QtCore.QRect(170*sz_w, 135*sz_h, 50*sz_w, 31*sz_h))
        self.text_MP_Sl.setObjectName("text_MP_Sl")
        self.text_MP_K = QtWidgets.QPlainTextEdit(self.Gb_MP_targets)
        self.text_MP_K.setGeometry(QtCore.QRect(170*sz_w, 170*sz_h, 50*sz_w, 31*sz_h))
        self.text_MP_K.setObjectName("text_MP_K")
        self.cB_tsp_T = QtWidgets.QCheckBox(self.Gb_MP_targets)
        self.cB_tsp_T.setGeometry(QtCore.QRect(10*sz_w, 70*sz_h, 181*sz_w, 22*sz_h))
        self.cB_tsp_T.setObjectName("cB_tsp_T")
        self.cB_tsp_T.setChecked(d_T_check)
        self.cB_tsp_igt = QtWidgets.QCheckBox(self.Gb_MP_targets)
        self.cB_tsp_igt.setGeometry(QtCore.QRect(10*sz_w, 105*sz_h, 181*sz_w, 22*sz_h))
        self.cB_tsp_igt.setObjectName("cB_tsp_igt")
        self.cB_tsp_igt.setChecked(d_ig_check)
        self.cB_tsp_Sl = QtWidgets.QCheckBox(self.Gb_MP_targets)
        self.cB_tsp_Sl.setGeometry(QtCore.QRect(10*sz_w, 140*sz_h, 181*sz_w, 22*sz_h))
        self.cB_tsp_Sl.setObjectName("cB_tsp_Sl")
        self.cB_tsp_Sl.setChecked(d_Sl_check)
        self.cB_tsp_K = QtWidgets.QCheckBox(self.Gb_MP_targets)
        self.cB_tsp_K.setGeometry(QtCore.QRect(10*sz_w, 175*sz_h, 181*sz_w, 22*sz_h))
        self.cB_tsp_K.setObjectName("cB_tsp_K")
        self.cB_tsp_K.setChecked(d_K_check)
        self.list_spec = QtWidgets.QListWidget(self.Gb_MP_targets)
        self.list_spec.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)
        self.list_spec.setGeometry(QtCore.QRect(250*sz_w, 30*sz_h, 142*sz_w, 170*sz_h))
        self.list_spec.setObjectName("list_spec")
        self.list_target = QtWidgets.QListWidget(self.Gb_MP_targets)
        self.list_target.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)
        self.list_target.setGeometry(QtCore.QRect(560*sz_w, 30*sz_h, 142*sz_w, 170*sz_h))
        self.list_target.setObjectName("list_target")

        self.pB_MP_add_spec = QtWidgets.QPushButton(self.Gb_MP_targets)
        self.pB_MP_add_spec.setGeometry(QtCore.QRect(415*sz_w, 60*sz_h, 120*sz_w, 30*sz_h))
        self.pB_MP_add_spec.setObjectName("pB_MP_add_spec")
        self.pB_MP_rem_spec = QtWidgets.QPushButton(self.Gb_MP_targets)
        self.pB_MP_rem_spec.setGeometry(QtCore.QRect(415*sz_w, 110*sz_h, 120*sz_w, 30*sz_h))
        self.pB_MP_rem_spec.setObjectName("pB_MP_rem_spec")

        self.label_verbose = QtWidgets.QLabel(self.Main_param)
        self.label_verbose.setGeometry(QtCore.QRect(230*sz_w, 536*sz_h, 71*sz_w, 18*sz_h))
        self.label_verbose.setObjectName("label_verbose")
        self.MP_verbose = QtWidgets.QDoubleSpinBox(self.Main_param)
        self.MP_verbose.setGeometry(QtCore.QRect(280*sz_w, 530*sz_h, 51*sz_w, 32*sz_h))
        self.MP_verbose.setDecimals(0)
        self.MP_verbose.setMaximum(10.0)
        self.MP_verbose.setProperty("value", d_verbose)
        self.MP_verbose.setObjectName("doubleSpinBox")


        self.cB_show_plots = QtWidgets.QCheckBox(self.Main_param)
        self.cB_show_plots.setGeometry(QtCore.QRect(360*sz_w, 525*sz_h, 250*sz_w, 22*sz_h))
        self.cB_show_plots.setObjectName("cB_show_plots")
        self.cB_show_plots.setChecked(d_show_plots)

        self.cB_chemkin_f = QtWidgets.QCheckBox(self.Main_param)
        self.cB_chemkin_f.setGeometry(QtCore.QRect(360*sz_w, 550*sz_h, 250*sz_w, 22*sz_h))
        self.cB_chemkin_f.setObjectName("cB_chemkin_f")
        self.cB_chemkin_f.setChecked(d_chemkin_f)


        self.pB_run = QtWidgets.QPushButton(self.Main_param)
        self.pB_run.setGeometry(QtCore.QRect(560*sz_w, 520*sz_h, 200*sz_w, 50*sz_h))
        self.pB_run.setObjectName("pB_run")

        self.pB_working_dir = QtWidgets.QPushButton(self.Main_param)
        self.pB_working_dir.setGeometry(QtCore.QRect(20*sz_w, 520*sz_h, 200*sz_w, 50*sz_h))
        self.pB_working_dir.setObjectName("pB_working_dir")




# =============================================================================
#       Conditions
# =============================================================================

        self.Conditions = QtWidgets.QWidget()
        self.Conditions.setObjectName("Conditions")
        self.scrollArea = QtWidgets.QScrollArea(self.Conditions)
        self.scrollArea.setGeometry(QtCore.QRect(0, 0, 761*sz_w, 541*sz_h))
        self.scrollArea.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
        self.scrollArea.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
        self.scrollArea.setWidgetResizable(False)
        self.scrollArea.setObjectName("scrollArea")
        self.scrollAreaWidgetContents = QtWidgets.QWidget()
        self.scrollAreaWidgetContents.setGeometry(QtCore.QRect(0, 0, 2700*sz_w, 500*sz_h))
        self.scrollAreaWidgetContents.setObjectName("scrollAreaWidgetContents")


        self.GCases         = []
        self.num_case       = []
        self.rB_reactor_UV  = []
        self.rB_reactor_HP  = []
        self.rB_JSR         = []
        self.rB_PFR         = []
        self.rB_fflame      = []
        self.rB_cfflame     = []
        self.pushButton     = []

        self.add_cases()

        self.scrollArea.setWidget(self.scrollAreaWidgetContents)
        self.tablet.addTab(self.Conditions, "")




        #  Conditions
        self.Condition_Gb   = []
        self.condition_activated = []
        self.txt8       = []
        self.Tmax       = []
        self.Pmax       = []
        self.txt2       = []
        self.txt7       = []
        self.Tincr      = []
        self.fuel_1     = []
        self.Diluent_1  = []
        self.eqmax      = []
        self.txt9       = []
        self.txt6       = []
        self.txt10      = []
        self.Tmin       = []
        self.Pmin       = []
        self.txt1       = []
        self.eqincr     = []
        self.txt5       = []
        self.Diluent_2  = []
        self.txt4       = []
        self.txt3       = []
        self.eqmin      = []
        self.Pincr      = []
        self.oxidant_1  = []

        # reactor options
        self.Gr             = []
        self.Gr1            = []
        self.txtr3          = []
        self.tol_ts_abs_r   = []
        self.txtr2          = []
        self.tol_ts_rel_r   = []
        self.txtr1          = []
        self.Gr2            = []
        self.delta_n_pts    = []
        self.txtr4          = []
        self.pts_num        = []
        self.txtr5          = []
        self.t_max_coeff    = []
        self.txtr8          = []
        self.txtr6          = []
        self.grad_curv      = []
        self.ref_spec       = []
        self.txtr7          = []
        self.Gr3            = []
        self.txtr9          = []
        self.txtr10         = []
        self.txt_Gr_ipn     = []
        self.txt_Gr_its     = []


        # flame options
        self.Gf             = []
        self.Gf1            = []
        self.txtf1          = []
        self.txtf2          = []
        self.txtf3          = []
        self.txtf4          = []
        self.txtf5          = []
        self.txtf5b         = []
        self.txtf6          = []
        self.txtf6b         = []
        self.txtf7          = []
        self.txtf8          = []
        self.tol_ts_rel_f   = []
        self.tol_ts_abs_f   = []
        self.tol_ss_rel_f   = []
        self.tol_ss_abs_f   = []
        self.Gf2            = []
        self.pts_scatter    = []
        self.slope_ff       = []
        self.curve_ff       = []
        self.ratio_ff       = []
        self.prune_ff       = []
        self.Gf3            = []
        self.mult           = []
        self.mix            = []
        self.txt_Gf_xmax    = []
        self.Gf4                 = []
        self.rB_flame_no_restore = []
        self.rB_flame_restore    = []
        self.label_flame_folder  = []

        # counterflow flames
        self.Conditions_diff_flame  = []
        self.Conditions_diff_flame_select = []
        self.rB_cff_diff         = []
        self.rB_cff_pp         = []
        self.rB_cff_tp         = []
        # Pressure
        self.df_P_Gb    = []
        self.txt6_2     = []
        self.df_Pmin     = []
        self.df_Pmax     = []
        self.df_Pincr    = []
        self.txt10_5    = []
        self.txt9_4     = []
        self.txt8_4     = []

        # Lines
        self.line       = []
        self.line_2     = []
        self.line_5     = []
        # Burner 1  -1
        self.df_G11         = []
        self.df_txt_burn1   = []
        self.df_txt         = []
        self.df_Diluent_1   = []
        self.df_oxidant_1   = []
        self.df_txt_2       = []
        self.df_txt_4       = []
        self.df_fuel_1      = []
        self.df_txt_3       = []
        self.df_Diluent_r_1 = []
        # Burner 1  -2
        self.df_G12         = []
        self.txt10_2        = []
        self.df_eqmin_1        = []
        self.txt9_2         = []
        self.txt5_2         = []
        self.df_eqmax_1        = []
        self.df_eqincr_1       = []
        self.txt8_2         = []
        self.df_T_1         = []
        self.txt7_2         = []
        self.df_mdot2_1        = []
        self.df_mdot1_1        = []
        self.txt7_3         = []
        self.df_mdot3_1        = []
        self.df_mdot4_1        = []
        # Burner 2  -1
        self.df_G21         = []
        self.df_txt_burn2         = []
        self.df_txt1_2      = []
        self.df_Diluent_2   = []
        self.df_oxidant_2   = []
        self.df_txt2_2      = []
        self.df_txt4_2      = []
        self.df_fuel_2       = []
        self.df_txt3_2      = []
        self.df_Diluent_r_2 = []
        # Burner 2  -2
        self.df_G22         = []
        self.txt10_3        = []
        self.df_eqmin_2        = []
        self.txt9_3         = []
        self.txt5_3         = []
        self.df_eqmax_2        = []
        self.df_eqincr_2       = []
        self.txt8_3         = []
        self.df_T_2         = []
        self.txt7_4         = []
        self.df_mdot2_2        = []
        self.df_mdot1_2        = []
        self.txt7_5         = []
        self.df_mdot3_2        = []
        self.df_mdot4_2        = []
        # jsr options
        self.GJSR           = []
        self.GJSR1          = []
        self.txtjsr_3       = []
        self.tol_ts_abs_jsr = []
        self.txtjsr_2       = []
        self.tol_ts_rel_jsr = []
        self.txtjsr_1       = []
        self.GJSR2          = []
        self.t_max_jsr      = []
        # pfr options
        self.Gpfr           = []
        self.Gpfr1          = []
        self.Gpfr2          = []
        self.Gpfr3          = []
        self.txtpfr1        = []
        self.txtpfr2        = []
        self.txtpfr3        = []
        self.txtpfr4        = []
        self.txtpfr5        = []
        self.txtpfr6        = []
        self.txtpfr7        = []
        self.txtpfr_8       = []
        self.txtpfr_9       = []
        self.txtpfr_10      = []
        self.txtpfr_11      = []
        self.txtpfr_12      = []
        self.tol_ts_abs_pfr = []
        self.tol_ts_rel_pfr = []
        self.pts_num_pfr    = []
        self.pfr_u_0        = []
        self.pfr_length     = []
        self.pfr_area       = []
        self.pfr_delta_n_pts= []
        self.pfr_pts_num    = []
        self.pfr_t_max_coeff= []
        self.pfr_grad_curv  = []
        self.pfr_ref_spec   = []
        self.rB_pfr_1       = []
        self.rB_pfr_2       = []


# =============================================================================
#        DRG
# =============================================================================

        self.list_operator = []

        self.DRG                    = []

        # DRG
        self.pB_remove_DRG          = []
        self.frame_3                = []
        self.label_DRG_espi         = []
        self.num_DRG_eps            = []
        self.label_DRG_deps         = []
        self.num_DRG_deps           = []
        self.label_Red_method_sp_4  = []
        self.num_DRG_pt_num         = []
        self.groupBox_3             = []
        self.num_DRG_tgt_error   = []
        self.pB_DRG_apply2all       = []
        self.pB_GA_drg              = []
        self.pB_PSO_drg              = []
        self.cB_ISI_drg             = []
        self.tableWidget_DRG        = []



# =============================================================================
#         SA
# =============================================================================

        self.SA                     = []
        self.frame_4                = []
        self.label_SA_espi_2        = []
        self.num_SA_eps             = []
        self.label_SA_deps_2        = []
        self.num_SA_deps            = []
        self.label_SA               = []
        self.num_SA_pt_num          = []
        self.pB_remove_SA           = []
        self.Gb_SA2                 = []
        self.tableWidget_SA         = []
        self.num_SA_tgt_error       = []
        self.pB_SA_apply2all        = []
        self.Gb_SA3                 = []
        self.num_SA_tol_1           = []
        self.num_SA_tol_2           = []
        self.rB_SA_tol_1            = []
        self.rB_SA_tol_2            = []
        self.label_Red_method_sp_8  = []
        self.label_Red_method_sp_9  = []
        self.pB_GA_sa               = []
        self.pB_PSO_sa              = []
        self.cB_ISI_sa              = []




# =============================================================================
#         GA
# =============================================================================

        self.GA                     = []
        self.opt_num                = []
        self.AG_group_gen_ind       = []
        self.label_SA_espi_3        = []
        self.num_GA_gen             = []
        self.label_SA_deps_3        = []
        self.num_GA_ind             = []
        self.AG_group_Arrh          = []
        self.label_scroll_subM      = []
        self.num_GA_A               = []
        self.label_SA_Arrh_2        = []
        self.num_GA_n               = []
        self.label_SA_Arrh_1        = []
        self.label_SA_Arrh_3        = []
        self.num_GA_Ea              = []
        self.GA_Box_Selection       = []
        self.Box_GA_Selection       = []
        self.GA_group5              = []
        self.scrollArea_subM        = []
        self.scrollArea_subM_contents = []
        self.cB_GA_meth             = []
        self.cB_GA_meth_DRG         = []
        self.cB_GA_meth_SA          = []
        self.cB_GA_meth_pts         = []
        self.GA_label_meth_pts      = []
        self.GA_label_51            = []
        self.GA_label_subM          = []
        self.cB_sub_C               = []
        self.cB_sub_CO              = []
        self.cB_sub_H               = []
        self.cB_sub_N               = []
        self.cB_sub_S               = []
        self.cB_sub_Si              = []
        self.opt_meth               = []
        self.num_GA_meth_fract      = []
        self.txt_GA_sel_opt         = []
        self.label_GA_sel_opt       = []
        self.GA_Box_Xover           = []
        self.txt_GA_Xover_opt_1     = []
        self.label_GA_Xover_opt     = []
        self.cB_Xover_op_1          = []
        self.cB_Xover_op_2          = []
        self.cB_Xover_op_3          = []
        self.cB_Xover_op_4          = []
        self.label_GA_Xover_int     = []
        self.txt_GA_Xover_int_1     = []
        self.txt_GA_Xover_int_2     = []
        self.txt_GA_Xover_int_3     = []
        self.txt_GA_Xover_int_4     = []
        self.txt_GA_Xover_opt_2     = []
        self.txt_GA_Xover_opt_3     = []
        self.txt_GA_Xover_opt_4     = []
        self.GA_Box_mut             = []
        self.txt_GA_mut_opt_1       = []
        self.label_GA_mut_int       = []
        self.cB_mut_op_1            = []
        self.cB_mut_op_2            = []
        self.cB_mut_op_3            = []
#        self.cB_mut_op_4            = []
        self.label_GA_mut_opt       = []
        self.txt_GA_mut_int_1       = []
        self.txt_GA_mut_int_2       = []
        self.txt_GA_mut_int_3       = []
#        self.txt_GA_mut_int_4       = []
        self.txt_GA_mut_opt_2       = []
        self.txt_GA_mut_opt_3       = []
#        self.txt_GA_mut_opt_4       = []
        self.GA_mut_frame           = []
        self.label_GA_mut_prob      = []
        self.num_GA_mut_prob        = []
        self.pB_GA_remove           = []
        self.Gb_GA_fitness          = []
        self.rB_GA_fit_1            = []
        self.rB_GA_fit_2            = []
        self.PSO_Box_Param          = []
        self.label_Inertia          = []
        self.Inertia_score          = []
        self.Inertia_min            = []
        self.Inertia_i              = []
        self.Inertia_end            = []
        self.label_cogn             = []
        self.Cogn_i                 = []
        self.Cogn_end               = []
        self.label_soc              = []
        self.Social_i               = []
        self.Social_end             = []
        self.label_pso_param_i      = []
        self.label_pso_param_e      = []
        self.label_pso_param_min    = []
        self.label_pso_param_im     = []
        self.label_pso_param_fm     = []


#        self.progress = QtWidgets.QWidget()
#        self.progress.setObjectName("progress")
#        self.tablet.addTab(self.progress, "")

        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 790*sz_w, 30*sz_h))
        self.menubar.setObjectName("menubar")
        self.menuKinetic_mechanism_optimization_tool = QtWidgets.QMenu(self.menubar)
        self.menuKinetic_mechanism_optimization_tool.setObjectName("menuKinetic_mechanism_optimization_tool")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)
        self.menubar.addAction(self.menuKinetic_mechanism_optimization_tool.menuAction())

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)



    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))





# =============================================================================
#         Main parameters
# =============================================================================

        self.label_folder_name.setText(_translate("MainWindow", "Folder_name"))
        self.Text_folder_name.setPlainText(_translate("MainWindow", "Folder name"))
        self.Text_file_name.setPlainText(_translate("MainWindow", "File name"))
#        self.label_file_name.setText(_translate("MainWindow", "File name"))
        self.pB_load.setText(_translate("MainWindow", "Load previous conditions"))
        self.label_verbose.setText(_translate("MainWindow", "Verbose"))
        self.pB_save.setText(_translate("MainWindow", "Save current conditions"))
        self.pB_Reference_mechanism.setText(_translate("MainWindow", "Reference mechanism"))
        self.pB_Reduced_mechanism.setText(_translate("MainWindow", "Reduced mechanism"))
        self.pB_External_data.setText(_translate("MainWindow", "Import external results"))
        self.pB_MP_add_spec.setText(_translate("MainWindow", "add       =>"))
        self.pB_MP_rem_spec.setText(_translate("MainWindow", "<=    remove"))
        self.label_reference_mechanism.setText(_translate("MainWindow", "no reference mechanism selected"))
        self.label_reduced_mechanism.setText(_translate("MainWindow", "no reduced mechanism selected"))
        self.label_External_data.setText(_translate("MainWindow", "no external results selected"))

        self.Gb_MP_targets.setTitle(_translate("MainWindow", "Targets"))

        self.label_MP_target.setText(_translate("MainWindow", "associated\n  species"))
        self.text_MP_T.setPlainText(_translate("MainWindow", d_sp_T))
        self.text_MP_ig.setPlainText(_translate("MainWindow", d_sp_ig))
        self.text_MP_Sl.setPlainText(_translate("MainWindow", d_sp_Sl))
        self.text_MP_K.setPlainText(_translate("MainWindow", d_sp_K))
        self.cB_tsp_T.setText(_translate("MainWindow", "Temperature"))
        self.cB_tsp_igt.setText(_translate("MainWindow", "Ignition time"))
        self.cB_tsp_Sl.setText(_translate("MainWindow", "Laminar flame speed"))
        self.cB_tsp_K.setText(_translate("MainWindow", "Extinction strain rate"))
        self.Gb_MP_reduction.setTitle(_translate("MainWindow", "Reduction methods"))
        self.label_Ref_method_r.setText(_translate("MainWindow", "Reaction oriented"))
        self.pB_DRG_sp.setText(_translate("MainWindow", "DRG_sp"))
        self.pB_DRGEP_sp.setText(_translate("MainWindow", "DRGEP_sp"))
        self.pB_SAR_sp.setText(_translate("MainWindow", "SAR_sp"))
        self.pB_SARGEP_sp.setText(_translate("MainWindow", "SARGEP_sp"))
        self.label_Red_method_sp.setText(_translate("MainWindow", "Species oriented"))
        self.pB_DRG_r.setText(_translate("MainWindow", "DRG_r"))
        self.pB_SAR_r.setText(_translate("MainWindow", "SAR_r"))

        self.Gb_MP_error.setTitle(_translate("MainWindow", "Errors calculation"))
        self.rB_MP_qoi_1.setText(_translate("MainWindow", "points"))
        self.rB_MP_qoi_2.setText(_translate("MainWindow", "QoI"))
        self.rB_MP_errmean_1.setText(_translate("MainWindow", "mean"))
        self.rB_MP_errmean_2.setText(_translate("MainWindow", "max"))
        if   d_error_calculation == 'points':  self.rB_MP_qoi_1.setChecked(True)
        elif d_error_calculation == 'QoI':     self.rB_MP_qoi_2.setChecked(True)
        if   d_error_coupling == 'mean':   self.rB_MP_errmean_1.setChecked(True)
        elif d_error_coupling == 'max':    self.rB_MP_errmean_2.setChecked(True)


        self.cB_show_plots.setText(_translate("MainWindow", "Show plots during reduction"))
        self.cB_chemkin_f.setText(_translate("MainWindow", "Write Chemkin file"))
        self.pB_run.setText(_translate("MainWindow", "RUN"))
        self.pB_working_dir.setText(_translate("MainWindow", "Working dir: "+self.WD_name))
        self.tablet.setTabText(self.tablet.indexOf(self.Main_param), _translate("MainWindow", "Main parameters"))


# =============================================================================
#         Conditions
# =============================================================================

        self.tablet.setTabText(self.tablet.indexOf(self.Conditions), _translate("MainWindow", "Conditions"))
        self.main_dir = os.getcwd()
        self.menuKinetic_mechanism_optimization_tool.setTitle(_translate("MainWindow", "&Brookesia GUI"))


# =============================================================================
#%%         Actions
# =============================================================================

        self.pB_DRG_sp.clicked.connect(lambda: self.DRG_clic('DRG_sp'))
        self.pB_DRGEP_sp.clicked.connect(lambda: self.DRG_clic('DRGEP_sp'))
        self.pB_DRG_r.clicked.connect(lambda: self.DRG_clic('DRG_r'))
        self.pB_SAR_sp.clicked.connect(lambda: self.SA_clic('SAR_sp'))
        self.pB_SARGEP_sp.clicked.connect(lambda: self.SA_clic('SARGEP_sp'))
        self.pB_SAR_r.clicked.connect(lambda: self.SA_clic('SAR_r'))

        self.pB_save.clicked.connect(lambda: self.write_parameters(True))
        self.pB_load.clicked.connect(self.load_parameters)
        self.pB_run.clicked.connect(self.run_reduction)
        self.pB_working_dir.clicked.connect(self.select_wd_bis)


        self.pB_Reference_mechanism.clicked.connect(lambda: self.get_mech('ref'))
        self.pB_Reduced_mechanism.clicked.connect(lambda: self.get_mech('red'))
        self.pB_External_data.clicked.connect(self.get_ext_res)

        self.pB_MP_add_spec.clicked.connect(lambda: self.add_target_spec(selected_items = self.list_spec.selectedItems()))
        self.pB_MP_rem_spec.clicked.connect(lambda: self.rem_target_spec(self.list_target.selectedItems()))

        self.cB_tsp_T.stateChanged.connect(self.update_error_table)
        self.cB_tsp_Sl.stateChanged.connect(self.update_error_table)
        self.cB_tsp_igt.stateChanged.connect(self.update_error_table)
        self.cB_tsp_K.stateChanged.connect(self.update_error_table)
        self.reference_mechanism = False
        self.reduced_mechanism   = False


    def get_mech(self, mech_type):
        global sz_w
        global sz_h
        os.chdir(self.WD_path)
        A=Files_windows()
        options = QtWidgets.QFileDialog.Options()
        options |= QtWidgets.QFileDialog.DontUseNativeDialog
        fileName, _ = QtWidgets.QFileDialog.getOpenFileName(A, "Select kinetic mechanism", "_kinetic_mech","cti Files (*.cti);;All Files (*);;Python Files (*.py)", options=options)
        _translate = QtCore.QCoreApplication.translate
        font = QtGui.QFont();font.setBold(True);font.setItalic(False);font.setWeight(75)
        if mech_type=='ref':
            self.label_reference_mechanism.setText(_translate("MainWindow", fileName.split('/')[-1]))
            self.label_reference_mechanism.setFont(font)
            self.reference_mechanism = fileName
            self.gas_ref = ct.Solution(self.reference_mechanism)
            self.mech_data = cdef.Mech_data(self.reference_mechanism,self.gas_ref)
            self.ns_ref = self.gas_ref.n_species
            self.nr_ref = self.gas_ref.n_reactions
            _translate = QtCore.QCoreApplication.translate
            for sp in range(self.ns_ref):
                item = QtWidgets.QListWidgetItem()
                self.list_spec.addItem(item)
                item = self.list_spec.item(sp)
                item.setText(_translate("MainWindow", self.gas_ref.species_name(sp)))

        elif mech_type=='red':
            self.label_reduced_mechanism.setText(_translate("MainWindow", fileName.split('/')[-1]))
            self.label_reduced_mechanism.setFont(font)
            self.reduced_mechanism = fileName

        # appearance of "import external data" option
        self.pB_External_data.setGeometry(QtCore.QRect(10*sz_w, 90*sz_h, 151*sz_w, 34*sz_h))
        self.label_External_data.setGeometry(QtCore.QRect(180*sz_w, 100*sz_h, 350*sz_w, 18*sz_h))


    def get_ext_res(self):
        global sz_w
        global sz_h
        A=Files_windows()
        options = QtWidgets.QFileDialog.Options()
        options |= QtWidgets.QFileDialog.DontUseNativeDialog
        fileName, _ = QtWidgets.QFileDialog.getOpenFileName(A, "Select result file", "_results_input","csv Files (*.csv);;txt Files (*.txt);; All Files (*)", options=options)
        _translate = QtCore.QCoreApplication.translate
        font = QtGui.QFont();font.setBold(True);font.setItalic(False);font.setWeight(75)

        self.label_External_data.setText(_translate("MainWindow", fileName.split('/')[-1]))
        self.label_External_data.setFont(font)
        self.external_results = fileName


        class Exp_or_Sim(QtWidgets.QWidget):
            def __init__(self):
                super().__init__()
                self.title = 'PyQt5 input dialogs - pythonspot.com'
                self.left = 10
                self.top = 10
                self.width = 640*sz_w
                self.height = 480*sz_h
                self.initUI()
            def initUI(self):
                self.setWindowTitle(self.title)
                self.setGeometry(self.left, self.top, self.width, self.height)
#                self.getChoice()
                self.getChoice_2()
                self.OptChoice()

                self.show()
            def getChoice(self):
                items = ("Experimental_data","Simulation_results")
                item, okPressed = QtWidgets.QInputDialog.getItem(self, "Get item","Type of data (Exp / simul):", items, 0, False)
                if okPressed and item:
                    self.ext_res_file_type = item
            def getChoice_2(self):
                items = ("Molar_fraction","mol/m3")
                item, okPressed = QtWidgets.QInputDialog.getItem(self, "Get item","Concentration unit:", items, 0, False)
                if okPressed and item:
                    self.ext_res_conc_unit = item
            def OptChoice(self):
                items = ("Genetic Algorithm","No optimization")
                item, okPressed = QtWidgets.QInputDialog.getItem(self, "Get item","Do you want to start with an optimization ?", items, 0, False)
                if okPressed and item:
                    self.opt_choice = item



        ex = Exp_or_Sim()
        ex.close()

        #self.ext_res_file_type = ex.ext_res_file_type
        self.ext_res_file_type = "Experimental_data"
        if ex.ext_res_conc_unit == "mol/m3": ex.ext_res_conc_unit = "mol_m3"
        self.ext_res_conc_unit = ex.ext_res_conc_unit
        # read external data
        conditions_list, ref_results_list = genf.read_ref_data\
        (fileName,self.gas_ref,self.ext_res_conc_unit,self.ext_res_file_type,'tspc',False)

        # if optimization:
        if ex.opt_choice != 'No optimization':
            self.tablet.setCurrentIndex(0)
            if ex.opt_choice == 'Genetic Algorithm':
                self.opt_clic('GA')
            elif ex.opt_choice == 'Particle Swarm Optimization':
                self.opt_clic('PSO')
        self.tablet.setCurrentIndex(0)
        global condition_tab_removed
        if condition_tab_removed == False:
            self.tablet.removeTab(1)
            condition_tab_removed = True


    def get_restore_flame_folder(self, appear_opt,num_case):
        global sz_w
        global sz_h
        global load_prev_cond

        if load_prev_cond:
            appear_opt =False ; load_prev_cond = False  # momentarily prevent the opening of the interface

        if appear_opt:
            A=Files_windows()
            options = QtWidgets.QFileDialog.Options()
            options |= QtWidgets.QFileDialog.DontUseNativeDialog
            folderName = QtWidgets.QFileDialog.getExistingDirectory(A, "Select former flame results folder (must contain .xml files) ", "_results_input/_flame_results", options=options)
            _translate = QtCore.QCoreApplication.translate
            font = QtGui.QFont();font.setBold(True);font.setItalic(False);font.setWeight(75)

            self.label_flame_folder[num_case].setText(_translate("MainWindow", folderName.split(os.getcwd()+'/_results_input/_flame_results/')[1]))
            self.label_flame_folder[num_case].setGeometry(QtCore.QRect(170*sz_w, 65*sz_h, 400*sz_w, 18*sz_h))

        else:
            self.label_flame_folder[num_case].setGeometry(QtCore.QRect(600*sz_w, 65*sz_h, 400*sz_w, 18*sz_h))

    def select_wd_bis(self):

        if os.name == 'nt': #different python path on windows
            python_path = self.root_path.split('lib')[0] + 'python.exe'
        else:
            python_path = self.root_path.split('/lib/')[0] + '/bin/python'

        os.system(python_path + ' ' + self.root_path + "wd_select_gui.py " + self.root_path)

        os.chdir(self.root_path)

        fs = open('wd_activ.txt', 'r')
        txt = fs.readline()
        self.WD_name = txt.split(';')[0]
        self.WD_path = txt.split(';')[1]

        os.chdir(self.WD_path)
        _translate = QtCore.QCoreApplication.translate
        self.pB_working_dir.setText(_translate("MainWindow", "Working dir: "+self.WD_name))


    def remove_tab(self,_option):
        global sz_w
        global sz_h
        idx = self.tablet.currentIndex()
        if condition_tab_removed:   _idx = idx + 1
        else:                       _idx = idx

        idx_opt, idx_DRG, idx_SA = -1, -1, -1
        for i in range(_idx-1):
            if self.list_operator[i]=='opt':
                idx_opt+=1
            elif 'DRG' in self.list_operator[i]:
                idx_DRG+=1
            elif 'SA' in self.list_operator[i]:
                idx_SA+=1


        # remove current operator
        current_operator = self.list_operator[_idx-2]
        if 'DRG' in current_operator:
            self.remove_DRG(idx_DRG)
        elif 'SA' in current_operator:
            self.remove_SA(idx_SA)
        elif 'opt' in current_operator:
            self.remove_opt(idx_opt)



        # remove optimization, if needed
        if len(self.list_operator)>_idx-1:
            next_operator = self.list_operator[_idx-1]
        else:
            next_operator = False
        if next_operator == 'opt':
            self.tablet.removeTab(idx+1)
            del self.list_operator[_idx-1]
            self.remove_opt(idx_opt)


        self.tablet.removeTab(idx)
        del self.list_operator[_idx-2]



    def add_target_spec(self, selected_items):
        global sz_w
        global sz_h
        # get the current target species
        c_tg = []
        for tg in range(self.list_target.count()):
            c_tg.append(self.list_target.item(tg).text())

        _translate = QtCore.QCoreApplication.translate
        for new_tg in selected_items:
            if type(new_tg) is not str: new_tg = new_tg.text()
            if new_tg not in c_tg:
                item = QtWidgets.QListWidgetItem()
                self.list_target.addItem(item)
                item = self.list_target.item(self.list_target.count()-1)
                item.setText(_translate("MainWindow", new_tg))
        self.update_error_table()

    def rem_target_spec(self, selected_items):
        global sz_w
        global sz_h
        for sp2rem in range(len(selected_items)):
            idx = self.list_target.indexFromItem(selected_items[-(sp2rem+1)]).row()
            self.list_target.takeItem(idx)
        self.update_error_table()

    def remove_DRG(self,idx_DRG):
        del self.DRG[idx_DRG]
        del self.pB_remove_DRG[idx_DRG]
        del self.frame_3[idx_DRG]
        del self.label_DRG_espi[idx_DRG]
        del self.num_DRG_eps[idx_DRG]
        del self.label_DRG_deps[idx_DRG]
        del self.num_DRG_deps[idx_DRG]
        del self.label_Red_method_sp_4[idx_DRG]
        del self.num_DRG_pt_num[idx_DRG]
        del self.groupBox_3[idx_DRG]
        del self.num_DRG_tgt_error[idx_DRG]
        del self.pB_DRG_apply2all[idx_DRG]
        del self.pB_GA_drg[idx_DRG]
        del self.pB_PSO_drg[idx_DRG]
        del self.cB_ISI_drg[idx_DRG]
        del self.tableWidget_DRG[idx_DRG]

    def remove_SA(self,idx_SA):
        del self.SA[idx_SA]
        del self.frame_4[idx_SA]
        del self.label_SA_espi_2[idx_SA]
        del self.num_SA_eps[idx_SA]
        del self.label_SA_deps_2[idx_SA]
        del self.num_SA_deps[idx_SA]
        del self.label_SA[idx_SA]
        del self.num_SA_pt_num[idx_SA]
        del self.pB_remove_SA[idx_SA]
        del self.Gb_SA2[idx_SA]
        del self.tableWidget_SA[idx_SA]
        del self.num_SA_tgt_error[idx_SA]
        del self.pB_SA_apply2all[idx_SA]
        del self.Gb_SA3[idx_SA]
        del self.num_SA_tol_1[idx_SA]
        del self.num_SA_tol_2[idx_SA]
        del self.rB_SA_tol_1[idx_SA]
        del self.rB_SA_tol_2[idx_SA]
        del self.label_Red_method_sp_8[idx_SA]
        del self.label_Red_method_sp_9[idx_SA]
        del self.pB_GA_sa[idx_SA]
        del self.pB_PSO_sa[idx_SA]
        del self.cB_ISI_sa[idx_SA]


    def remove_opt(self,idx_opt):
        del self.GA[idx_opt]
        del self.AG_group_gen_ind[idx_opt]
        del self.label_SA_espi_3[idx_opt]
        del self.num_GA_gen[idx_opt]
        del self.label_SA_deps_3[idx_opt]
        del self.num_GA_ind[idx_opt]
        del self.AG_group_Arrh[idx_opt]
        del self.label_scroll_subM[idx_opt]
        del self.num_GA_A[idx_opt]
        del self.label_SA_Arrh_2[idx_opt]
        del self.num_GA_n[idx_opt]
        del self.label_SA_Arrh_1[idx_opt]
        del self.label_SA_Arrh_3[idx_opt]
        del self.num_GA_Ea[idx_opt]
        del self.GA_Box_Selection[idx_opt]
        del self.Box_GA_Selection[idx_opt]
        del self.GA_group5[idx_opt]
        del self.scrollArea_subM[idx_opt]
        del self.scrollArea_subM_contents[idx_opt]
        del self.cB_GA_meth[idx_opt]
        del self.cB_GA_meth_DRG[idx_opt]
        del self.cB_GA_meth_SA[idx_opt]
        del self.cB_GA_meth_pts[idx_opt]
        del self.GA_label_meth_pts[idx_opt]
        del self.GA_label_51[idx_opt]
        del self.GA_label_subM[idx_opt]
        del self.cB_sub_C[idx_opt]
        del self.cB_sub_CO[idx_opt]
        del self.cB_sub_H[idx_opt]
        del self.cB_sub_N[idx_opt]
        del self.cB_sub_S[idx_opt]
        del self.cB_sub_Si[idx_opt]
        del self.num_GA_meth_fract[idx_opt]
        del self.txt_GA_sel_opt[idx_opt]
        del self.label_GA_sel_opt[idx_opt]
        del self.GA_Box_Xover[idx_opt]
        del self.txt_GA_Xover_opt_1[idx_opt]
        del self.label_GA_Xover_opt[idx_opt]
        del self.cB_Xover_op_1[idx_opt]
        del self.cB_Xover_op_2[idx_opt]
        del self.cB_Xover_op_3[idx_opt]
        del self.cB_Xover_op_4[idx_opt]
        del self.label_GA_Xover_int[idx_opt]
        del self.txt_GA_Xover_int_1[idx_opt]
        del self.txt_GA_Xover_int_2[idx_opt]
        del self.txt_GA_Xover_int_3[idx_opt]
        del self.txt_GA_Xover_int_4[idx_opt]
        del self.txt_GA_Xover_opt_2[idx_opt]
        del self.txt_GA_Xover_opt_3[idx_opt]
        del self.txt_GA_Xover_opt_4[idx_opt]
        del self.GA_Box_mut[idx_opt]
        del self.txt_GA_mut_opt_1[idx_opt]
        del self.label_GA_mut_int[idx_opt]
        del self.cB_mut_op_1[idx_opt]
        del self.cB_mut_op_2[idx_opt]
        del self.cB_mut_op_3[idx_opt]
        del self.label_GA_mut_opt[idx_opt]
        del self.txt_GA_mut_int_1[idx_opt]
        del self.txt_GA_mut_int_2[idx_opt]
        del self.txt_GA_mut_int_3[idx_opt]
#        del self.txt_GA_mut_int_4[idx_opt]
        del self.txt_GA_mut_opt_2[idx_opt]
        del self.txt_GA_mut_opt_3[idx_opt]
#        del self.txt_GA_mut_opt_4[idx_opt]
        del self.GA_mut_frame[idx_opt]
        del self.label_GA_mut_prob[idx_opt]
        del self.num_GA_mut_prob[idx_opt]
        del self.pB_GA_remove[idx_opt]
        del self.Gb_GA_fitness[idx_opt]
        del self.rB_GA_fit_1[idx_opt]
        del self.rB_GA_fit_2[idx_opt]
        del self.PSO_Box_Param[idx_opt]
        del self.label_Inertia[idx_opt]
        del self.Inertia_score[idx_opt]
        del self.Inertia_min[idx_opt]
        del self.Inertia_i[idx_opt]
        del self.Inertia_end[idx_opt]
        del self.label_cogn[idx_opt]
        del self.Cogn_i[idx_opt]
        del self.Cogn_end[idx_opt]
        del self.label_soc[idx_opt]
        del self.Social_i[idx_opt]
        del self.Social_end[idx_opt]
        del self.label_pso_param_i[idx_opt]
        del self.label_pso_param_e[idx_opt]
        del self.label_pso_param_min[idx_opt]
        del self.label_pso_param_im[idx_opt]
        del self.label_pso_param_fm[idx_opt]




    def DRG_clic(self,_option):
        global sz_w
        global sz_h
        self.DRG.append(QtWidgets.QWidget())

        self.DRG[-1].setObjectName("DRG")
        self.pB_remove_DRG.append(QtWidgets.QPushButton(self.DRG[-1]))
        self.pB_remove_DRG[-1].setGeometry(QtCore.QRect(670*sz_w, 20*sz_h, 80*sz_w, 70*sz_h))
        self.pB_remove_DRG[-1].setObjectName("pB_remove_DRG")
        self.frame_3.append(QtWidgets.QFrame(self.DRG[-1]))
        self.frame_3[-1].setGeometry(QtCore.QRect(40*sz_w, 30*sz_h, 201*sz_w, 91*sz_h))
        self.frame_3[-1].setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_3[-1].setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_3[-1].setObjectName("frame_3")
        self.label_DRG_espi.append(QtWidgets.QLabel(self.frame_3[-1]))
        self.label_DRG_espi[-1].setGeometry(QtCore.QRect(10*sz_w, 20*sz_h, 131*sz_w, 18*sz_h))
        self.label_DRG_espi[-1].setObjectName("label_DRG_espi")
        self.num_DRG_eps.append(QtWidgets.QDoubleSpinBox(self.frame_3[-1]))
        self.num_DRG_eps[-1].setGeometry(QtCore.QRect(110*sz_w, 10*sz_h, 81*sz_w, 32*sz_h))
        self.num_DRG_eps[-1].setDecimals(4)
        self.num_DRG_eps[-1].setMaximum(10.0)
        self.num_DRG_eps[-1].setSingleStep(0.01)
        self.num_DRG_eps[-1].setObjectName("num_DRG_eps")
        self.label_DRG_deps.append(QtWidgets.QLabel(self.frame_3[-1]))
        self.label_DRG_deps[-1].setGeometry(QtCore.QRect(10*sz_w, 60*sz_h, 131*sz_w, 18*sz_h))
        self.label_DRG_deps[-1].setObjectName("label_DRG_deps")
        self.num_DRG_deps.append(QtWidgets.QDoubleSpinBox(self.frame_3[-1]))
        self.num_DRG_deps[-1].setGeometry(QtCore.QRect(110*sz_w, 50*sz_h, 81*sz_w, 32*sz_h))
        self.num_DRG_deps[-1].setDecimals(4)
        self.num_DRG_deps[-1].setMaximum(2.0)
        self.num_DRG_deps[-1].setObjectName("num_DRG_deps")
        self.label_Red_method_sp_4.append(QtWidgets.QLabel(self.DRG[-1]))
        self.label_Red_method_sp_4[-1].setGeometry(QtCore.QRect(350*sz_w, 60*sz_h, 191*sz_w, 18*sz_h))
        self.label_Red_method_sp_4[-1].setObjectName("label_Red_method_sp_4")
        self.num_DRG_pt_num.append(QtWidgets.QDoubleSpinBox(self.DRG[-1]))
        self.num_DRG_pt_num[-1].setGeometry(QtCore.QRect(550*sz_w, 50*sz_h, 61*sz_w, 32*sz_h))
        self.num_DRG_pt_num[-1].setDecimals(0)
        self.num_DRG_pt_num[-1].setMaximum(500.0)
        self.num_DRG_pt_num[-1].setSingleStep(5.0)
        self.num_DRG_pt_num[-1].setObjectName("num_DRG_pt_num")
        self.groupBox_3.append(QtWidgets.QGroupBox(self.DRG[-1]))
        self.groupBox_3[-1].setGeometry(QtCore.QRect(40*sz_w, 130*sz_h, 261*sz_w, 301*sz_h))
        self.groupBox_3[-1].setAlignment(QtCore.Qt.AlignCenter)
        self.groupBox_3[-1].setObjectName("groupBox_3")
        self.tableWidget_DRG.append(QtWidgets.QTableWidget(self.groupBox_3[-1]))
        self.tableWidget_DRG[-1].setGeometry(QtCore.QRect(10*sz_w, 30*sz_h, 241*sz_w, 221*sz_h))
        self.tableWidget_DRG[-1].setObjectName("tableWidget_DRG")

        self.num_DRG_eps[-1].setProperty("value", d_DRG_eps)
        self.num_DRG_deps[-1].setProperty("value", d_DRG_deps)
        self.num_DRG_pt_num[-1].setProperty("value", d_DRG_pt_num)

        self.update_error_table()

        self.num_DRG_tgt_error.append(QtWidgets.QDoubleSpinBox(self.groupBox_3[-1]))
        self.num_DRG_tgt_error[-1].setGeometry(QtCore.QRect(40*sz_w, 260*sz_h, 61*sz_w, 32*sz_h))
        self.num_DRG_tgt_error[-1].setDecimals(0)
        self.num_DRG_tgt_error[-1].setMaximum(100.0)
        self.num_DRG_tgt_error[-1].setSingleStep(5.0)
        self.num_DRG_tgt_error[-1].setObjectName("num_DRG_tgt_error")
        self.num_DRG_tgt_error[-1].setProperty("value", d_DRG_tgt_error)

        self.pB_DRG_apply2all.append(QtWidgets.QPushButton(self.groupBox_3[-1]))
        self.pB_DRG_apply2all[-1].setGeometry(QtCore.QRect(140*sz_w, 260*sz_h, 88*sz_w, 34*sz_h))
        self.pB_DRG_apply2all[-1].setObjectName("pB_DRG_apply2all")

        self.pB_GA_drg.append(QtWidgets.QPushButton(self.DRG[-1]))
        self.pB_GA_drg[-1].setGeometry(QtCore.QRect(415*sz_w, 330*sz_h, 190*sz_w, 34*sz_h))
        self.pB_GA_drg[-1].setObjectName("pB_GA_drg")
        self.pB_GA_drg[-1].clicked.connect(lambda: self.opt_clic('GA'))

        self.pB_PSO_drg.append(QtWidgets.QPushButton(self.DRG[-1]))
        self.pB_PSO_drg[-1].setGeometry(QtCore.QRect(41500*sz_w, 370*sz_h, 190*sz_w, 34*sz_h))
        self.pB_PSO_drg[-1].setObjectName("pB_PSO_drg")
        self.pB_PSO_drg[-1].clicked.connect(lambda: self.opt_clic('PSO'))

        # inter-species interaction
        self.cB_ISI_drg.append(QtWidgets.QCheckBox(self.DRG[-1]))
        self.cB_ISI_drg[-1].setGeometry(QtCore.QRect(430*sz_w, 250*sz_h, 181*sz_w, 22*sz_h))
        self.cB_ISI_drg[-1].setObjectName("cB_ISI_drg")

        self.tablet.addTab(self.DRG[-1], "")

        _translate = QtCore.QCoreApplication.translate
        self.pB_remove_DRG[-1].setText(_translate("MainWindow", "Remove"))
        self.label_DRG_espi[-1].setText(_translate("MainWindow", "Initial epsilon"))
        self.label_DRG_deps[-1].setText(_translate("MainWindow", "Delta epsilon"))
        self.label_Red_method_sp_4[-1].setText(_translate("MainWindow", "DRG calculation point number"))
        self.groupBox_3[-1].setTitle(_translate("MainWindow", "Errors"))
        self.pB_DRG_apply2all[-1].setText(_translate("MainWindow", "Apply to all"))
        self.cB_ISI_drg[-1].setText(_translate("MainWindow", "Inter-species interactions"))
        self.pB_GA_drg[-1].setText(_translate("MainWindow", "Genetic Algorithm Optimization"))
        self.pB_PSO_drg[-1].setText(_translate("MainWindow", "Particle Swarm Optimization"))


        if _option == 'DRG_sp':
            self.tablet.setTabText(self.tablet.indexOf(self.DRG[-1]), _translate("MainWindow", "DRG_sp"))
            self.list_operator.append('DRG_sp')
        elif _option == 'DRGEP_sp':
            self.tablet.setTabText(self.tablet.indexOf(self.DRG[-1]), _translate("MainWindow", "DRGEP_sp"))
            self.list_operator.append('DRGEP_sp')
        if _option == 'DRG_r':
            self.tablet.setTabText(self.tablet.indexOf(self.DRG[-1]), _translate("MainWindow", "DRG_r"))
            self.list_operator.append('DRG_r')

        self.pB_remove_DRG[-1].clicked.connect(lambda: self.remove_tab(_option))
        self.pB_DRG_apply2all[-1].clicked.connect(self.change_error_table)
        self.cB_ISI_drg[-1].setChecked(d_DRG_ISI)

    def SA_clic(self,_option):
        global sz_w
        global sz_h
        self.SA.append(QtWidgets.QWidget())

        _translate = QtCore.QCoreApplication.translate

        self.SA[-1].setObjectName("SA")
        self.frame_4.append(QtWidgets.QFrame(self.SA[-1]))
        self.frame_4[-1].setGeometry(QtCore.QRect(40*sz_w, 30*sz_h, 201*sz_w, 91*sz_h))
        self.frame_4[-1].setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_4[-1].setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_4[-1].setObjectName("frame_4")
        self.label_SA_espi_2.append(QtWidgets.QLabel(self.frame_4[-1]))
        self.label_SA_espi_2[-1].setGeometry(QtCore.QRect(10*sz_w, 20*sz_h, 131*sz_w, 18*sz_h))
        self.label_SA_espi_2[-1].setObjectName("label_SA_espi_2")
        self.label_SA_espi_2[-1].setText(_translate("MainWindow", "Initial epsilon"))
        self.num_SA_eps.append(QtWidgets.QDoubleSpinBox(self.frame_4[-1]))
        self.num_SA_eps[-1].setGeometry(QtCore.QRect(110*sz_w, 10*sz_h, 81*sz_w, 32*sz_h))
        self.num_SA_eps[-1].setDecimals(4)
        self.num_SA_eps[-1].setMaximum(10.0)
        self.num_SA_eps[-1].setSingleStep(0.01)
        self.num_SA_eps[-1].setObjectName("num_SA_eps")
        self.label_SA_deps_2.append(QtWidgets.QLabel(self.frame_4[-1]))
        self.label_SA_deps_2[-1].setGeometry(QtCore.QRect(10*sz_w, 60*sz_h, 131*sz_w, 18*sz_h))
        self.label_SA_deps_2[-1].setObjectName("label_SA_deps_2")
        self.label_SA_deps_2[-1].setText(_translate("MainWindow", "Delta epsilon"))
        self.num_SA_deps.append(QtWidgets.QDoubleSpinBox(self.frame_4[-1]))
        self.num_SA_deps[-1].setGeometry(QtCore.QRect(110*sz_w, 50*sz_h, 81*sz_w, 32*sz_h))
        self.num_SA_deps[-1].setDecimals(4)
        self.num_SA_deps[-1].setMaximum(2.0)
        self.num_SA_deps[-1].setObjectName("num_SA_deps")
        self.label_SA.append(QtWidgets.QLabel(self.SA[-1]))
        self.label_SA[-1].setGeometry(QtCore.QRect(350*sz_w, 60*sz_h, 191*sz_w, 18*sz_h))
        self.label_SA[-1].setObjectName("label_SA")
        self.label_SA[-1].setText(_translate("MainWindow", "SA calculation point number"))
        self.num_SA_pt_num.append(QtWidgets.QDoubleSpinBox(self.SA[-1]))
        self.num_SA_pt_num[-1].setGeometry(QtCore.QRect(550*sz_w, 50*sz_h, 61*sz_w, 32*sz_h))
        self.num_SA_pt_num[-1].setDecimals(0)
        self.num_SA_pt_num[-1].setMaximum(500.0)
        self.num_SA_pt_num[-1].setSingleStep(5.0)
        self.num_SA_pt_num[-1].setObjectName("num_SA_pt_num")
        self.pB_remove_SA.append(QtWidgets.QPushButton(self.SA[-1]))
        self.pB_remove_SA[-1].setGeometry(QtCore.QRect(670*sz_w, 20*sz_h, 80*sz_w, 70*sz_h))
        self.pB_remove_SA[-1].setObjectName("pB_remove_SA")
        self.pB_remove_SA[-1].setText(_translate("MainWindow", "Remove"))

        self.num_SA_eps[-1].setProperty("value", d_SA_eps)
        self.num_SA_deps[-1].setProperty("value", d_SA_deps)
        self.num_SA_pt_num[-1].setProperty("value", d_SA_pt_num)

        self.Gb_SA2.append(QtWidgets.QGroupBox(self.SA[-1]))
        self.Gb_SA2[-1].setGeometry(QtCore.QRect(40*sz_w, 130*sz_h, 241*sz_w, 301*sz_h))
        self.Gb_SA2[-1].setAlignment(QtCore.Qt.AlignCenter)
        self.Gb_SA2[-1].setObjectName("Gb_SA2")
        self.Gb_SA2[-1].setTitle(_translate("MainWindow", "Errors"))
        self.tableWidget_SA.append(QtWidgets.QTableWidget(self.Gb_SA2[-1]))
        self.tableWidget_SA[-1].setGeometry(QtCore.QRect(10*sz_w, 30*sz_h, 241*sz_w, 221*sz_h))
        self.tableWidget_SA[-1].setObjectName("tableWidget_SA")
        self.update_error_table()
        self.num_SA_tgt_error.append(QtWidgets.QDoubleSpinBox(self.Gb_SA2[-1]))
        self.num_SA_tgt_error[-1].setGeometry(QtCore.QRect(40*sz_w, 260*sz_h, 61*sz_w, 32*sz_h))
        self.num_SA_tgt_error[-1].setDecimals(0)
        self.num_SA_tgt_error[-1].setMaximum(100.0)
        self.num_SA_tgt_error[-1].setSingleStep(5.0)
        self.num_SA_tgt_error[-1].setObjectName("num_SA_tgt_error")
        self.pB_SA_apply2all.append(QtWidgets.QPushButton(self.Gb_SA2[-1]))
        self.pB_SA_apply2all[-1].setGeometry(QtCore.QRect(140*sz_w, 260*sz_h, 88*sz_w, 34*sz_h))
        self.pB_SA_apply2all[-1].setObjectName("pB_SA_apply2all")
        self.pB_SA_apply2all[-1].setText(_translate("MainWindow", "Apply to all"))

        self.num_SA_tgt_error[-1].setProperty("value", d_SA_tgt_error)

        self.Gb_SA3.append(QtWidgets.QGroupBox(self.SA[-1]))
        self.Gb_SA3[-1].setGeometry(QtCore.QRect(310*sz_w, 130*sz_h, 431*sz_w, 171*sz_h))
        self.Gb_SA3[-1].setAlignment(QtCore.Qt.AlignCenter)
        self.Gb_SA3[-1].setObjectName("Gb_SA3")
        self.Gb_SA3[-1].setTitle(_translate("MainWindow", "Sensitivity tolerances"))

        self.rB_SA_tol_1.append(QtWidgets.QRadioButton(self.Gb_SA3[-1]))
        self.rB_SA_tol_1[-1].setGeometry(QtCore.QRect(20*sz_w, 70*sz_h, 231*sz_w, 18*sz_h))
        self.rB_SA_tol_1[-1].setObjectName("rB_SA_tol_1")
        self.rB_SA_tol_1[-1].setText(_translate("MainWindow", "same as simulation"))
        self.rB_SA_tol_2.append(QtWidgets.QRadioButton(self.Gb_SA3[-1]))
        self.rB_SA_tol_2[-1].setGeometry(QtCore.QRect(20*sz_w, 110*sz_h, 231*sz_w, 18*sz_h))
        self.rB_SA_tol_2[-1].setObjectName("rB_SA_tol_2")
        self.rB_SA_tol_2[-1].setText(_translate("MainWindow", "define specific tolerances"))

        self.label_Red_method_sp_8.append(QtWidgets.QLabel(self.Gb_SA3[-1]))
        self.label_Red_method_sp_8[-1].setGeometry(QtCore.QRect(2300*sz_w, 70*sz_h, 131*sz_w, 18*sz_h))
        self.label_Red_method_sp_8[-1].setObjectName("label_Red_method_sp_8")
        self.label_Red_method_sp_8[-1].setText(_translate("MainWindow", "relative tolerance"))
        self.num_SA_tol_1.append(QtWidgets.QPlainTextEdit(self.Gb_SA3[-1]))
        self.num_SA_tol_1[-1].setGeometry(QtCore.QRect(3500*sz_w, 60*sz_h, 61*sz_w, 32*sz_h))
        self.num_SA_tol_1[-1].setObjectName("num_SA_tol_1")
        self.label_Red_method_sp_9.append(QtWidgets.QLabel(self.Gb_SA3[-1]))
        self.label_Red_method_sp_9[-1].setGeometry(QtCore.QRect(2300*sz_w, 110*sz_h, 131*sz_w, 18*sz_h))
        self.label_Red_method_sp_9[-1].setObjectName("label_Red_method_sp_9")
        self.label_Red_method_sp_9[-1].setText(_translate("MainWindow", "absolute tolerance"))
        self.num_SA_tol_2.append(QtWidgets.QPlainTextEdit(self.Gb_SA3[-1]))
        self.num_SA_tol_2[-1].setGeometry(QtCore.QRect(3500*sz_w, 100*sz_h, 61*sz_w, 32*sz_h))
        self.num_SA_tol_2[-1].setObjectName("num_SA_tol_2")

        self.num_SA_tol_1[-1].setPlainText(_translate("MainWindow", d_SA_tol_1))
        self.num_SA_tol_2[-1].setPlainText(_translate("MainWindow", d_SA_tol_2))

        self.cB_ISI_sa.append(QtWidgets.QCheckBox(self.SA[-1]))
        self.cB_ISI_sa[-1].setGeometry(QtCore.QRect(430*sz_w, 330*sz_h, 181*sz_w, 22*sz_h))
        self.cB_ISI_sa[-1].setObjectName("cB_ISI_sa")
        self.cB_ISI_sa[-1].setText(_translate("MainWindow", "Inter-species interactions"))
        self.pB_GA_sa.append(QtWidgets.QPushButton(self.SA[-1]))
        self.pB_GA_sa[-1].setGeometry(QtCore.QRect(415*sz_w, 380*sz_h, 190*sz_w, 34*sz_h))
        self.pB_GA_sa[-1].setObjectName("pB_GA_sa")
        self.pB_GA_sa[-1].clicked.connect(lambda: self.opt_clic('GA'))
        self.pB_GA_sa[-1].setText(_translate("MainWindow", "Genetic Algorithm Optimization"))

        self.pB_PSO_sa.append(QtWidgets.QPushButton(self.SA[-1]))
        self.pB_PSO_sa[-1].setGeometry(QtCore.QRect(41500*sz_w, 420*sz_h, 190*sz_w, 34*sz_h))
        self.pB_PSO_sa[-1].setObjectName("pB_PSO_sa")
        self.pB_PSO_sa[-1].clicked.connect(lambda: self.opt_clic('PSO'))
        self.pB_PSO_sa[-1].setText(_translate("MainWindow", "Particle Swarm Optimization"))



        self.tablet.addTab(self.SA[-1], "")


        if _option == 'SAR_sp':
            self.tablet.setTabText(self.tablet.indexOf(self.SA[-1]), _translate("MainWindow", "SAR_sp"))
            self.list_operator.append('SAR_sp')
        elif _option == 'SARGEP_sp':
            self.tablet.setTabText(self.tablet.indexOf(self.SA[-1]), _translate("MainWindow", "SARGEP_sp"))
            self.list_operator.append('SARGEP_sp')
        if _option == 'SAR_r':
            self.tablet.setTabText(self.tablet.indexOf(self.SA[-1]), _translate("MainWindow", "SAR_r"))
            self.list_operator.append('SAR_r')

        self.pB_remove_SA[-1].clicked.connect(lambda: self.remove_tab(_option))
        self.pB_SA_apply2all[-1].clicked.connect(self.change_error_table)
        self.rB_SA_tol_1[-1].toggled.connect(lambda: self.SA_tol_appear('out'))
        self.rB_SA_tol_2[-1].toggled.connect(lambda: self.SA_tol_appear('in'))
        if d_SA_tol:    self.rB_SA_tol_1[-1].setChecked(True)
        else:           self.rB_SA_tol_2[-1].setChecked(True)
        self.cB_ISI_sa[-1].setChecked(d_SA_ISI)

    def SA_tol_appear(self,option):
        global sz_w
        global sz_h
        idx = self.tablet.currentIndex()
        idx_op=0
        for i in range(idx-2):
            if 'SA' in self.list_operator[i] or 'DSRG' in self.list_operator[i]:
                idx_op+=1
        if option == 'in':
            self.label_Red_method_sp_8[idx_op].setGeometry(QtCore.QRect(230*sz_w, 70*sz_h, 131*sz_w, 18*sz_h))
            self.num_SA_tol_1[idx_op].setGeometry(QtCore.QRect(350*sz_w, 60*sz_h, 61*sz_w, 32*sz_h))
            self.label_Red_method_sp_9[idx_op].setGeometry(QtCore.QRect(230*sz_w, 110*sz_h, 131*sz_w, 18*sz_h))
            self.num_SA_tol_2[idx_op].setGeometry(QtCore.QRect(350*sz_w, 100*sz_h, 61*sz_w, 32*sz_h))
        else:
            self.label_Red_method_sp_8[idx_op].setGeometry(QtCore.QRect(2300*sz_w, 70*sz_h, 131*sz_w, 18*sz_h))
            self.num_SA_tol_1[idx_op].setGeometry(QtCore.QRect(3500*sz_w, 60*sz_h, 61*sz_w, 32*sz_h))
            self.label_Red_method_sp_9[idx_op].setGeometry(QtCore.QRect(2300*sz_w, 110*sz_h, 131*sz_w, 18*sz_h))
            self.num_SA_tol_2[idx_op].setGeometry(QtCore.QRect(3500*sz_w, 100*sz_h, 61*sz_w, 32*sz_h))

    def PSO_score_appear(self):
        global sz_w
        global sz_h

        idx = self.tablet.currentIndex()
        if condition_tab_removed:   _idx = idx + 1
        else:                       _idx = idx
        idx_op=-1
        for i in range(_idx-1):
            if 'opt' in self.list_operator[i]:
                idx_op+=1

        if self.Inertia_score[idx_op].isChecked():
            self.Inertia_min[idx_op].setGeometry(QtCore.QRect(260*sz_w, 30*sz_h, 50*sz_w, 30*sz_h))
            self.Inertia_i[idx_op].setGeometry(QtCore.QRect(  260*sz_w, 60*sz_h, 50*sz_w, 30*sz_h))
            self.Inertia_end[idx_op].setGeometry(QtCore.QRect(260*sz_w, 90*sz_h, 50*sz_w, 30*sz_h))
            self.label_pso_param_i[idx_op].setGeometry(QtCore.QRect(170*sz_w, 130*sz_h, 131*sz_w, 22*sz_h))
            self.label_pso_param_e[idx_op].setGeometry(QtCore.QRect(260*sz_w, 130*sz_h, 131*sz_w, 22*sz_h))
            self.label_pso_param_min[idx_op].setGeometry(QtCore.QRect(170*sz_w, 35*sz_h, 131*sz_w, 22*sz_h))
            self.label_pso_param_im[idx_op].setGeometry(QtCore.QRect(170*sz_w, 65*sz_h, 131*sz_w, 22*sz_h))
            self.label_pso_param_fm[idx_op].setGeometry(QtCore.QRect(170*sz_w, 95*sz_h, 131*sz_w, 22*sz_h))
        else:
            self.Inertia_min[idx_op].setGeometry(QtCore.QRect(1700*sz_w, 60*sz_h, 50*sz_w, 32*sz_h))
            self.Inertia_i[idx_op].setGeometry(QtCore.QRect(170*sz_w, 60*sz_h, 50*sz_w, 32*sz_h))
            self.Inertia_end[idx_op].setGeometry(QtCore.QRect(260*sz_w, 60*sz_h, 50*sz_w, 32*sz_h))
            self.label_pso_param_i[idx_op].setGeometry(QtCore.QRect(170*sz_w, 30*sz_h, 131*sz_w, 22*sz_h))
            self.label_pso_param_e[idx_op].setGeometry(QtCore.QRect(260*sz_w, 30*sz_h, 131*sz_w, 22*sz_h))
            self.label_pso_param_min[idx_op].setGeometry(QtCore.QRect(1700*sz_w, 40*sz_h, 131*sz_w, 22*sz_h))
            self.label_pso_param_im[idx_op].setGeometry(QtCore.QRect(1700*sz_w, 60*sz_h, 131*sz_w, 22*sz_h))
            self.label_pso_param_fm[idx_op].setGeometry(QtCore.QRect(1700*sz_w, 80*sz_h, 131*sz_w, 22*sz_h))


    def update_error_table(self):
        global sz_w
        global sz_h
        tspc = []
        for tg in range(self.list_target.count()):
            tspc.append(self.list_target.item(tg).text())

        num_row = len(tspc)
        if self.cB_tsp_T.isChecked():
            num_row+=1 ; tspc.insert(0,'T')
        if self.cB_tsp_Sl.isChecked():
            num_row+=1 ; tspc.insert(0,'Sl')
        if self.cB_tsp_igt.isChecked():
            num_row+=1 ; tspc.insert(0,'Ig_t')
        if self.cB_tsp_K.isChecked():
            num_row+=1 ; tspc.insert(0,'K')


        # update DRG table
        for tw in range(len(self.tableWidget_DRG)):
            self.tableWidget_DRG[tw].setColumnCount(2)
            self.tableWidget_DRG[tw].setRowCount(len(tspc))

            _translate = QtCore.QCoreApplication.translate
            # row titles :
            for r in range(len(tspc)):
                item = QtWidgets.QTableWidgetItem()
                self.tableWidget_DRG[tw].setVerticalHeaderItem(r, item)
                item = self.tableWidget_DRG[tw].verticalHeaderItem(r)
                item.setText(_translate("MainWindow", str(r+1)))
            # column titles :
            col_title = ["Target","Error (%)"]
            for col in range(2):
                item = QtWidgets.QTableWidgetItem()
                self.tableWidget_DRG[tw].setHorizontalHeaderItem(col, item)
                item = self.tableWidget_DRG[tw].horizontalHeaderItem(col)
                item.setText(_translate("MainWindow", col_title[col]))


            # fill the table
            for col in range(2):
                if col==0: txt = tspc
                else:      txt = [30]*len(tspc)
                for r in range(len(tspc)):
                    item = QtWidgets.QTableWidgetItem()
                    self.tableWidget_DRG[tw].setItem(r, col, item)
                    item = self.tableWidget_DRG[tw].item(r, col)
                    item.setText(_translate("MainWindow", str(txt[r])))

        # update SA table
        for tw in range(len(self.tableWidget_SA)):
            self.tableWidget_SA[tw].setColumnCount(2)
            self.tableWidget_SA[tw].setRowCount(len(tspc))

            _translate = QtCore.QCoreApplication.translate
            # row titles :
            for r in range(len(tspc)):
                item = QtWidgets.QTableWidgetItem()
                self.tableWidget_SA[tw].setVerticalHeaderItem(r, item)
                item = self.tableWidget_SA[tw].verticalHeaderItem(r)
                item.setText(_translate("MainWindow", str(r+1)))
            # column titles :
            col_title = ["Target","Error (%)"]
            for col in range(2):
                item = QtWidgets.QTableWidgetItem()
                self.tableWidget_SA[tw].setHorizontalHeaderItem(col, item)
                item = self.tableWidget_SA[tw].horizontalHeaderItem(col)
                item.setText(_translate("MainWindow", col_title[col]))


            # fill the table
            for col in range(2):
                if col==0: txt = tspc
                else:      txt = ['30']*len(tspc)
                for r in range(len(tspc)):
                    item = QtWidgets.QTableWidgetItem()
                    self.tableWidget_SA[tw].setItem(r, col, item)
                    item = self.tableWidget_SA[tw].item(r, col)
                    item.setText(_translate("MainWindow", txt[r]))

    def change_error_table(self):
        global sz_w
        global sz_h
        _translate = QtCore.QCoreApplication.translate
        idx = self.tablet.currentIndex()
        operator = self.list_operator[idx-2]
        idx_op=0
        if 'DRG' in operator :
            for i in range(idx-2):
                if 'DRG' in self.list_operator[i]:
                    idx_op+=1
#            tsp_nb = self.list_target.count()
#            txt = [self.num_DRG_tgt_error[idx_op].text()]*self.list_target.count()
            txt = self.num_DRG_tgt_error[idx_op].text()
            for r in range(self.tableWidget_DRG[idx_op].rowCount()):
                item = QtWidgets.QTableWidgetItem()
                self.tableWidget_DRG[idx_op].setItem(r, 1, item)
                item = self.tableWidget_DRG[idx_op].item(r, 1)
                item.setText(_translate("MainWindow", txt))

        if 'SA' in operator or 'DSRG' in operator:
            for i in range(idx-2):
                if 'SA' in self.list_operator[i] or 'DSRG' in self.list_operator[i]:
                    idx_op+=1
#            txt = [self.num_SA_tgt_error[idx_op].text()]*self.list_target.count()
            txt = self.num_SA_tgt_error[idx_op].text()
            for r in range(self.tableWidget_SA[idx_op].rowCount()):
                item = QtWidgets.QTableWidgetItem()
                self.tableWidget_SA[idx_op].setItem(r, 1, item)
                item = self.tableWidget_SA[idx_op].item(r, 1)
                item.setText(_translate("MainWindow", txt))

    def mdot_incr_calc(self,num_case):
        global sz_w
        global sz_h
        _translate = QtCore.QCoreApplication.translate

        try:
            phi1_min  = float(self.df_eqmin_1[num_case].document().toPlainText())
            phi1_max  = float(self.df_eqmax_1[num_case].document().toPlainText())
            phi1_incr = float(self.df_eqincr_1[num_case].document().toPlainText())

            nb_it = np.floor((phi1_max-phi1_min)/phi1_incr)

            eqmax_2    = float(self.df_eqmax_2[num_case].document().toPlainText())
            eqmin_2    = float(self.df_eqmin_2[num_case].document().toPlainText())
            mdot_min_1 = float(self.df_mdot1_1[num_case].document().toPlainText())
            mdot_max_1 = float(self.df_mdot2_1[num_case].document().toPlainText())
            mdot_min_2 = float(self.df_mdot1_2[num_case].document().toPlainText())
            mdot_max_2 = float(self.df_mdot2_2[num_case].document().toPlainText())

            if nb_it != 0:
                incr_1 = (eqmax_2-eqmin_2)/nb_it
                incr_2 = (mdot_max_1-mdot_min_1)/nb_it
                incr_3 = (mdot_max_2-mdot_min_2)/nb_it
            else:
                incr_1 = 0
                incr_2 = 0
                incr_3 = 0
        except:
            incr_1='error' ; incr_2='error' ; incr_3='error'

        self.df_eqincr_2[num_case].setText(_translate("MainWindow", str(incr_1)))
        self.df_mdot3_1[num_case].setText(_translate("MainWindow", str(incr_2)))
        self.df_mdot3_2[num_case].setText(_translate("MainWindow", str(incr_3)))



    def opt_clic(self, opt_meth, op_null=False):
        global sz_w
        global sz_h
        global opt_num

        idx = self.tablet.currentIndex()
        add_GA = True
        if condition_tab_removed:   _idx = idx + 1
        else:                       _idx = idx
        if idx == 0: idx = 1

        idx_opt = 0
        for i in range(_idx-1):
            if self.list_operator[i]=='opt':
                idx_opt+=1

        if not self.reference_mechanism:
            add_GA = False
            print('Warning: Select a reference mechanism first')

        if len(self.list_operator)>=max(1,_idx):
            if self.list_operator[_idx-1]=='opt':
                add_GA = False

        if add_GA:
            _translate = QtCore.QCoreApplication.translate

            self.GA.insert(idx_opt,QtWidgets.QWidget())
            self.GA[idx_opt].setObjectName("GA")

            self.tablet.insertTab(idx+1,self.GA[idx_opt],"")
            opt_num+=1
            self.opt_num.insert(idx_opt,opt_num)

            # gen / ind number
            self.AG_group_gen_ind.insert(idx_opt,QtWidgets.QFrame(self.GA[idx_opt]))
            self.AG_group_gen_ind[idx_opt].setGeometry(QtCore.QRect(20*sz_w, 20*sz_h, 290*sz_w, 91*sz_h))
            self.AG_group_gen_ind[idx_opt].setFrameShape(QtWidgets.QFrame.StyledPanel)
            self.AG_group_gen_ind[idx_opt].setFrameShadow(QtWidgets.QFrame.Raised)
            self.AG_group_gen_ind[idx_opt].setObjectName("AG_group_gen_ind")
            self.label_SA_espi_3.insert(idx_opt,QtWidgets.QLabel(self.AG_group_gen_ind[idx_opt]))
            self.label_SA_espi_3[idx_opt].setGeometry(QtCore.QRect(20*sz_w, 20*sz_h, 131*sz_w, 18*sz_h))
            self.label_SA_espi_3[idx_opt].setObjectName("label_SA_espi_3")
            self.num_GA_gen.insert(idx_opt,QtWidgets.QDoubleSpinBox(self.AG_group_gen_ind[idx_opt]))
            self.num_GA_gen[idx_opt].setGeometry(QtCore.QRect(180*sz_w, 10*sz_h, 61*sz_w, 32*sz_h))
            self.num_GA_gen[idx_opt].setDecimals(0)
            self.num_GA_gen[idx_opt].setMaximum(500.0)
            self.num_GA_gen[idx_opt].setSingleStep(5.0)
            self.num_GA_gen[idx_opt].setObjectName("num_GA_gen")
            self.label_SA_deps_3.insert(idx_opt,QtWidgets.QLabel(self.AG_group_gen_ind[idx_opt]))
            self.label_SA_deps_3[idx_opt].setGeometry(QtCore.QRect(20*sz_w, 60*sz_h, 131*sz_w, 18*sz_h))
            self.label_SA_deps_3[idx_opt].setObjectName("label_SA_deps_3")
            self.num_GA_ind.insert(idx_opt,QtWidgets.QDoubleSpinBox(self.AG_group_gen_ind[idx_opt]))
            self.num_GA_ind[idx_opt].setGeometry(QtCore.QRect(180*sz_w, 50*sz_h, 61*sz_w, 32*sz_h))
            self.num_GA_ind[idx_opt].setDecimals(0)
            self.num_GA_ind[idx_opt].setMaximum(500.0)
            self.num_GA_ind[idx_opt].setSingleStep(5.0)
            self.num_GA_ind[idx_opt].setObjectName("num_GA_ind")

            self.num_GA_gen[idx_opt].setProperty("value", d_GA_gen)
            self.num_GA_ind[idx_opt].setProperty("value", d_GA_ind)



            # -------------------
            # Sub mech Selection
            self.label_scroll_subM.insert(idx_opt,QtWidgets.QLabel(self.GA[idx_opt]))
            self.label_scroll_subM[idx_opt].setGeometry(QtCore.QRect(80*sz_w, 120*sz_h, 200*sz_w, 18*sz_h))
            self.label_scroll_subM[idx_opt].setObjectName("label_scroll_subM")
            self.label_scroll_subM[idx_opt].setText(_translate("MainWindow", "Sub-mechanism selection"))


            self.scrollArea_subM.insert(idx_opt,QtWidgets.QScrollArea(self.GA[idx_opt]))
            self.scrollArea_subM[idx_opt].setGeometry(QtCore.QRect(20*sz_w,140*sz_h, 290*sz_w, 101*sz_h))
            self.scrollArea_subM[idx_opt].setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
            self.scrollArea_subM[idx_opt].setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
            self.scrollArea_subM[idx_opt].setWidgetResizable(False)
            self.scrollArea_subM[idx_opt].setObjectName("scrollAreasubM")
            self.scrollArea_subM_contents.insert(idx_opt,QtWidgets.QWidget())
            self.scrollArea_subM_contents[idx_opt].setGeometry(QtCore.QRect(0, 0, 265*sz_w, 101*sz_h))
            self.scrollArea_subM_contents[idx_opt].setObjectName("scrollArea_subM_contents")
            self.scrollArea_subM[idx_opt].setWidget(self.scrollArea_subM_contents[idx_opt])

    #        self.scrollArea_subM.setWidget(self.GA_groupsubM)
            self.GA_label_subM.insert(idx_opt,QtWidgets.QLabel(self.scrollArea_subM_contents[idx_opt]))
            self.GA_label_subM[idx_opt].setGeometry(QtCore.QRect(10*sz_w, 65*sz_h, 191*sz_w, 18*sz_h))
            self.GA_label_subM[idx_opt].setObjectName("label_Red_method_sp_5")




            # column 1: H2 CO
            self.cB_sub_H.insert(idx_opt,QtWidgets.QCheckBox(self.scrollArea_subM_contents[idx_opt]))
            self.cB_sub_H[idx_opt].setGeometry(QtCore.QRect(10*sz_w, 10*sz_h, 50*sz_w, 22*sz_h))
            self.cB_sub_H[idx_opt].setObjectName("cB_sub_H")
            self.cB_sub_H[idx_opt].setText(_translate("MainWindow", "H2"))
            self.cB_sub_H[idx_opt].setChecked(True)

            self.cB_sub_CO.insert(idx_opt,QtWidgets.QCheckBox(self.scrollArea_subM_contents[idx_opt]))
            self.cB_sub_CO[idx_opt].setObjectName("cB_sub_CO")
            self.cB_sub_CO[idx_opt].setText(_translate("MainWindow", "CO"))
            if True in self.mech_data.react.subm_CO:
                self.cB_sub_CO[idx_opt].setGeometry(QtCore.QRect(10*sz_w, 35*sz_h, 50*sz_w, 22*sz_h))
                self.cB_sub_CO[idx_opt].setChecked(True)

            else:
                self.cB_sub_CO[idx_opt].setGeometry(QtCore.QRect(-100*sz_w, 30*sz_h, 50*sz_w, 22*sz_h))
                self.cB_sub_CO[idx_opt].setChecked(True)




            # column 2: Nitrogen sulfur silane chemistry
            mv_y = 0

            # N
            self.cB_sub_N.insert(idx_opt,QtWidgets.QCheckBox(self.scrollArea_subM_contents[idx_opt]))
            self.cB_sub_N[idx_opt].setObjectName("cB_sub_N")
            self.cB_sub_N[idx_opt].setText(_translate("MainWindow", "N"))
            if max(self.mech_data.react.subm_N)>0:
                self.cB_sub_N[idx_opt].setGeometry(QtCore.QRect(80*sz_w, 10*sz_h, 50*sz_w, 22*sz_h))
                self.cB_sub_N[idx_opt].setChecked(True)
                mv_y +=2
            else:
                self.cB_sub_N[idx_opt].setGeometry(QtCore.QRect(-80*sz_w, 10*sz_h, 50*sz_w, 22*sz_h))
                self.cB_sub_N[idx_opt].setChecked(False)

            # S
            self.cB_sub_S.insert(idx_opt,QtWidgets.QCheckBox(self.scrollArea_subM_contents[idx_opt]))
            self.cB_sub_S[idx_opt].setObjectName("cB_sub_S")
            self.cB_sub_S[idx_opt].setText(_translate("MainWindow", "S"))
            if max(self.mech_data.react.subm_S)>0:
                self.cB_sub_S[idx_opt].setGeometry(QtCore.QRect(80*sz_w, (10+mv_y)*sz_h, 50*sz_w, 22*sz_h))
                self.cB_sub_S[idx_opt].setChecked(True)
                mv_y +=25
            else :
                self.cB_sub_S[idx_opt].setGeometry(QtCore.QRect(-80*sz_w, (10+mv_y)*sz_h, 50*sz_w, 22*sz_h))
                self.cB_sub_S[idx_opt].setChecked(False)

            # Si
            self.cB_sub_Si.insert(idx_opt,QtWidgets.QCheckBox(self.scrollArea_subM_contents[idx_opt]))
            self.cB_sub_Si[idx_opt].setObjectName("cB_sub_Si")
            self.cB_sub_Si[idx_opt].setText(_translate("MainWindow", "Si"))
            if max(self.mech_data.react.subm_Si)>0:
                self.cB_sub_Si[idx_opt].setGeometry(QtCore.QRect(80*sz_w, (10+mv_y)*sz_h, 50*sz_w, 22*sz_h))
                self.cB_sub_Si[idx_opt].setChecked(True)
            else:
                self.cB_sub_Si[idx_opt].setGeometry(QtCore.QRect(-80*sz_w, (10+mv_y)*sz_h, 50*sz_w, 22*sz_h))
                self.cB_sub_Si[idx_opt].setChecked(False)


            # column 3: C1-Cx
            self.cB_sub_C.insert(idx_opt,[])
            max_subm_C = max(self.mech_data.react.subm_C)
            for sC in range(max_subm_C):
                self.cB_sub_C[idx_opt].append(QtWidgets.QCheckBox(self.scrollArea_subM_contents[idx_opt]))
                self.cB_sub_C[idx_opt][-1].setGeometry(QtCore.QRect(150*sz_w, (10+25*sC)*sz_h, 50*sz_w, 22*sz_h))
                self.cB_sub_C[idx_opt][-1].setObjectName("cB_sub_C"+str(sC))
                self.cB_sub_C[idx_opt][-1].setText(_translate("MainWindow", "C"+str(sC+1)))
                self.cB_sub_C[idx_opt][-1].setChecked(True)
            if max_subm_C>3:
                self.scrollArea_subM[idx_opt].setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
                self.scrollArea_subM_contents[idx_opt].setGeometry(QtCore.QRect(0, 0, 250*sz_w, (max_subm_C+1)*28*sz_h))#(sC+1)*28))

            # -------------------
            # Arrhenius
            self.AG_group_Arrh.insert(idx_opt,QtWidgets.QGroupBox(self.GA[idx_opt]))
            self.AG_group_Arrh[idx_opt].setGeometry(QtCore.QRect(20*sz_w, 250*sz_h, 290*sz_w, 91*sz_h))
            self.AG_group_Arrh[idx_opt].setAlignment(QtCore.Qt.AlignCenter)
            self.AG_group_Arrh[idx_opt].setObjectName("AG_group_Arrh")
            self.AG_group_Arrh[idx_opt].setTitle(_translate("MainWindow", "Arrhenius parameters maximal variation (%)"))
            self.num_GA_A.insert(idx_opt,QtWidgets.QDoubleSpinBox(self.AG_group_Arrh[idx_opt]))
            self.num_GA_A[idx_opt].setGeometry(QtCore.QRect(30*sz_w, 40*sz_h, 50*sz_w, 32*sz_h))
            self.num_GA_A[idx_opt].setDecimals(0)
            self.num_GA_A[idx_opt].setMaximum(500.0)
            self.num_GA_A[idx_opt].setSingleStep(5.0)
            self.num_GA_A[idx_opt].setObjectName("num_GA_A")
            self.num_GA_n.insert(idx_opt,QtWidgets.QDoubleSpinBox(self.AG_group_Arrh[idx_opt]))
            self.num_GA_n[idx_opt].setGeometry(QtCore.QRect(115*sz_w, 40*sz_h, 50*sz_w, 32*sz_h))
            self.num_GA_n[idx_opt].setDecimals(0)
            self.num_GA_n[idx_opt].setMaximum(500.0)
            self.num_GA_n[idx_opt].setSingleStep(5.0)
            self.num_GA_n[idx_opt].setObjectName("num_GA_n")
            self.num_GA_Ea.insert(idx_opt,QtWidgets.QDoubleSpinBox(self.AG_group_Arrh[idx_opt]))
            self.num_GA_Ea[idx_opt].setGeometry(QtCore.QRect(205*sz_w, 40*sz_h, 50*sz_w, 32*sz_h))
            self.num_GA_Ea[idx_opt].setDecimals(0)
            self.num_GA_Ea[idx_opt].setMaximum(500.0)
            self.num_GA_Ea[idx_opt].setSingleStep(5.0)
            self.num_GA_Ea[idx_opt].setObjectName("num_GA_Ea")
            self.label_SA_Arrh_1.insert(idx_opt,QtWidgets.QLabel(self.AG_group_Arrh[idx_opt]))
            self.label_SA_Arrh_1[idx_opt].setGeometry(QtCore.QRect(15*sz_w, 50*sz_h, 21*sz_w, 18*sz_h))
            self.label_SA_Arrh_1[idx_opt].setObjectName("label_SA_Arrh_1")
            self.label_SA_Arrh_1[idx_opt].setText(_translate("MainWindow", "A"))
            self.label_SA_Arrh_2.insert(idx_opt,QtWidgets.QLabel(self.AG_group_Arrh[idx_opt]))
            self.label_SA_Arrh_2[idx_opt].setGeometry(QtCore.QRect(100*sz_w, 50*sz_h, 21*sz_w, 18*sz_h))
            self.label_SA_Arrh_2[idx_opt].setObjectName("label_SA_Arrh_2")
            self.label_SA_Arrh_2[idx_opt].setText(_translate("MainWindow", "n"))
            self.label_SA_Arrh_3.insert(idx_opt,QtWidgets.QLabel(self.AG_group_Arrh[idx_opt]))
            self.label_SA_Arrh_3[idx_opt].setGeometry(QtCore.QRect(180*sz_w, 50*sz_h, 21*sz_w, 18*sz_h))
            self.label_SA_Arrh_3[idx_opt].setObjectName("label_SA_Arrh_3")
            self.label_SA_Arrh_3[idx_opt].setText(_translate("MainWindow", "Ea"))

            self.num_GA_A[idx_opt].setProperty("value", d_GA_A)
            self.num_GA_n[idx_opt].setProperty("value", d_GA_n)
            self.num_GA_Ea[idx_opt].setProperty("value", d_GA_Ea)


            # -------------------
            # Optimization on method
            self.GA_group5.insert(idx_opt,QtWidgets.QGroupBox(self.GA[idx_opt]))
            self.GA_group5[idx_opt].setGeometry(QtCore.QRect(20*sz_w, 355*sz_h, 290*sz_w, 101*sz_h))
            self.GA_group5[idx_opt].setAlignment(QtCore.Qt.AlignCenter)
            self.GA_group5[idx_opt].setObjectName("GA_group5")
            self.cB_GA_meth.insert(idx_opt,QtWidgets.QCheckBox(self.GA_group5[idx_opt]))
            self.cB_GA_meth[idx_opt].setGeometry(QtCore.QRect(10*sz_w, 30*sz_h, 231*sz_w, 22*sz_h))
            self.cB_GA_meth[idx_opt].setObjectName("checkBox_2")
            self.GA_label_51.insert(idx_opt,QtWidgets.QLabel(self.GA_group5[idx_opt]))
            self.GA_label_51[idx_opt].setGeometry(QtCore.QRect(10*sz_w, 65*sz_h, 200*sz_w, 18*sz_h))
            self.GA_label_51[idx_opt].setObjectName("label_Red_method_sp_5")
            self.num_GA_meth_fract.insert(idx_opt,QtWidgets.QDoubleSpinBox(self.GA_group5[idx_opt]))
            self.num_GA_meth_fract[idx_opt].setGeometry(QtCore.QRect(215*sz_w, 60*sz_h, 61*sz_w, 32*sz_h))
            self.num_GA_meth_fract[idx_opt].setDecimals(0)
            self.num_GA_meth_fract[idx_opt].setMaximum(100.0)
            self.num_GA_meth_fract[idx_opt].setSingleStep(5.0)
            self.num_GA_meth_fract[idx_opt].setObjectName("num_GA_meth_fract")
            self.num_GA_meth_fract[idx_opt].setProperty("value", d_GA_meth_fract)

            # -------------------
            # Selection
            self.GA_Box_Selection.insert(idx_opt,QtWidgets.QGroupBox(self.GA[idx_opt]))
            self.GA_Box_Selection[idx_opt].setGeometry(QtCore.QRect(3300*sz_w, 20*sz_h, 321*sz_w, 100*sz_h))
            self.GA_Box_Selection[idx_opt].setAlignment(QtCore.Qt.AlignCenter)
            self.GA_Box_Selection[idx_opt].setObjectName("GA_Box_Selection")
            self.GA_Box_Selection[idx_opt].setTitle(_translate("MainWindow", "Selection operator"))
            self.Box_GA_Selection.insert(idx_opt,QtWidgets.QComboBox(self.GA_Box_Selection[idx_opt]))
            self.Box_GA_Selection[idx_opt].setGeometry(QtCore.QRect(25*sz_w, 30*sz_h, 221*sz_w, 25*sz_h))
            self.Box_GA_Selection[idx_opt].setObjectName("Box_GA_Selection")
            self.Box_GA_Selection[idx_opt].addItem("")
            self.Box_GA_Selection[idx_opt].addItem("")
            self.Box_GA_Selection[idx_opt].addItem("")
            self.Box_GA_Selection[idx_opt].addItem("")
            self.Box_GA_Selection[idx_opt].setItemText(0, _translate("MainWindow", "Roulette"))
            self.Box_GA_Selection[idx_opt].setItemText(1, _translate("MainWindow", "Rank"))
            self.Box_GA_Selection[idx_opt].setItemText(2, _translate("MainWindow", "Geometric_norm"))
            self.Box_GA_Selection[idx_opt].setItemText(3, _translate("MainWindow", "Elitism"))
#            self.Box_GA_Selection[idx_opt].setItemText(4, _translate("MainWindow", "Tournament"))
            self.txt_GA_sel_opt.insert(idx_opt,QtWidgets.QPlainTextEdit(self.GA_Box_Selection[idx_opt]))
            self.txt_GA_sel_opt[idx_opt].setGeometry(QtCore.QRect(145*sz_w, 60*sz_h, 101*sz_w, 31*sz_h))
            self.txt_GA_sel_opt[idx_opt].setObjectName("txt_GA_sel_opt")
            self.label_GA_sel_opt.insert(idx_opt,QtWidgets.QLabel(self.GA_Box_Selection[idx_opt]))
            self.label_GA_sel_opt[idx_opt].setGeometry(QtCore.QRect(25*sz_w, 70*sz_h, 131*sz_w, 18*sz_h))
            self.label_GA_sel_opt[idx_opt].setObjectName("label_GA_sel_opt")
            self.label_GA_sel_opt[idx_opt].setText(_translate("MainWindow", "Selection option"))

            self.txt_GA_sel_opt[idx_opt].setPlainText(_translate("MainWindow", d_GA_sel_opt))

            # -------------------
            # Xover
            self.GA_Box_Xover.insert(idx_opt,QtWidgets.QGroupBox(self.GA[idx_opt]))
            self.GA_Box_Xover[idx_opt].setGeometry(QtCore.QRect(3300*sz_w, 130*sz_h, 420*sz_w, 190*sz_h)) #321
            self.GA_Box_Xover[idx_opt].setAlignment(QtCore.Qt.AlignCenter)
            self.GA_Box_Xover[idx_opt].setObjectName("GA_Box_Xover")
            self.GA_Box_Xover[idx_opt].setTitle(_translate("MainWindow", "Cross-over operator"))

            self.cB_Xover_op_1.insert(idx_opt,QtWidgets.QCheckBox(self.GA_Box_Xover[idx_opt]))
            self.cB_Xover_op_1[idx_opt].setGeometry(QtCore.QRect(20*sz_w, 60*sz_h, 131*sz_w, 22*sz_h))
            self.cB_Xover_op_1[idx_opt].setObjectName("cB_Xover_op_1")
            self.cB_Xover_op_1[idx_opt].setText(_translate("MainWindow", "Simple Xover"))
            self.cB_Xover_op_2.insert(idx_opt,QtWidgets.QCheckBox(self.GA_Box_Xover[idx_opt]))
            self.cB_Xover_op_2[idx_opt].setGeometry(QtCore.QRect(20*sz_w, 90*sz_h, 121*sz_w, 22*sz_h))
            self.cB_Xover_op_2[idx_opt].setObjectName("cB_Xover_op_2")
            self.cB_Xover_op_2[idx_opt].setText(_translate("MainWindow", "Multiple Xover"))
            self.cB_Xover_op_3.insert(idx_opt,QtWidgets.QCheckBox(self.GA_Box_Xover[idx_opt]))
            self.cB_Xover_op_3[idx_opt].setGeometry(QtCore.QRect(20*sz_w, 120*sz_h, 121*sz_w, 22*sz_h))
            self.cB_Xover_op_3[idx_opt].setObjectName("cB_Xover_op_3")
            self.cB_Xover_op_3[idx_opt].setText(_translate("MainWindow", "Arith Xover"))
            self.cB_Xover_op_4.insert(idx_opt,QtWidgets.QCheckBox(self.GA_Box_Xover[idx_opt]))
            self.cB_Xover_op_4[idx_opt].setGeometry(QtCore.QRect(20*sz_w, 150*sz_h, 151*sz_w, 22*sz_h))
            self.cB_Xover_op_4[idx_opt].setObjectName("cB_Xover_op_4")
            self.cB_Xover_op_4[idx_opt].setText(_translate("MainWindow", "Heuristic Xover"))
            self.txt_GA_Xover_int_1.insert(idx_opt,QtWidgets.QPlainTextEdit(self.GA_Box_Xover[idx_opt]))
            self.txt_GA_Xover_int_1[idx_opt].setGeometry(QtCore.QRect(220*sz_w, 50*sz_h, 71*sz_w, 31*sz_h))
            self.txt_GA_Xover_int_1[idx_opt].setObjectName("txt_GA_Xover_int_1")
            self.txt_GA_Xover_int_1[idx_opt].setPlainText(_translate("MainWindow", d_GA_Xover_int_1))
            self.txt_GA_Xover_int_2.insert(idx_opt,QtWidgets.QPlainTextEdit(self.GA_Box_Xover[idx_opt]))
            self.txt_GA_Xover_int_2[idx_opt].setGeometry(QtCore.QRect(220*sz_w, 80*sz_h, 71*sz_w, 31*sz_h))
            self.txt_GA_Xover_int_2[idx_opt].setObjectName("txt_GA_Xover_int_2")
            self.txt_GA_Xover_int_2[idx_opt].setPlainText(_translate("MainWindow", d_GA_Xover_int_2))
            self.txt_GA_Xover_int_3.insert(idx_opt,QtWidgets.QPlainTextEdit(self.GA_Box_Xover[idx_opt]))
            self.txt_GA_Xover_int_3[idx_opt].setGeometry(QtCore.QRect(220*sz_w, 110*sz_h, 71*sz_w, 31*sz_h))
            self.txt_GA_Xover_int_3[idx_opt].setObjectName("txt_GA_Xover_int_3")
            self.txt_GA_Xover_int_3[idx_opt].setPlainText(_translate("MainWindow", d_GA_Xover_int_3))
            self.txt_GA_Xover_int_4.insert(idx_opt,QtWidgets.QPlainTextEdit(self.GA_Box_Xover[idx_opt]))
            self.txt_GA_Xover_int_4[idx_opt].setGeometry(QtCore.QRect(220*sz_w, 140*sz_h, 71*sz_w, 31*sz_h))
            self.txt_GA_Xover_int_4[idx_opt].setObjectName("txt_GA_Xover_int_4")
            self.txt_GA_Xover_int_4[idx_opt].setPlainText(_translate("MainWindow", d_GA_Xover_int_4))
            self.txt_GA_Xover_opt_1.insert(idx_opt,QtWidgets.QPlainTextEdit(self.GA_Box_Xover[idx_opt]))
            self.txt_GA_Xover_opt_1[idx_opt].setGeometry(QtCore.QRect(310*sz_w, 50*sz_h, 71*sz_w, 31*sz_h))
            self.txt_GA_Xover_opt_1[idx_opt].setPlainText(_translate("MainWindow", d_GA_Xover_opt_1))
            self.txt_GA_Xover_opt_1[idx_opt].setObjectName("txt_GA_Xover_opt_1")
            self.txt_GA_Xover_opt_2.insert(idx_opt,QtWidgets.QPlainTextEdit(self.GA_Box_Xover[idx_opt]))
            self.txt_GA_Xover_opt_2[idx_opt].setGeometry(QtCore.QRect(310*sz_w, 80*sz_h, 71*sz_w, 31*sz_h))
            self.txt_GA_Xover_opt_2[idx_opt].setPlainText(_translate("MainWindow", d_GA_Xover_opt_2))
            self.txt_GA_Xover_opt_2[idx_opt].setObjectName("txt_GA_Xover_opt_2")
            self.txt_GA_Xover_opt_3.insert(idx_opt,QtWidgets.QPlainTextEdit(self.GA_Box_Xover[idx_opt]))
            self.txt_GA_Xover_opt_3[idx_opt].setGeometry(QtCore.QRect(310*sz_w, 110*sz_h, 71*sz_w, 31*sz_h))
            self.txt_GA_Xover_opt_3[idx_opt].setPlainText(_translate("MainWindow", d_GA_Xover_opt_3))
            self.txt_GA_Xover_opt_3[idx_opt].setObjectName("txt_GA_Xover_opt_3")
            self.txt_GA_Xover_opt_4.insert(idx_opt,QtWidgets.QPlainTextEdit(self.GA_Box_Xover[idx_opt]))
            self.txt_GA_Xover_opt_4[idx_opt].setGeometry(QtCore.QRect(310*sz_w, 140*sz_h, 71*sz_w, 31*sz_h))
            self.txt_GA_Xover_opt_4[idx_opt].setPlainText(_translate("MainWindow", d_GA_Xover_opt_4))
            self.txt_GA_Xover_opt_4[idx_opt].setObjectName("txt_GA_Xover_opt_4")
            self.label_GA_Xover_int.insert(idx_opt,QtWidgets.QLabel(self.GA_Box_Xover[idx_opt]))
            self.label_GA_Xover_int[idx_opt].setGeometry(QtCore.QRect(220*sz_w, 30*sz_h, 81*sz_w, 18*sz_h))
            self.label_GA_Xover_int[idx_opt].setObjectName("label_GA_Xover_int")
            self.label_GA_Xover_int[idx_opt].setText(_translate("MainWindow", "Intensity(%)"))
            self.label_GA_Xover_opt.insert(idx_opt,QtWidgets.QLabel(self.GA_Box_Xover[idx_opt]))
            self.label_GA_Xover_opt[idx_opt].setGeometry(QtCore.QRect(320*sz_w, 30*sz_h, 71*sz_w, 18*sz_h))
            self.label_GA_Xover_opt[idx_opt].setObjectName("label_GA_Xover_opt")
            self.label_GA_Xover_opt[idx_opt].setText(_translate("MainWindow", "Options"))
            self.pB_GA_remove.insert(idx_opt,QtWidgets.QPushButton(self.GA[idx_opt]))
            self.pB_GA_remove[idx_opt].setGeometry(QtCore.QRect(670*sz_w, 20*sz_h, 80*sz_w, 70*sz_h))
            self.pB_GA_remove[idx_opt].setObjectName("pB_GA_remove")

            # -------------------
            # Mutation
            self.GA_Box_mut.insert(idx_opt,QtWidgets.QGroupBox(self.GA[idx_opt]))
            self.GA_Box_mut[idx_opt].setGeometry(QtCore.QRect(3300*sz_w, 330*sz_h, 420*sz_w, 220*sz_h))
            self.GA_Box_mut[idx_opt].setAlignment(QtCore.Qt.AlignCenter)
            self.GA_Box_mut[idx_opt].setObjectName("GA_Box_mut")
            self.GA_Box_mut[idx_opt].setTitle(_translate("MainWindow", "Mutation operator"))

            self.cB_mut_op_1.insert(idx_opt,QtWidgets.QCheckBox(self.GA_Box_mut[idx_opt]))
            self.cB_mut_op_1[idx_opt].setGeometry(QtCore.QRect(20*sz_w, 60*sz_h, 141*sz_w, 22*sz_h))
            self.cB_mut_op_1[idx_opt].setObjectName("cB_mut_op_1")
            self.cB_mut_op_1[idx_opt].setText(_translate("MainWindow", "Uniform mutation"))
            self.cB_mut_op_2.insert(idx_opt,QtWidgets.QCheckBox(self.GA_Box_mut[idx_opt]))
            self.cB_mut_op_2[idx_opt].setGeometry(QtCore.QRect(20*sz_w, 90*sz_h, 171*sz_w, 22*sz_h))
            self.cB_mut_op_2[idx_opt].setObjectName("cB_mut_op_2")
            self.cB_mut_op_2[idx_opt].setText(_translate("MainWindow", "Non-uniform mutation"))
            self.cB_mut_op_3.insert(idx_opt,QtWidgets.QCheckBox(self.GA_Box_mut[idx_opt]))
            self.cB_mut_op_3[idx_opt].setGeometry(QtCore.QRect(20*sz_w, 120*sz_h, 171*sz_w, 22*sz_h))
            self.cB_mut_op_3[idx_opt].setObjectName("cB_mut_op_3")
            self.cB_mut_op_3[idx_opt].setText(_translate("MainWindow", "Boundary mutation"))
#            self.cB_mut_op_4.insert(idx_opt,QtWidgets.QCheckBox(self.GA_Box_mut[idx_opt]))
#            self.cB_mut_op_4[idx_opt].setGeometry(QtCore.QRect(20, 160, 201, 22))
#            self.cB_mut_op_4[idx_opt].setObjectName("cB_mut_op_4")
#            self.cB_mut_op_4[idx_opt].setText(_translate("MainWindow", "Multi non-uniform mutation"))


            self.txt_GA_mut_int_1.insert(idx_opt,QtWidgets.QPlainTextEdit(self.GA_Box_mut[idx_opt]))
            self.txt_GA_mut_int_1[idx_opt].setGeometry(QtCore.QRect(220*sz_w, 50*sz_h, 71*sz_w, 31*sz_h))
            self.txt_GA_mut_int_1[idx_opt].setObjectName("txt_GA_mut_int_1")
            self.txt_GA_mut_int_1[idx_opt].setPlainText(_translate("MainWindow", d_GA_mut_int_1))
            self.txt_GA_mut_int_2.insert(idx_opt,QtWidgets.QPlainTextEdit(self.GA_Box_mut[idx_opt]))
            self.txt_GA_mut_int_2[idx_opt].setGeometry(QtCore.QRect(220*sz_w, 80*sz_h, 71*sz_w, 31*sz_h))
            self.txt_GA_mut_int_2[idx_opt].setObjectName("txt_GA_mut_int_2")
            self.txt_GA_mut_int_2[idx_opt].setPlainText(_translate("MainWindow", d_GA_mut_int_2))
            self.txt_GA_mut_int_3.insert(idx_opt,QtWidgets.QPlainTextEdit(self.GA_Box_mut[idx_opt]))
            self.txt_GA_mut_int_3[idx_opt].setGeometry(QtCore.QRect(220*sz_w, 110*sz_h, 71*sz_w, 31*sz_h))
            self.txt_GA_mut_int_3[idx_opt].setObjectName("txt_GA_mut_int_3")
            self.txt_GA_mut_int_3[idx_opt].setPlainText(_translate("MainWindow", d_GA_mut_int_3))
#            self.txt_GA_mut_int_4[idx_opt].setPlainText(_translate("MainWindow", "25"))
#            self.txt_GA_mut_int_4.insert(idx_opt,QtWidgets.QPlainTextEdit(self.GA_Box_mut[idx_opt]))
#            self.txt_GA_mut_int_4[idx_opt].setGeometry(QtCore.QRect(250, 150, 71, 31))
#            self.txt_GA_mut_int_4[idx_opt].setObjectName("txt_GA_mut_int_4")
            self.txt_GA_mut_opt_1.insert(idx_opt,QtWidgets.QPlainTextEdit(self.GA_Box_mut[idx_opt]))
            self.txt_GA_mut_opt_1[idx_opt].setGeometry(QtCore.QRect(310*sz_w, 50*sz_h, 71*sz_w, 31*sz_h))
            self.txt_GA_mut_opt_1[idx_opt].setObjectName("txt_GA_mut_opt_1")
            self.txt_GA_mut_opt_2.insert(idx_opt,QtWidgets.QPlainTextEdit(self.GA_Box_mut[idx_opt]))
            self.txt_GA_mut_opt_2[idx_opt].setGeometry(QtCore.QRect(310*sz_w, 80*sz_h, 71*sz_w, 31*sz_h))
            self.txt_GA_mut_opt_2[idx_opt].setObjectName("txt_GA_mut_opt_2")
            self.txt_GA_mut_opt_3.insert(idx_opt,QtWidgets.QPlainTextEdit(self.GA_Box_mut[idx_opt]))
            self.txt_GA_mut_opt_3[idx_opt].setGeometry(QtCore.QRect(310*sz_w, 110*sz_h, 71*sz_w, 31*sz_h))
            self.txt_GA_mut_opt_3[idx_opt].setObjectName("txt_GA_mut_opt_3")
#            self.txt_GA_mut_opt_4.insert(idx_opt,QtWidgets.QPlainTextEdit(self.GA_Box_mut[idx_opt]))
#            self.txt_GA_mut_opt_4[idx_opt].setGeometry(QtCore.QRect(350, 150, 71, 31))
#            self.txt_GA_mut_opt_4[idx_opt].setObjectName("txt_GA_mut_opt_4")

            self.txt_GA_mut_opt_1[idx_opt].setPlainText(_translate("MainWindow", d_GA_mut_opt_1))
            self.txt_GA_mut_opt_2[idx_opt].setPlainText(_translate("MainWindow", d_GA_mut_opt_2))
            self.txt_GA_mut_opt_3[idx_opt].setPlainText(_translate("MainWindow", d_GA_mut_opt_3))
#            self.txt_GA_mut_opt_4[idx_opt].setPlainText("")

            self.label_GA_mut_int.insert(idx_opt,QtWidgets.QLabel(self.GA_Box_mut[idx_opt]))
            self.label_GA_mut_int[idx_opt].setGeometry(QtCore.QRect(320*sz_w, 30*sz_h, 71*sz_w, 18*sz_h))
            self.label_GA_mut_int[idx_opt].setObjectName("label_GA_mut_int")
            self.label_GA_mut_int[idx_opt].setText(_translate("MainWindow", "Options"))
            self.label_GA_mut_opt.insert(idx_opt,QtWidgets.QLabel(self.GA_Box_mut[idx_opt]))
            self.label_GA_mut_opt[idx_opt].setGeometry(QtCore.QRect(210*sz_w, 30*sz_h, 81*sz_w, 18*sz_h))
            self.label_GA_mut_opt[idx_opt].setObjectName("label_GA_mut_opt")
            self.label_GA_mut_opt[idx_opt].setText(_translate("MainWindow", "Intensity(%)"))

            self.GA_mut_frame.insert(idx_opt,QtWidgets.QFrame(self.GA_Box_mut[idx_opt]))
            self.GA_mut_frame[idx_opt].setGeometry(QtCore.QRect(80*sz_w, 155*sz_h, 250*sz_w, 50*sz_h))
            self.GA_mut_frame[idx_opt].setFrameShape(QtWidgets.QFrame.StyledPanel)
            self.GA_mut_frame[idx_opt].setFrameShadow(QtWidgets.QFrame.Raised)
            self.GA_mut_frame[idx_opt].setObjectName("GA_mut_frame")
            self.label_GA_mut_prob.insert(idx_opt,QtWidgets.QLabel(self.GA_mut_frame[idx_opt]))
            self.label_GA_mut_prob[idx_opt].setGeometry(QtCore.QRect(30*sz_w, 10*sz_h, 111*sz_w, 30*sz_h))
            self.label_GA_mut_prob[idx_opt].setObjectName("label_GA_mut_prob")
            self.label_GA_mut_prob[idx_opt].setText(_translate("MainWindow", "Gene Mutation\n   probability"))
            self.num_GA_mut_prob.insert(idx_opt,QtWidgets.QDoubleSpinBox(self.GA_mut_frame[idx_opt]))
            self.num_GA_mut_prob[idx_opt].setGeometry(QtCore.QRect(180*sz_w, 10*sz_h, 61*sz_w, 32*sz_h))
            self.num_GA_mut_prob[idx_opt].setDecimals(0)
            self.num_GA_mut_prob[idx_opt].setMaximum(100.0)
            self.num_GA_mut_prob[idx_opt].setSingleStep(5.0)
            self.num_GA_mut_prob[idx_opt].setProperty("value", d_GA_mut_prob)
            self.num_GA_mut_prob[idx_opt].setObjectName("num_GA_mut_prob")


            # -------------------
            # PSO options
            self.PSO_Box_Param.insert(idx_opt,QtWidgets.QGroupBox(self.GA[idx_opt]))
#            self.PSO_Box_Param[idx_opt].setGeometry(QtCore.QRect(3300*sz_w, 130*sz_h, 420*sz_w, 250*sz_h)) #321
            self.PSO_Box_Param[idx_opt].setAlignment(QtCore.Qt.AlignCenter)
            self.PSO_Box_Param[idx_opt].setObjectName("PSO_Box_Param")
            self.PSO_Box_Param[idx_opt].setTitle(_translate("MainWindow", "PSO parameters"))

            self.label_pso_param_i.insert(idx_opt,QtWidgets.QLabel(self.PSO_Box_Param[idx_opt]))
            self.label_pso_param_i[idx_opt].setObjectName("label_pso_param_i1")
            self.label_pso_param_i[idx_opt].setText(_translate("MainWindow", "Initial value"))

            self.label_pso_param_e.insert(idx_opt,QtWidgets.QLabel(self.PSO_Box_Param[idx_opt]))
            self.label_pso_param_e[idx_opt].setObjectName("label_pso_param_e1")
            self.label_pso_param_e[idx_opt].setText(_translate("MainWindow", "Final value"))

            self.label_pso_param_min.insert(idx_opt,QtWidgets.QLabel(self.PSO_Box_Param[idx_opt]))
            self.label_pso_param_min[idx_opt].setObjectName("label_pso_param_im")
            self.label_pso_param_min[idx_opt].setText(_translate("MainWindow", "Min value"))

            self.label_pso_param_im.insert(idx_opt,QtWidgets.QLabel(self.PSO_Box_Param[idx_opt]))
            self.label_pso_param_im[idx_opt].setObjectName("label_pso_param_im")
            self.label_pso_param_im[idx_opt].setText(_translate("MainWindow", "Initial max value"))

            self.label_pso_param_fm.insert(idx_opt,QtWidgets.QLabel(self.PSO_Box_Param[idx_opt]))
            self.label_pso_param_fm[idx_opt].setObjectName("label_pso_param_im")
            self.label_pso_param_fm[idx_opt].setText(_translate("MainWindow", "Final max value"))


            self.label_Inertia.insert(idx_opt,QtWidgets.QLabel(self.PSO_Box_Param[idx_opt]))
            self.label_Inertia[idx_opt].setGeometry(QtCore.QRect(50*sz_w, 55*sz_h, 131*sz_w, 22*sz_h))
            self.label_Inertia[idx_opt].setObjectName("label_Inertia")
            self.label_Inertia[idx_opt].setText(_translate("MainWindow", "Inertia"))

            self.Inertia_score.insert(idx_opt,QtWidgets.QCheckBox(self.PSO_Box_Param[idx_opt]))
            self.Inertia_score[idx_opt].setGeometry(QtCore.QRect(20*sz_w, 75*sz_h, 121*sz_w, 22*sz_h))
            self.Inertia_score[idx_opt].setObjectName("Inertia_score")
            self.Inertia_score[idx_opt].setText(_translate("MainWindow", "Based on score"))

            self.Inertia_min.insert(idx_opt,QtWidgets.QDoubleSpinBox(self.PSO_Box_Param[idx_opt]))
            self.Inertia_min[idx_opt].setDecimals(2)
            self.Inertia_min[idx_opt].setMaximum(10.0)
            self.Inertia_min[idx_opt].setSingleStep(.1)
            self.Inertia_min[idx_opt].setObjectName("PSO_inertia_min")

            self.Inertia_i.insert(idx_opt,QtWidgets.QDoubleSpinBox(self.PSO_Box_Param[idx_opt]))
            self.Inertia_i[idx_opt].setDecimals(2)
            self.Inertia_i[idx_opt].setMaximum(10.0)
            self.Inertia_i[idx_opt].setSingleStep(.1)
            self.Inertia_i[idx_opt].setObjectName("PSO_inertia_i")

            self.Inertia_end.insert(idx_opt,QtWidgets.QDoubleSpinBox(self.PSO_Box_Param[idx_opt]))
            self.Inertia_end[idx_opt].setDecimals(2)
            self.Inertia_end[idx_opt].setMaximum(10.0)
            self.Inertia_end[idx_opt].setSingleStep(.1)
            self.Inertia_end[idx_opt].setObjectName("PSO_inertia_end")


            self.label_cogn.insert(idx_opt,QtWidgets.QLabel(self.PSO_Box_Param[idx_opt]))
            self.label_cogn[idx_opt].setGeometry(QtCore.QRect(20*sz_w, 165*sz_h, 131*sz_w, 22*sz_h))
            self.label_cogn[idx_opt].setObjectName("label_cogn")
            self.label_cogn[idx_opt].setText(_translate("MainWindow", "Cognitive component"))

            self.Cogn_i.insert(idx_opt,QtWidgets.QDoubleSpinBox(self.PSO_Box_Param[idx_opt]))
            self.Cogn_i[idx_opt].setGeometry(QtCore.QRect(170*sz_w, 160*sz_h, 50*sz_w, 32*sz_h))
            self.Cogn_i[idx_opt].setDecimals(2)
            self.Cogn_i[idx_opt].setMaximum(10.0)
            self.Cogn_i[idx_opt].setSingleStep(.1)
            self.Cogn_i[idx_opt].setObjectName("PSO_cogn_i")

            self.Cogn_end.insert(idx_opt,QtWidgets.QDoubleSpinBox(self.PSO_Box_Param[idx_opt]))
            self.Cogn_end[idx_opt].setGeometry(QtCore.QRect(260*sz_w, 160*sz_h, 50*sz_w, 32*sz_h))
            self.Cogn_end[idx_opt].setDecimals(2)
            self.Cogn_end[idx_opt].setMaximum(10.0)
            self.Cogn_end[idx_opt].setSingleStep(.1)
            self.Cogn_end[idx_opt].setObjectName("PSO_cogn_end")


            self.label_soc.insert(idx_opt,QtWidgets.QLabel(self.PSO_Box_Param[idx_opt]))
            self.label_soc[idx_opt].setGeometry(QtCore.QRect(20*sz_w, 205*sz_h, 131*sz_w, 22*sz_h))
            self.label_soc[idx_opt].setObjectName("label_cogn")
            self.label_soc[idx_opt].setText(_translate("MainWindow", "Social component"))

            self.Social_i.insert(idx_opt,QtWidgets.QDoubleSpinBox(self.PSO_Box_Param[idx_opt]))
            self.Social_i[idx_opt].setGeometry(QtCore.QRect(170*sz_w, 200*sz_h, 50*sz_w, 32*sz_h))
            self.Social_i[idx_opt].setDecimals(2)
            self.Social_i[idx_opt].setMaximum(10.0)
            self.Social_i[idx_opt].setSingleStep(.1)
            self.Social_i[idx_opt].setObjectName("PSO_social_i")

            self.Social_end.insert(idx_opt,QtWidgets.QDoubleSpinBox(self.PSO_Box_Param[idx_opt]))
            self.Social_end[idx_opt].setGeometry(QtCore.QRect(260*sz_w, 200*sz_h, 50*sz_w, 32*sz_h))
            self.Social_end[idx_opt].setDecimals(2)
            self.Social_end[idx_opt].setMaximum(10.0)
            self.Social_end[idx_opt].setSingleStep(.1)
            self.Social_end[idx_opt].setObjectName("PSO_social_end")

            # PSO default values
            self.Inertia_score[idx_opt].setChecked(d_PSO_inert_score)
            self.Inertia_min[idx_opt].setProperty("value", d_PSO_inert_min)
            self.Inertia_i[idx_opt].setProperty("value", d_PSO_inert_i)
            self.Inertia_end[idx_opt].setProperty("value", d_PSO_inert_end)
            if not d_PSO_inert_score:
                self.Inertia_i[idx_opt].setProperty("value", d_PSO_inert_i_ws)
                self.Inertia_end[idx_opt].setProperty("value", d_PSO_inert_end_ws)
            self.Cogn_i[idx_opt].setProperty("value", d_PSO_cogn_i)
            self.Cogn_end[idx_opt].setProperty("value", d_PSO_cogn_end)
            self.Social_i[idx_opt].setProperty("value", d_PSO_social_i)
            self.Social_end[idx_opt].setProperty("value", d_PSO_social_end)


            self.Inertia_score[idx_opt].toggled.connect(self.PSO_score_appear)



            self.opt_meth.insert(idx_opt,opt_meth)

            if opt_meth == 'GA':
                # show GA
                self.GA_Box_Selection[idx_opt].setGeometry(QtCore.QRect(330*sz_w, 20*sz_h, 321*sz_w, 100*sz_h))
                self.GA_Box_Xover[idx_opt].setGeometry(QtCore.QRect(330*sz_w, 130*sz_h, 420*sz_w, 190*sz_h))
                self.GA_Box_mut[idx_opt].setGeometry(QtCore.QRect(330*sz_w, 330*sz_h, 420*sz_w, 220*sz_h))
                self.tablet.setTabText(self.tablet.indexOf(self.GA[idx_opt]), _translate("MainWindow", "GA"))
                self.label_SA_espi_3[idx_opt].setText(_translate("MainWindow", "Generation number"))
                # remove PSO
                self.PSO_Box_Param[idx_opt].setGeometry(QtCore.QRect(33000*sz_w, 130*sz_h, 420*sz_w, 250*sz_h)) #321
            elif opt_meth == 'PSO':
                # show PSO
                self.PSO_Box_Param[idx_opt].setGeometry(QtCore.QRect(330*sz_w, 130*sz_h, 420*sz_w, 250*sz_h)) #321
                self.tablet.setTabText(self.tablet.indexOf(self.GA[idx_opt]), _translate("MainWindow", "PSO"))
                self.label_SA_espi_3[idx_opt].setText(_translate("MainWindow", "Iteration number"))
                # remove GA
                self.GA_Box_Selection[idx_opt].setGeometry(QtCore.QRect(3300*sz_w, 20*sz_h, 321*sz_w, 100*sz_h))
                self.GA_Box_Xover[idx_opt].setGeometry(QtCore.QRect(3300*sz_w, 130*sz_h, 420*sz_w, 190*sz_h))
                self.GA_Box_mut[idx_opt].setGeometry(QtCore.QRect(3300*sz_w, 330*sz_h, 420*sz_w, 220*sz_h))



            if self.Inertia_score[idx_opt].isChecked():
                self.Inertia_min[idx_opt].setGeometry(QtCore.QRect(260*sz_w, 30*sz_h, 50*sz_w, 30*sz_h))
                self.Inertia_i[idx_opt].setGeometry(QtCore.QRect(  260*sz_w, 60*sz_h, 50*sz_w, 30*sz_h))
                self.Inertia_end[idx_opt].setGeometry(QtCore.QRect(260*sz_w, 90*sz_h, 50*sz_w, 30*sz_h))
                self.label_pso_param_i[idx_opt].setGeometry(QtCore.QRect(170*sz_w, 130*sz_h, 131*sz_w, 22*sz_h))
                self.label_pso_param_e[idx_opt].setGeometry(QtCore.QRect(260*sz_w, 130*sz_h, 131*sz_w, 22*sz_h))
                self.label_pso_param_min[idx_opt].setGeometry(QtCore.QRect(170*sz_w, 35*sz_h, 131*sz_w, 22*sz_h))
                self.label_pso_param_im[idx_opt].setGeometry(QtCore.QRect(170*sz_w, 65*sz_h, 131*sz_w, 22*sz_h))
                self.label_pso_param_fm[idx_opt].setGeometry(QtCore.QRect(170*sz_w, 95*sz_h, 131*sz_w, 22*sz_h))
            else:
                self.Inertia_min[idx_opt].setGeometry(QtCore.QRect(1700*sz_w, 60*sz_h, 50*sz_w, 32*sz_h))
                self.Inertia_i[idx_opt].setGeometry(QtCore.QRect(170*sz_w, 60*sz_h, 50*sz_w, 32*sz_h))
                self.Inertia_end[idx_opt].setGeometry(QtCore.QRect(260*sz_w, 60*sz_h, 50*sz_w, 32*sz_h))
                self.label_pso_param_i[idx_opt].setGeometry(QtCore.QRect(170*sz_w, 30*sz_h, 131*sz_w, 22*sz_h))
                self.label_pso_param_e[idx_opt].setGeometry(QtCore.QRect(260*sz_w, 30*sz_h, 131*sz_w, 22*sz_h))
                self.label_pso_param_min[idx_opt].setGeometry(QtCore.QRect(1700*sz_w, 40*sz_h, 131*sz_w, 22*sz_h))
                self.label_pso_param_im[idx_opt].setGeometry(QtCore.QRect(1700*sz_w, 60*sz_h, 131*sz_w, 22*sz_h))
                self.label_pso_param_fm[idx_opt].setGeometry(QtCore.QRect(1700*sz_w, 80*sz_h, 131*sz_w, 22*sz_h))





            # -------------------
            # Fitness
            self.Gb_GA_fitness.insert(idx_opt,QtWidgets.QGroupBox(self.GA[idx_opt]))
            self.Gb_GA_fitness[idx_opt].setGeometry(QtCore.QRect(20*sz_w, 470*sz_h, 290*sz_w, 80*sz_h))
            self.Gb_GA_fitness[idx_opt].setAlignment(QtCore.Qt.AlignCenter)
            self.Gb_GA_fitness[idx_opt].setObjectName("Gb_GA_fitness")
            self.Gb_GA_fitness[idx_opt].setTitle(_translate("MainWindow", "Fitness calculation"))

            self.rB_GA_fit_1.insert(idx_opt,QtWidgets.QRadioButton(self.Gb_GA_fitness[idx_opt]))
            self.rB_GA_fit_1[idx_opt].setGeometry(QtCore.QRect(10*sz_w, 0, 230*sz_w, 80*sz_h))
            self.rB_GA_fit_1[idx_opt].setObjectName("rB_GA_fit_1")
            self.rB_GA_fit_1[idx_opt].setText(_translate("MainWindow", "average among all conditions"))
            self.rB_GA_fit_2.insert(idx_opt,QtWidgets.QRadioButton(self.Gb_GA_fitness[idx_opt]))
            self.rB_GA_fit_2[idx_opt].setGeometry(QtCore.QRect(10, 20*sz_h, 230, 80*sz_h))
            self.rB_GA_fit_2[idx_opt].setObjectName("rB_GA_fit_2")
            self.rB_GA_fit_2[idx_opt].setText(_translate("MainWindow", "maximum among all conditions"))
            if d_GA_fit == 'mean':  self.rB_GA_fit_1[idx_opt].setChecked(True)
            else:                   self.rB_GA_fit_2[idx_opt].setChecked(True)

            self.label_SA_deps_3[idx_opt].setText(_translate("MainWindow", "Individual number"))


            # Optimization on method
            self.GA_group5[idx_opt].setTitle(_translate("MainWindow", "Optimization on method"))
#            self.scrollArea_subM[idx_opt].setTitle(_translate("MainWindow", "Sub-Mechanism optimization "))
            if _idx != 0 and not op_null:
                self.cB_GA_meth[idx_opt].setText(_translate("MainWindow", "Optimization based on "+self.list_operator[_idx-2]))
                self.cB_GA_meth_DRG.insert(idx_opt,False)
                self.cB_GA_meth_SA.insert(idx_opt,False)
                self.cB_GA_meth_pts.insert(idx_opt,False)
                self.GA_label_meth_pts.insert(idx_opt,False)
            else:
                self.cB_GA_meth[idx_opt].setGeometry(QtCore.QRect(-110*sz_w, 30*sz_h, 231*sz_w, 22*sz_h))

                self.cB_GA_meth_DRG.insert(idx_opt,QtWidgets.QRadioButton(self.GA_group5[idx_opt]))
                self.cB_GA_meth_DRG[idx_opt].setGeometry(QtCore.QRect(10*sz_w, 30*sz_h, 231*sz_w, 22*sz_h))
                self.cB_GA_meth_DRG[idx_opt].setObjectName("cB_GA_meth_DRG")
                self.cB_GA_meth_DRG[idx_opt].setText(_translate("MainWindow", "DRG"))

                self.cB_GA_meth_SA.insert(idx_opt,QtWidgets.QRadioButton(self.GA_group5[idx_opt]))
                self.cB_GA_meth_SA[idx_opt].setGeometry(QtCore.QRect(70*sz_w, 30*sz_h, 231*sz_w, 22*sz_h))
                self.cB_GA_meth_SA[idx_opt].setObjectName("cB_GA_meth_SA")
                self.cB_GA_meth_SA[idx_opt].setText(_translate("MainWindow", "SA "))

                self.GA_label_meth_pts.insert(idx_opt,QtWidgets.QLabel(self.GA_group5[idx_opt]))
                self.GA_label_meth_pts[idx_opt].setGeometry(QtCore.QRect(150*sz_w, 30*sz_h, 90*sz_w, 18*sz_h))
                self.GA_label_meth_pts[idx_opt].setObjectName("label_Red_method_sp_5")
                self.GA_label_meth_pts[idx_opt].setText(_translate("MainWindow", "points: "))

                self.cB_GA_meth_pts.insert(idx_opt,QtWidgets.QDoubleSpinBox(self.GA_group5[idx_opt]))
                self.cB_GA_meth_pts[idx_opt].setGeometry(QtCore.QRect(215*sz_w, 25*sz_h, 61*sz_w, 32*sz_h))
                self.cB_GA_meth_pts[idx_opt].setDecimals(0)
                self.cB_GA_meth_pts[idx_opt].setMaximum(100.0)
                self.cB_GA_meth_pts[idx_opt].setSingleStep(1.0)
                self.cB_GA_meth_pts[idx_opt].setObjectName("cB_GA_meth_pts")
                self.cB_GA_meth_pts[idx_opt].setProperty("value", d_GA_meth_pts)


            self.GA_label_51[idx_opt].setText(_translate("MainWindow", "Number of reactions to optimize"))

            self.pB_GA_remove[idx_opt].setText(_translate("MainWindow", "Remove"))
            self.pB_GA_remove[idx_opt].clicked.connect(lambda: self.remove_tab('GA'))

            self.cB_Xover_op_1[idx_opt].setChecked(d_GA_Xover_op_1)
            self.cB_Xover_op_2[idx_opt].setChecked(d_GA_Xover_op_2)
            self.cB_Xover_op_3[idx_opt].setChecked(d_GA_Xover_op_3)
            self.cB_Xover_op_4[idx_opt].setChecked(d_GA_Xover_op_4)
            self.cB_mut_op_1[idx_opt].setChecked(d_GA_mut_op_1)
            self.cB_mut_op_2[idx_opt].setChecked(d_GA_mut_op_2)
            self.cB_mut_op_3[idx_opt].setChecked(d_GA_mut_op_3)
            self.cB_GA_meth[idx_opt].setChecked(d_GA_meth)

            if condition_tab_removed:   _idx = idx+1
            else: _idx = idx
            self.list_operator.insert(_idx-1,'opt')





    def add_cases(self):
        global sz_w
        global sz_h
        self.GCases.append(QtWidgets.QGroupBox(self.scrollAreaWidgetContents))
        self.GCases[-1].setGeometry(QtCore.QRect(10*sz_w, (10+(len(self.GCases)-1)*181)*sz_h, 171*sz_w, 171*sz_h))
        self.GCases[-1].setObjectName("GCases")
        self.num_case.append(len(self.GCases))
        self.rB_reactor_UV.append(QtWidgets.QRadioButton(self.GCases[-1]))
        self.rB_reactor_UV[-1].setGeometry(QtCore.QRect(10*sz_w, 25*sz_h, 105*sz_w, 22*sz_h))
        self.rB_reactor_UV[-1].setObjectName("rB_reactor_UV")
        self.rB_reactor_HP.append(QtWidgets.QRadioButton(self.GCases[-1]))
        self.rB_reactor_HP[-1].setGeometry(QtCore.QRect(10*sz_w, 43*sz_h, 105*sz_w, 22*sz_h))
        self.rB_reactor_HP[-1].setObjectName("rB_reactor_HP")
        self.rB_JSR.append(QtWidgets.QRadioButton(self.GCases[-1]))
        self.rB_JSR[-1].setGeometry(QtCore.QRect(10*sz_w, 61*sz_h, 105*sz_w, 22*sz_h))
        self.rB_JSR[-1].setObjectName("rB_JSR")
        self.rB_PFR.append(QtWidgets.QRadioButton(self.GCases[-1]))
        self.rB_PFR[-1].setGeometry(QtCore.QRect(10*sz_w, 79*sz_h, 151*sz_w, 22*sz_h))
        self.rB_PFR[-1].setObjectName("rB_PFR")
        self.rB_fflame.append(QtWidgets.QRadioButton(self.GCases[-1]))
        self.rB_fflame[-1].setGeometry(QtCore.QRect(10*sz_w, 97*sz_h, 105*sz_w, 22*sz_h))
        self.rB_fflame[-1].setObjectName("rB_fflame")
        self.rB_cfflame.append(QtWidgets.QRadioButton(self.GCases[-1]))
        self.rB_cfflame[-1].setGeometry(QtCore.QRect(10*sz_w, 115*sz_h, 151*sz_w, 22*sz_h))
        self.rB_cfflame[-1].setObjectName("rB_cfflame")

        self.pushButton.append(QtWidgets.QPushButton(self.GCases[-1]))
        self.pushButton[-1].setGeometry(QtCore.QRect(40*sz_w, 140*sz_h, 88*sz_w, 21*sz_h))
        self.pushButton[-1].setObjectName("pushButton")

        _translate = QtCore.QCoreApplication.translate
        self.GCases[-1].setTitle(_translate("MainWindow", "Cases"))
        self.rB_reactor_UV[-1].setText(_translate("MainWindow", "Reactor (&U,V)"))
        self.rB_reactor_HP[-1].setText(_translate("MainWindow", "Reactor (&H,p)"))
        self.rB_JSR[-1].setText(_translate("MainWindow", "J&SR"))
        self.rB_PFR[-1].setText(_translate("MainWindow", "&PFR"))
        self.rB_fflame[-1].setText(_translate("MainWindow", "&Flat flame"))
        self.rB_cfflame[-1].setText(_translate("MainWindow", "&Counter-flow flame"))
        self.pushButton[-1].setText(_translate("MainWindow", "add"))
        self.pushButton[-1].clicked.connect(lambda: self.add_conditions(self.num_case[-1]))



    def add_conditions(self,num_case):
        global sz_w
        global sz_h
        _translate = QtCore.QCoreApplication.translate

        self.tablet.setCurrentIndex(0)

        add_new_condition = True
        if not  self.rB_reactor_UV[num_case-1].isChecked() \
        and not self.rB_reactor_HP[num_case-1].isChecked() \
        and not self.rB_JSR[num_case-1].isChecked() \
        and not self.rB_fflame[num_case-1].isChecked()\
        and not self.rB_cfflame[num_case-1].isChecked()\
        and not self.rB_PFR[num_case-1].isChecked():
            print('Please, select a simulation configuration')
            add_new_condition = False

        if add_new_condition:
            if len(self.Condition_Gb)<num_case:
                # condition grp
                self.Condition_Gb.append(QtWidgets.QGroupBox(self.scrollAreaWidgetContents))
                self.condition_activated.append(True)
                self.txt8.append(QtWidgets.QLabel(self.Condition_Gb[-1]))

                self.eqmin.append(QtWidgets.QPlainTextEdit(self.Condition_Gb[-1]))
                self.eqmin[-1].setGeometry(QtCore.QRect(320*sz_w, 60*sz_h, 60*sz_w, 31*sz_h))
                self.eqmin[-1].setObjectName("eqmin")
                self.eqmin[-1].setPlainText(_translate("MainWindow", d_phi_min))
                self.eqmax.append(QtWidgets.QPlainTextEdit(self.Condition_Gb[-1]))
                self.eqmax[-1].setGeometry(QtCore.QRect(380*sz_w, 60*sz_h, 60*sz_w, 31*sz_h))
                self.eqmax[-1].setObjectName("eqmax")
                self.eqmax[-1].setPlainText(_translate("MainWindow", d_phi_max))
                self.eqincr.append(QtWidgets.QPlainTextEdit(self.Condition_Gb[-1]))
                self.eqincr[-1].setGeometry(QtCore.QRect(450*sz_w, 60*sz_h, 60*sz_w, 31*sz_h))
                self.eqincr[-1].setObjectName("eqincr")
                self.eqincr[-1].setPlainText(_translate("MainWindow", d_phi_incr))

                self.Tmin.append(QtWidgets.QPlainTextEdit(self.Condition_Gb[-1]))
                self.Tmin[-1].setGeometry(QtCore.QRect(320*sz_w, 120*sz_h, 60*sz_w, 31*sz_h))
                self.Tmin[-1].setObjectName("Tmin")
                self.Tmin[-1].setPlainText(_translate("MainWindow", "1500"))
                self.Tmax.append(QtWidgets.QPlainTextEdit(self.Condition_Gb[-1]))
                self.Tmax[-1].setGeometry(QtCore.QRect(380*sz_w, 120*sz_h, 60*sz_w, 31*sz_h))
                self.Tmax[-1].setObjectName("Tmax")
                self.Tmax[-1].setPlainText(_translate("MainWindow", "1700"))
                self.Tincr.append(QtWidgets.QPlainTextEdit(self.Condition_Gb[-1]))
                self.Tincr[-1].setGeometry(QtCore.QRect(450*sz_w, 120*sz_h, 60*sz_w, 31*sz_h))
                self.Tincr[-1].setObjectName("Tincr")
                self.Tincr[-1].setPlainText(_translate("MainWindow", "200"))

                font = QtGui.QFont()
                font.setPointSize(8)
                self.Pmin.append(QtWidgets.QPlainTextEdit(self.Condition_Gb[-1]))
                self.Pmin[-1].setGeometry(QtCore.QRect(320*sz_w, 90*sz_h, 60*sz_w, 31*sz_h))
                self.Pmin[-1].setObjectName("Pmin")
                self.Pmin[-1].setPlainText(_translate("MainWindow", d_P_min))
                self.Pmin[-1].setFont(font)
                self.Pmax.append(QtWidgets.QPlainTextEdit(self.Condition_Gb[-1]))
                self.Pmax[-1].setGeometry(QtCore.QRect(380*sz_w, 90*sz_h, 60*sz_w, 31*sz_h))
                self.Pmax[-1].setObjectName("Pmax")
                self.Pmax[-1].setFont(font)
                self.Pmax[-1].setPlainText(_translate("MainWindow", d_P_max))
                self.Pincr.append(QtWidgets.QPlainTextEdit(self.Condition_Gb[-1]))
                self.Pincr[-1].setGeometry(QtCore.QRect(450*sz_w, 90*sz_h, 60*sz_w, 31*sz_h))
                self.Pincr[-1].setObjectName("Pincr")
                self.Pincr[-1].setFont(font)
                self.Pincr[-1].setPlainText(_translate("MainWindow", d_P_incr))

                self.txt2.append(QtWidgets.QLabel(self.Condition_Gb[-1]))
                self.txt7.append(QtWidgets.QLabel(self.Condition_Gb[-1]))
                self.fuel_1.append(QtWidgets.QPlainTextEdit(self.Condition_Gb[-1]))
                self.Diluent_1.append(QtWidgets.QPlainTextEdit(self.Condition_Gb[-1]))
                self.txt9.append(QtWidgets.QLabel(self.Condition_Gb[-1]))
                self.txt6.append(QtWidgets.QLabel(self.Condition_Gb[-1]))
                self.txt10.append(QtWidgets.QLabel(self.Condition_Gb[-1]))
                self.txt1.append(QtWidgets.QLabel(self.Condition_Gb[-1]))
                self.txt5.append(QtWidgets.QLabel(self.Condition_Gb[-1]))
                self.Diluent_2.append(QtWidgets.QPlainTextEdit(self.Condition_Gb[-1]))
                self.txt4.append(QtWidgets.QLabel(self.Condition_Gb[-1]))
                self.txt3.append(QtWidgets.QLabel(self.Condition_Gb[-1]))
                self.oxidant_1.append(QtWidgets.QPlainTextEdit(self.Condition_Gb[-1]))

                # grp counterflow
                self.Conditions_diff_flame.append(QtWidgets.QGroupBox(self.scrollAreaWidgetContents))
                self.Conditions_diff_flame[-1].setObjectName("Conditions_diff_flame")
                self.Conditions_diff_flame[-1].setTitle(_translate("MainWindow", "Conditions"))

                self.Conditions_diff_flame_select.append(QtWidgets.QGroupBox(self.scrollAreaWidgetContents))
                self.Conditions_diff_flame_select[-1].setObjectName("Conditions_diff_flame_select")
                self.Conditions_diff_flame_select[-1].setTitle(_translate("MainWindow", "Configuration"))

                self.rB_cff_diff.append(QtWidgets.QRadioButton(self.Conditions_diff_flame_select[-1]))
                self.rB_cff_pp.append(QtWidgets.QRadioButton(self.Conditions_diff_flame_select[-1]))
                self.rB_cff_tp.append(QtWidgets.QRadioButton(self.Conditions_diff_flame_select[-1]))
                self.rB_cff_diff[-1].setObjectName("rB_cff_diff")
                self.rB_cff_pp[-1].setObjectName("rB_cff_pp")
                self.rB_cff_tp[-1].setObjectName("rB_cff_tp")
                self.rB_cff_diff[-1].setText(_translate("MainWindow", "Diffusion flame"))
                self.rB_cff_pp[-1].setText(_translate("MainWindow", "Partiall&y-premixed\nflame"))
                self.rB_cff_tp[-1].setText(_translate("MainWindow", "Premixed flame"))
                self.rB_cff_diff[-1].setGeometry(QtCore.QRect(10*sz_w, 20*sz_h, 141*sz_w, 22*sz_h))
                self.rB_cff_pp[-1].setGeometry(QtCore.QRect(10*sz_w, 70*sz_h, 141*sz_w, 41*sz_h))
                self.rB_cff_tp[-1].setGeometry(QtCore.QRect(10*sz_w, 120*sz_h, 141*sz_w, 22*sz_h))
                self.rB_cff_diff[-1].toggled.connect(lambda: self.add_cf_diff_options(num_case))
                self.rB_cff_pp[-1].toggled.connect(lambda: self.add_cf_pprem_options(num_case))
                self.rB_cff_tp[-1].toggled.connect(lambda: self.add_cf_tprem_options(num_case))
                # Pressure
                self.df_P_Gb.append(QtWidgets.QWidget(self.Conditions_diff_flame[-1]))
                self.txt6_2.append(QtWidgets.QLabel(self.df_P_Gb[-1]))
                self.df_Pmin.append(QtWidgets.QPlainTextEdit(self.df_P_Gb[-1]))
                self.df_Pmax.append(QtWidgets.QPlainTextEdit(self.df_P_Gb[-1]))
                self.df_Pincr.append(QtWidgets.QPlainTextEdit(self.df_P_Gb[-1]))
                self.txt10_5.append(QtWidgets.QLabel(self.df_P_Gb[-1]))
                self.txt9_4.append(QtWidgets.QLabel(self.df_P_Gb[-1]))
                self.txt8_4.append(QtWidgets.QLabel(self.df_P_Gb[-1]))
                self.df_P_Gb[-1].setObjectName("df_P_Gb")
                self.txt6_2[-1].setObjectName("txt6_2")
                self.df_Pmin[-1].setObjectName("df_Pmin")
                self.df_Pmax[-1].setObjectName("df_Pmax")
                self.df_Pincr[-1].setObjectName("df_Pincr")
                self.txt10_5[-1].setObjectName("txt10_5")
                self.txt9_4[-1].setObjectName("txt9_4")
                self.txt8_4[-1].setObjectName("txt8_4")
                self.df_P_Gb[-1].setGeometry(QtCore.QRect(5*sz_w, 20*sz_h, 111*sz_w, 141*sz_h))
                self.txt6_2[-1].setGeometry(QtCore.QRect(30*sz_w, 10*sz_h, 61*sz_w, 18*sz_h))
                self.df_Pmin[-1].setGeometry(QtCore.QRect(50*sz_w, 30*sz_h, 51*sz_w, 31*sz_h))
                self.df_Pmax[-1].setGeometry(QtCore.QRect(50*sz_w, 60*sz_h, 51*sz_w, 31*sz_h))
                self.df_Pincr[-1].setGeometry(QtCore.QRect(50*sz_w, 100*sz_h, 51*sz_w, 31*sz_h))
                self.txt10_5[-1].setGeometry(QtCore.QRect(10*sz_w, 110*sz_h, 31*sz_w, 18*sz_h))
                self.txt9_4[-1].setGeometry(QtCore.QRect(10*sz_w, 70*sz_h, 31*sz_w, 18*sz_h))
                self.txt8_4[-1].setGeometry(QtCore.QRect(10*sz_w, 40*sz_h, 31*sz_w, 18*sz_h))
                self.txt6_2[-1].setText(_translate("MainWindow", "Pressure (Pa)"))
                self.txt10_5[-1].setText(_translate("MainWindow", "incr"))
                self.txt9_4[-1].setText(_translate("MainWindow", "max"))
                self.txt8_4[-1].setText(_translate("MainWindow", "min"))
                self.df_Pmin[-1].setPlainText(_translate("MainWindow", "1e5"))
                self.df_Pmax[-1].setPlainText(_translate("MainWindow", "1e5"))
                self.df_Pincr[-1].setPlainText(_translate("MainWindow", "1e5"))

                # Lines
                self.line.append(QtWidgets.QFrame(self.Conditions_diff_flame[-1]))
                self.line[-1].setGeometry(QtCore.QRect(110*sz_w, 23*sz_h, 20*sz_w, 141*sz_h))
                self.line[-1].setFrameShape(QtWidgets.QFrame.VLine)
                self.line[-1].setFrameShadow(QtWidgets.QFrame.Sunken)
                self.line[-1].setObjectName("line")
                self.line_2.append(QtWidgets.QFrame(self.Conditions_diff_flame[-1]))
#                self.line_2[-1].setGeometry(QtCore.QRect(630, 23, 20, 141))
                self.line_2[-1].setFrameShape(QtWidgets.QFrame.VLine)
                self.line_2[-1].setFrameShadow(QtWidgets.QFrame.Sunken)
                self.line_2[-1].setObjectName("line_2")
                # Burner 1  -1
                self.df_G11.append(QtWidgets.QWidget(self.Conditions_diff_flame[-1]))
                self.df_txt_burn1.append(QtWidgets.QLabel(self.Conditions_diff_flame[-1]))
                self.df_txt_burn1[-1].setText(_translate("MainWindow", "Burner 1"))

                self.df_fuel_1.append(QtWidgets.QPlainTextEdit(self.df_G11[-1]))
                self.df_oxidant_1.append(QtWidgets.QPlainTextEdit(self.df_G11[-1]))
                self.df_Diluent_1.append(QtWidgets.QPlainTextEdit(self.df_G11[-1]))
                self.df_Diluent_r_1.append(QtWidgets.QPlainTextEdit(self.df_G11[-1]))
                self.df_txt.append(QtWidgets.QLabel(self.df_G11[-1]))
                self.df_txt_2.append(QtWidgets.QLabel(self.df_G11[-1]))
                self.df_txt_4.append(QtWidgets.QLabel(self.df_G11[-1]))
                self.df_txt_3.append(QtWidgets.QLabel(self.df_G11[-1]))
                self.df_G11[-1].setObjectName("df_G11")
                self.df_txt_burn1[-1].setObjectName("df_txt_burn1")
                self.df_txt[-1].setObjectName("df_txt")
                self.df_txt_2[-1].setObjectName("df_txt_2")
                self.df_txt_3[-1].setObjectName("df_txt_3")
                self.df_txt_4[-1].setObjectName("df_txt_4")
                self.df_fuel_1[-1].setObjectName("df_fuel_1")
                self.df_oxidant_1[-1].setObjectName("df_oxidant_1")
                self.df_Diluent_1[-1].setObjectName("df_Diluent_1")
                self.df_Diluent_r_1[-1].setObjectName("df_Diluent_r_1")
                self.df_txt[-1].setText(_translate("MainWindow", "Fuel"))
                self.df_txt_2[-1].setText(_translate("MainWindow", "Oxidant"))
                self.df_txt_4[-1].setText(_translate("MainWindow", "% diluent"))
                self.df_txt_3[-1].setText(_translate("MainWindow", "Diluent"))
                # Burner 1  -2
                self.df_G12.append(QtWidgets.QWidget(self.Conditions_diff_flame[-1]))
                self.txt5_2.append(QtWidgets.QLabel(self.df_G12[-1]))
                self.txt7_2.append(QtWidgets.QLabel(self.df_G12[-1]))
                self.txt7_3.append(QtWidgets.QLabel(self.df_G12[-1]))
                self.txt8_2.append(QtWidgets.QLabel(self.df_G12[-1]))
                self.txt9_2.append(QtWidgets.QLabel(self.df_G12[-1]))
                self.txt10_2.append(QtWidgets.QLabel(self.df_G12[-1]))
                self.df_eqmin_1.append(QtWidgets.QPlainTextEdit(self.df_G12[-1]))
                self.df_eqmax_1.append(QtWidgets.QPlainTextEdit(self.df_G12[-1]))
                self.df_eqincr_1.append(QtWidgets.QPlainTextEdit(self.df_G12[-1]))
                self.df_T_1.append(QtWidgets.QPlainTextEdit(self.df_G12[-1]))
                self.df_mdot1_1.append(QtWidgets.QPlainTextEdit(self.df_G12[-1]))
                self.df_mdot2_1.append(QtWidgets.QPlainTextEdit(self.df_G12[-1]))
                self.df_mdot3_1.append(QtWidgets.QLabel(self.df_G12[-1]))
                self.df_mdot4_1.append(QtWidgets.QPlainTextEdit(self.df_G12[-1]))
                self.df_G12[-1].setObjectName("df_G12")
                self.txt5_2[-1].setObjectName("txt5_2")
                self.txt7_2[-1].setObjectName("txt7_2")
                self.txt7_3[-1].setObjectName("txt7_3")
                self.txt8_2[-1].setObjectName("txt8_2")
                self.txt9_2[-1].setObjectName("txt9_2")
                self.txt10_2[-1].setObjectName("txt10_2")
                self.df_eqmin_1[-1].setObjectName("df_eqmin_1")
                self.df_eqmax_1[-1].setObjectName("df_eqmax_1")
                self.df_eqincr_1[-1].setObjectName("df_eqincr_1")
                self.df_T_1[-1].setObjectName("df_T_1")
                self.df_mdot1_1[-1].setObjectName("df_mdot1_1")
                self.df_mdot2_1[-1].setObjectName("df_mdot2_1")
                self.df_mdot3_1[-1].setAlignment(QtCore.Qt.AlignCenter)
                self.df_mdot3_1[-1].setObjectName("df_mdot3_1")
                self.df_mdot4_1[-1].setObjectName("df_mdot4_1")
                self.txt5_2[-1].setText(_translate("MainWindow", "phi"))
                self.txt7_2[-1].setText(_translate("MainWindow", "Tin"))
                self.txt7_3[-1].setText(_translate("MainWindow", "mdot"))
                self.txt8_2[-1].setText(_translate("MainWindow", "min"))
                self.txt9_2[-1].setText(_translate("MainWindow", "max"))
                self.txt10_2[-1].setText(_translate("MainWindow", "incr"))
                # Burner 2  -1
                self.df_G21.append(QtWidgets.QWidget(self.Conditions_diff_flame[-1]))
                self.df_txt_burn2.append(QtWidgets.QLabel(self.Conditions_diff_flame[-1]))
                self.df_txt_burn2[-1].setObjectName("df_txt_burn2")
                self.df_txt_burn2[-1].setText(_translate("MainWindow", "Burner 2"))
                self.df_fuel_2.append(QtWidgets.QPlainTextEdit(self.df_G21[-1]))
                self.df_oxidant_2.append(QtWidgets.QPlainTextEdit(self.df_G21[-1]))
                self.df_Diluent_2.append(QtWidgets.QPlainTextEdit(self.df_G21[-1]))
                self.df_Diluent_r_2.append(QtWidgets.QPlainTextEdit(self.df_G21[-1]))
                self.df_txt1_2.append(QtWidgets.QLabel(self.df_G21[-1]))
                self.df_txt2_2.append(QtWidgets.QLabel(self.df_G21[-1]))
                self.df_txt3_2.append(QtWidgets.QLabel(self.df_G21[-1]))
                self.df_txt4_2.append(QtWidgets.QLabel(self.df_G21[-1]))
                self.df_G21[-1].setObjectName("df_G21")
                self.df_txt1_2[-1].setObjectName("df_txt1_2")
                self.df_fuel_2[-1].setObjectName("df_fuel_2")
                self.df_oxidant_2[-1].setObjectName("df_oxidant_2")
                self.df_Diluent_2[-1].setObjectName("df_Diluent_2")
                self.df_Diluent_r_2[-1].setObjectName("df_Diluent_r_2")
                self.df_txt2_2[-1].setObjectName("df_txt2_2")
                self.df_txt3_2[-1].setObjectName("df_txt3_2")
                self.df_txt4_2[-1].setObjectName("df_txt4_2")
                self.df_txt1_2[-1].setText(_translate("MainWindow", "Fuel"))
                self.df_txt2_2[-1].setText(_translate("MainWindow", "Oxidant"))
                self.df_txt4_2[-1].setText(_translate("MainWindow", "Dilution ratio"))
                self.df_txt3_2[-1].setText(_translate("MainWindow", "Diluent"))
                # Burner 2  -2
                self.df_G22.append(QtWidgets.QWidget(self.Conditions_diff_flame[-1]))
                self.txt5_3.append(QtWidgets.QLabel(self.df_G22[-1]))
                self.txt7_4.append(QtWidgets.QLabel(self.df_G22[-1]))
                self.txt7_5.append(QtWidgets.QLabel(self.df_G22[-1]))
                self.txt8_3.append(QtWidgets.QLabel(self.df_G22[-1]))
                self.txt9_3.append(QtWidgets.QLabel(self.df_G22[-1]))
                self.txt10_3.append(QtWidgets.QLabel(self.df_G22[-1]))
                self.df_T_2.append(QtWidgets.QPlainTextEdit(self.df_G22[-1]))
                self.df_eqmin_2.append(QtWidgets.QPlainTextEdit(self.df_G22[-1]))
                self.df_eqmax_2.append(QtWidgets.QPlainTextEdit(self.df_G22[-1]))
                self.df_eqincr_2.append(QtWidgets.QLabel(self.df_G22[-1]))
                self.df_mdot1_2.append(QtWidgets.QPlainTextEdit(self.df_G22[-1]))
                self.df_mdot2_2.append(QtWidgets.QPlainTextEdit(self.df_G22[-1]))
                self.df_mdot3_2.append(QtWidgets.QLabel(self.df_G22[-1]))
                self.df_mdot4_2.append(QtWidgets.QPlainTextEdit(self.df_G22[-1]))
                self.df_G22[-1].setObjectName("df_G22")
                self.txt5_3[-1].setObjectName("txt5_3")
                self.txt7_4[-1].setObjectName("txt7_4")
                self.txt7_5[-1].setObjectName("txt7_5")
                self.txt8_3[-1].setObjectName("txt8_3")
                self.txt9_3[-1].setObjectName("txt9_3")
                self.txt10_3[-1].setObjectName("txt10_3")
                self.df_T_2[-1].setObjectName("df_T_2")
                self.df_eqmin_2[-1].setObjectName("df_eqmin_2")
                self.df_eqmax_2[-1].setObjectName("df_eqmax_2")
                self.df_eqincr_2[-1].setObjectName("df_eqincr_2")
                self.df_mdot1_2[-1].setObjectName("df_mdot1_2")
                self.df_mdot2_2[-1].setObjectName("df_mdot2_2")
                self.df_mdot3_2[-1].setObjectName("df_mdot3_2")
                self.df_mdot4_2[-1].setObjectName("df_mdot4_2")
                self.txt5_3[-1].setText(_translate("MainWindow", "phi"))
                self.txt7_4[-1].setText(_translate("MainWindow", "Tin"))
                self.txt7_5[-1].setText(_translate("MainWindow", "mdot"))
                self.txt8_3[-1].setText(_translate("MainWindow", "min"))
                self.txt9_3[-1].setText(_translate("MainWindow", "max"))
                self.txt10_3[-1].setText(_translate("MainWindow", "incr"))

                self.df_mdot1_1[-1].textChanged.connect(lambda:  self.mdot_incr_calc(num_case-1))
                self.df_mdot2_1[-1].textChanged.connect(lambda:  self.mdot_incr_calc(num_case-1))
                self.df_mdot1_2[-1].textChanged.connect(lambda:  self.mdot_incr_calc(num_case-1))
                self.df_mdot2_2[-1].textChanged.connect(lambda:  self.mdot_incr_calc(num_case-1))
                self.df_eqmin_1[-1].textChanged.connect(lambda:  self.mdot_incr_calc(num_case-1))
                self.df_eqmax_1[-1].textChanged.connect(lambda:  self.mdot_incr_calc(num_case-1))
                self.df_eqincr_1[-1].textChanged.connect(lambda: self.mdot_incr_calc(num_case-1))
                self.df_eqmin_2[-1].textChanged.connect(lambda:  self.mdot_incr_calc(num_case-1))
                self.df_eqmax_2[-1].textChanged.connect(lambda:  self.mdot_incr_calc(num_case-1))

                # add options
                self.add_reactor_options(num_case)
                self.add_jsr_options(num_case)
                self.add_flame_options(num_case)
                self.add_pfr_options(num_case)



            self.Condition_Gb[num_case-1].setGeometry(QtCore.QRect(180*sz_w, (10+(num_case-1)*181)*sz_h, 521*sz_w, 171*sz_h))
            self.Condition_Gb[num_case-1].setObjectName("groupBox")
            self.Conditions_diff_flame[num_case-1].setGeometry(QtCore.QRect(18000*sz_w, (10+(num_case-1)*181)*sz_h, 1311*sz_w, 171*sz_h))
            self.Conditions_diff_flame_select[num_case-1].setGeometry(QtCore.QRect(18000*sz_w, (10+(num_case-1)*181)*sz_h, 171*sz_w, 171*sz_h))


            if self.rB_reactor_UV[num_case-1].isChecked(): # reactor(U,V)
                self.Gr[num_case-1].setGeometry(QtCore.QRect(700*sz_w, (10+(num_case-1)*181)*sz_h, 731*sz_w, 171*sz_h))
                self.Tmin[num_case-1].setPlainText(_translate("MainWindow", d_T_min_r))
                self.Tmax[num_case-1].setPlainText(_translate("MainWindow", d_T_max_r))
                self.Tincr[num_case-1].setPlainText(_translate("MainWindow", d_T_incr_r))
                self.tol_ts_abs_r[num_case-1].setPlainText(_translate("MainWindow", d_tol_ts_r.split(',')[1]))
                self.tol_ts_rel_r[num_case-1].setPlainText(_translate("MainWindow", d_tol_ts_r.split(',')[0]))

            elif self.rB_reactor_HP[num_case-1].isChecked(): # reactor(H,P)
                self.Gr[num_case-1].setGeometry(QtCore.QRect(700*sz_w, (10+(num_case-1)*181)*sz_h, 731*sz_w, 171*sz_h))
                self.Tmin[num_case-1].setPlainText(_translate("MainWindow", d_T_min_r))
                self.Tmax[num_case-1].setPlainText(_translate("MainWindow", d_T_max_r))
                self.Tincr[num_case-1].setPlainText(_translate("MainWindow", d_T_incr_r))
                self.tol_ts_abs_r[num_case-1].setPlainText(_translate("MainWindow", d_tol_ts_r.split(',')[1]))
                self.tol_ts_rel_r[num_case-1].setPlainText(_translate("MainWindow", d_tol_ts_r.split(',')[0]))

            elif self.rB_JSR[num_case-1].isChecked(): # jsr
                self.GJSR[num_case-1].setGeometry(QtCore.QRect(700*sz_w, (10+(num_case-1)*181)*sz_h, 311*sz_w, 171*sz_h))
                self.Tmin[num_case-1].setPlainText(_translate("MainWindow", d_T_min_jsr))
                self.Tmax[num_case-1].setPlainText(_translate("MainWindow", d_T_max_jsr))
                self.Tincr[num_case-1].setPlainText(_translate("MainWindow",  d_T_incr_jsr))
                self.tol_ts_abs_jsr[num_case-1].setPlainText(_translate("MainWindow", d_tol_ts_jsr.split(',')[1]))
                self.tol_ts_rel_jsr[num_case-1].setPlainText(_translate("MainWindow", d_tol_ts_jsr.split(',')[0]))

            elif self.rB_PFR[num_case-1].isChecked(): # jsr
                self.Gpfr[num_case-1].setGeometry(QtCore.QRect(700*sz_w, (10+(num_case-1)*181)*sz_h, 810*sz_w, 171*sz_h))
                self.Tmin[num_case-1].setPlainText(_translate("MainWindow", d_T_min_pfr))
                self.Tmax[num_case-1].setPlainText(_translate("MainWindow", d_T_max_pfr))
                self.Tincr[num_case-1].setPlainText(_translate("MainWindow",  d_T_incr_pfr))
                self.tol_ts_abs_pfr[num_case-1].setPlainText(_translate("MainWindow", d_tol_ts_pfr.split(',')[1]))
                self.tol_ts_rel_pfr[num_case-1].setPlainText(_translate("MainWindow", d_tol_ts_pfr.split(',')[0]))

            elif self.rB_fflame[num_case-1].isChecked(): # flat flame
                self.Gf[num_case-1].setGeometry(QtCore.QRect(700*sz_w, (10+(num_case-1)*181)*sz_h, 1200*sz_w, 171*sz_h))
                self.Tmin[num_case-1].setPlainText(_translate("MainWindow", d_T_min_ff))
                self.Tmax[num_case-1].setPlainText(_translate("MainWindow", d_T_max_ff))
                self.Tincr[num_case-1].setPlainText(_translate("MainWindow", d_T_incr_ff))
                self.tol_ts_abs_f[num_case-1].setPlainText(_translate("MainWindow", d_tol_ts_ff.split(',')[1]))
                self.tol_ts_rel_f[num_case-1].setPlainText(_translate("MainWindow", d_tol_ts_ff.split(',')[0]))
                self.tol_ss_abs_f[num_case-1].setPlainText(_translate("MainWindow", d_tol_ss_ff.split(',')[1]))
                self.tol_ss_rel_f[num_case-1].setPlainText(_translate("MainWindow", d_tol_ss_ff.split(',')[0]))
                self.txtf8[num_case-1].setText(_translate("MainWindow", "xmax (m)"))


            elif self.rB_cfflame[num_case-1].isChecked(): # counterflow flame
                self.Condition_Gb[num_case-1].setGeometry(QtCore.QRect(10000*sz_w, (10+(num_case-1)*181)*sz_h, 521*sz_w, 171*sz_h))
                self.Conditions_diff_flame_select[num_case-1].setGeometry(QtCore.QRect(180*sz_w, (10+(num_case-1)*181)*sz_h, 171*sz_w, 171*sz_h))
                self.Tmin[num_case-1].setPlainText(_translate("MainWindow", d_T_min_ff))
                self.Tmax[num_case-1].setPlainText(_translate("MainWindow", d_T_max_ff))
                self.Tincr[num_case-1].setPlainText(_translate("MainWindow", d_T_incr_ff))
                self.tol_ts_abs_f[num_case-1].setPlainText(_translate("MainWindow", d_tol_ts_ff.split(',')[1]))
                self.tol_ts_rel_f[num_case-1].setPlainText(_translate("MainWindow", d_tol_ts_ff.split(',')[0]))
                self.tol_ss_abs_f[num_case-1].setPlainText(_translate("MainWindow", d_tol_ss_ff.split(',')[1]))
                self.tol_ss_rel_f[num_case-1].setPlainText(_translate("MainWindow", d_tol_ss_ff.split(',')[0]))
                self.txtf8[num_case-1].setText(_translate("MainWindow", "width (m)"))

            if num_case>=2:
                self.scrollAreaWidgetContents.setGeometry(QtCore.QRect(0, 0, 2500*sz_w, (15+181*(num_case+1))*sz_h))

            # properties of main conditions parameters
            self.condition_activated[num_case-1]=True
            self.txt8[num_case-1].setGeometry(QtCore.QRect(325*sz_w, 40*sz_h, 31*sz_w, 18*sz_h))
            self.txt8[num_case-1].setObjectName("txt8")
            self.txt2[num_case-1].setGeometry(QtCore.QRect(10*sz_w, 80*sz_h, 61*sz_w, 18*sz_h))
            self.txt2[num_case-1].setObjectName("txt2")
            self.txt7[num_case-1].setGeometry(QtCore.QRect(200*sz_w, 130*sz_h, 131*sz_w, 18*sz_h))
            self.txt7[num_case-1].setObjectName("txt7")
            self.fuel_1[num_case-1].setGeometry(QtCore.QRect(100*sz_w, 40*sz_h, 91*sz_w, 31*sz_h))
            self.fuel_1[num_case-1].setObjectName("fuel_1")
            self.Diluent_1[num_case-1].setGeometry(QtCore.QRect(100*sz_w, 100*sz_h, 91*sz_w, 31*sz_h))
            self.Diluent_1[num_case-1].setObjectName("Diluent_1")
            self.txt9[num_case-1].setGeometry(QtCore.QRect(390*sz_w, 40*sz_h, 31*sz_w, 18*sz_h))
            self.txt9[num_case-1].setObjectName("txt9")
            self.txt6[num_case-1].setGeometry(QtCore.QRect(200*sz_w, 100*sz_h, 131*sz_w, 18*sz_h))
            self.txt6[num_case-1].setObjectName("txt6")
            self.txt10[num_case-1].setGeometry(QtCore.QRect(465*sz_w, 40*sz_h, 31*sz_w, 18*sz_h))
            self.txt10[num_case-1].setObjectName("txt10")
            self.txt1[num_case-1].setGeometry(QtCore.QRect(10*sz_w, 50*sz_h, 41*sz_w, 18*sz_h))
            self.txt1[num_case-1].setObjectName("txt1")
            self.txt5[num_case-1].setGeometry(QtCore.QRect(200*sz_w, 70*sz_h, 131*sz_w, 18*sz_h))
            self.txt5[num_case-1].setObjectName("txt5")
            self.Diluent_2[num_case-1].setGeometry(QtCore.QRect(100*sz_w, 130*sz_h, 91*sz_w, 31*sz_h))
            self.Diluent_2[num_case-1].setObjectName("Diluent_2")
            self.txt4[num_case-1].setGeometry(QtCore.QRect(10*sz_w, 140*sz_h, 91*sz_w, 18*sz_h))
            self.txt4[num_case-1].setObjectName("txt4")
            self.txt3[num_case-1].setGeometry(QtCore.QRect(10*sz_w, 110*sz_h, 61*sz_w, 18*sz_h))
            self.txt3[num_case-1].setObjectName("txt3")
            self.oxidant_1[num_case-1].setGeometry(QtCore.QRect(100*sz_w, 70*sz_h, 91*sz_w, 31*sz_h))
            self.oxidant_1[num_case-1].setObjectName("oxidant_1")

            self.Condition_Gb[num_case-1].setTitle(_translate("MainWindow", "Conditions"))
            self.txt8[num_case-1].setText(_translate("MainWindow", "min"))
            self.txt2[num_case-1].setText(_translate("MainWindow", "Oxidant"))
            self.txt7[num_case-1].setText(_translate("MainWindow", "Temperature (K)"))
            self.fuel_1[num_case-1].setPlainText(_translate("MainWindow", d_fuel))
            self.Diluent_1[num_case-1].setPlainText(_translate("MainWindow", d_diluent))
            self.txt9[num_case-1].setText(_translate("MainWindow", "max"))
            self.txt6[num_case-1].setText(_translate("MainWindow", "Pressure (Pa)"))
            self.txt10[num_case-1].setText(_translate("MainWindow", "incr"))

            self.txt1[num_case-1].setText(_translate("MainWindow", "Fuel"))
            self.txt5[num_case-1].setText(_translate("MainWindow", "Equivalence ratio"))
            d_dil_ratio = str(d_diluent_ratio).replace('[','').replace(']','').replace(',','').replace("'","")
            self.Diluent_2[num_case-1].setPlainText(_translate("MainWindow", d_dil_ratio));
            self.txt4[num_case-1].setText(_translate("MainWindow", "Dilution ratio"))
            self.txt3[num_case-1].setText(_translate("MainWindow", "Diluent"))
            self.oxidant_1[num_case-1].setPlainText(_translate("MainWindow", d_oxidant))

            self.pushButton.append(QtWidgets.QPushButton(self.GCases[num_case-1]))
            self.pushButton[-1].setGeometry(QtCore.QRect(40*sz_w, 140*sz_h, 88*sz_w, 21*sz_h))
            self.pushButton[-1].setObjectName("pushButton")
            self.pushButton[-1].setText(_translate("MainWindow", "remove"))
            self.pushButton[-1].clicked.connect(lambda: self.remove_conditions(num_case))

            if num_case == len(self.GCases):
                self.add_cases()

        self.tablet.setCurrentIndex(1)


    def remove_conditions(self, num_case):
        global sz_w
        global sz_h
        self.tablet.setCurrentIndex(0)
        self.Condition_Gb[num_case-1].setGeometry(QtCore.QRect(180*sz_w, 10000*sz_h, 521*sz_w, 171*sz_h))
        self.Conditions_diff_flame[num_case-1].setGeometry(QtCore.QRect(2000*sz_w, 10000*sz_h, 670*sz_w, 171*sz_h))
        self.Conditions_diff_flame_select[num_case-1].setGeometry(QtCore.QRect(2000*sz_w, 10000*sz_h, 670*sz_w, 171*sz_h))
        self.Gr[num_case-1].setGeometry(QtCore.QRect(20000*sz_w, 10*sz_h, 731*sz_w, 171*sz_h))
        self.GJSR[num_case-1].setGeometry(QtCore.QRect(20000*sz_w, 10*sz_h, 311*sz_w, 171*sz_h))
        self.Gpfr[num_case-1].setGeometry(QtCore.QRect(20000*sz_w, 10*sz_h, 311*sz_w, 171*sz_h))
        self.Gf[num_case-1].setGeometry(QtCore.QRect(20000*sz_w, 10*sz_h, 1200*sz_w, 171*sz_h))
        self.condition_activated[num_case-1] = False

        self.pushButton.append(QtWidgets.QPushButton(self.GCases[num_case-1]))
        self.pushButton[-1].setGeometry(QtCore.QRect(40*sz_w, 140*sz_h, 88*sz_w, 21*sz_h))
        self.pushButton[-1].setObjectName("pushButton")

        _translate = QtCore.QCoreApplication.translate
        self.pushButton[-1].setText(_translate("MainWindow", "add"))
        self.pushButton[-1].clicked.connect(lambda: self.add_conditions(num_case))
        self.tablet.setCurrentIndex(1)


    def add_flame_options(self, num_case):
        global sz_w
        global sz_h
        _translate = QtCore.QCoreApplication.translate

        # Options for flat flames
        self.Gf.append(QtWidgets.QGroupBox(self.scrollAreaWidgetContents))
        self.Gf[num_case-1].setGeometry(QtCore.QRect(20000*sz_w, 10*sz_h, 1200*sz_w, 171*sz_h))
        self.Gf[num_case-1].setObjectName("Gf")
        self.Gf1.append(QtWidgets.QGroupBox(self.Gf[num_case-1]))
        self.Gf1[num_case-1].setGeometry(QtCore.QRect(10*sz_w, 30*sz_h, 161*sz_w, 131*sz_h))
        self.Gf1[num_case-1].setObjectName("Gf1")
        self.txtf3.append(QtWidgets.QLabel(self.Gf1[num_case-1]))
        self.txtf3[num_case-1].setGeometry(QtCore.QRect(10*sz_w, 70*sz_h, 71*sz_w, 18*sz_h))
        self.txtf3[num_case-1].setObjectName("txtf3")
        self.txtf1.append(QtWidgets.QLabel(self.Gf1[num_case-1]))
        self.txtf1[num_case-1].setGeometry(QtCore.QRect(60*sz_w, 40*sz_h, 31*sz_w, 18*sz_h))
        self.txtf1[num_case-1].setObjectName("txtf1")
        self.tol_ts_rel_f.append(QtWidgets.QPlainTextEdit(self.Gf1[num_case-1]))
        self.tol_ts_rel_f[num_case-1].setGeometry(QtCore.QRect(100*sz_w, 60*sz_h, 51*sz_w, 31*sz_h))
        self.tol_ts_rel_f[num_case-1].setObjectName("tol_ts_rel_f")
        self.tol_ts_abs_f.append(QtWidgets.QPlainTextEdit(self.Gf1[num_case-1]))
        self.tol_ts_abs_f[num_case-1].setGeometry(QtCore.QRect(50*sz_w, 60*sz_h, 51*sz_w, 31*sz_h))
        self.tol_ts_abs_f[num_case-1].setObjectName("tol_ts_abs_f")
        self.txtf4.append(QtWidgets.QLabel(self.Gf1[num_case-1]))
        self.txtf4[num_case-1].setGeometry(QtCore.QRect(10*sz_w, 100*sz_h, 51*sz_w, 18*sz_h))
        self.txtf4[num_case-1].setObjectName("txtf4")
        self.tol_ss_rel_f.append(QtWidgets.QPlainTextEdit(self.Gf1[num_case-1]))
        self.tol_ss_rel_f[num_case-1].setGeometry(QtCore.QRect(100*sz_w, 90*sz_h, 51*sz_w, 31*sz_h))
        self.tol_ss_rel_f[num_case-1].setObjectName("tol_ss_rel_f")
        self.txtf2.append(QtWidgets.QLabel(self.Gf1[num_case-1]))
        self.txtf2[num_case-1].setGeometry(QtCore.QRect(110*sz_w, 40*sz_h, 31*sz_w, 18*sz_h))
        self.txtf2[num_case-1].setObjectName("txtf2")
        self.tol_ss_abs_f.append(QtWidgets.QPlainTextEdit(self.Gf1[num_case-1]))
        self.tol_ss_abs_f[num_case-1].setGeometry(QtCore.QRect(50*sz_w, 90*sz_h, 51*sz_w, 31*sz_h))
        self.tol_ss_abs_f[num_case-1].setObjectName("tol_ss_abs_f")

        self.Gf2.append(QtWidgets.QGroupBox(self.Gf[num_case-1]))
        self.Gf2[num_case-1].setGeometry(QtCore.QRect(170*sz_w, 30*sz_h, 370*sz_w, 131*sz_h))
        self.Gf2[num_case-1].setObjectName("Gf2")
        self.txtf7.append(QtWidgets.QLabel(self.Gf2[num_case-1]))
        self.txtf7[num_case-1].setGeometry(QtCore.QRect(10*sz_w, 100*sz_h, 71*sz_w, 18*sz_h))
        self.txtf7[num_case-1].setObjectName("txtf7")
        self.txtf5.append(QtWidgets.QLabel(self.Gf2[num_case-1]))
        self.txtf5[num_case-1].setGeometry(QtCore.QRect(10*sz_w, 40*sz_h, 61*sz_w, 18*sz_h))
        self.txtf5[num_case-1].setObjectName("txtf5")
        self.txtf6.append(QtWidgets.QLabel(self.Gf2[num_case-1]))
        self.txtf6[num_case-1].setGeometry(QtCore.QRect(10*sz_w, 70*sz_h, 81*sz_w, 18*sz_h))
        self.txtf6[num_case-1].setObjectName("txtf6")

        self.txtf5b.append(QtWidgets.QLabel(self.Gf2[num_case-1]))
        self.txtf5b[num_case-1].setGeometry(QtCore.QRect(160*sz_w, 40*sz_h, 61*sz_w, 18*sz_h))
        self.txtf5b[num_case-1].setObjectName("txtf5b")
        self.txtf6b.append(QtWidgets.QLabel(self.Gf2[num_case-1]))
        self.txtf6b[num_case-1].setGeometry(QtCore.QRect(160*sz_w, 70*sz_h, 81*sz_w, 18*sz_h))
        self.txtf6b[num_case-1].setObjectName("txtf6b")

        self.txtf8.append(QtWidgets.QLabel(self.Gf2[num_case-1]))
        self.txtf8[num_case-1].setGeometry(QtCore.QRect(310*sz_w, 30*sz_h, 61*sz_w, 18*sz_h))
        self.txtf8[num_case-1].setObjectName("txtf8")
        self.txt_Gf_xmax.append(QtWidgets.QPlainTextEdit(self.Gf2[num_case-1]))
        self.txt_Gf_xmax[num_case-1].setGeometry(QtCore.QRect(300*sz_w, 50*sz_h, 51*sz_w, 31*sz_h))
        self.txt_Gf_xmax[num_case-1].setObjectName("txt_Gf_xmax")
        self.pts_scatter.append(QtWidgets.QPlainTextEdit(self.Gf2[num_case-1]))
        self.pts_scatter[num_case-1].setGeometry(QtCore.QRect(80*sz_w, 90*sz_h, 221*sz_w, 31*sz_h))
        self.pts_scatter[num_case-1].setObjectName("pts_scatter")

        self.slope_ff.append(QtWidgets.QDoubleSpinBox(self.Gf2[num_case-1]))
        self.slope_ff[num_case-1].setGeometry(QtCore.QRect(80*sz_w, 30*sz_h, 61*sz_w, 32*sz_h))
        self.slope_ff[num_case-1].setDecimals(2)
        self.slope_ff[num_case-1].setMaximum(1.0)
        self.slope_ff[num_case-1].setSingleStep(0.01)
        self.slope_ff[num_case-1].setProperty("value", d_slope_ff)
        self.slope_ff[num_case-1].setObjectName("slope_ff")
        self.curve_ff.append(QtWidgets.QDoubleSpinBox(self.Gf2[num_case-1]))
        self.curve_ff[num_case-1].setGeometry(QtCore.QRect(80*sz_w, 60*sz_h, 61*sz_w, 32*sz_h))
        self.curve_ff[num_case-1].setDecimals(2)
        self.curve_ff[num_case-1].setMaximum(1.0)
        self.curve_ff[num_case-1].setSingleStep(0.01)
        self.curve_ff[num_case-1].setProperty("value", d_curve_ff)
        self.curve_ff[num_case-1].setObjectName("curve_ff")

        self.ratio_ff.append(QtWidgets.QDoubleSpinBox(self.Gf2[num_case-1]))
        self.ratio_ff[num_case-1].setGeometry(QtCore.QRect(210*sz_w, 30*sz_h, 61*sz_w, 32*sz_h))
        self.ratio_ff[num_case-1].setDecimals(1)
        self.ratio_ff[num_case-1].setMaximum(10.0)
        self.ratio_ff[num_case-1].setSingleStep(.5)
        self.ratio_ff[num_case-1].setProperty("value", d_ratio_ff)
        self.ratio_ff[num_case-1].setObjectName("slope_ff")
        self.prune_ff.append(QtWidgets.QDoubleSpinBox(self.Gf2[num_case-1]))
        self.prune_ff[num_case-1].setGeometry(QtCore.QRect(210*sz_w, 60*sz_h, 61*sz_w, 32*sz_h))
        self.prune_ff[num_case-1].setDecimals(3)
        self.prune_ff[num_case-1].setMaximum(1.0)
        self.prune_ff[num_case-1].setSingleStep(0.01)
        self.prune_ff[num_case-1].setProperty("value", d_prune_ff)
        self.prune_ff[num_case-1].setObjectName("curve_ff")


        self.Gf3.append(QtWidgets.QGroupBox(self.Gf[num_case-1]))
        self.Gf3[num_case-1].setGeometry(QtCore.QRect(540*sz_w, 30*sz_h, 120*sz_w, 131*sz_h))
        self.Gf3[num_case-1].setObjectName("Gf3")
        self.mult.append(QtWidgets.QRadioButton(self.Gf3[num_case-1]))
        self.mult[num_case-1].setGeometry(QtCore.QRect(30*sz_w, 80*sz_h, 61*sz_w, 22*sz_h))
        self.mult[num_case-1].setObjectName("mult")
        self.mix.append(QtWidgets.QRadioButton(self.Gf3[num_case-1]))
        self.mix[num_case-1].setGeometry(QtCore.QRect(30*sz_w, 50*sz_h, 61*sz_w, 22*sz_h))
        self.mix[num_case-1].setObjectName("mix")

        # restore flame option
        self.Gf4.append(QtWidgets.QGroupBox(self.Gf[num_case-1]))
        self.Gf4[num_case-1].setGeometry(QtCore.QRect(670*sz_w, 30*sz_h, 500*sz_w, 131*sz_h))
        self.Gf4[num_case-1].setObjectName("Gf3")
        self.rB_flame_no_restore.append(QtWidgets.QRadioButton(self.Gf4[num_case-1]))
        self.rB_flame_no_restore[num_case-1].setGeometry(QtCore.QRect(30*sz_w, 80*sz_h, 120*sz_w, 22*sz_h))
        self.rB_flame_no_restore[num_case-1].setObjectName("no_restore")
        self.rB_flame_restore.append(QtWidgets.QRadioButton(self.Gf4[num_case-1]))
        self.rB_flame_restore[num_case-1].setGeometry(QtCore.QRect(30*sz_w, 50*sz_h, 120*sz_w, 22*sz_h))
        self.rB_flame_restore[num_case-1].setObjectName("restore")
        self.label_flame_folder.append(QtWidgets.QLabel(self.Gf4[num_case-1]))
        self.label_flame_folder[num_case-1].setGeometry(QtCore.QRect(170*sz_w, 65*sz_h, 230*sz_w, 18*sz_h))
        self.label_flame_folder[num_case-1].setObjectName("label_flame_folder")


        self.Gf[num_case-1].setTitle(_translate("MainWindow", "Options for flames"))
        self.Gf1[num_case-1].setTitle(_translate("MainWindow", "Tolerances"))
        self.txtf3[num_case-1].setText(_translate("MainWindow", "tol ts"))
        self.txtf1[num_case-1].setText(_translate("MainWindow", "abs"))
        self.txtf2[num_case-1].setText(_translate("MainWindow", "rel"))
        self.tol_ts_abs_f[num_case-1].setPlainText(_translate("MainWindow", d_tol_ts_ff.split(',')[0]))
        self.tol_ts_rel_f[num_case-1].setPlainText(_translate("MainWindow", d_tol_ts_ff.split(',')[1]))
        self.txtf4[num_case-1].setText(_translate("MainWindow", "tol ss"))
        self.tol_ss_abs_f[num_case-1].setPlainText(_translate("MainWindow", d_tol_ss_ff.split(',')[0]))
        self.tol_ss_rel_f[num_case-1].setPlainText(_translate("MainWindow", d_tol_ss_ff.split(',')[1]))
        self.Gf2[num_case-1].setTitle(_translate("MainWindow", "Meshing"))
        self.txtf5[num_case-1].setText(_translate("MainWindow", "slope"))
        self.txtf6[num_case-1].setText(_translate("MainWindow", "curve"))
        self.txtf5b[num_case-1].setText(_translate("MainWindow", "ratio"))
        self.txtf6b[num_case-1].setText(_translate("MainWindow", "prune"))
        self.txtf7[num_case-1].setText(_translate("MainWindow", "pts scatter"))


        self.txt_Gf_xmax[num_case-1].setPlainText(_translate("MainWindow", d_xmax))
        self.pts_scatter[num_case-1].setPlainText(_translate("MainWindow", d_pts_scatter))



        self.Gf3[num_case-1].setTitle(_translate("MainWindow", "Transport model"))
        self.mult[num_case-1].setText(_translate("MainWindow", "Mult"))
        self.mix[num_case-1].setText(_translate("MainWindow", "&Mix"))
        if   d_transport_model=='Mix':
            self.mix[-1].setChecked(True)
        elif d_transport_model=='Mult':
            self.mult[-1].setChecked(True)

        self.Gf4[num_case-1].setTitle(_translate("MainWindow", "Restore previous flame results"))
        self.rB_flame_restore[num_case-1].setText(_translate("MainWindow", "Restore flame "))
        self.rB_flame_no_restore[num_case-1].setText(_translate("MainWindow", "Start from scratch"))
        self.rB_flame_no_restore[num_case-1].setChecked(True)
#        self.rB_flame_no_restore[num_case-1].toggled.connect(lambda: self.get_restore_flame_folder(self.rB_flame_no_restore[num_case-1].isChecked(), num_case-1))
        self.rB_flame_restore[num_case-1].toggled.connect(lambda: self.get_restore_flame_folder(self.rB_flame_restore[num_case-1].isChecked(), num_case-1))



    def add_reactor_options(self,num_case):
        global sz_w
        global sz_h
        # Options for reactors
        _translate = QtCore.QCoreApplication.translate

        self.Gr.append(QtWidgets.QGroupBox(self.scrollAreaWidgetContents))
        self.Gr[num_case-1].setGeometry(QtCore.QRect(20000*sz_w, 10*sz_h, 862*sz_w, 171*sz_h))
        self.Gr[num_case-1].setObjectName("Gr")
        self.Gr[num_case-1].setTitle(_translate("MainWindow", "Options for reactor"))

        #tolerances
        self.Gr1.append(QtWidgets.QGroupBox(self.Gr[num_case-1]))
        self.Gr1[num_case-1].setGeometry(QtCore.QRect(10*sz_w, 30*sz_h, 161*sz_w, 131*sz_h))
        self.Gr1[num_case-1].setObjectName("Gr1")
        self.Gr1[num_case-1].setTitle(_translate("MainWindow", "Tolerances"))
        self.txtr3.append(QtWidgets.QLabel(self.Gr1[num_case-1]))
        self.txtr3[num_case-1].setGeometry(QtCore.QRect(10*sz_w, 70*sz_h, 71*sz_w, 18*sz_h))
        self.txtr3[num_case-1].setObjectName("txtr3")
        self.txtr3[num_case-1].setText(_translate("MainWindow", "tol ts"))
        self.tol_ts_abs_r.append(QtWidgets.QPlainTextEdit(self.Gr1[num_case-1]))
        self.tol_ts_abs_r[num_case-1].setGeometry(QtCore.QRect(50*sz_w, 60*sz_h, 51*sz_w, 31*sz_h))
        self.tol_ts_abs_r[num_case-1].setObjectName("tol_ts_abs_r")
        self.txtr2.append(QtWidgets.QLabel(self.Gr1[num_case-1]))
        self.txtr2[num_case-1].setGeometry(QtCore.QRect(110*sz_w, 40*sz_h, 31*sz_w, 18*sz_h))
        self.txtr2[num_case-1].setObjectName("txtr2")
        self.txtr2[num_case-1].setText(_translate("MainWindow", "rel"))
        self.tol_ts_rel_r.append(QtWidgets.QPlainTextEdit(self.Gr1[num_case-1]))
        self.tol_ts_rel_r[num_case-1].setGeometry(QtCore.QRect(100*sz_w, 60*sz_h, 51*sz_w, 31*sz_h))
        self.tol_ts_rel_r[num_case-1].setObjectName("tol_ts_rel_r")
        self.txtr1.append(QtWidgets.QLabel(self.Gr1[num_case-1]))
        self.txtr1[num_case-1].setGeometry(QtCore.QRect(60*sz_w, 40*sz_h, 31*sz_w, 18*sz_h))
        self.txtr1[num_case-1].setObjectName("txtr1")
        self.txtr1[num_case-1].setText(_translate("MainWindow", "abs"))

        self.tol_ts_abs_r[num_case-1].setPlainText(_translate("MainWindow", d_tol_ts_r.split(',')[0]))
        self.tol_ts_rel_r[num_case-1].setPlainText(_translate("MainWindow", d_tol_ts_r.split(',')[1]))

        # Time discretization
        self.Gr2.append(QtWidgets.QGroupBox(self.Gr[num_case-1]))
        self.Gr2[num_case-1].setGeometry(QtCore.QRect(170*sz_w, 30*sz_h, 351*sz_w, 131*sz_h))
        self.Gr2[num_case-1].setObjectName("Gr2")
        self.Gr2[num_case-1].setTitle(_translate("MainWindow", "Time discretization"))
        self.delta_n_pts.append(QtWidgets.QDoubleSpinBox(self.Gr2[num_case-1]))
        self.delta_n_pts[num_case-1].setGeometry(QtCore.QRect(80*sz_w, 60*sz_h, 61*sz_w, 32*sz_h))
        self.delta_n_pts[num_case-1].setDecimals(0)
        self.delta_n_pts[num_case-1].setMaximum(100.0)
        self.delta_n_pts[num_case-1].setSingleStep(5.0)
        self.delta_n_pts[num_case-1].setObjectName("delta_n_pts")
        self.txtr4.append(QtWidgets.QLabel(self.Gr2[num_case-1]))
        self.txtr4[num_case-1].setGeometry(QtCore.QRect(10*sz_w, 40*sz_h, 61*sz_w, 18*sz_h))
        self.txtr4[num_case-1].setObjectName("txtr4")
        self.txtr4[num_case-1].setText(_translate("MainWindow", "pt num"))
        self.pts_num.append(QtWidgets.QDoubleSpinBox(self.Gr2[num_case-1]))
        self.pts_num[num_case-1].setGeometry(QtCore.QRect(80*sz_w, 30*sz_h, 61*sz_w, 32*sz_h))
        self.pts_num[num_case-1].setDecimals(0)
        self.pts_num[num_case-1].setMaximum(1000.0)
        self.pts_num[num_case-1].setSingleStep(50.0)
        self.pts_num[num_case-1].setObjectName("pts_num")
        self.txtr5.append(QtWidgets.QLabel(self.Gr2[num_case-1]))
        self.txtr5[num_case-1].setGeometry(QtCore.QRect(10*sz_w, 70*sz_h, 81*sz_w, 18*sz_h))
        self.txtr5[num_case-1].setObjectName("txtr5")
        self.txtr5[num_case-1].setText(_translate("MainWindow", "delta n pts"))
        self.t_max_coeff.append(QtWidgets.QDoubleSpinBox(self.Gr2[num_case-1]))
        self.t_max_coeff[num_case-1].setGeometry(QtCore.QRect(280*sz_w, 30*sz_h, 61*sz_w, 32*sz_h))
        self.t_max_coeff[num_case-1].setDecimals(0)
        self.t_max_coeff[num_case-1].setMaximum(10.0)
        self.t_max_coeff[num_case-1].setSingleStep(1.0)
        self.t_max_coeff[num_case-1].setObjectName("t_max_coeff")
        self.txtr8.append(QtWidgets.QLabel(self.Gr2[num_case-1]))
        self.txtr8[num_case-1].setGeometry(QtCore.QRect(170*sz_w, 100*sz_h, 81*sz_w, 18*sz_h))
        self.txtr8[num_case-1].setObjectName("txtr8")
        self.txtr8[num_case-1].setText(_translate("MainWindow", "grad/curv"))
        self.txtr6.append(QtWidgets.QLabel(self.Gr2[num_case-1]))
        self.txtr6[num_case-1].setGeometry(QtCore.QRect(170*sz_w, 40*sz_h, 101*sz_w, 18*sz_h))
        self.txtr6[num_case-1].setObjectName("txtr6")
        self.txtr6[num_case-1].setText(_translate("MainWindow", "time max coeff"))
        self.grad_curv.append(QtWidgets.QDoubleSpinBox(self.Gr2[num_case-1]))
        self.grad_curv[num_case-1].setGeometry(QtCore.QRect(280*sz_w, 90*sz_h, 61*sz_w, 32*sz_h))
        self.grad_curv[num_case-1].setDecimals(2)
        self.grad_curv[num_case-1].setMaximum(1.0)
        self.grad_curv[num_case-1].setSingleStep(0.1)
        self.grad_curv[num_case-1].setObjectName("grad_curv")
        self.ref_spec.append(QtWidgets.QPlainTextEdit(self.Gr2[num_case-1]))
        self.ref_spec[num_case-1].setGeometry(QtCore.QRect(280*sz_w, 60*sz_h, 61*sz_w, 31*sz_h))
        self.ref_spec[num_case-1].setObjectName("ref_spec")
        self.txtr7.append(QtWidgets.QLabel(self.Gr2[num_case-1]))
        self.txtr7[num_case-1].setGeometry(QtCore.QRect(170*sz_w, 70*sz_h, 111*sz_w, 18*sz_h))
        self.txtr7[num_case-1].setObjectName("txtr7")
        self.txtr7[num_case-1].setText(_translate("MainWindow", "reference species"))

        self.pts_num[num_case-1].setProperty("value", d_n_pts_r)
        self.delta_n_pts[num_case-1].setProperty("value", d_delta_npts)
        self.t_max_coeff[num_case-1].setProperty("value", d_t_max_coeff)
        self.ref_spec[num_case-1].setPlainText(_translate("MainWindow", d_Scal_ref))
        self.grad_curv[num_case-1].setProperty("value", d_grad_curv_ratio)


        # Ignition detection parameters # depredicated
        self.Gr3.append(QtWidgets.QGroupBox(self.Gr[num_case-1]))
        self.Gr3[num_case-1].setGeometry(QtCore.QRect(2000*sz_w, 30*sz_h, 202*sz_w, 131*sz_h))
        self.Gr3[num_case-1].setObjectName("Gr3")
        self.Gr3[num_case-1].setTitle(_translate("MainWindow", "Ignition detection parameters"))
        self.txtr9.append(QtWidgets.QLabel(self.Gr3[num_case-1]))
        self.txtr9[num_case-1].setGeometry(QtCore.QRect(10*sz_w, 45*sz_h, 150*sz_w, 18*sz_h))
        self.txtr9[num_case-1].setObjectName("txtr9")
        self.txtr9[num_case-1].setText(_translate("MainWindow", "initial point number"))
        self.txtr10.append(QtWidgets.QLabel(self.Gr3[num_case-1]))
        self.txtr10[num_case-1].setGeometry(QtCore.QRect(10*sz_w, 90*sz_h, 150*sz_w, 18*sz_h))
        self.txtr10[num_case-1].setObjectName("txtr10")
        self.txtr10[num_case-1].setText(_translate("MainWindow", "initial time step"))
        self.txt_Gr_ipn.append(QtWidgets.QPlainTextEdit(self.Gr3[num_case-1]))
        self.txt_Gr_ipn[num_case-1].setGeometry(QtCore.QRect(140*sz_w, 35*sz_h, 51*sz_w, 31*sz_h))
        self.txt_Gr_ipn[num_case-1].setObjectName("txt_Gr_ipn")
        self.txt_Gr_its.append(QtWidgets.QPlainTextEdit(self.Gr3[num_case-1]))
        self.txt_Gr_its[num_case-1].setGeometry(QtCore.QRect(140*sz_w, 80*sz_h, 51*sz_w, 31*sz_h))
        self.txt_Gr_its[num_case-1].setObjectName("txt_Gr_its")

        self.txt_Gr_ipn[num_case-1].setPlainText(_translate("MainWindow", d_tign_nPoints))
        self.txt_Gr_its[num_case-1].setPlainText(_translate("MainWindow", d_tign_dt))

    def add_pfr_options(self,num_case):
        global sz_w
        global sz_h
        # Options for reactors
        _translate = QtCore.QCoreApplication.translate

        self.Gpfr.append(QtWidgets.QGroupBox(self.scrollAreaWidgetContents))
        self.Gpfr[num_case-1].setGeometry(QtCore.QRect(20000*sz_w, 10*sz_h, 810*sz_w, 171*sz_h))
        self.Gpfr[num_case-1].setObjectName("Gr")
        self.Gpfr[num_case-1].setTitle(_translate("MainWindow", "Options for PFR"))

        #tolerances
        self.Gpfr1.append(QtWidgets.QGroupBox(self.Gpfr[num_case-1]))
        self.Gpfr1[num_case-1].setGeometry(QtCore.QRect(10*sz_w, 30*sz_h, 161*sz_w, 131*sz_h))
        self.Gpfr1[num_case-1].setObjectName("Gpfr1")
        self.Gpfr1[num_case-1].setTitle(_translate("MainWindow", "Tolerances"))
        self.txtpfr3.append(QtWidgets.QLabel(self.Gpfr1[num_case-1]))
        self.txtpfr3[num_case-1].setGeometry(QtCore.QRect(10*sz_w, 70*sz_h, 71*sz_w, 18*sz_h))
        self.txtpfr3[num_case-1].setObjectName("txtpfr3")
        self.txtpfr3[num_case-1].setText(_translate("MainWindow", "tol ts"))
        self.tol_ts_abs_pfr.append(QtWidgets.QPlainTextEdit(self.Gpfr1[num_case-1]))
        self.tol_ts_abs_pfr[num_case-1].setGeometry(QtCore.QRect(50*sz_w, 60*sz_h, 51*sz_w, 31*sz_h))
        self.tol_ts_abs_pfr[num_case-1].setObjectName("tol_ts_abs_pfr")
        self.txtpfr2.append(QtWidgets.QLabel(self.Gpfr1[num_case-1]))
        self.txtpfr2[num_case-1].setGeometry(QtCore.QRect(110*sz_w, 40*sz_h, 31*sz_w, 18*sz_h))
        self.txtpfr2[num_case-1].setObjectName("txtpfr2")
        self.txtpfr2[num_case-1].setText(_translate("MainWindow", "rel"))
        self.tol_ts_rel_pfr.append(QtWidgets.QPlainTextEdit(self.Gpfr1[num_case-1]))
        self.tol_ts_rel_pfr[num_case-1].setGeometry(QtCore.QRect(100*sz_w, 60*sz_h, 51*sz_w, 31*sz_h))
        self.tol_ts_rel_pfr[num_case-1].setObjectName("tol_ts_rel_pfr")
        self.txtpfr1.append(QtWidgets.QLabel(self.Gpfr1[num_case-1]))
        self.txtpfr1[num_case-1].setGeometry(QtCore.QRect(60*sz_w, 40*sz_h, 31*sz_w, 18*sz_h))
        self.txtpfr1[num_case-1].setObjectName("txtpfr1")
        self.txtpfr1[num_case-1].setText(_translate("MainWindow", "abs"))

        self.tol_ts_abs_pfr[num_case-1].setPlainText(_translate("MainWindow", d_tol_ts_pfr.split(',')[0]))
        self.tol_ts_rel_pfr[num_case-1].setPlainText(_translate("MainWindow", d_tol_ts_pfr.split(',')[1]))

        # Plug Flow Reactor properties
        self.Gpfr2.append(QtWidgets.QGroupBox(self.Gpfr[num_case-1]))
        self.Gpfr2[num_case-1].setGeometry(QtCore.QRect(170*sz_w, 30*sz_h, 280*sz_w, 131*sz_h))
        self.Gpfr2[num_case-1].setObjectName("Gpfr2")
        self.Gpfr2[num_case-1].setTitle(_translate("MainWindow", "Plug Flow Reactor properties"))

        self.rB_pfr_1.append(QtWidgets.QRadioButton(self.Gpfr2[-1]))
        self.rB_pfr_1[-1].setGeometry(QtCore.QRect(10*sz_w, 20*sz_h, 180*sz_w, 18*sz_h))
        self.rB_pfr_1[-1].setObjectName("rB_pfr_1")
        self.rB_pfr_1[-1].setText(_translate("MainWindow", "Automatic time discretization"))
        self.rB_pfr_2.append(QtWidgets.QRadioButton(self.Gpfr2[-1]))
        self.rB_pfr_2[-1].setGeometry(QtCore.QRect(10*sz_w, 40*sz_h, 180*sz_w, 18*sz_h))
        self.rB_pfr_2[-1].setObjectName("rB_pfr_2")
        self.rB_pfr_2[-1].setText(_translate("MainWindow", "define length max"))


        self.txtpfr4.append(QtWidgets.QLabel(self.Gpfr2[num_case-1]))
        self.txtpfr4[num_case-1].setGeometry(QtCore.QRect(10*sz_w, 40*sz_h, 61*sz_w, 18*sz_h))
        self.txtpfr4[num_case-1].setObjectName("txtpfr4")
        self.txtpfr4[num_case-1].setText(_translate("MainWindow", "pt num"))
        self.pts_num_pfr.append(QtWidgets.QDoubleSpinBox(self.Gpfr2[num_case-1]))
        self.pts_num_pfr[num_case-1].setGeometry(QtCore.QRect(60*sz_w, 30*sz_h, 61*sz_w, 32*sz_h))
        self.pts_num_pfr[num_case-1].setDecimals(0)
        self.pts_num_pfr[num_case-1].setMaximum(10000.0)
        self.pts_num_pfr[num_case-1].setSingleStep(100.0)
        self.pts_num_pfr[num_case-1].setObjectName("pts_num_pfr")

        self.pfr_u_0.append(QtWidgets.QPlainTextEdit(self.Gpfr2[num_case-1]))
        self.pfr_u_0[num_case-1].setGeometry(QtCore.QRect(60*sz_w, 95*sz_h, 61*sz_w, 32*sz_h))
        self.pfr_u_0[num_case-1].setObjectName("pfr_length")
        self.txtpfr5.append(QtWidgets.QLabel(self.Gpfr2[num_case-1]))
        self.txtpfr5[num_case-1].setGeometry(QtCore.QRect(10*sz_w, 105*sz_h, 81*sz_w, 18*sz_h))
        self.txtpfr5[num_case-1].setObjectName("txtpfr5")
        self.txtpfr5[num_case-1].setText(_translate("MainWindow", "u_0"))

        self.pfr_length.append(QtWidgets.QPlainTextEdit(self.Gpfr2[num_case-1]))
        self.pfr_length[num_case-1].setGeometry(QtCore.QRect(210*sz_w, 70*sz_h, 61*sz_w, 32*sz_h))
        self.pfr_length[num_case-1].setObjectName("pfr_length")
        self.txtpfr6.append(QtWidgets.QLabel(self.Gpfr2[num_case-1]))
        self.txtpfr6[num_case-1].setGeometry(QtCore.QRect(150*sz_w, 80*sz_h, 101*sz_w, 18*sz_h))
        self.txtpfr6[num_case-1].setObjectName("txtpfr6")
        self.txtpfr6[num_case-1].setText(_translate("MainWindow", "length"))

        self.pfr_area.append(QtWidgets.QPlainTextEdit(self.Gpfr2[num_case-1]))
        self.pfr_area[num_case-1].setGeometry(QtCore.QRect(210*sz_w, 95*sz_h, 61*sz_w, 31*sz_h))
        self.pfr_area[num_case-1].setObjectName("pfr_area")
        self.txtpfr7.append(QtWidgets.QLabel(self.Gpfr2[num_case-1]))
        self.txtpfr7[num_case-1].setGeometry(QtCore.QRect(150*sz_w, 105*sz_h, 111*sz_w, 18*sz_h))
        self.txtpfr7[num_case-1].setObjectName("txtpfr7")
        self.txtpfr7[num_case-1].setText(_translate("MainWindow", "area"))

        self.pts_num_pfr[num_case-1].setProperty("value", d_n_pts_pfr)
        self.pfr_u_0[num_case-1].setPlainText(_translate("MainWindow", str(d_u_0_pfr)))
        self.pfr_length[num_case-1].setPlainText(_translate("MainWindow", str(d_length_pfr)))
        self.pfr_area[num_case-1].setPlainText(_translate("MainWindow", str(d_area_pfr)))


        # Time discretization
        self.Gpfr3.append(QtWidgets.QGroupBox(self.Gpfr[num_case-1]))
        self.Gpfr3[num_case-1].setGeometry(QtCore.QRect(460*sz_w, 30*sz_h, 351*sz_w, 131*sz_h))
        self.Gpfr3[num_case-1].setObjectName("Gpfr3")
        self.Gpfr3[num_case-1].setTitle(_translate("MainWindow", "Time discretization"))
        self.pfr_delta_n_pts.append(QtWidgets.QDoubleSpinBox(self.Gpfr3[num_case-1]))
        self.pfr_delta_n_pts[num_case-1].setGeometry(QtCore.QRect(80*sz_w, 60*sz_h, 61*sz_w, 32*sz_h))
        self.pfr_delta_n_pts[num_case-1].setDecimals(0)
        self.pfr_delta_n_pts[num_case-1].setMaximum(100.0)
        self.pfr_delta_n_pts[num_case-1].setSingleStep(5.0)
        self.pfr_delta_n_pts[num_case-1].setObjectName("pfr_delta_n_pts")
        self.txtpfr_8.append(QtWidgets.QLabel(self.Gpfr3[num_case-1]))
        self.txtpfr_8[num_case-1].setGeometry(QtCore.QRect(10*sz_w, 40*sz_h, 61*sz_w, 18*sz_h))
        self.txtpfr_8[num_case-1].setObjectName("txtpfr_8")
        self.txtpfr_8[num_case-1].setText(_translate("MainWindow", "pt num"))
        self.pfr_pts_num.append(QtWidgets.QDoubleSpinBox(self.Gpfr3[num_case-1]))
        self.pfr_pts_num[num_case-1].setGeometry(QtCore.QRect(80*sz_w, 30*sz_h, 61*sz_w, 32*sz_h))
        self.pfr_pts_num[num_case-1].setDecimals(0)
        self.pfr_pts_num[num_case-1].setMaximum(1000.0)
        self.pfr_pts_num[num_case-1].setSingleStep(50.0)
        self.pfr_pts_num[num_case-1].setObjectName("pfr_pts_num")
        self.txtpfr_9.append(QtWidgets.QLabel(self.Gpfr3[num_case-1]))
        self.txtpfr_9[num_case-1].setGeometry(QtCore.QRect(10*sz_w, 70*sz_h, 81*sz_w, 18*sz_h))
        self.txtpfr_9[num_case-1].setObjectName("txtpfr_9")
        self.txtpfr_9[num_case-1].setText(_translate("MainWindow", "delta n pts"))
        self.pfr_t_max_coeff.append(QtWidgets.QDoubleSpinBox(self.Gpfr3[num_case-1]))
        self.pfr_t_max_coeff[num_case-1].setGeometry(QtCore.QRect(280*sz_w, 30*sz_h, 61*sz_w, 32*sz_h))
        self.pfr_t_max_coeff[num_case-1].setDecimals(0)
        self.pfr_t_max_coeff[num_case-1].setMaximum(10.0)
        self.pfr_t_max_coeff[num_case-1].setSingleStep(1.0)
        self.pfr_t_max_coeff[num_case-1].setObjectName("pfr_t_max_coeff")
        self.txtpfr_10.append(QtWidgets.QLabel(self.Gpfr3[num_case-1]))
        self.txtpfr_10[num_case-1].setGeometry(QtCore.QRect(170*sz_w, 100*sz_h, 81*sz_w, 18*sz_h))
        self.txtpfr_10[num_case-1].setObjectName("txtpfr_10")
        self.txtpfr_10[num_case-1].setText(_translate("MainWindow", "grad/curv"))
        self.txtpfr_11.append(QtWidgets.QLabel(self.Gpfr3[num_case-1]))
        self.txtpfr_11[num_case-1].setGeometry(QtCore.QRect(170*sz_w, 40*sz_h, 101*sz_w, 18*sz_h))
        self.txtpfr_11[num_case-1].setObjectName("txtpfr_11")
        self.txtpfr_11[num_case-1].setText(_translate("MainWindow", "time max coeff"))
        self.pfr_grad_curv.append(QtWidgets.QDoubleSpinBox(self.Gpfr3[num_case-1]))
        self.pfr_grad_curv[num_case-1].setGeometry(QtCore.QRect(280*sz_w, 90*sz_h, 61*sz_w, 32*sz_h))
        self.pfr_grad_curv[num_case-1].setDecimals(2)
        self.pfr_grad_curv[num_case-1].setMaximum(1.0)
        self.pfr_grad_curv[num_case-1].setSingleStep(0.1)
        self.pfr_grad_curv[num_case-1].setObjectName("pfr_grad_curv")
        self.pfr_ref_spec.append(QtWidgets.QPlainTextEdit(self.Gpfr3[num_case-1]))
        self.pfr_ref_spec[num_case-1].setGeometry(QtCore.QRect(280*sz_w, 60*sz_h, 61*sz_w, 31*sz_h))
        self.pfr_ref_spec[num_case-1].setObjectName("pfr_ref_spec")
        self.txtpfr_12.append(QtWidgets.QLabel(self.Gpfr3[num_case-1]))
        self.txtpfr_12[num_case-1].setGeometry(QtCore.QRect(170*sz_w, 70*sz_h, 111*sz_w, 18*sz_h))
        self.txtpfr_12[num_case-1].setObjectName("txtpfr_12")
        self.txtpfr_12[num_case-1].setText(_translate("MainWindow", "reference species"))

        self.pfr_pts_num[num_case-1].setProperty("value", d_n_pts_r)
        self.pfr_delta_n_pts[num_case-1].setProperty("value", d_delta_npts)
        self.pfr_t_max_coeff[num_case-1].setProperty("value", d_t_max_coeff)
        self.pfr_ref_spec[num_case-1].setPlainText(_translate("MainWindow", d_Scal_ref))
        self.pfr_grad_curv[num_case-1].setProperty("value", d_grad_curv_ratio)




        self.rB_pfr_1[-1].toggled.connect(lambda: self.PFR_discr_appear(['in',num_case]))
        self.rB_pfr_2[-1].toggled.connect(lambda: self.PFR_discr_appear(['out',num_case]))
        if d_pfr_discr:    self.rB_pfr_1[-1].setChecked(True)
        else:           self.rB_pfr_2[-1].setChecked(True)


    def PFR_discr_appear(self,option):
        global sz_w
        global sz_h
#        idx_op=0
        appear   = option[0]
        num_case = option[1]
        if appear == 'in':
            self.Gpfr3[num_case-1].setGeometry(QtCore.QRect(449*sz_w, 30*sz_h, 351*sz_w, 131*sz_h))
            self.pfr_length[num_case-1].setGeometry(QtCore.QRect(2100*sz_w, 70*sz_h, 61*sz_w, 32*sz_h))
            self.txtpfr6[num_case-1].setGeometry(QtCore.QRect(1500*sz_w, 80*sz_h, 101*sz_w, 18*sz_h))
            self.txtpfr4[num_case-1].setGeometry(QtCore.QRect(1000*sz_w, 80*sz_h, 61*sz_w, 18*sz_h))
            self.pts_num_pfr[num_case-1].setGeometry(QtCore.QRect(6000*sz_w, 70*sz_h, 61*sz_w, 32*sz_h))
        else:
            self.Gpfr3[num_case-1].setGeometry(QtCore.QRect(4600*sz_w, 30*sz_h, 351*sz_w, 131*sz_h))
            self.pfr_length[num_case-1].setGeometry(QtCore.QRect(210*sz_w, 60*sz_h, 61*sz_w, 32*sz_h))
            self.txtpfr6[num_case-1].setGeometry(QtCore.QRect(150*sz_w, 70*sz_h, 101*sz_w, 18*sz_h))
            self.txtpfr4[num_case-1].setGeometry(QtCore.QRect(10*sz_w, 70*sz_h, 61*sz_w, 18*sz_h))
            self.pts_num_pfr[num_case-1].setGeometry(QtCore.QRect(60*sz_w, 60*sz_h, 61*sz_w, 32*sz_h))


    def add_jsr_options(self, num_case):
        global sz_w
        global sz_h
        # Options for JSR
        _translate = QtCore.QCoreApplication.translate

        self.GJSR.append(QtWidgets.QGroupBox(self.scrollAreaWidgetContents))
        self.GJSR[num_case-1].setGeometry(QtCore.QRect(20000*sz_w, 10*sz_h, 311*sz_w, 171*sz_h))
        self.GJSR[num_case-1].setObjectName("GJSR")
        self.GJSR[num_case-1].setTitle(_translate("MainWindow", "Options for JSR"))


        # tolerances
        self.GJSR1.append(QtWidgets.QGroupBox(self.GJSR[num_case-1]))
        self.GJSR1[num_case-1].setGeometry(QtCore.QRect(10*sz_w, 30*sz_h, 161*sz_w, 131*sz_h))
        self.GJSR1[num_case-1].setObjectName("GJSR1")
        self.GJSR1[num_case-1].setTitle(_translate("MainWindow", "Tolerances"))
        self.txtjsr_3.append(QtWidgets.QLabel(self.GJSR1[num_case-1]))
        self.txtjsr_3[num_case-1].setGeometry(QtCore.QRect(10*sz_w, 70*sz_h, 71*sz_w, 18*sz_h))
        self.txtjsr_3[num_case-1].setObjectName("txtjsr_3")
        self.tol_ts_abs_jsr.append(QtWidgets.QPlainTextEdit(self.GJSR1[num_case-1]))
        self.tol_ts_abs_jsr[num_case-1].setGeometry(QtCore.QRect(50*sz_w, 60*sz_h, 51*sz_w, 31*sz_h))
        self.tol_ts_abs_jsr[num_case-1].setObjectName("tol_ts_abs_jsr")
        self.txtjsr_2.append(QtWidgets.QLabel(self.GJSR1[num_case-1]))
        self.txtjsr_2[num_case-1].setGeometry(QtCore.QRect(110*sz_w, 40*sz_h, 31*sz_w, 18*sz_h))
        self.txtjsr_2[num_case-1].setObjectName("txtjsr_2")
        self.tol_ts_rel_jsr.append(QtWidgets.QPlainTextEdit(self.GJSR1[num_case-1]))
        self.tol_ts_rel_jsr[num_case-1].setGeometry(QtCore.QRect(100*sz_w, 60*sz_h, 51*sz_w, 31*sz_h))
        self.tol_ts_rel_jsr[num_case-1].setObjectName("tol_ts_rel_jsr")
        self.txtjsr_1.append(QtWidgets.QLabel(self.GJSR1[num_case-1]))
        self.txtjsr_1[num_case-1].setGeometry(QtCore.QRect(60*sz_w, 40*sz_h, 31*sz_w, 18*sz_h))
        self.txtjsr_1[num_case-1].setObjectName("txtjsr_1")
        self.txtjsr_3[num_case-1].setText(_translate("MainWindow", "tol ts"))
        self.txtjsr_2[num_case-1].setText(_translate("MainWindow", "rel"))
        self.txtjsr_1[num_case-1].setText(_translate("MainWindow", "abs"))
        self.tol_ts_abs_jsr[num_case-1].setPlainText(_translate("MainWindow", d_tol_ts_jsr.split(',')[0]))
        self.tol_ts_rel_jsr[num_case-1].setPlainText(_translate("MainWindow", d_tol_ts_jsr.split(',')[1]))

        # resident time
        self.GJSR2.append(QtWidgets.QGroupBox(self.GJSR[num_case-1]))
        self.GJSR2[num_case-1].setGeometry(QtCore.QRect(170*sz_w, 30*sz_h, 131*sz_w, 131*sz_h))
        self.GJSR2[num_case-1].setObjectName("GJSR2")
        self.GJSR2[num_case-1].setTitle(_translate("MainWindow", "Resident time"))
        self.t_max_jsr.append(QtWidgets.QDoubleSpinBox(self.GJSR2[num_case-1]))
        self.t_max_jsr[num_case-1].setGeometry(QtCore.QRect(40*sz_w, 50*sz_h, 61*sz_w, 32*sz_h))
        self.t_max_jsr[num_case-1].setDecimals(2)
        self.t_max_jsr[num_case-1].setMaximum(10.0)
        self.t_max_jsr[num_case-1].setSingleStep(0.1)
        self.t_max_jsr[num_case-1].setObjectName("t_max_jsr")
        self.t_max_jsr[num_case-1].setProperty("value", d_t_max)
        self.Diluent_2[num_case-1].setPlainText(_translate("MainWindow", d_diluent_ratio_jsr));


    def add_cf_pprem_options(self,num_case): # partially premixed flame conditions
        global sz_w
        global sz_h
        _translate = QtCore.QCoreApplication.translate

        self.Conditions_diff_flame[num_case-1].setGeometry(QtCore.QRect(351*sz_w, (10+(num_case-1)*181)*sz_h, 1161*sz_w, 171*sz_h))
        self.Gf[num_case-1].setGeometry(QtCore.QRect(1512*sz_w, (10+(num_case-1)*181)*sz_h, 1200*sz_w, 171*sz_h))
        # Burner 1  -1
        self.df_G11[num_case-1].setGeometry(QtCore.QRect(200*sz_w, 20*sz_h, 201*sz_w, 141*sz_h))
        self.df_txt_burn1[num_case-1].setGeometry(QtCore.QRect(130*sz_w, 30*sz_h, 61*sz_w, 18*sz_h))
        self.df_txt[num_case-1].setGeometry(QtCore.QRect(10*sz_w, 20*sz_h, 41*sz_w, 18*sz_h))
        self.df_txt_2[num_case-1].setGeometry(QtCore.QRect(10*sz_w, 50*sz_h, 61*sz_w, 18*sz_h))
        self.df_txt_3[num_case-1].setGeometry(QtCore.QRect(10*sz_w, 80*sz_h, 61*sz_w, 18*sz_h))
        self.df_txt_4[num_case-1].setGeometry(QtCore.QRect(10*sz_w, 110*sz_h, 91*sz_w, 18*sz_h))
        self.df_fuel_1[num_case-1].setGeometry(QtCore.QRect(100*sz_w, 10*sz_h, 91*sz_w, 31*sz_h))
        self.df_oxidant_1[num_case-1].setGeometry(QtCore.QRect(100*sz_w, 40*sz_h, 91*sz_w, 31*sz_h))
        self.df_Diluent_1[num_case-1].setGeometry(QtCore.QRect(100*sz_w, 70*sz_h, 91*sz_w, 31*sz_h))
        self.df_Diluent_r_1[num_case-1].setGeometry(QtCore.QRect(100*sz_w, 100*sz_h, 91*sz_w, 31*sz_h))
        self.df_fuel_1[num_case-1].setPlainText(_translate("MainWindow", "CH4"))
        self.df_Diluent_1[num_case-1].setPlainText(_translate("MainWindow", "N2"))
        self.df_oxidant_1[num_case-1].setPlainText(_translate("MainWindow", "O2"))
        d_dil_ratio = str(d_diluent_ratio).replace('[','').replace(']','').replace(',','').replace("'","")
        self.df_Diluent_r_1[num_case-1].setPlainText(_translate("MainWindow", d_dil_ratio))
        # Burner 1  -2
        self.df_G12[num_case-1].setGeometry(QtCore.QRect(400*sz_w, 20*sz_h, 231, 141*sz_h))
        self.txt10_2[num_case-1].setGeometry(QtCore.QRect(80*sz_w, 110*sz_h, 31, 18*sz_h))
        self.df_eqmin_1[num_case-1].setGeometry(QtCore.QRect(120*sz_w, 30*sz_h, 51, 31*sz_h))
        self.txt9_2[num_case-1].setGeometry(QtCore.QRect(80*sz_w, 70*sz_h, 31, 18*sz_h))
        self.txt5_2[num_case-1].setGeometry(QtCore.QRect(140*sz_w, 10*sz_h, 31, 18*sz_h))
        self.df_eqmax_1[num_case-1].setGeometry(QtCore.QRect(120*sz_w, 60*sz_h, 51, 31*sz_h))
        self.df_eqincr_1[num_case-1].setGeometry(QtCore.QRect(120*sz_w, 100*sz_h, 51, 31*sz_h))
        self.txt8_2[num_case-1].setGeometry(QtCore.QRect(80*sz_w, 40*sz_h, 31*sz_w, 18*sz_h))
        self.df_T_1[num_case-1].setGeometry(QtCore.QRect(10*sz_w, 30*sz_h, 51*sz_w, 31*sz_h))
        self.txt7_2[num_case-1].setGeometry(QtCore.QRect(20*sz_w, 10*sz_h, 21*sz_w, 18*sz_h))
        self.df_mdot2_1[num_case-1].setGeometry(QtCore.QRect(170*sz_w, 60*sz_h, 51*sz_w, 31*sz_h))
        self.df_mdot1_1[num_case-1].setGeometry(QtCore.QRect(170*sz_w, 30*sz_h, 51*sz_w, 31*sz_h))
        self.txt7_3[num_case-1].setGeometry(QtCore.QRect(180*sz_w, 10*sz_h, 41*sz_w, 20*sz_h))
        self.df_mdot3_1[num_case-1].setGeometry(QtCore.QRect(170*sz_w, 97*sz_h, 51*sz_w, 31*sz_h))
        self.df_mdot4_1[num_case-1].setGeometry(QtCore.QRect(170*sz_w, -97*sz_h, 51*sz_w, 31*sz_h))
        self.df_eqmin_1[num_case-1].setPlainText(_translate("MainWindow", "0.5"))
        self.df_eqmax_1[num_case-1].setPlainText(_translate("MainWindow", "1.5"))
        self.df_eqincr_1[num_case-1].setPlainText(_translate("MainWindow", ".5"))
        self.df_T_1[num_case-1].setPlainText(_translate("MainWindow", "300"))
        self.df_mdot1_1[num_case-1].setPlainText(_translate("MainWindow", "1"))
        self.df_mdot2_1[num_case-1].setPlainText(_translate("MainWindow", "3"))
        self.df_mdot3_1[num_case-1].setText(_translate("MainWindow", "1"))
        # Lines
        self.line_2[num_case-1].setGeometry(QtCore.QRect(630*sz_w, 23*sz_h, 20*sz_w, 141*sz_h))
        # Burner 2  -1
        self.df_G21[num_case-1].setGeometry(QtCore.QRect(720*sz_w, 20*sz_h, 201*sz_w, 141*sz_h))
        self.df_txt_burn2[num_case-1].setGeometry(QtCore.QRect(650*sz_w, 30*sz_h, 61*sz_w, 18*sz_h))
        self.df_fuel_2[num_case-1].setGeometry(QtCore.QRect(100*sz_w, 10*sz_h, 91*sz_w, 31*sz_h))
        self.df_oxidant_2[num_case-1].setGeometry(QtCore.QRect(100*sz_w, 40*sz_h, 91*sz_w, 31*sz_h))
        self.df_Diluent_2[num_case-1].setGeometry(QtCore.QRect(100*sz_w, 70*sz_h, 91*sz_w, 31*sz_h))
        self.df_Diluent_r_2[num_case-1].setGeometry(QtCore.QRect(100*sz_w, 100*sz_h, 91*sz_w, 31*sz_h))
        self.df_txt1_2[num_case-1].setGeometry(QtCore.QRect(10*sz_w, 20*sz_h, 41*sz_w, 18*sz_h))
        self.df_txt2_2[num_case-1].setGeometry(QtCore.QRect(10*sz_w, 50*sz_h, 61*sz_w, 18*sz_h))
        self.df_txt4_2[num_case-1].setGeometry(QtCore.QRect(10*sz_w, 110*sz_h, 91*sz_w, 18*sz_h))
        self.df_txt3_2[num_case-1].setGeometry(QtCore.QRect(10*sz_w, 80*sz_h, 61*sz_w, 18*sz_h))
        self.df_fuel_2[num_case-1].setPlainText(_translate("MainWindow", "CH4"))
        self.df_oxidant_2[num_case-1].setPlainText(_translate("MainWindow", "O2"))
        self.df_Diluent_2[num_case-1].setPlainText(_translate("MainWindow", "N2"))
        d_dil_ratio = str(d_diluent_ratio).replace('[','').replace(']','').replace(',','').replace("'","")
        self.df_Diluent_r_2[num_case-1].setPlainText(_translate("MainWindow", d_dil_ratio))
        # Burner 2  -2
        self.df_G22[num_case-1].setGeometry(QtCore.QRect(920*sz_w, 20*sz_h, 231*sz_w, 141*sz_h))
        self.df_T_2[num_case-1].setGeometry(QtCore.QRect(10*sz_w, 30*sz_h, 51*sz_w, 31*sz_h))
        self.df_eqmin_2[num_case-1].setGeometry(QtCore.QRect(120*sz_w, 30*sz_h, 51*sz_w, 31*sz_h))
        self.df_eqmax_2[num_case-1].setGeometry(QtCore.QRect(120*sz_w, 60*sz_h, 51*sz_w, 31*sz_h))
        self.df_eqincr_2[num_case-1].setGeometry(QtCore.QRect(120*sz_w, 100*sz_h, 51*sz_w, 31*sz_h))
        self.df_eqincr_2[num_case-1].setAlignment(QtCore.Qt.AlignCenter)
        self.df_mdot1_2[num_case-1].setGeometry(QtCore.QRect(170*sz_w, 30*sz_h, 51*sz_w, 31*sz_h))
        self.df_mdot2_2[num_case-1].setGeometry(QtCore.QRect(170*sz_w, 60*sz_h, 51*sz_w, 31*sz_h))
        self.df_mdot3_2[num_case-1].setGeometry(QtCore.QRect(170*sz_w, 100*sz_h, 51*sz_w, 31*sz_h))
        self.df_mdot3_2[num_case-1].setAlignment(QtCore.Qt.AlignCenter)
        self.df_mdot4_2[num_case-1].setGeometry(QtCore.QRect(120*sz_w, -100*sz_h, 51*sz_w, 31*sz_h))
        self.txt5_3[num_case-1].setGeometry(QtCore.QRect(140*sz_w, 10*sz_h, 31*sz_w, 18*sz_h))
        self.txt7_4[num_case-1].setGeometry(QtCore.QRect(20*sz_w, 10*sz_h, 21*sz_w, 18*sz_h))
        self.txt7_5[num_case-1].setGeometry(QtCore.QRect(180*sz_w, 10*sz_h, 41*sz_w, 20*sz_h))
        self.txt8_3[num_case-1].setGeometry(QtCore.QRect(80*sz_w, 40*sz_h, 31*sz_w, 18*sz_h))
        self.txt9_3[num_case-1].setGeometry(QtCore.QRect(80*sz_w, 70*sz_h, 31*sz_w, 18*sz_h))
        self.txt10_3[num_case-1].setGeometry(QtCore.QRect(80*sz_w, 110*sz_h, 31*sz_w, 18*sz_h))
        self.df_eqmin_2[num_case-1].setPlainText(_translate("MainWindow", "0.5"))
        self.df_eqmax_2[num_case-1].setPlainText(_translate("MainWindow", "1.5"))
        self.df_eqincr_2[num_case-1].setText(_translate("MainWindow", ".5"))
        self.df_T_2[num_case-1].setPlainText(_translate("MainWindow", "300"))
        self.df_mdot1_2[num_case-1].setPlainText(_translate("MainWindow", "1"))
        self.df_mdot2_2[num_case-1].setPlainText(_translate("MainWindow", "3"))
        self.df_mdot3_2[num_case-1].setText(_translate("MainWindow", "1"))

        self.pts_scatter[num_case-1].setPlainText(_translate("MainWindow", d_pts_scatter_cf))


    def add_cf_diff_options(self,num_case): # diffusion flame conditions
        global sz_w
        global sz_h
        _translate = QtCore.QCoreApplication.translate

        self.Conditions_diff_flame[num_case-1].setGeometry(QtCore.QRect(351*sz_w, (10+(num_case-1)*181)*sz_h, 1050*sz_w, 171*sz_h))
        self.Gf[num_case-1].setGeometry(QtCore.QRect(1401*sz_w, (10+(num_case-1)*181)*sz_h, 1200*sz_w, 171*sz_h))
        # Burner 1  -1
        self.df_G11[num_case-1].setGeometry(QtCore.QRect(200*sz_w, 20*sz_h, 201*sz_w, 141*sz_h))
        self.df_txt_burn1[num_case-1].setGeometry(QtCore.QRect(130*sz_w, 30*sz_h, 61*sz_w, 18*sz_h))
        self.df_fuel_1[num_case-1].setGeometry(QtCore.QRect(100*sz_w, 30*sz_h, 91*sz_w, 31*sz_h))
        self.df_Diluent_1[num_case-1].setGeometry(QtCore.QRect(100*sz_w, 60*sz_h, 91*sz_w, 31*sz_h))
        self.df_oxidant_1[num_case-1].setGeometry(QtCore.QRect(-100*sz_w, 40*sz_h, 91*sz_w, 31*sz_h))
        self.df_Diluent_r_1[num_case-1].setGeometry(QtCore.QRect(100*sz_w, 90*sz_h, 91*sz_w, 31*sz_h))
        self.df_txt[num_case-1].setGeometry(QtCore.QRect(10*sz_w, 40*sz_h, 41*sz_w, 18*sz_h))
        self.df_txt_2[num_case-1].setGeometry(QtCore.QRect(10*sz_w, -70*sz_h, 61*sz_w, 18*sz_h))
        self.df_txt_4[num_case-1].setGeometry(QtCore.QRect(10*sz_w, 100*sz_h, 91*sz_w, 18*sz_h))
        self.df_txt_3[num_case-1].setGeometry(QtCore.QRect(10*sz_w, 70*sz_h, 61*sz_w, 18*sz_h))
        self.df_fuel_1[num_case-1].setPlainText(_translate("MainWindow", "CH4"))
        self.df_Diluent_1[num_case-1].setPlainText(_translate("MainWindow", "N2"))
        self.df_oxidant_1[num_case-1].setPlainText(_translate("MainWindow", "O2"))
        self.df_Diluent_r_1[num_case-1].setPlainText(_translate("MainWindow", str(d_diluent_ratio_diff1)))
        # Burner 1  -2
        self.df_G12[num_case-1].setGeometry(QtCore.QRect(400*sz_w, 20*sz_h, 171*sz_w, 141*sz_h))
        self.df_eqmin_1[num_case-1].setGeometry(QtCore.QRect(120*sz_w, -300*sz_h, 51*sz_w, 31*sz_h))
        self.df_eqmax_1[num_case-1].setGeometry(QtCore.QRect(120*sz_w, -300*sz_h, 51*sz_w, 31*sz_h))
        self.df_eqincr_1[num_case-1].setGeometry(QtCore.QRect(120*sz_w, -300*sz_h, 51*sz_w, 31*sz_h))
        self.txt5_2[num_case-1].setGeometry(QtCore.QRect(140*sz_w, -100*sz_h, 31*sz_w, 18*sz_h))
        self.txt7_2[num_case-1].setGeometry(QtCore.QRect(20*sz_w, 10*sz_h, 21*sz_w, 18*sz_h))
        self.txt7_3[num_case-1].setGeometry(QtCore.QRect(120*sz_w, 10*sz_h, 41*sz_w, 20*sz_h))
        self.txt8_2[num_case-1].setGeometry(QtCore.QRect(80*sz_w, 40*sz_h, 31*sz_w, 18*sz_h))
        self.txt9_2[num_case-1].setGeometry(QtCore.QRect(80*sz_w, 70*sz_h, 31*sz_w, 18*sz_h))
        self.txt10_2[num_case-1].setGeometry(QtCore.QRect(80*sz_w, 110*sz_h, 31*sz_w, 18*sz_h))
        self.df_T_1[num_case-1].setGeometry(QtCore.QRect(10*sz_w, 30*sz_h, 51*sz_w, 31*sz_h))
        self.df_mdot1_1[num_case-1].setGeometry(QtCore.QRect(120*sz_w, 30*sz_h, 51*sz_w, 31*sz_h))
        self.df_mdot2_1[num_case-1].setGeometry(QtCore.QRect(120*sz_w, 60*sz_h, 51*sz_w, 31*sz_h))
        self.df_mdot3_1[num_case-1].setGeometry(QtCore.QRect(120*sz_w, -97*sz_h, 51*sz_w, 31*sz_h))
        self.df_mdot4_1[num_case-1].setGeometry(QtCore.QRect(120*sz_w, 100*sz_h, 51*sz_w, 31*sz_h))
        self.df_eqmin_1[num_case-1].setPlainText(_translate("MainWindow", "0.5"))
        self.df_eqmax_1[num_case-1].setPlainText(_translate("MainWindow", "1.5"))
        self.df_eqincr_1[num_case-1].setPlainText(_translate("MainWindow", ".5"))
        self.df_T_1[num_case-1].setPlainText(_translate("MainWindow", "300"))
        self.df_mdot1_1[num_case-1].setPlainText(_translate("MainWindow", "1"))
        self.df_mdot2_1[num_case-1].setPlainText(_translate("MainWindow", "3"))
        self.df_mdot3_1[num_case-1].setText(_translate("MainWindow", "1"))
        self.df_mdot4_1[num_case-1].setPlainText(_translate("MainWindow", "1"))
        # line
        self.line_2[num_case-1].setGeometry(QtCore.QRect(580*sz_w, 23*sz_h, 20*sz_w, 141*sz_h))
        # Burner 2  -1
        self.df_txt_burn2[num_case-1].setGeometry(QtCore.QRect(600*sz_w, 30*sz_h, 61*sz_w, 18*sz_h))
        self.df_G21[num_case-1].setGeometry(QtCore.QRect(660*sz_w, 20*sz_h, 201*sz_w, 141*sz_h))
        self.df_fuel_2[num_case-1].setGeometry(QtCore.QRect(-100*sz_w, 10*sz_h, 91*sz_w, 31*sz_h))
        self.df_Diluent_2[num_case-1].setGeometry(QtCore.QRect(100*sz_w, 60*sz_h, 91*sz_w, 31*sz_h))
        self.df_oxidant_2[num_case-1].setGeometry(QtCore.QRect(100*sz_w, 30*sz_h, 91*sz_w, 31*sz_h))
        self.df_Diluent_r_2[num_case-1].setGeometry(QtCore.QRect(100*sz_w, 90*sz_h, 91*sz_w, 31*sz_h))
        self.df_txt1_2[num_case-1].setGeometry(QtCore.QRect(10*sz_w, -200*sz_h, 41*sz_w, 18*sz_h))
        self.df_txt2_2[num_case-1].setGeometry(QtCore.QRect(10*sz_w, 40*sz_h, 61*sz_w, 18*sz_h))
        self.df_txt4_2[num_case-1].setGeometry(QtCore.QRect(10*sz_w, 100*sz_h, 91*sz_w, 18*sz_h))
        self.df_txt3_2[num_case-1].setGeometry(QtCore.QRect(10*sz_w, 70*sz_h, 61*sz_w, 18*sz_h))
        self.df_oxidant_2[num_case-1].setPlainText(_translate("MainWindow", "O2"))
        self.df_Diluent_2[num_case-1].setPlainText(_translate("MainWindow", "N2"))
        self.df_Diluent_r_2[num_case-1].setPlainText(_translate("MainWindow", str(d_diluent_ratio_diff2)))
        # Burner 2  -2
        self.df_G22[num_case-1].setGeometry(QtCore.QRect(860*sz_w, 20*sz_h, 171*sz_w, 141*sz_h))
        self.df_eqmin_2[num_case-1].setGeometry(QtCore.QRect(120*sz_w, -50*sz_h, 51*sz_w, 31*sz_h))
        self.df_eqmax_2[num_case-1].setGeometry(QtCore.QRect(120*sz_w, -60*sz_h, 51*sz_w, 31*sz_h))
        self.df_eqincr_2[num_case-1].setGeometry(QtCore.QRect(120*sz_w, -100*sz_h, 51*sz_w, 31*sz_h))
        self.df_T_2[num_case-1].setGeometry(QtCore.QRect(10*sz_w, 30*sz_h, 51*sz_w, 31*sz_h))
        self.txt5_3[num_case-1].setGeometry(QtCore.QRect(140*sz_w, -100*sz_h, 31*sz_w, 18*sz_h))
        self.txt7_4[num_case-1].setGeometry(QtCore.QRect(20*sz_w, 10*sz_h, 21*sz_w, 18*sz_h))
        self.txt7_5[num_case-1].setGeometry(QtCore.QRect(130*sz_w, 10*sz_h, 41*sz_w, 20*sz_h))
        self.txt8_3[num_case-1].setGeometry(QtCore.QRect(80*sz_w, 40*sz_h, 31*sz_w, 18*sz_h))
        self.txt9_3[num_case-1].setGeometry(QtCore.QRect(80*sz_w, 70*sz_h, 31*sz_w, 18*sz_h))
        self.txt10_3[num_case-1].setGeometry(QtCore.QRect(80*sz_w, 110*sz_h, 31*sz_w, 18*sz_h))
        self.df_mdot1_2[num_case-1].setGeometry(QtCore.QRect(120*sz_w, 30*sz_h, 51*sz_w, 31*sz_h))
        self.df_mdot2_2[num_case-1].setGeometry(QtCore.QRect(120*sz_w, 60*sz_h, 51*sz_w, 31*sz_h))
        self.df_mdot3_2[num_case-1].setGeometry(QtCore.QRect(120*sz_w, -100*sz_h, 51*sz_w, 31*sz_h))
        self.df_mdot4_2[num_case-1].setGeometry(QtCore.QRect(120*sz_w, 100*sz_h, 51*sz_w, 31*sz_h))
        self.df_mdot3_2[num_case-1].setAlignment(QtCore.Qt.AlignCenter)
        self.df_eqmin_2[num_case-1].setPlainText(_translate("MainWindow", "0.5"))
        self.df_eqmax_2[num_case-1].setPlainText(_translate("MainWindow", "1.5"))
        self.df_eqincr_2[num_case-1].setText(_translate("MainWindow", ".5"))
        self.df_T_2[num_case-1].setPlainText(_translate("MainWindow", "300"))
        self.df_mdot1_2[num_case-1].setPlainText(_translate("MainWindow", "1"))
        self.df_mdot2_2[num_case-1].setPlainText(_translate("MainWindow", "3"))
        self.df_mdot3_2[num_case-1].setText(_translate("MainWindow", "1"))
        self.df_mdot4_2[num_case-1].setPlainText(_translate("MainWindow", "1"))

        self.pts_scatter[num_case-1].setPlainText(_translate("MainWindow", d_pts_scatter_cf))


    def add_cf_tprem_options(self,num_case):
        global sz_w
        global sz_h
        _translate = QtCore.QCoreApplication.translate

        self.Conditions_diff_flame[num_case-1].setGeometry(QtCore.QRect(351*sz_w, (10+(num_case-1)*181)*sz_h, 570*sz_w, 171*sz_h))
        self.Gf[num_case-1].setGeometry(QtCore.QRect(921*sz_w, (10+(num_case-1)*181)*sz_h, 1200*sz_w, 171*sz_h))
        # Burner 1  -1
        self.df_txt_burn1[num_case-1].setGeometry(QtCore.QRect(10*sz_w, -50*sz_h, 61*sz_w, 18*sz_h))
        self.df_G11[num_case-1].setGeometry(QtCore.QRect(130*sz_w, 20*sz_h, 201*sz_w, 141*sz_h))
        self.df_fuel_1[num_case-1].setGeometry(QtCore.QRect(100*sz_w, 10*sz_h, 91*sz_w, 31*sz_h))
        self.df_Diluent_1[num_case-1].setGeometry(QtCore.QRect(100*sz_w, 70*sz_h, 91*sz_w, 31*sz_h))
        self.df_oxidant_1[num_case-1].setGeometry(QtCore.QRect(100*sz_w, 40*sz_h, 91*sz_w, 31*sz_h))
        self.df_Diluent_r_1[num_case-1].setGeometry(QtCore.QRect(100*sz_w, 100*sz_h, 91*sz_w, 31*sz_h))
        self.df_txt[num_case-1].setGeometry(QtCore.QRect(10*sz_w, 20*sz_h, 41*sz_w, 18*sz_h))
        self.df_txt_2[num_case-1].setGeometry(QtCore.QRect(10*sz_w, 50*sz_h, 61*sz_w, 18*sz_h))
        self.df_txt_4[num_case-1].setGeometry(QtCore.QRect(10*sz_w, 110*sz_h, 91*sz_w, 18*sz_h))
        self.df_txt_3[num_case-1].setGeometry(QtCore.QRect(10*sz_w, 80*sz_h, 61*sz_w, 18*sz_h))
        self.df_fuel_1[num_case-1].setPlainText(_translate("MainWindow", "CH4"))
        self.df_Diluent_1[num_case-1].setPlainText(_translate("MainWindow", "N2"))
        self.df_oxidant_1[num_case-1].setPlainText(_translate("MainWindow", "O2"))
        self.df_Diluent_r_1[num_case-1].setPlainText(_translate("MainWindow", "O2/N2  3.76"))
        # Burner 1  -2
        self.df_G12[num_case-1].setGeometry(QtCore.QRect(330*sz_w, 20*sz_h, 231*sz_w, 141*sz_h))
        self.df_T_1[num_case-1].setGeometry(QtCore.QRect(10*sz_w, 30*sz_h, 51*sz_w, 31*sz_h))
        self.df_eqmin_1[num_case-1].setGeometry(QtCore.QRect(120*sz_w, 30*sz_h, 51*sz_w, 31*sz_h))
        self.df_eqmax_1[num_case-1].setGeometry(QtCore.QRect(120*sz_w, 60*sz_h, 51*sz_w, 31*sz_h))
        self.df_eqincr_1[num_case-1].setGeometry(QtCore.QRect(120*sz_w, 100*sz_h, 51*sz_w, 31*sz_h))
        self.df_mdot1_1[num_case-1].setGeometry(QtCore.QRect(170*sz_w, 30*sz_h, 51*sz_w, 31*sz_h))
        self.df_mdot2_1[num_case-1].setGeometry(QtCore.QRect(170*sz_w, 60*sz_h, 51*sz_w, 31*sz_h))
        self.df_mdot3_1[num_case-1].setGeometry(QtCore.QRect(170*sz_w, 97*sz_h, 51*sz_w, 31*sz_h))
        self.df_mdot4_1[num_case-1].setGeometry(QtCore.QRect(170*sz_w, -97*sz_h, 51*sz_w, 31*sz_h))
        self.txt5_2[num_case-1].setGeometry(QtCore.QRect(140*sz_w, 10*sz_h, 31*sz_w, 18*sz_h))
        self.txt7_2[num_case-1].setGeometry(QtCore.QRect(20*sz_w, 10*sz_h, 21*sz_w, 18*sz_h))
        self.txt7_3[num_case-1].setGeometry(QtCore.QRect(180*sz_w, 10*sz_h, 41*sz_w, 20*sz_h))
        self.txt8_2[num_case-1].setGeometry(QtCore.QRect(80*sz_w, 40*sz_h, 31*sz_w, 18*sz_h))
        self.txt9_2[num_case-1].setGeometry(QtCore.QRect(80*sz_w, 70*sz_h, 31*sz_w, 18*sz_h))
        self.txt10_2[num_case-1].setGeometry(QtCore.QRect(80*sz_w, 110*sz_h, 31*sz_w, 18*sz_h))
        self.df_T_1[num_case-1].setPlainText(_translate("MainWindow", "300"))
        self.df_eqmin_1[num_case-1].setPlainText(_translate("MainWindow", "0.5"))
        self.df_eqmax_1[num_case-1].setPlainText(_translate("MainWindow", "1.5"))
        self.df_eqincr_1[num_case-1].setPlainText(_translate("MainWindow", ".5"))
        self.df_mdot1_1[num_case-1].setPlainText(_translate("MainWindow", "1"))
        self.df_mdot2_1[num_case-1].setPlainText(_translate("MainWindow", "3"))
        self.df_mdot3_1[num_case-1].setText(_translate("MainWindow", "1"))
        # Lines
        self.line_2[num_case-1].setGeometry(QtCore.QRect(-570*sz_w, 23*sz_h, 20*sz_w, 141*sz_h))
        # Burner 2  -1
        self.df_G21[num_case-1].setGeometry(QtCore.QRect(-870*sz_w, 20*sz_h, 201*sz_w, 141*sz_h))
        self.df_txt_burn2[num_case-1].setGeometry(QtCore.QRect(800*sz_w, 30*sz_h, 61*sz_w, 18*sz_h))
        # Burner 2  -2
        self.df_G22[num_case-1].setGeometry(QtCore.QRect(-1070*sz_w, 20*sz_h, 231*sz_w, 141*sz_h))

        self.pts_scatter[num_case-1].setPlainText(_translate("MainWindow", d_pts_scatter_cf))



    def write_parameters(self, get_fn):
        global sz_w
        global sz_h
        os.chdir(self.main_dir)


# =============================================================================
#         Get Values
# =============================================================================

# ==================   Main parameters   ======================================

        main_path       = self.Text_folder_name.document().toPlainText()
        mech            = self.reference_mechanism.split('/')[-1]
        mech_prev_red   = self.reduced_mechanism
        verbose         = int(self.MP_verbose.text())
        show_plots      = self.cB_show_plots.isChecked()
        write_ck        = self.cB_chemkin_f.isChecked()

        # targets
        tspc = []
        for tg in range(self.list_target.count()):
            tspc.append(self.list_target.item(tg).text())
        T_check  = self.cB_tsp_T.isChecked()   ; sp_T  = self.text_MP_T.document().toPlainText().split(' ')
        Sl_check = self.cB_tsp_Sl.isChecked()  ; sp_Sl = self.text_MP_Sl.document().toPlainText().split(' ')
        ig_check = self.cB_tsp_igt.isChecked() ; sp_ig = self.text_MP_ig.document().toPlainText().split(' ')
        K_check  = self.cB_tsp_K.isChecked()   ; sp_K  = self.text_MP_K.document().toPlainText().split(' ')
        # error calculation
        if self.rB_MP_qoi_1.isChecked() :   error_calculation = 'points'
        else:                               error_calculation = 'QoI'
        if self.rB_MP_errmean_1.isChecked() :  error_coupling = 'mean'
        else:                                  error_coupling = 'max'




# ====================   Operators   =========================================

        eps                 = []
        delta_eps           = []
        npoints             = []
        max_error_sp        = []
        max_error_T         = []
        max_error_ig        = []
        max_error_Sl        = []
        max_error_K         = []
        optim               = []
        ttol_sensi          = []
        isi                 = []

        n_gen               = []
        n_indiv             = []
        Arrh_max_variation  = []
        optim_on_meth       = []
        nb_r2opt            = []
        sub_mech_sel        = []

        opt_meth            = []

        # genetic Algorithm
        selection_operator  = []
        selection_options   = []
        Xover_operator      = []
        Xover_pct           = []
        Xover_opt           = []
        mut_operator        = []
        mut_pct             = []
        mut_opt             = []
        mut_intensity       = []

        # PSO
        inertia             = []
        cogn_comp           = []
        social_comp         = []

        error_fitness       = []

        drg = -1 ; sa = -1 ; ga = -1


        for red_i in range(len(self.list_operator)):
            if self.list_operator[red_i] != 'opt':

                if 'DRG' in self.list_operator[red_i]:
                    drg+=1
                    # eps
                    eps.append(float(self.num_DRG_eps[drg].text().replace(',','.')))
                    delta_eps.append(float(self.num_DRG_deps[drg].text().replace(',','.')))
                    # npts
                    npoints.append(float(self.num_DRG_pt_num[drg].text().replace(',','.')))
                    # errors
                    max_error_sp.append([])
                    for r in range(self.tableWidget_DRG[drg].rowCount()):
                        if self.tableWidget_DRG[drg].item(r, 0).text()=='T':
                            max_error_T.append(self.tableWidget_DRG[drg].item(r, 1).text())
                        elif self.tableWidget_DRG[drg].item(r, 0).text()=='Ig_t':
                            max_error_ig.append(self.tableWidget_DRG[drg].item(r, 1).text())
                        elif self.tableWidget_DRG[drg].item(r, 0).text()=='Sl':
                            max_error_Sl.append(self.tableWidget_DRG[drg].item(r, 1).text())
                        elif self.tableWidget_DRG[drg].item(r, 0).text()=='K':
                            max_error_K.append(self.tableWidget_DRG[drg].item(r, 1).text())
                        else:
                            max_error_sp[-1].append(self.tableWidget_DRG[drg].item(r, 1).text())
                    # Inter-species interations
                    isi.append(self.cB_ISI_drg[sa].isChecked())
                    # optimization
                    if red_i+1<len(self.list_operator):
                        if self.list_operator[red_i+1]=='opt': optim.append(True)
                        else:                                 optim.append(False)
                    else:                                     optim.append(False)
                    ttol_sensi.append([])

                if 'SA' in self.list_operator[red_i] or 'DSRG' in self.list_operator[red_i] :
                    sa+=1
                    # eps
                    eps.append(float(self.num_SA_eps[sa].text().replace(',','.')))
                    delta_eps.append(float(self.num_SA_deps[sa].text().replace(',','.')))
                    # npts
                    npoints.append(float(self.num_SA_pt_num[sa].text().replace(',','.')))
                    # errors
                    max_error_sp.append([])
                    for r in range(self.tableWidget_SA[sa].rowCount()):
                        if self.tableWidget_SA[sa].item(r, 0).text()=='T':
                            max_error_T.append(self.tableWidget_SA[sa].item(r, 1).text())
                        elif self.tableWidget_SA[sa].item(r, 0).text()=='Ig_t':
                            max_error_ig.append(self.tableWidget_SA[sa].item(r, 1).text())
                        elif self.tableWidget_SA[sa].item(r, 0).text()=='Sl':
                            max_error_Sl.append(self.tableWidget_SA[sa].item(r, 1).text())
                        elif self.tableWidget_SA[sa].item(r, 0).text()=='K':
                            max_error_K.append(self.tableWidget_SA[sa].item(r, 1).text())
                        else:
                            max_error_sp[-1].append(self.tableWidget_SA[sa].item(r, 1).text())
                    # Inter-species interations
                    isi.append(self.cB_ISI_sa[sa].isChecked())
                    # optimization
                    if red_i+1<len(self.list_operator):
                        if self.list_operator[red_i+1]=='opt': optim.append(True)
                        else:                                 optim.append(False)
                    else:                                     optim.append(False)
                    # tolerances
                    if self.rB_SA_tol_1[sa].isChecked():
                        ttol_sensi.append([False,False])
                    else:
                        ttol_sensi.append([float(self.num_SA_tol_1[sa].document().toPlainText().replace(',','.')),\
                                           float(self.num_SA_tol_2[sa].document().toPlainText().replace(',','.'))])

            else: #GA
                ga+=1
                # fitness calculation
                if self.rB_GA_fit_1[ga].isChecked(): error_fitness.append('mean')
                else:                                error_fitness.append('max')

                n_gen.append(int(self.num_GA_gen[ga].text()))
                n_indiv.append(int(self.num_GA_ind[ga].text()))
                Arrh_max_variation.append([int(self.num_GA_A[ga].text())])
                Arrh_max_variation[-1].append(int(self.num_GA_n[ga].text()))
                Arrh_max_variation[-1].append(int(self.num_GA_Ea[ga].text()))

                optim_on_meth.append(self.cB_GA_meth[ga].isChecked())
                nb_r2opt.append(int(self.num_GA_meth_fract[ga].text()))

                sub_mech_sel.append([])
                if self.cB_sub_H[ga].isChecked():
                    sub_mech_sel[-1].append('H2')
                if self.cB_sub_CO[ga].isChecked():
                    sub_mech_sel[-1].append('CO')
                for sC in range(len(self.cB_sub_C[ga])):
                    if self.cB_sub_C[ga][sC].isChecked():
                        sub_mech_sel[-1].append('C'+str(sC+1))
                if self.cB_sub_N[ga].isChecked():
                    sub_mech_sel[-1].append('N')
                if self.cB_sub_S[ga].isChecked():
                    sub_mech_sel[-1].append('S')
                if self.cB_sub_Si[ga].isChecked():
                    sub_mech_sel[-1].append('Si')

                opt_meth.append(self.opt_meth[ga])

                # Genetic Algorithm options
                selection_operator.append(self.Box_GA_Selection[ga].currentText())
                selection_options.append([float(self.txt_GA_sel_opt[ga].document().toPlainText())])
                Xover_operator.append([]) ; Xover_pct.append([]) ; Xover_opt.append([])
                if self.cB_Xover_op_1[ga].isChecked(): # 1- simpleXover
                    Xover_operator[-1].append('simple_Xover')
                    Xover_pct[-1].append(int(self.txt_GA_Xover_int_1[ga].document().toPlainText()))
                    Xover_opt[-1].append(self.txt_GA_Xover_opt_1[ga].document().toPlainText())
                if self.cB_Xover_op_2[ga].isChecked(): # 2- multipleXover
                    Xover_operator[-1].append('multiple_Xover')
                    Xover_pct[-1].append(int(self.txt_GA_Xover_int_2[ga].document().toPlainText()))
                    Xover_opt[-1].append(self.txt_GA_Xover_opt_2[ga].document().toPlainText())
                if self.cB_Xover_op_3[ga].isChecked(): # 3- arithXover
                    Xover_operator[-1].append('arith_Xover')
                    Xover_pct[-1].append(int(self.txt_GA_Xover_int_3[ga].document().toPlainText()))
                    Xover_opt[-1].append(self.txt_GA_Xover_opt_3[ga].document().toPlainText())
                if self.cB_Xover_op_4[ga].isChecked(): # 4- heuristicXover
                    Xover_operator[-1].append('heuristic_Xover')
                    Xover_pct[-1].append(int(self.txt_GA_Xover_int_4[ga].document().toPlainText()))
                    Xover_opt[-1].append(self.txt_GA_Xover_opt_4[ga].document().toPlainText())
                mut_operator.append([]) ; mut_pct.append([]) ; mut_opt.append([])
                if self.cB_mut_op_1[ga].isChecked():
                    mut_operator[-1].append('uniform_mutation')
                    mut_pct[-1].append(int(self.txt_GA_mut_int_1[ga].document().toPlainText()))
                    mut_opt[-1].append(self.txt_GA_mut_opt_1[ga].document().toPlainText())
                if self.cB_mut_op_2[ga].isChecked():
                    mut_operator[-1].append('non_uniform_mutation')
                    mut_pct[-1].append(int(self.txt_GA_mut_int_2[ga].document().toPlainText()))
                    mut_opt[-1].append(self.txt_GA_mut_opt_2[ga].document().toPlainText())
                if self.cB_mut_op_3[ga].isChecked():
                    mut_operator[-1].append('boundary_mutation')
                    mut_pct[-1].append(int(self.txt_GA_mut_int_3[ga].document().toPlainText()))
                    mut_opt[-1].append(self.txt_GA_mut_opt_3[ga].document().toPlainText())
#                if self.cB_mut_op_4[ga].isChecked(): # 4- heuristicXover
#                    mut_operator[-1].append(4)
#                    mut_pct[-1].append(round(n_indiv[-1]*0.01*float(self.txt_GA_mut_int_4[ga].document().toPlainText())))
#                    mut_opt[-1].append(self.txt_GA_mut_opt_4[ga].document().toPlainText())
                mut_intensity.append(int(self.num_GA_mut_prob[ga].text()))

                # PSO options

                inertia.append([])
                inertia[-1].append(str(self.Inertia_score[ga].isChecked()))
                if self.Inertia_score[ga].isChecked():
                    inertia[-1].append(self.Inertia_min[ga].text())
                inertia[-1].append(self.Inertia_i[ga].text())
                inertia[-1].append(self.Inertia_end[ga].text())
                cogn_comp.append([])
                cogn_comp[-1].append(self.Cogn_i[ga].text())
                cogn_comp[-1].append(self.Cogn_end[ga].text())
                social_comp.append([])
                social_comp[-1].append(self.Social_i[ga].text())
                social_comp[-1].append(self.Social_end[ga].text())



# ================   Simulations cases   ======================================

        configs          = []
        fuel             = []
        oxidant          = []
        diluent          = []
        diluent_ratio    = []
        Ts               = []
        Ps               = []
        phis             = []
        tol_ts           = []
        n_pts            = []
        delta_npts       = []
        t_max_coeff      = []
        Scal_ref         = []
        grad_curv_ratio  = []
#        timeVector_param = [n_pts, delta_npts, t_max_coeff, Scal_ref, grad_curv_ratio]
#        timeVector_param = []
        pfr_autom        = []
        pfr_n_pts        = []
        pfr_delta_npts   = []
        pfr_t_max_coeff  = []
        pfr_Scal_ref     = []
        pfr_grad_curv  = []
#        pfr_timeVector_param = []
        tign_nPoints     = []
        tign_dt          = []
        t_max            = []
        transport_model  = []
        xmax             = []
        pts_scatter      = []
        tol_ss           = []
        slope_ff         = []
        curve_ff         = []
        ratio_ff         = []
        prune_ff         = []
        mdots_1          = []
        restore_flame_folder = []
        # Burner 2
        fuel_2           = []
        oxidant_2        = []
        diluent_2        = []
        diluent_ratio_2  = []
        Ts_2             = []
        phis_2           = []
        mdots_2          = []
        #pfr
        n_pts_pfr        = []
        area             = []
        length           = []
        u_0              = []



        for case in range(len(self.condition_activated)):
            if self.condition_activated[case]:
                if self.rB_reactor_UV[case].isChecked(): configs.append('reactor_UV')
                elif self.rB_reactor_HP[case].isChecked(): configs.append('reactor_HP')
                elif self.rB_JSR[case].isChecked(): configs.append('JSR')
                elif self.rB_PFR[case].isChecked(): configs.append('PFR')
                elif self.rB_fflame[case].isChecked(): configs.append('free_flame')
                elif  self.rB_cfflame[case].isChecked():
                    if self.rB_cff_diff[case].isChecked(): configs.append('diff_flame')
                    elif self.rB_cff_pp[case].isChecked(): configs.append('pp_flame')
                    elif self.rB_cff_tp[case].isChecked(): configs.append('tp_flame')


                if not self.rB_cfflame[case].isChecked():
                    fuel.append(self.fuel_1[case].document().toPlainText())
                    oxidant.append(self.oxidant_1[case].document().toPlainText())
                    diluent.append(self.Diluent_1[case].document().toPlainText())
                    diluent_ratio.append(self.Diluent_2[case].document().toPlainText())
                    T_min  = float(self.Tmin[case].document().toPlainText())
                    T_max  = float(self.Tmax[case].document().toPlainText())
                    T_incr = float(self.Tincr[case].document().toPlainText())
                    if T_incr==0: T_incr=T_max-T_min+T_max*5
                    Ts.append(list(np.arange(T_min, T_max+T_incr/2, T_incr)))
                    P_min  = float(self.Pmin[case].document().toPlainText())
                    P_max  = float(self.Pmax[case].document().toPlainText())
                    P_incr = float(self.Pincr[case].document().toPlainText())
                    if P_incr==0: P_incr=P_max-P_min+P_max*5
                    Ps.append(list(np.arange(P_min, P_max+P_incr/2, P_incr)))
                    phi_min  = float(self.eqmin[case].document().toPlainText())
                    phi_max  = float(self.eqmax[case].document().toPlainText())
                    phi_incr = float(self.eqincr[case].document().toPlainText())
                    if phi_incr==0: phi_incr=phi_max-phi_min+phi_max*5
                    phis.append(list(np.arange(phi_min, phi_max+phi_incr/2, phi_incr)))

                    mdots_1.append(False)
                    fuel_2.append(False)
                    oxidant_2.append(False)
                    diluent_2.append(False)
                    diluent_ratio_2.append(False)
                    Ts_2.append(False)
                    phis_2.append(False)
                    mdots_2.append(False)
                else: # counterflow flame configuration
                    P_min  = float(self.df_Pmin[case].document().toPlainText())
                    P_max  = float(self.df_Pmax[case].document().toPlainText())
                    P_incr = float(self.df_Pincr[case].document().toPlainText())
                    Ps.append(list(np.arange(P_min, P_max+P_incr/2, P_incr)))
                    # Burner 1
                    fuel.append(self.df_fuel_1[case].document().toPlainText())
                    oxidant.append(self.df_oxidant_1[case].document().toPlainText())
                    diluent.append(self.df_Diluent_1[case].document().toPlainText())
                    diluent_ratio.append(self.df_Diluent_r_1[case].document().toPlainText())
                    Ts.append(float(self.df_T_1[case].document().toPlainText()))
                    phi_min  = float(self.df_eqmin_1[case].document().toPlainText())
                    phi_max  = float(self.df_eqmax_1[case].document().toPlainText())
                    phi_incr = float(self.df_eqincr_1[case].document().toPlainText())
                    if phi_incr==0: phis.append([phi_min])
                    else: phis.append(list(np.arange(phi_min, phi_max+phi_incr/2, phi_incr)))
                    mdot_min  = float(self.df_mdot1_1[case].document().toPlainText())
                    mdot_max  = float(self.df_mdot2_1[case].document().toPlainText())
                    if 'diff' in configs[-1]:
                        mdot_incr  = float(self.df_mdot4_1[case].document().toPlainText())
                    elif self.df_mdot3_1[case].text() == 'error':
                        mdot_incr  = 0
                    else:
                        mdot_incr  = float(self.df_mdot3_1[case].text())
                    if mdot_incr==0: mdot_incr=mdot_max-mdot_min+mdot_max*5
                    mdots_1.append(list(np.arange(mdot_min, mdot_max+mdot_min/2, mdot_incr)))
                    # Burner 2
                    if 'tp_' not in configs[-1]:
                        fuel_2.append(self.df_fuel_2[case].document().toPlainText())
                        oxidant_2.append(self.df_oxidant_2[case].document().toPlainText())
                        diluent_2.append(self.df_Diluent_2[case].document().toPlainText())
                        diluent_ratio_2.append(self.df_Diluent_r_2[case].document().toPlainText())
                        Ts_2.append(float(self.df_T_2[case].document().toPlainText()))
                        phi_min  = float(self.df_eqmin_2[case].document().toPlainText())
                        phi_max  = float(self.df_eqmax_2[case].document().toPlainText())
                        phi_incr = float(self.df_eqincr_2[case].text())
                        if phi_incr!=0:
                            phis_2.append(list(np.arange(phi_min, phi_max+phi_incr/2, phi_incr)))
                        else:
                            phis_2.append([phi_min]*len(phis[-1]))
                        mdot_min  = float(self.df_mdot1_2[case].document().toPlainText())
                        mdot_max  = float(self.df_mdot2_2[case].document().toPlainText())
                        if 'diff' in configs[-1]:
                            mdot_incr  = float(self.df_mdot4_1[case].document().toPlainText())
                        else:
                            mdot_incr  = float(self.df_mdot3_1[case].text())
                        if mdot_incr==0: mdot_incr=mdot_max-mdot_min+mdot_max*5
                        mdots_2.append(list(np.arange(mdot_min, mdot_max+mdot_min/2, mdot_incr)))
                    else:
                        fuel_2.append(False)
                        oxidant_2.append(False)
                        diluent_2.append(False)
                        diluent_ratio_2.append(False)
                        Ts_2.append(False)
                        phis_2.append(False)
                        mdots_2.append(False)

                # JSR + reactor + flame options
                if 'reactor' in configs[-1]:
                    tol_ts.append([float(self.tol_ts_rel_r[case].document().toPlainText()),float(self.tol_ts_abs_r[case].document().toPlainText())])
                elif configs[-1] == 'JSR':
                    tol_ts.append([float(self.tol_ts_rel_jsr[case].document().toPlainText()),float(self.tol_ts_abs_jsr[case].document().toPlainText())])
                elif configs[-1] == 'PFR':
                    tol_ts.append([float(self.tol_ts_rel_pfr[case].document().toPlainText()),float(self.tol_ts_abs_pfr[case].document().toPlainText())])
                elif 'flame' in configs[-1]:
                    tol_ts.append([float(self.tol_ts_rel_f[case].document().toPlainText()),float(self.tol_ts_abs_f[case].document().toPlainText())])

                # reactor options
                n_pts.append(float(self.pts_num[case].text()))
                delta_npts.append(float(self.delta_n_pts[case].text()))
                t_max_coeff.append(float(self.t_max_coeff[case].text()))
                Scal_ref.append(self.ref_spec[case].document().toPlainText())
                grad_curv_ratio.append(float(self.grad_curv[case].text().replace(',','.')))
#                timeVector_param.append([n_pts, delta_npts, t_max_coeff, Scal_ref, grad_curv_ratio])
                tign_nPoints.append(float(self.txt_Gr_ipn[case].document().toPlainText()))       # number of time step for auto-ignition detection
                tign_dt.append(float(self.txt_Gr_its[case].document().toPlainText()))            # initial time step for auto-ignition detection

                # JSR options
                t_max.append(float(self.t_max_jsr[case].text().replace(',','.')))

                # PFR options
                n_pts_pfr.append(float(self.pts_num_pfr[case].text()))
                area.append(self.pfr_area[case].document().toPlainText())
                length.append(self.pfr_length[case].document().toPlainText())
                u_0.append(self.pfr_u_0[case].document().toPlainText())
                if self.rB_pfr_1[case].isChecked(): pfr_autom.append(True)
                else:                               pfr_autom.append(False)
                pfr_n_pts.append(float(self.pfr_pts_num[case].text()))
                pfr_delta_npts.append(float(self.pfr_delta_n_pts[case].text()))
                pfr_t_max_coeff.append(float(self.pfr_t_max_coeff[case].text()))
                pfr_Scal_ref.append(self.pfr_ref_spec[case].document().toPlainText())
                pfr_grad_curv.append(float(self.pfr_grad_curv[case].text().replace(',','.')))
#                pfr_timeVector_param.append([pfr_n_pts, pfr_delta_npts, pfr_t_max_coeff, pfr_Scal_ref, pfr_grad_curv])


                # flame options
                tol_ss.append([float(self.tol_ss_rel_f[case].document().toPlainText()),float(self.tol_ss_abs_f[case].document().toPlainText())])
                if self.mix[case].isChecked():   transport_model.append('Mix')
                else:                            transport_model.append('Mult')
                xmax.append(float(self.txt_Gf_xmax[case].document().toPlainText()))
                pts_scatter.append('['+self.pts_scatter[case].document().toPlainText()+']')
                slope_ff.append(float(self.slope_ff[case].text().replace(',','.')))
                curve_ff.append(float(self.curve_ff[case].text().replace(',','.')))
                ratio_ff.append(float(self.ratio_ff[case].text().replace(',','.')))
                prune_ff.append(float(self.prune_ff[case].text().replace(',','.')))
                if self.rB_flame_restore[case].isChecked():   restore_flame_folder.append(self.label_flame_folder[case].text())
                else:                                         restore_flame_folder.append(False)


# =============================================================================
#         Write data
# =============================================================================

        os.chdir(self.WD_path)

        if get_fn:
            filename = '_conditions_input/'+self.Text_file_name.document().toPlainText()+'.inp'
        else:
            filename = '_conditions_input/last_condition.inp'
        fd = open(filename, 'w')

        fd.write('#=============================================\n')
        fd.write('#           Main parameters\n')
        fd.write('#=============================================\n')
        fd.write('main_path         = ' + main_path          + '\n')
        fd.write('mech              = ' + mech               + '\n')
        if str(mech_prev_red)!='False':
            if mech_prev_red:
                fd.write('mech_prev_red     = ' + str(mech_prev_red) + '\n')
        if self.external_results:
            fd.write('ext_results_file  = ' + self.external_results.split('/')[-1] + '\n')
            fd.write('conc_units        = ' + self.ext_res_conc_unit    + '\n')
#            fd.write('ext_data_type     = ' + self.ext_res_file_type    + '\n')
        fd.write('verbose           = ' + str(verbose)       + '\n')
        fd.write('show_plots        = ' + str(show_plots)    + '\n')
        fd.write('write_ck          = ' + str(write_ck)    + '\n')
        fd.write('tspc              = ' + list2txt(tspc)     + '\n')
        fd.write('T_check           = ' + str(T_check)       + '\n')
        fd.write('sp_T              = ' + list2txt(sp_T)     + '\n')
        fd.write('Sl_check          = ' + str(Sl_check)      + '\n')
        fd.write('sp_Sl             = ' + list2txt(sp_Sl)    + '\n')
        fd.write('ig_check          = ' + str(ig_check)      + '\n')
        fd.write('sp_ig             = ' + list2txt(sp_ig)    + '\n')
        fd.write('K_check           = ' + str(K_check)       + '\n')
        fd.write('sp_K              = ' + list2txt(sp_K)     + '\n')
        fd.write('error_calculation = ' + error_calculation  + '\n')
        fd.write('error_coupling    = ' + error_coupling     + '\n')


        fd.write('\n\n#=============================================\n')
        fd.write('#           Simulation cases\n')
        fd.write('#=============================================\n')
        case_act = 0
        for case in range(len(self.condition_activated)):
            if self.condition_activated[case]:
                fd.write('\n#======> Case '     + str(case_act+1)              + '\n')
                fd.write('config            = ' + configs[case_act]            + '\n')
                fd.write('Ps                = ' + list2txt(Ps[case_act])       + '\n')
                if 'diff_' not in configs[case_act] and 'pp_' not in configs[case_act]:
                    if ':' in fuel[case_act]:
                        fd.write('mixt              = ' + fuel[case_act]               + '\n')
                    else:
                        fd.write('fuel              = ' + fuel[case_act]               + '\n')
                        fd.write('oxidant           = ' + oxidant[case_act]            + '\n')
                        fd.write('diluent           = ' + diluent[case_act]            + '\n')
                        fd.write('diluent_ratio     = ' + str(diluent_ratio[case_act]) + '\n')
                    fd.write('Ts                = ' + list2txt(Ts[case_act])       + '\n')
#                    fd.write('Ps                = ' + list2txt(Ps[case_act])       + '\n')
                    fd.write('phis              = ' + list2txt(phis[case_act])     + '\n')
                else:
                    if ':' in fuel[case_act]:
                        fd.write('mixt              = ' + fuel[case_act]               + '\n')
                    else:
                        fd.write('fuel_1            = ' + fuel[case_act]               + '\n')
                        if 'diff_' not in configs[case_act]:
                            fd.write('oxidant_1         = ' + oxidant[case_act]            + '\n')
                        fd.write('diluent_1         = ' + diluent[case_act]            + '\n')
                        fd.write('diluent_ratio_1   = ' + str(diluent_ratio[case_act]) + '\n')
                    fd.write('Ts_1              = ' + list2txt(Ts[case_act])       + '\n')
                    fd.write('phis_1            = ' + list2txt(phis[case_act])     + '\n')
                if 'diff_' in configs[case_act] or 'pp_' in configs[case_act]:
                    fd.write('mdots_1           = ' + list2txt(mdots_1[case_act])      + '\n')
                if 'reactor' in configs[case_act]:
                    fd.write('n_pts             = ' + str(n_pts[case_act])             + '\n')
                    fd.write('delta_npts        = ' + str(delta_npts[case_act])        + '\n')
                    fd.write('t_max_coeff       = ' + str(t_max_coeff[case_act])       + '\n')
                    fd.write('Scal_ref          = ' + Scal_ref[case_act]               + '\n')
                    fd.write('grad_curv_ratio   = ' + str(grad_curv_ratio[case_act])   + '\n')
                    fd.write('tign_det_nPoints  = ' + str(tign_nPoints[case_act])      + '\n')
                    fd.write('tign_det_tmax_sec = ' + str(tign_dt[case_act])           + '\n')
                elif 'JSR' in configs[case_act]:
                    fd.write('t_max             = ' + str(t_max[case_act])             + '\n')
                elif 'PFR' in configs[case_act]:
                    fd.write('pfr_autom         = ' + str(pfr_autom[case_act])             + '\n')
                    if pfr_autom[case_act]:
                        fd.write('n_pts             = ' + str(pfr_n_pts[case_act])             + '\n')
                        fd.write('delta_npts        = ' + str(pfr_delta_npts[case_act])        + '\n')
                        fd.write('t_max_coeff       = ' + str(pfr_t_max_coeff[case_act])       + '\n')
                        fd.write('Scal_ref          = ' + str(pfr_Scal_ref[case_act])          + '\n')
                        fd.write('grad_curv_ratio   = ' + str(pfr_grad_curv[case_act])   + '\n')
                    else:
                        fd.write('n_pts             = ' + str(n_pts_pfr[case_act])             + '\n')
                        fd.write('length (m)        = ' + str(length[case_act])                + '\n')
                    fd.write('area (m2)             = ' + str(area[case_act])                  + '\n')
                    fd.write('u_0 (m/s)             = ' + str(u_0[case_act])                   + '\n')
                elif 'free_' in configs[case_act]:
                    fd.write('xmax              = ' + str(xmax[case_act])              + '\n')
                elif 'diff_' in configs[case_act] or 'pp_' in configs[case_act]:
                    if 'diff_' not in configs[case_act]:
                        if ':' in fuel_2[case_act]:
                            fd.write('mixt_2            = ' + fuel_2[case_act]             + '\n')
                        else:
                            fd.write('fuel_2            = ' + fuel_2[case_act]             + '\n')
                    if ':' in oxidant_2[case_act]:
                        fd.write('mixt_2            = ' + oxidant_2[case_act]             + '\n')
                    fd.write('oxidant_2         = ' + oxidant_2[case_act]              + '\n')
                    fd.write('diluent_2         = ' + diluent_2[case_act]              + '\n')
                    fd.write('diluent_ratio_2   = ' + diluent_ratio_2[case_act]        + '\n')
                    fd.write('Ts_2              = ' + list2txt(Ts_2[case_act])         + '\n')
                    if 'diff_' not in configs[case_act]:
                        fd.write('phis_2            = ' + list2txt(phis_2[case_act])       + '\n')
                    fd.write('mdots_2           = ' + list2txt(mdots_2[case_act])      + '\n')
                    fd.write('width             = ' + str(xmax[case_act])              + '\n')
                elif 'tp_' in configs[case_act]:
                    fd.write('mdots             = ' + list2txt(mdots_1[case_act])      + '\n')
                    fd.write('width             = ' + str(xmax[case_act])              + '\n')
                fd.write('tol_ts            = ' + list2txt(tol_ts[case_act])   + '\n')
                if 'flame' in configs[case_act]:
                    fd.write('tol_ss            = ' + list2txt(tol_ss[case_act])            + '\n')
                    fd.write('transport_model   = ' + transport_model[case_act]        + '\n')
                    fd.write('pts_scatter       = ' + str(pts_scatter[case_act])       + '\n')
                    fd.write('slope             = ' + str(slope_ff[case_act])       + '\n')
                    fd.write('curve             = ' + str(curve_ff[case_act])       + '\n')
                    fd.write('ratio             = ' + str(ratio_ff[case_act])       + '\n')
                    fd.write('prune             = ' + str(prune_ff[case_act])       + '\n')
                    if restore_flame_folder[case_act]:
                        fd.write('restore_flame_folder = ' + str(restore_flame_folder[case_act])       + '\n')
                    else:
                        fd.write('restore_flame_folder = False'  + '\n')
                case_act+=1


        fd.write('\n\n\n#=============================================\n')
        fd.write('#           Operators\n')
        fd.write('#=============================================\n\n')

        op = -1 ; ag = -1

        if len(self.list_operator)>0:
            if self.list_operator[0] == 'opt':
                ag+=1
                fd.write('#===========> Op: Optimization without reduction\n')
                fd.write('operator        = NULL \n')
                fd.write('optim           = '  + str(opt_meth[ag])         +  '\n')
                fd.write('#====> Optimization\n')
                fd.write('n_gen              = ' + str(n_gen[ag])              + '\n')
                fd.write('n_indiv            = ' + str(n_indiv[ag])            + '\n')
                fd.write('error_fitness      = ' + str(error_fitness[ag])      + '\n')
                fd.write('Arrh_max_variation = ' + list2txt(Arrh_max_variation[ag]) + '\n')
                if type(self.cB_GA_meth_DRG[0]) is bool:
                    if self.cB_GA_meth_DRG[0]:
                        fd.write('optim_on_meth      = DRG \n')
                        fd.write('optim_on_meth_pts  = ' +  str(self.cB_GA_meth_pts[ag].text()) + '\n')
                        fd.write('nb_r2opt            = ' + str(nb_r2opt[ag])         + '\n')
                else:
                    if self.cB_GA_meth_DRG[0].isChecked():
                        fd.write('optim_on_meth      = DRG \n')
                        fd.write('optim_on_meth_pts  = ' +  str(self.cB_GA_meth_pts[ag].text()) + '\n')
                        fd.write('nb_r2opt            = ' + str(nb_r2opt[ag])         + '\n')
                if type(self.cB_GA_meth_SA[0]) is bool:
                    if self.cB_GA_meth_SA[0]:
                        fd.write('optim_on_meth      = SA \n')
                        fd.write('optim_on_meth_pts  = ' +  str(self.cB_GA_meth_pts[ag].text()) + '\n')
                        fd.write('nb_r2opt            = ' + str(nb_r2opt[ag])         + '\n')
                elif self.cB_GA_meth_SA[0].isChecked():
                    fd.write('optim_on_meth      = SA \n')
                    fd.write('optim_on_meth_pts  = ' +  str(self.cB_GA_meth_pts[ag].text()) + '\n')
                    fd.write('nb_r2opt            = ' + str(nb_r2opt[ag])         + '\n')
                if opt_meth[ag] == 'GA':
                    fd.write('sub_mech_sel        = ' + list2txt(sub_mech_sel[ag])       + '\n\n')
                    # GA options
                    fd.write('selection_operator  = ' + str(selection_operator[ag]) + '\n')
                    fd.write('selection_options   = ' + list2txt(selection_options[ag])  + '\n')
                    fd.write('Xover_operator      = ' + list2txt(Xover_operator[ag])     + '\n')
                    fd.write('Xover_pct           = ' + list2txt(Xover_pct[ag])          + '\n')
#                    fd.write('Xover_opt           = ' + list2txt(Xover_opt[ag])          + '\n')
                    fd.write('mut_operator        = ' + list2txt(mut_operator[ag])       + '\n')
                    fd.write('mut_pct             = ' + list2txt(mut_pct[ag])            + '\n')
                    fd.write('mut_opt             = ' + list2txt(mut_opt[ag])            + '\n')
                    fd.write('mut_intensity       = ' + str(mut_intensity[ag])           + '\n')
                if opt_meth[ag] == 'PSO':
                    # PSO options
                    if 'True' in inertia[ag][0]:
                        fd.write('inertia_score       = True \n')
                        fd.write('inertia_min         = ' + inertia[ag][1].replace(',','.')     + '\n')
                        fd.write('inertia_max_i       = ' + inertia[ag][2].replace(',','.')     + '\n')
                        fd.write('inertia_max_end     = ' + inertia[ag][3].replace(',','.')     + '\n')
                    else:
                        fd.write('inertia_score       = False \n')
                        fd.write('inertia_i           = ' + inertia[ag][1].replace(',','.')     + '\n')
                        fd.write('inertia_end         = ' + inertia[ag][2].replace(',','.')     + '\n')
                    fd.write('cognitive_accel_i   = ' + cogn_comp[ag][0].replace(',','.')       + '\n')
                    fd.write('cognitive_accel_end = ' + cogn_comp[ag][0].replace(',','.')       + '\n')
                    fd.write('social_accel_i      = ' + social_comp[ag][0].replace(',','.')     + '\n')
                    fd.write('social_accel_end    = ' + social_comp[ag][1].replace(',','.')     + '\n\n')
                self.cB_GA_meth_DRG.append(False)
                self.cB_GA_meth_SA.append(False)
                self.cB_GA_meth_pts.append(False)
                self.GA_label_meth_pts.append(False)
        else:
            fd.write('#===========> Op: no reduction\n')
            fd.write('operator        = NULL \n')


        for red_i in range(len(self.list_operator)):
            if self.list_operator[red_i] != 'opt':
                op+=1
#                if 'DRG' in self.list_operator[red_i]:
                fd.write('\n#===========> Op: ' + self.list_operator[red_i] + '\n')
                fd.write('operator        = ' + self.list_operator[red_i] + ' \n')
#                else:
#                    fd.write('\n#===========> Op: ' + self.list_operator[red_i] + '\n')
#                    fd.write('operator        = ' + self.list_operator[red_i] + ' \n')
                fd.write('eps             = ' + str(eps[op])           +  '\n')
                fd.write('delta_eps       = ' + str(delta_eps[op])     +  '\n')
                fd.write('n_points        = ' + str(npoints[op])       +  '\n')
                fd.write('max_error_sp    = ' + str(max_error_sp[op]).replace('[','').replace(']','').replace("'","")  +  '\n')
                if T_check:  fd.write('max_error_T     = ' + str(max_error_T[op])   +  '\n')
                if ig_check: fd.write('max_error_ig    = ' + str(max_error_ig[op])  +  '\n')
                if Sl_check: fd.write('max_error_Sl    = ' + str(max_error_Sl[op])  +  '\n')
                if K_check:  fd.write('max_error_K     = ' + str(max_error_K[op])   +  '\n')
                fd.write('inter_sp_inter  = '  + str(isi[op])         +  '\n')
                if 'DRG' not in self.list_operator[red_i]:
                    fd.write('ttol_sensi      = '    + list2txt(ttol_sensi[op])    +  '\n')
                if optim[op]:
                    ag+=1
                    fd.write('optim           = '  + str(opt_meth[ag])         +  '\n')
                    if opt_meth[ag] == 'GA':
                        fd.write('#====> Genetic Algorithm Optimization\n')
                        fd.write('n_gen               = ' + str(n_gen[ag])              + '\n')
                    elif opt_meth[ag] == 'PSO':
                        fd.write('#====> Particle Swarm Optimization\n')
                        fd.write('n_it                = ' + str(n_gen[ag])              + '\n')
                    fd.write('n_indiv             = ' + str(n_indiv[ag])            + '\n')
                    fd.write('error_fitness       = ' + str(error_fitness[ag])      + '\n')
                    fd.write('Arrh_max_variation  = ' + list2txt(Arrh_max_variation[ag]) + '\n')
                    fd.write('optim_on_meth       = ' + str(optim_on_meth[ag])      + '\n')
                    if str(optim_on_meth[ag]) != 'False':
                        fd.write('nb_r2opt            = ' + str(nb_r2opt[ag])         + '\n')
                    fd.write('sub_mech_sel        = ' + list2txt(sub_mech_sel[ag])       + '\n')
                    if opt_meth[ag] == 'GA':
                        # GA options
                        fd.write('selection_operator  = ' + str(selection_operator[ag]) + '\n')
                        fd.write('selection_options   = ' + list2txt(selection_options[ag])  + '\n')
                        fd.write('Xover_operator      = ' + list2txt(Xover_operator[ag])     + '\n')
                        fd.write('Xover_pct           = ' + list2txt(Xover_pct[ag])          + '\n')
    #                    fd.write('Xover_opt           = ' + list2txt(Xover_opt[ag])          + '\n')
                        fd.write('mut_operator        = ' + list2txt(mut_operator[ag])       + '\n')
                        fd.write('mut_pct             = ' + list2txt(mut_pct[ag])            + '\n')
                        fd.write('mut_opt             = ' + list2txt(mut_opt[ag])            + '\n')
                        fd.write('mut_intensity       = ' + str(mut_intensity[ag])           + '\n\n')
                    elif opt_meth[ag] == 'PSO':
                        # PSO options
                        if 'True' in inertia[ag][0]:
                            fd.write('inertia_score       = True \n')
                            fd.write('inertia_min         = ' + inertia[ag][1].replace(',','.')     + '\n')
                            fd.write('inertia_max_i       = ' + inertia[ag][2].replace(',','.')     + '\n')
                            fd.write('inertia_max_end     = ' + inertia[ag][3].replace(',','.')     + '\n')
                        else:
                            fd.write('inertia_score       = False \n')
                            fd.write('inertia_i           = ' + inertia[ag][1].replace(',','.')     + '\n')
                            fd.write('inertia_end         = ' + inertia[ag][2].replace(',','.')     + '\n')
                        fd.write('cognitive_accel_i   = ' + cogn_comp[ag][0].replace(',','.')       + '\n')
                        fd.write('cognitive_accel_end = ' + cogn_comp[ag][0].replace(',','.')       + '\n')
                        fd.write('social_accel_i      = ' + social_comp[ag][0].replace(',','.')     + '\n')
                        fd.write('social_accel_end    = ' + social_comp[ag][1].replace(',','.')     + '\n\n')


        fd.close()


    def load_parameters(self):
        global sz_w
        global sz_h
        global load_prev_cond
        os.chdir(self.WD_path)
        load_prev_cond = True

        #remove conditions and operators
        for case in range(len(self.condition_activated)):
            self.remove_conditions(case)
        n_GA = self.list_operator.count('opt')
        for red_i in range(len(self.list_operator)-n_GA):
            self.tablet.setCurrentIndex(2)
            self.remove_tab('opt')
        self.tablet.setCurrentIndex(0)
        for idx in range(len(self.list_target)):
            self.list_target.takeItem(0)


        # get_file
        A=Files_windows()
        options = QtWidgets.QFileDialog.Options()
        options |= QtWidgets.QFileDialog.DontUseNativeDialog
        filename, _ = QtWidgets.QFileDialog.getOpenFileName(A, "Select the input file", "_conditions_input","Input Files (*.inp);;All Files (*)", options=options)

        _translate = QtCore.QCoreApplication.translate

        fs = open(filename, 'r')
        txt = fs.readline().split('=')

        # Main parameters
        while 'Case' not in txt[-1] and 'Operator' not in txt[-1]  and txt[0] != '':
            txt = fs.readline().split('=')
            txt[0]=txt[0].replace(' ','')

            if txt[0] == 'main_path':         main_path_ext = genf.clean_txt(txt[1])
            if txt[0] == 'mech':              mech          = genf.clean_txt(txt[1])
            if txt[0] == 'mech_prev_red':     mech_prev_red = genf.clean_txt(txt[1])
            if txt[0] == 'ext_results_file':  ext_results_file = genf.clean_txt(txt[1])
            if txt[0] == 'ext_data_type':     ext_data_type    = genf.clean_txt(txt[1])
            if txt[0] == 'conc_units':        conc_units       = genf.clean_txt(txt[1])
            if txt[0] == 'verbose':           verbose       = int(txt[1])
            if txt[0] == 'show_plots':        show_plots    = genf.str2bool(txt[1])
            if txt[0] == 'write_ck':          write_ck      = genf.str2bool(txt[1])
            if txt[0] == 'tspc':
                tspc          = genf.txt2list_string(txt[1])
                n_tspc = len(tspc)
            if txt[0] == 'T_check':           T_check       = genf.str2bool(txt[1])
            if txt[0] == 'sp_T':              sp_T          = genf.clean_txt(txt[1])
            if txt[0] == 'Sl_check':          Sl_check      = genf.str2bool(txt[1])
            if txt[0] == 'sp_Sl':             sp_Sl         = genf.clean_txt(txt[1])
            if txt[0] == 'ig_check':          ig_check      = genf.str2bool(txt[1])
            if txt[0] == 'sp_ig':             sp_ig         = genf.clean_txt(txt[1])
            if txt[0] == 'K_check':           K_check       = genf.str2bool(txt[1])
            if txt[0] == 'sp_K':              sp_K          = genf.clean_txt(txt[1])
            if txt[0] == 'error_calculation': error_calculation = genf.clean_txt(txt[1])
            if txt[0] == 'error_coupling':    error_coupling    = genf.clean_txt(txt[1])


        if 'main_path_ext' in locals():
            self.Text_folder_name.setPlainText(_translate("MainWindow", main_path_ext))
        if 'mech' in locals():
            self.label_reference_mechanism.setText(_translate("MainWindow", mech))
            font = QtGui.QFont();font.setBold(True);font.setItalic(False);font.setWeight(75)
            self.label_reference_mechanism.setFont(font)
            self.reference_mechanism = mech
            ct.suppress_thermo_warnings()
            os.chdir(self.WD_path)
            try:
                self.gas_ref = ct.Solution(self.reference_mechanism)
                self.mech_data = cdef.Mech_data(self.reference_mechanism,self.gas_ref)
            except:
                self.gas_ref = ct.Solution('_kinetic_mech/'+self.reference_mechanism)
                self.mech_data = cdef.Mech_data('_kinetic_mech/'+self.reference_mechanism,self.gas_ref)
            self.ns_ref = self.gas_ref.n_species
            self.nr_ref = self.gas_ref.n_reactions
            for sp in range(self.ns_ref):
                item = QtWidgets.QListWidgetItem()
                self.list_spec.addItem(item)
                item = self.list_spec.item(sp)
                item.setText(_translate("MainWindow", self.gas_ref.species_name(sp)))
            # appearance of "import external data" option
            self.pB_External_data.setGeometry(QtCore.QRect(10*sz_w, 90*sz_h, 151*sz_w, 34*sz_h))
            self.label_External_data.setGeometry(QtCore.QRect(180*sz_w, 100*sz_h, 350*sz_w, 18*sz_h))
        filename_save = filename.split('/')[-1].split('.')[0]
        self.Text_file_name.setPlainText(_translate("MainWindow", filename_save))
        if 'mech_prev_red' in locals():
            self.label_reduced_mechanism.setText(_translate("MainWindow", mech_prev_red))
            font = QtGui.QFont();font.setBold(True);font.setItalic(False);font.setWeight(75)
            self.label_reduced_mechanism.setFont(font)
            self.reduced_mechanism = mech_prev_red

        if 'ext_results_file' in locals():
            self.label_External_data.setText(_translate("MainWindow", ext_results_file.split('/')[-1]))
            self.label_External_data.setFont(font)
            self.external_results = ext_results_file
            if 'ext_data_type' in locals():
                self.ext_res_file_type = ext_data_type
            if 'conc_units' in locals():
                self.ext_res_conc_unit = conc_units

        if 'verbose' in locals():     self.MP_verbose.setProperty("value", verbose)
        if 'show_plots' in locals():
            if show_plots:  self.cB_show_plots.setChecked(True)
            else:           self.cB_show_plots.setChecked(False)
        if 'write_ck' in locals():
            if write_ck:    self.cB_chemkin_f.setChecked(True)
            else:           self.cB_chemkin_f.setChecked(False)
        if 'tspc' in locals():        self.add_target_spec(tspc)
        if 'T_check' in locals():
            if T_check:     self.cB_tsp_T.setChecked(True)
            else:           self.cB_tsp_T.setChecked(False)
        if 'Sl_check' in locals():
            if Sl_check:    self.cB_tsp_Sl.setChecked(True)
            else:           self.cB_tsp_Sl.setChecked(False)
        if 'ig_check' in locals():
            if ig_check:    self.cB_tsp_igt.setChecked(True)
            else:           self.cB_tsp_igt.setChecked(False)
        if 'K_check' in locals():
            if K_check:     self.cB_tsp_K.setChecked(True)
            else:           self.cB_tsp_K.setChecked(False)

        if 'sp_T'  in locals(): self.text_MP_T.setPlainText(_translate("MainWindow", sp_T))
        if 'sp_Sl' in locals(): self.text_MP_Sl.setPlainText(_translate("MainWindow", sp_Sl))
        if 'sp_ig' in locals(): self.text_MP_ig.setPlainText(_translate("MainWindow", sp_ig))
        if 'sp_K' in locals(): self.text_MP_K.setPlainText(_translate("MainWindow", sp_K))

        if 'T_check' not in locals():   T_check = False
        if 'Sl_check' not in locals():  Sl_check = False
        if 'ig_check' not in locals():  ig_check = False
        if 'K_check' not in locals():   K_check = False

        if 'error_calculation' in locals():
            if   error_calculation == 'points':  self.rB_MP_qoi_1.setChecked(True)
            elif error_calculation == 'QoI':     self.rB_MP_qoi_2.setChecked(True)
        if 'error_coupling' in locals():
            if   error_coupling == 'mean':   self.rB_MP_errmean_1.setChecked(True)
            elif error_coupling == 'max':    self.rB_MP_errmean_2.setChecked(True)


        # Simulation cases
        conditions_list = [] ; load_conditions=False ; case_n=-1
        while 'Operators' not in txt[-1] and txt[0] != '':
            txt = fs.readline().split('=')
            # get data
            while 'Case' not in txt[-1] and 'Operators' not in txt[-1]:
                txt[0]=txt[0].replace(' ','')
                if txt[0]== 'config':
                    config = genf.clean_txt(txt[1]) ; case_n+=1 ; load_conditions = True
                if txt[0] == 'fuel':              fuel              = genf.clean_txt(txt[1])
                if txt[0] == 'mixt':              fuel              = genf.clean_txt(txt[1])
                if txt[0] == 'oxidant':           oxidant           = genf.clean_txt(txt[1])
                if txt[0] == 'diluent':           diluent           = genf.clean_txt(txt[1])
                if txt[0] == 'diluent_ratio':     diluent_ratio     = genf.diluent_fct(txt[1])
                if txt[0] == 'T_min':             T_min             = genf.clean_txt(txt[1])
                if txt[0] == 'T_max':             T_max             = genf.clean_txt(txt[1])
                if txt[0] == 'T_incr':            T_incr            = genf.clean_txt(txt[1])
                if txt[0] == 'Ts':                Ts                = genf.txt2list_float(txt[1])
                if txt[0] == 'phi_min':           phi_min           = genf.clean_txt(txt[1])
                if txt[0] == 'phi_max':           phi_max           = genf.clean_txt(txt[1])
                if txt[0] == 'phi_incr':          phi_incr          = genf.clean_txt(txt[1])
                if txt[0] == 'phis':              phis              = genf.txt2list_float(txt[1])

                if txt[0] == 'fuel_1':              fuel              = genf.clean_txt(txt[1])
                if txt[0] == 'mixt_1':              fuel              = genf.clean_txt(txt[1])
                if txt[0] == 'oxidant_1':           oxidant           = genf.clean_txt(txt[1])
                if txt[0] == 'diluent_1':           diluent           = genf.clean_txt(txt[1])
                if txt[0] == 'diluent_ratio_1':     diluent_ratio     = genf.diluent_fct(txt[1])
                if txt[0] == 'T_min_1':             T_min             = genf.clean_txt(txt[1])
                if txt[0] == 'T_max_1':             T_max             = genf.clean_txt(txt[1])
                if txt[0] == 'T_incr_1':            T_incr            = genf.clean_txt(txt[1])
                if txt[0] == 'Ts_1':                Ts                = genf.txt2list_float(txt[1])
                if txt[0] == 'phi_min_1':           phi_min           = genf.clean_txt(txt[1])
                if txt[0] == 'phi_max_1':           phi_max           = genf.clean_txt(txt[1])
                if txt[0] == 'phi_incr_1':          phi_incr          = genf.clean_txt(txt[1])
                if txt[0] == 'phis_1':              phis              = genf.txt2list_float(txt[1])

                if txt[0] == 'P_min':             P_min             = genf.clean_txt(txt[1])
                if txt[0] == 'P_max':             P_max             = genf.clean_txt(txt[1])
                if txt[0] == 'P_incr':            P_incr            = genf.clean_txt(txt[1])
                if txt[0] == 'Ps':                Ps                = genf.txt2list_float(txt[1])

                if txt[0] == 'tol_ts':            tol_ts            = genf.clean_txt(txt[1])
                # options for reactor and prf
                if txt[0] == 'n_pts':             n_pts             = float(txt[1])
                if txt[0] == 'delta_npts':        delta_npts        = float(txt[1])
                if txt[0] == 't_max_coeff':       t_max_coeff       = float(txt[1])
                if txt[0] == 'Scal_ref':          Scal_ref          = genf.clean_txt(txt[1])
                if txt[0] == 'grad_curv_ratio':   grad_curv_ratio   = float(txt[1])
                if txt[0] == 'tign_nPoints':      tign_nPoints      = genf.clean_txt(txt[1])
                if txt[0] == 'tign_dt':           tign_dt           = genf.clean_txt(txt[1])
                # options for pfr
                if txt[0] == 'area':              area              = float(txt[1])
                if txt[0] == 'length':            length            = float(txt[1])
                if txt[0] == 'u_0':               u_0               = float(txt[1])
                if txt[0] == 'pfr_autom':         pfr_autom         = genf.str2bool(txt[1])
                # options for jsr
                if txt[0] == 't_max':             t_max             = float(txt[1]); caution_opt_jsr=False
                # options for flame
                if txt[0] == 'tol_ss':            tol_ss          = genf.clean_txt(txt[1])
                if txt[0] == 'transport_model':   transport_model = genf.clean_txt(txt[1])
                if txt[0] == 'pts_scatter':       pts_scatter     = genf.clean_txt(txt[1])
                if txt[0] == 'slope':             slope_ff        = float(txt[1])
                if txt[0] == 'curve':             curve_ff        = float(txt[1])
                if txt[0] == 'ratio':             ratio_ff        = float(txt[1])
                if txt[0] == 'prune':             prune_ff        = float(txt[1])
                if txt[0] == 'xmax':              xmax            = genf.clean_txt(txt[1]);  caution_opt_flame=False
                if txt[0] == 'restore_flame_folder': restore_flame_folder = genf.clean_txt(txt[1])
                # options for counterflow flames
                if txt[0] == 'mdot_min':          mdot_min = genf.txt2list_float(txt[1])
                if txt[0] == 'mdot_max':          mdot_max = genf.txt2list_float(txt[1])
                if txt[0] == 'mdot_incr':         mdot_incr = genf.txt2list_float(txt[1])
                if txt[0] == 'mdots':             mdots = genf.txt2list_float(txt[1])
                if txt[0] == 'mdot_min_1':        mdot_min = genf.txt2list_float(txt[1])
                if txt[0] == 'mdot_max_1':        mdot_max = genf.txt2list_float(txt[1])
                if txt[0] == 'mdot_incr_1':       mdot_incr = genf.txt2list_float(txt[1])
                if txt[0] == 'mdots_1':           mdots = genf.txt2list_float(txt[1])
                if txt[0] == 'mdot_2_min':        mdot_2_min = genf.txt2list_float(txt[1])
                if txt[0] == 'mdot_2_max':        mdot_2_max = genf.txt2list_float(txt[1])
                if txt[0] == 'mdot_2_incr':       mdot_2_incr = genf.txt2list_float(txt[1])
                if txt[0] == 'mdots_2':           mdots_2 = genf.txt2list_float(txt[1])
                if txt[0] == 'mixt_2':            fuel_2            = genf.clean_txt(txt[1])
                if txt[0] == 'fuel_2':            fuel_2 = genf.clean_txt(txt[1])
                if txt[0] == 'oxidant_2':         oxidant_2 = genf.clean_txt(txt[1])
                if txt[0] == 'diluent_2':         diluent_2 = genf.clean_txt(txt[1])
                if txt[0] == 'diluent_ratio_2':   diluent_ratio_2 = genf.diluent_fct(txt[1])
                if txt[0] == 'phi_min_2':         phi2_min           = genf.clean_txt(txt[1])
                if txt[0] == 'phi_max_2':         phi2_max           = genf.clean_txt(txt[1])
                if txt[0] == 'phi_incr_2':        phi2_incr          = genf.clean_txt(txt[1])
                if txt[0] == 'phis_2':            phis2              = genf.txt2list_float(txt[1])



                if txt[0] == 'width':             width = float(txt[1])

                txt = fs.readline().split('=')


            if load_conditions:
                if 'config' in locals():
                    if config=='reactor_UV' or config=='reactor_HP':
                        if config=='reactor_UV':
                            self.rB_reactor_UV[case_n].setChecked(True)
                        elif config=='reactor_HP':
                            self.rB_reactor_HP[case_n].setChecked(True)
                        self.add_conditions(case_n+1)
                        self.Gr[case_n].setGeometry(QtCore.QRect(700*sz_w, (10+(case_n)*181)*sz_h, 731*sz_w, 171*sz_h))
                        # options for reactor
                        if 'tol_ts' in locals():
                            self.tol_ts_rel_r[case_n].setPlainText(_translate("MainWindow", tol_ts.split(',')[0]))
                            self.tol_ts_abs_r[case_n].setPlainText(_translate("MainWindow", tol_ts.split(',')[1]))
                            del tol_ts
                        if 'n_pts' in locals():
                            self.pts_num[case_n].setProperty("value", int(n_pts));            del n_pts
                        if 'delta_npts' in locals():
                            self.delta_n_pts[case_n].setProperty("value", int(delta_npts));   del delta_npts
                        if 't_max_coeff' in locals():
                            self.t_max_coeff[case_n].setProperty("value", int(t_max_coeff));  del t_max_coeff
                        if 'Scal_ref' in locals():
                            self.ref_spec[case_n].setPlainText(_translate("MainWindow", Scal_ref));   del Scal_ref
                        if 'grad_curv_ratio' in locals():
                            self.grad_curv[case_n].setProperty("value", float(grad_curv_ratio));      del grad_curv_ratio
                        if 'tign_nPoints' in locals():
                            self.txt_Gr_ipn[case_n].setPlainText(_translate("MainWindow", tign_nPoints)); del tign_nPoints
                        if 'tign_dt' in locals():
                            self.txt_Gr_its[case_n].setPlainText(_translate("MainWindow", tign_dt));  del tign_dt

                    if config=='JSR':
                        self.rB_JSR[case_n].setChecked(True)
                        self.add_conditions(case_n+1)
                        self.GJSR[case_n].setGeometry(QtCore.QRect(700*sz_w, (10+(case_n)*181)*sz_h, 311*sz_w, 171*sz_h))
                        # options for JSR
                        if 'tol_ts' in locals():
                            self.tol_ts_rel_jsr[case_n].setPlainText(_translate("MainWindow", tol_ts.split(',')[0]))
                            self.tol_ts_abs_jsr[case_n].setPlainText(_translate("MainWindow", tol_ts.split(',')[1]))
                            del tol_ts
                        if 't_max' in locals():
                            self.t_max_jsr[case_n].setProperty("value", float(t_max)); del t_max

                    if config=='PFR':
                        self.rB_PFR[case_n].setChecked(True)
                        self.add_conditions(case_n+1)
                        self.Gpfr[case_n].setGeometry(QtCore.QRect(700*sz_w, (10+(case_n)*181)*sz_h, 810*sz_w, 171*sz_h))
                        # options for PFR
                        if 'tol_ts' in locals():
                            self.tol_ts_rel_pfr[case_n].setPlainText(_translate("MainWindow", tol_ts.split(',')[0]))
                            self.tol_ts_abs_pfr[case_n].setPlainText(_translate("MainWindow", tol_ts.split(',')[1]))
                            del tol_ts
                        if 'pfr_autom' in locals():
                            if pfr_autom: self.rB_pfr_1[case_n].setChecked(True)
                            else:         self.rB_pfr_2[case_n].setChecked(True)
                            del pfr_autom
                        if 'n_pts' in locals():
                            self.pts_num_pfr[case_n].setProperty("value", int(n_pts))
                            self.pfr_pts_num[case_n].setProperty("value", int(n_pts));                      del n_pts
                        if 'delta_npts' in locals():
                            self.pfr_delta_n_pts[case_n].setProperty("value", int(delta_npts));                 del delta_npts
                        if 'area' in locals():
                            self.pfr_area[case_n].setPlainText(_translate("MainWindow", str(area)));        del area
                        if 'length' in locals():
                            self.pfr_length[case_n].setPlainText(_translate("MainWindow", str(length)));    del length
                        if 'u_0' in locals():
                            self.pfr_u_0[case_n].setPlainText(_translate("MainWindow", str(u_0)));          del u_0
                        if 't_max_coeff' in locals():
                            self.pfr_t_max_coeff[case_n].setProperty("value", int(t_max_coeff));            del t_max_coeff
                        if 'Scal_ref' in locals():
                            self.pfr_ref_spec[case_n].setPlainText(_translate("MainWindow", Scal_ref));     del Scal_ref
                        if 'grad_curv_ratio' in locals():
                            self.pfr_grad_curv[case_n].setProperty("value", float(grad_curv_ratio));  del grad_curv_ratio

                    if 'flame' in config:
                        if config=='free_flame':
                            self.rB_fflame[case_n].setChecked(True)
                            self.add_conditions(case_n+1)
                            self.Gf[case_n].setGeometry(QtCore.QRect(700*sz_w, (10+(case_n)*181)*sz_h, 1000*sz_w, 171*sz_h))
                        if config=='diff_flame' or config=='pp_flame' or config=='tp_flame':
                            self.rB_cfflame[case_n].setChecked(True)
                            self.add_conditions(case_n+1)
                            self.Conditions_diff_flame[case_n].setGeometry(QtCore.QRect(351*sz_w, (10+(case_n)*181)*sz_h, 1050*sz_w, 171*sz_h))
                            self.Gf[case_n].setGeometry(QtCore.QRect(1401*sz_w, (10+(case_n)*181)*sz_h, 670*sz_w, 171*sz_h))
                        if config=='diff_flame': self.rB_cff_diff[case_n].setChecked(True)
                        if config=='pp_flame': self.rB_cff_pp[case_n].setChecked(True)
                        if config=='tp_flame': self.rB_cff_tp[case_n].setChecked(True)

                        # options for flame
                        if 'tol_ts' in locals():
                            self.tol_ts_rel_f[case_n].setPlainText(_translate("MainWindow", tol_ts.split(',')[0]))
                            self.tol_ts_abs_f[case_n].setPlainText(_translate("MainWindow", tol_ts.split(',')[1]))
                            del tol_ts
                        if 'tol_ss' in locals():
                            self.tol_ss_rel_f[case_n].setPlainText(_translate("MainWindow", tol_ss.split(',')[0]))
                            self.tol_ss_abs_f[case_n].setPlainText(_translate("MainWindow", tol_ss.split(',')[1]))
                            del tol_ss
                        if 'transport_model' in locals():
                            if transport_model=='Mix':
                                self.mix[case_n].setChecked(True)
                            elif transport_model=='Mult':
                                self.mult[case_n].setChecked(True)
                            del transport_model
                        if 'xmax' in locals():
                            self.txt_Gf_xmax[case_n].setPlainText(_translate("MainWindow", xmax)); del xmax
                        if 'pts_scatter' in locals():
                            self.pts_scatter[case_n].setPlainText(_translate("MainWindow", pts_scatter.replace('[','').replace(']','')))
                            del pts_scatter
                        if 'slope_ff' in locals():
                            self.slope_ff[case_n].setProperty("value", slope_ff);            del slope_ff
                        if 'curve_ff' in locals():
                            self.curve_ff[case_n].setProperty("value", curve_ff);            del curve_ff
                        if 'ratio_ff' in locals():
                            self.ratio_ff[case_n].setProperty("value", ratio_ff);            del ratio_ff
                        if 'prune_ff' in locals():
                            self.prune_ff[case_n].setProperty("value", prune_ff);            del prune_ff
                        if 'restore_flame_folder' in locals():
                            if restore_flame_folder != 'False':
                                self.rB_flame_restore[case_n].setChecked(True)
                                self.label_flame_folder[case_n].setText(_translate("MainWindow", restore_flame_folder))
                                self.label_flame_folder[case_n].setGeometry(QtCore.QRect(170*sz_w, 65*sz_h, 400*sz_w, 18*sz_h))
                            del restore_flame_folder

                        # options for counterflow flame
                        if 'fuel' in locals():
                            self.df_fuel_1[case_n].setPlainText(_translate("MainWindow", fuel));            del fuel
                        if 'oxidant' in locals():
                            self.df_oxidant_1[case_n].setPlainText(_translate("MainWindow", oxidant));      del oxidant
                        if 'diluent' in locals():
                            self.df_Diluent_1[case_n].setPlainText(_translate("MainWindow", diluent));      del diluent
#                        if 'diluent_ratio' in locals():
#                            self.df_Diluent_r_1[case_n].setPlainText(_translate("MainWindow", str(diluent_ratio))); del diluent_ratio
                        if 'mdot_min'         in locals(): self.df_mdot1_1[case_n].setPlainText(_translate("MainWindow",  mdot_min));    del mdot_min
                        if 'mdot_max'         in locals(): self.df_mdot2_1[case_n].setPlainText(_translate("MainWindow",  mdot_max));    del mdot_max
                        if 'mdot_incr'        in locals(): self.df_mdot3_1[case_n].setPlainText(_translate("MainWindow",  mdot_incr));   del mdot_incr
                        if 'mdots'            in locals():
                            self.df_mdot1_1[case_n].setPlainText(_translate("MainWindow",  format(mdots[0],'.0f')))
                            self.df_mdot2_1[case_n].setPlainText(_translate("MainWindow",  format(mdots[-1],'.0f')))
                            if len(mdots)==1:  div_ = 1
                            else:              div_ = len(mdots)-1
                            self.df_mdot4_1[case_n].setPlainText(_translate("MainWindow", format((mdots[-1]-mdots[0])/div_,'.0f')))
                            del mdots
                        if 'fuel_2' in locals():
                            self.df_fuel_2[case_n].setPlainText(_translate("MainWindow", fuel_2));           del fuel_2
                        if 'oxidant_2' in locals():
                            self.df_oxidant_2[case_n].setPlainText(_translate("MainWindow", oxidant_2));     del oxidant_2
                        if 'diluent_2' in locals():
                            self.df_Diluent_2[case_n].setPlainText(_translate("MainWindow", diluent_2));     del diluent_2
#                        if 'diluent_ratio_2' in locals():
#                            self.df_Diluent_r_2[case_n].setPlainText(_translate("MainWindow", str(diluent_ratio_2))); del diluent_ratio_2
                        if 'mdot_2_min'       in locals(): self.df_mdot1_1[case_n].setPlainText(_translate("MainWindow", mdot_2_min));    del mdot_2_min
                        if 'mdot_2_max'       in locals(): self.df_mdot2_1[case_n].setPlainText(_translate("MainWindow", mdot_2_max));    del mdot_2_max
                        if 'mdot_2_incr'      in locals(): self.df_mdot3_1[case_n].setPlainText(_translate("MainWindow", mdot_2_incr));   del mdot_2_incr
                        if 'mdots_2'            in locals():
                            self.df_mdot1_2[case_n].setPlainText(_translate("MainWindow",  format(mdots_2[0],'.0f')))
                            self.df_mdot2_2[case_n].setPlainText(_translate("MainWindow",  format(mdots_2[-1],'.0f')))
                            if len(mdots_2)==1:  div_ = 1
                            else:                div_ = len(mdots_2)-1
                            self.df_mdot4_2[case_n].setPlainText(_translate("MainWindow", format((mdots_2[-1]-mdots_2[0])/div_,'.0f')))
                            del mdots_2

#                        # options for flame
#                        if 'tol_ts' in locals():
#                            self.tol_ts_rel_f[case_n].setPlainText(_translate("MainWindow", tol_ts.split(',')[0]))
#                            self.tol_ts_abs_f[case_n].setPlainText(_translate("MainWindow", tol_ts.split(',')[1]))
#                            del tol_ts
#                        if 'tol_ss' in locals():
#                            self.tol_ss_rel_f[case_n].setPlainText(_translate("MainWindow", tol_ss.split(',')[0]))
#                            self.tol_ss_abs_f[case_n].setPlainText(_translate("MainWindow", tol_ss.split(',')[1]))
#                            del tol_ss
#                        if 'transport_model' in locals():
#                            if transport_model=='Mix':
#                                self.mix[case_n].setChecked(True)
#                            elif transport_model=='Mult':
#                                self.mult[case_n].setChecked(True)
#                            del transport_model
#                        if 'xmax' in locals():
#                            self.txt_Gf_xmax[case_n].setPlainText(_translate("MainWindow", xmax)); del xmax
#                        if 'pts_scatter' in locals():
#                            self.pts_scatter[case_n].setPlainText(_translate("MainWindow", pts_scatter.replace('[','').replace(']','')))
#                            del pts_scatter
#                        if 'slope_ff' in locals():
#                            self.slope_ff[case_n].setProperty("value", slope_ff);            del slope_ff
#                        if 'curve_ff' in locals():
#                            self.curve_ff[case_n].setProperty("value", curve_ff);            del curve_ff
#                        if 'ratio_ff' in locals():
#                            self.ratio_ff[case_n].setProperty("value", ratio_ff);            del ratio_ff
#                        if 'prune_ff' in locals():
#                            self.prune_ff[case_n].setProperty("value", prune_ff);            del prune_ff


                    # Main options
                    if 'fuel'          in locals(): self.fuel_1[case_n].setPlainText(_translate("MainWindow", fuel)); del fuel
                    if 'oxidant'       in locals(): self.oxidant_1[case_n].setPlainText(_translate("MainWindow", oxidant)); del oxidant
                    if 'diluent'       in locals(): self.Diluent_1[case_n].setPlainText(_translate("MainWindow", diluent)); del diluent

                    if 'diluent_ratio' in locals():
                        diluent_ratio = str(diluent_ratio).replace('[','').replace(']','').replace(',','').replace("'","")
                        self.Diluent_2[case_n].setPlainText(_translate("MainWindow", diluent_ratio))
                        self.df_Diluent_r_1[case_n].setPlainText(_translate("MainWindow", str(diluent_ratio)));    del diluent_ratio
                    if 'diluent_ratio_2' in locals():
                        diluent_ratio_2 = str(diluent_ratio_2).replace('[','').replace(']','').replace(',','').replace("'","")
                        self.df_Diluent_r_2[case_n].setPlainText(_translate("MainWindow", str(diluent_ratio_2)));    del diluent_ratio_2

                    if 'T_min'         in locals(): self.Tmin[case_n].setPlainText(_translate("MainWindow",  T_min));    del T_min
                    if 'T_max'         in locals(): self.Tmax[case_n].setPlainText(_translate("MainWindow",  T_max));    del T_max
                    if 'T_incr'        in locals(): self.Tincr[case_n].setPlainText(_translate("MainWindow", T_incr));   del T_incr
                    if 'Ts'            in locals():
                        self.Tmin[case_n].setPlainText(_translate("MainWindow",  format(Ts[0],'.0f')))
                        self.Tmax[case_n].setPlainText(_translate("MainWindow",  format(Ts[-1],'.0f')))
                        if len(Ts)==1:  div_ = 1
                        else:           div_ = len(Ts)-1
                        self.Tincr[case_n].setPlainText(_translate("MainWindow", format((Ts[-1]-Ts[0])/div_,'.0f')))
                        del Ts
                    if 'P_min'         in locals(): self.Pmin[case_n].setPlainText(_translate("MainWindow",  P_min));    del P_min
                    if 'P_max'         in locals(): self.Pmax[case_n].setPlainText(_translate("MainWindow",  P_max));    del P_max
                    if 'P_incr'        in locals(): self.Pincr[case_n].setPlainText(_translate("MainWindow", P_incr));   del P_incr
                    if 'Ps'            in locals():
                        self.Pmin[case_n].setPlainText(_translate("MainWindow",  format(Ps[0],'.1e')))
                        self.Pmax[case_n].setPlainText(_translate("MainWindow",  format(Ps[-1],'.1e')))
                        if len(Ps)==1:  div_ = 1
                        else:           div_ = len(Ps)-1
                        self.Pincr[case_n].setPlainText(_translate("MainWindow", format((Ps[-1]-Ps[0])/div_,'.1e')))
                        del Ps
                    if config=='diff_flame' or config=='pp_flame' or config=='tp_flame':
                        if 'phi_min'       in locals(): self.df_eqmin_1[case_n].setPlainText(_translate("MainWindow", phi_min));  del phi_min
                        if 'phi_max'       in locals(): self.df_eqmax_1[case_n].setPlainText(_translate("MainWindow", phi_max));  del phi_max
                        if 'phi_incr'      in locals(): self.df_eqincr_1[case_n].setPlainText(_translate("MainWindow",phi_incr)); del phi_incr
                        if 'phis'            in locals():
                            self.df_eqmin_1[case_n].setPlainText(_translate("MainWindow",  format(phis[0],'.2f')))
                            self.df_eqmax_1[case_n].setPlainText(_translate("MainWindow",  format(phis[-1],'.2f')))
                            if len(phis)==1:    div_ = 1
                            else:               div_ = len(phis)-1
                            self.df_eqincr_1[case_n].setPlainText(_translate("MainWindow", format((phis[-1]-phis[0])/div_,'.2f')))
                            del phis
                        if 'phi2_min'       in locals(): self.df_eqmin_2[case_n].setPlainText(_translate("MainWindow", phi_min));  del phi2_min
                        if 'phi2_max'       in locals(): self.df_eqmax_2[case_n].setPlainText(_translate("MainWindow", phi_max));  del phi2_max
                        if 'phi2_incr'      in locals(): self.df_eqincr_2[case_n].setPlainText(_translate("MainWindow",phi_incr)); del phi2_incr
                        if 'phis2'            in locals():
                            self.df_eqmin_2[case_n].setPlainText(_translate("MainWindow",  format(phis2[0],'.2f')))
                            self.df_eqmax_2[case_n].setPlainText(_translate("MainWindow",  format(phis2[-1],'.2f')))
                            if len(phis2)==1:    div_ = 1
                            else:               div_ = len(phis2)-1
#                            self.df_eqincr_2[case_n].setPlainText(_translate("MainWindow", format((phis2[-1]-phis2[0])/div_,'.2f')))
                            del phis2
                    else:
                        if 'phi_min'       in locals(): self.eqmin[case_n].setPlainText(_translate("MainWindow", phi_min));  del phi_min
                        if 'phi_max'       in locals(): self.eqmax[case_n].setPlainText(_translate("MainWindow", phi_max));  del phi_max
                        if 'phi_incr'      in locals(): self.eqincr[case_n].setPlainText(_translate("MainWindow",phi_incr)); del phi_incr
                        if 'phis'            in locals():
                            self.eqmin[case_n].setPlainText(_translate("MainWindow",  format(phis[0],'.2f')))
                            self.eqmax[case_n].setPlainText(_translate("MainWindow",  format(phis[-1],'.2f')))
                            if len(phis)==1:    div_ = 1
                            else:               div_ = len(phis)-1
                            self.eqincr[case_n].setPlainText(_translate("MainWindow", format((phis[-1]-phis[0])/div_,'.2f')))
                            del phis

        while '> Op:' not in txt[-1] and txt[0] != '':
            txt = fs.readline().split('=')

        global condition_tab_removed

        # get data
        red_data_list = [] ; save_op = False ; op_n=0
        while txt[0] != '':
            while '> Op:' not in txt[-1] and txt[0] != '':
                txt[0]=txt[0].replace(' ','')

                if txt[0]== 'operator':
                    reduction_operator = genf.clean_txt(txt[1])
                if txt[0] == 'eps':
                    eps = genf.txt2list_float(txt[1])
                    eps = eps[0]
                if txt[0] == 'delta_eps':
                    delta_eps = genf.txt2list_float(txt[1])
                    delta_eps = delta_eps[0]
                if txt[0] == 'n_points':
                    n_points           = float(txt[1])
                if txt[0] == 'max_error_sp':
                    max_error_sp = genf.txt2list_float(txt[1])
                    while len(max_error_sp)<n_tspc: max_error_sp.append(max_error_sp[-1])
                if txt[0] == 'max_error_T':       max_error_T        = float(txt[1])
                if txt[0] == 'max_error_ig':      max_error_ig       = float(txt[1])
                if txt[0] == 'max_error_Sl':      max_error_Sl       = float(txt[1])
                if txt[0] == 'max_error_K':       max_error_K        = float(txt[1])
                if txt[0] == 'inter_sp_inter':    inter_sp_inter     = genf.str2bool(txt[1])
                if txt[0] == 'optim':
                    optim              = genf.clean_txt(txt[1])
                    if 'False' in optim: optim = False

                if txt[0] == 'ttol_sensi':
                    try:    ttol_sensi = genf.txt2list_float(txt[1])
                    except: ttol_sensi = genf.txt2list_bool(txt[1])

                if 'optim' in locals() and optim != False:
                    if txt[0] == 'n_gen':               n_gen              = int(txt[1])
                    if txt[0] == 'n_it':                n_gen              = int(txt[1])
                    if txt[0] == 'n_indiv':             n_indiv            = int(txt[1])
                    if txt[0] == 'error_fitness':       error_fitness      = genf.clean_txt(txt[1])
                    if txt[0] == 'Arrh_max_variation':  Arrh_max_variation = genf.txt2list_float(txt[1])
                    if txt[0] == 'optim_on_meth':       optim_on_meth      = genf.clean_txt(txt[1])
                    if txt[0] == 'nb_r2opt':            nb_r2opt           = int(txt[1])
                    if txt[0] == 'selection_operator':  selection_operator = genf.clean_txt(txt[1])
                    if txt[0] == 'selection_options':   selection_options  = genf.clean_txt(txt[1])
                    if txt[0] == 'Xover_operator':      Xover_operator     = genf.txt2list_string(txt[1])
                    if txt[0] == 'Xover_pct':           Xover_pct          = genf.txt2list_string(txt[1])
                    if txt[0] == 'Xover_opt':           Xover_opt          = genf.txt2list_string(txt[1])
                    if txt[0] == 'mut_operator':        mut_operator       = genf.txt2list_string(txt[1])
                    if txt[0] == 'mut_pct':             mut_pct            = genf.txt2list_string(txt[1])
                    if txt[0] == 'mut_opt':             mut_opt            = genf.txt2list_string(txt[1])
                    if txt[0] == 'mut_intensity':       mut_intensity      = float(txt[1])
                    if txt[0] == 'sub_mech_sel':        sub_mech_sel       = genf.txt2list_string(txt[1])
                    if txt[0] == 'inertia_score':       inertia_score      = genf.str2bool(txt[1])
                    if txt[0] == 'inertia_min':         inertia_min        = float(txt[1])
                    if txt[0] == 'inertia_max_i'\
                    or txt[0] == 'inertia_i':           inertia_i          = float(txt[1])
                    if txt[0] == 'inertia_max_endi'\
                    or txt[0] == 'inertia_end':         inertia_end        = float(txt[1])
                    if txt[0] == 'cognitive_accel_i':   cognitive_accel_i   = float(txt[1])
                    if txt[0] == 'score':               score_integ         = str2bool(txt[1])
                    if txt[0] == 'cognitive_accel_end': cognitive_accel_end = float(txt[1])
                    if txt[0] == 'social_accel_i':      social_accel_i      = float(txt[1])
                    if txt[0] == 'social_accel_end':    social_accel_end    = float(txt[1])



                txt = fs.readline().split('=')
                save_op = True


            # store data
            if save_op:

                if    reduction_operator=='DRG_sp':   self.DRG_clic('DRG_sp')
                elif  reduction_operator=='DRGEP_sp': self.DRG_clic('DRGEP_sp')
                elif  reduction_operator=='DRG_r':    self.DRG_clic('DRG_r')
                elif  reduction_operator=='SAR_sp':    self.SA_clic('SAR_sp')
                elif  reduction_operator=='SARGEP_sp':  self.SA_clic('SARGEP_sp')
                elif  reduction_operator=='SAR_r':     self.SA_clic('SAR_r')
                if reduction_operator != 'NULL' and 'ext_results_file' in locals():
                    if condition_tab_removed == False:
                        self.tablet.removeTab(1)
                        condition_tab_removed = True

                tspc_red = copy.deepcopy(tspc)
                if 'max_error_sp' not in locals():  max_error_sp = []
                errors   = copy.deepcopy(max_error_sp)
                num_row = len(tspc_red)
                if T_check:
                    num_row+=1 ; tspc_red.insert(0,'T')
                if Sl_check:
                    num_row+=1 ; tspc_red.insert(0,'Sl')
                if ig_check:
                    num_row+=1 ; tspc_red.insert(0,'Ig_t')
                if K_check:
                    num_row+=1 ; tspc_red.insert(0,'K')


                if 'SAR' in reduction_operator:
                    if 'eps'            in locals():
                        self.num_SA_eps[-1].setProperty("value", eps);         del eps
                    if 'delta_eps'      in locals():
                        self.num_SA_deps[-1].setProperty("value", delta_eps);  del delta_eps
                    if 'n_points'       in locals():
                        self.num_SA_pt_num[-1].setProperty("value", n_points); del n_points
                    if 'inter_sp_inter'       in locals():
                        if inter_sp_inter:  self.cB_ISI_sa[-1].setChecked(True)
                        else:               self.cB_ISI_sa[-1].setChecked(False)
                    if 'ttol_sensi'       in locals():
                        if ttol_sensi[0]:
                            self.rB_SA_tol_1[-1].setChecked(True)
                            self.SA_tol_appear('in')
                            self.num_SA_tol_1[-1].setPlainText(_translate("MainWindow", str(ttol_sensi[0])))
                            self.num_SA_tol_2[-1].setPlainText(_translate("MainWindow", str(ttol_sensi[1])))
                        else: self.rB_SA_tol_2[-1].setChecked(True)
                        del ttol_sensi

                    if T_check:
                        if 'max_error_T' in locals():  errors.insert(0,max_error_T);  del max_error_T
                        else:                          errors.insert(0,30)
                    if Sl_check:
                        if 'max_error_Sl' in locals(): errors.insert(0,max_error_Sl); del max_error_Sl
                        else:                          errors.insert(0,30)
                    if ig_check:
                        if 'max_error_ig' in locals(): errors.insert(0,max_error_ig); del max_error_ig
                        else:                          errors.insert(0,30)
                    if K_check:
                        if 'max_error_K' in locals():  errors.insert(0,max_error_K); del max_error_K
                        else:                          errors.insert(0,30)

                    self.tableWidget_SA[-1].setColumnCount(2)
                    self.tableWidget_SA[-1].setRowCount(len(tspc_red))

                    _translate = QtCore.QCoreApplication.translate
                    # row titles :
                    for r in range(len(tspc_red)):
                        item = QtWidgets.QTableWidgetItem()
                        self.tableWidget_SA[-1].setVerticalHeaderItem(r, item)
                        item = self.tableWidget_SA[-1].verticalHeaderItem(r)
                        item.setText(_translate("MainWindow", str(r+1)))
                    # column titles :
                    col_title = ["Target","Error (%)"]
                    for col in range(2):
                        item = QtWidgets.QTableWidgetItem()
                        self.tableWidget_SA[-1].setHorizontalHeaderItem(col, item)
                        item = self.tableWidget_SA[-1].horizontalHeaderItem(col)
                        item.setText(_translate("MainWindow", col_title[col]))

                    # fill the table
                    for col in range(2):
                        if col==0: txt = copy.deepcopy(tspc_red)
                        else:      txt = copy.deepcopy(errors)
                        for r in range(len(tspc_red)):
                            item = QtWidgets.QTableWidgetItem()
                            self.tableWidget_SA[-1].setItem(r, col, item)
                            item = self.tableWidget_SA[-1].item(r, col)
                            item.setText(_translate("MainWindow", str(txt[r])))

                    idx_tab = self.tablet.indexOf(self.SA[-1])
                    self.tablet.setCurrentIndex(idx_tab)

                if 'DRG' in reduction_operator:
                    if 'eps'            in locals():
                        self.num_DRG_eps[-1].setProperty("value", eps);        del eps
                    if 'delta_eps'      in locals():
                        self.num_DRG_deps[-1].setProperty("value", delta_eps); del delta_eps
                    if 'n_points'       in locals():
                        self.num_DRG_pt_num[-1].setProperty("value", n_points);del n_points
                    if 'inter_sp_inter'       in locals():
                        if inter_sp_inter:  self.cB_ISI_drg[-1].setChecked(True)
                        else:               self.cB_ISI_drg[-1].setChecked(True)
                        del inter_sp_inter

                    if T_check:
                        if 'max_error_T' in locals():  errors.insert(0,max_error_T); del max_error_T
                        else:                          errors.insert(0,30)
                    if Sl_check:
                        if 'max_error_Sl' in locals(): errors.insert(0,max_error_Sl); del max_error_Sl
                        else:                          errors.insert(0,30)
                    if ig_check:
                        if 'max_error_ig' in locals(): errors.insert(0,max_error_ig); del max_error_ig
                        else:                          errors.insert(0,30)
                    if K_check:
                        if 'max_error_K' in locals():  errors.insert(0,max_error_K); del max_error_K
                        else:                          errors.insert(0,30)

                    self.tableWidget_DRG[-1].setColumnCount(2)
                    self.tableWidget_DRG[-1].setRowCount(len(tspc_red))

                    _translate = QtCore.QCoreApplication.translate
                    # row titles :
                    for r in range(len(tspc_red)):
                        item = QtWidgets.QTableWidgetItem()
                        self.tableWidget_DRG[-1].setVerticalHeaderItem(r, item)
                        item = self.tableWidget_DRG[-1].verticalHeaderItem(r)
                        item.setText(_translate("MainWindow", str(r+1)))
                    # column titles :
                    col_title = ["Target","Error (%)"]
                    for col in range(2):
                        item = QtWidgets.QTableWidgetItem()
                        self.tableWidget_DRG[-1].setHorizontalHeaderItem(col, item)
                        item = self.tableWidget_DRG[-1].horizontalHeaderItem(col)
                        item.setText(_translate("MainWindow", col_title[col]))

                    # fill the table
                    for col in range(2):
                        if col==0: txt = copy.deepcopy(tspc_red)
                        else:      txt = copy.deepcopy(errors)
                        for r in range(len(tspc_red)):
                            item = QtWidgets.QTableWidgetItem()
                            self.tableWidget_DRG[-1].setItem(r, col, item)
                            item = self.tableWidget_DRG[-1].item(r, col)
                            item.setText(_translate("MainWindow", str(txt[r])))

                    idx_tab = self.tablet.indexOf(self.DRG[-1])
                    self.tablet.setCurrentIndex(idx_tab)

            if 'optim' in locals():
                if optim:
                    if reduction_operator == 'NULL' and 'ext_results_file' in locals():
                        self.tablet.setCurrentIndex(0)
                        self.opt_clic(optim)
                        if condition_tab_removed == False:
                            self.tablet.removeTab(1)
                            condition_tab_removed = True
                        self.tablet.setCurrentIndex(0)
                    else:
                        if optim == True: optim = 'GA'
                        if reduction_operator == 'NULL':
                            self.opt_clic(optim, True)
                        else:
                            self.opt_clic(optim)
                    if 'n_gen'              in locals():
                        self.num_GA_gen[-1].setProperty("value", n_gen); del n_gen
                    if 'n_indiv'            in locals():
                        self.num_GA_ind[-1].setProperty("value", n_indiv); del n_indiv
                    if 'Arrh_max_variation' in locals():
                        self.num_GA_A[-1].setProperty("value",  Arrh_max_variation[0])
                        self.num_GA_n[-1].setProperty("value",  Arrh_max_variation[1])
                        self.num_GA_Ea[-1].setProperty("value", Arrh_max_variation[2])
                        del Arrh_max_variation
                    if 'optim_on_meth' in locals():
                        if reduction_operator != 'NULL':
                            if optim_on_meth == 'True':   self.cB_GA_meth[-1].setChecked(True)
                            else :                        self.cB_GA_meth[-1].setChecked(False)
                        else:
                            if   optim_on_meth == 'SA':   self.cB_GA_meth_SA[-1].setChecked(True)
                            elif optim_on_meth == 'DRG':  self.cB_GA_meth_DRG[-1].setChecked(True)
                        del optim_on_meth
                    if 'nb_r2opt'in locals():
                        self.num_GA_meth_fract[-1].setProperty("value", nb_r2opt)
                        del nb_r2opt
                    if 'sub_mech_sel' in locals():
                        if 'H2' in sub_mech_sel: self.cB_sub_H[-1].setChecked(True)
                        else:                    self.cB_sub_H[-1].setChecked(False)
                        if 'CO' in sub_mech_sel: self.cB_sub_CO[-1].setChecked(True)
                        else:                    self.cB_sub_CO[-1].setChecked(False)
                        if 'N'  in sub_mech_sel: self.cB_sub_N[-1].setChecked(True)
                        else:                    self.cB_sub_N[-1].setChecked(False)
                        if 'S'  in sub_mech_sel: self.cB_sub_S[-1].setChecked(True)
                        else:                    self.cB_sub_S[-1].setChecked(False)
                        if 'Si' in sub_mech_sel: self.cB_sub_Si[-1].setChecked(True)
                        else:                    self.cB_sub_Si[-1].setChecked(False)
                        for sC in range(max(self.mech_data.react.subm_C)):
                            if 'C'+str(sC+1) in sub_mech_sel: self.cB_sub_C[-1][sC].setChecked(True)
                            else:                             self.cB_sub_C[-1][sC].setChecked(False)
                        del sub_mech_sel
                    if 'error_fitness'      in locals():
                        if error_fitness=='mean': self.rB_GA_fit_1[-1].setChecked(True)
                        else:                     self.rB_GA_fit_2[-1].setChecked(True)
                        del error_fitness

                    # GA options
                    if 'selection_operator' in locals():
                        if selection_operator=='Roulette':
                            self.Box_GA_Selection[-1].setCurrentIndex(0)
                        if selection_operator=='Rank':
                            self.Box_GA_Selection[-1].setCurrentIndex(1)
                        if selection_operator=='Geometric_norm':
                            self.Box_GA_Selection[-1].setCurrentIndex(2)
                        if selection_operator=='Elitism':
                            self.Box_GA_Selection[-1].setCurrentIndex(3)
                        del selection_operator
                    if 'selection_options'  in locals():
                        selection_options=selection_options.replace('[','').replace(']','')
                        self.txt_GA_sel_opt[-1].setPlainText(_translate("MainWindow", selection_options))
                        del selection_options

                    X_nb = 0
                    if 'Xover_operator' in locals():
                        if 'simple_Xover' in Xover_operator:
                            self.cB_Xover_op_1[-1].setChecked(True)
                            if 'Xover_pct' in locals():
                                self.txt_GA_Xover_int_1[-1].setPlainText(_translate("MainWindow", Xover_pct[X_nb])) ; X_nb+=1
                        else: self.cB_Xover_op_1[-1].setChecked(False)
                        if 'multiple_Xover' in Xover_operator:
                            self.cB_Xover_op_2[-1].setChecked(True)
                            if 'Xover_pct' in locals():
                                self.txt_GA_Xover_int_2[-1].setPlainText(_translate("MainWindow", Xover_pct[X_nb])) ; X_nb+=1
                        else: self.cB_Xover_op_2[-1].setChecked(False)
                        if 'arith_Xover' in Xover_operator:
                            self.cB_Xover_op_3[-1].setChecked(True)
                            if 'Xover_pct' in locals():
                                self.txt_GA_Xover_int_3[-1].setPlainText(_translate("MainWindow", Xover_pct[X_nb])) ; X_nb+=1
                        else: self.cB_Xover_op_3[-1].setChecked(False)
                        if 'heuristic_Xover' in Xover_operator:
                            self.cB_Xover_op_4[-1].setChecked(True)
                            if 'Xover_pct' in locals():
                                self.txt_GA_Xover_int_4[-1].setPlainText(_translate("MainWindow", Xover_pct[X_nb])) ; X_nb+=1
                        else: self.cB_Xover_op_4[-1].setChecked(False)
                        del Xover_operator

                    mut_nb = 0
                    if 'mut_operator' in locals():
                        if 'uniform_mutation' in mut_operator:
                            self.cB_mut_op_1[-1].setChecked(True)
                            if 'mut_opt' in locals():
                                self.txt_GA_mut_opt_1[-1].setPlainText(mut_opt[mut_nb])
                            if 'mut_pct' in locals():
                                self.txt_GA_mut_int_1[-1].setPlainText(_translate("MainWindow", mut_pct[mut_nb])) ; mut_nb+=1
                        else: self.cB_mut_op_1[-1].setChecked(False)
                        if 'non_uniform_mutation' in mut_operator:
                            self.cB_mut_op_2[-1].setChecked(True)
                            if 'mut_opt' in locals():
                                self.txt_GA_mut_opt_3[-1].setPlainText(mut_opt[mut_nb])
                            if 'mut_pct' in locals():
                                self.txt_GA_mut_int_2[-1].setPlainText(_translate("MainWindow", mut_pct[mut_nb])) ; mut_nb+=1
                        else: self.cB_mut_op_2[-1].setChecked(False)
                        if 'boundary_mutation' in mut_operator:
                            self.cB_mut_op_3[-1].setChecked(True)
                            if 'mut_opt' in locals():
                                self.txt_GA_mut_opt_3[-1].setPlainText(mut_opt[mut_nb])
                            if 'mut_pct' in locals():
                                self.txt_GA_mut_int_3[-1].setPlainText(_translate("MainWindow", mut_pct[mut_nb])) ; mut_nb+=1
                        else: self.cB_mut_op_3[-1].setChecked(False)
                        del mut_operator
                    if 'mut_intensity' in locals():
                        self.num_GA_mut_prob[-1].setProperty("value", mut_intensity)
                        del mut_intensity

                    # PSO options
                    if 'inertia_score' in locals():
                        self.Inertia_score[-1].setChecked(inertia_score)
                    if 'inertia_min' in locals():
                        self.Inertia_min[-1].setProperty("value", inertia_min)
                    if 'inertia_max_i' in locals():
                        self.Inertia_i[-1].setProperty("value", inertia_max_i)
                    if 'inertia_i' in locals():
                        self.Inertia_i[-1].setProperty("value", inertia_i)
                    if 'inertia_max_end' in locals():
                        self.Inertia_end[-1].setProperty("value", inertia_max_end)
                    if 'inertia_end' in locals():
                        self.Inertia_end[-1].setProperty("value", inertia_end)
                    if 'cognitive_accel_i' in locals():
                        self.Cogn_i[-1].setProperty("value", cognitive_accel_i)
                    if 'cognitive_accel_end' in locals():
                        self.Cogn_end[-1].setProperty("value", cognitive_accel_end)
                    if 'social_accel_i' in locals():
                        self.Social_i[-1].setProperty("value", social_accel_i)
                    if 'social_accel_end' in locals():
                        self.Social_end[-1].setProperty("value", social_accel_end)



                del optim
                save_op = False

            txt = fs.readline().split('=')
        self.tablet.setCurrentIndex(0)


    def run_reduction(self):
        self.write_parameters(False)
        genf.main_redopt_algo('last_condition.inp',self.WD_path)



def get_wd_folder(WD_p=True, WD_n=True):
    global sz_w
    global sz_h
    global WD_path

    app = QtWidgets.QApplication(sys.argv)
    Window_f = QtWidgets.QMainWindow()

    _height = QtWidgets.QDesktopWidget().screenGeometry(-1).height()
    _width = QtWidgets.QDesktopWidget().screenGeometry(-1).width()
    sz_w = _width/1366
    sz_h = _height/768

    Window_f.setObjectName("Brookesia")
    Window_f.resize(790*sz_w, 671*sz_h)

    A=Files_windows()
    options = QtWidgets.QFileDialog.Options()
    options |= QtWidgets.QFileDialog.DontUseNativeDialog

    WD_path = False ; WD_name = False
    if WD_p:
        WD_path = QtWidgets.QFileDialog.getExistingDirectory(A, "Select the new working directory", options=options)
    if WD_n:
        WD_name, ok = QtWidgets.QInputDialog.getText(A, 'Text Input Dialog', 'Enter the new name:')

    return WD_path, WD_name





def list2txt(_list):
    _txt = str(_list)
    _txt  = _txt.replace('[','')
    _txt  = _txt.replace(']','')
    _txt  = _txt.replace('\n','')
    _txt  = _txt.replace("'","")
    _txt  = _txt.replace('"','')
    return _txt

class Files_windows(QtWidgets.QWidget):

    def __init__(self):
        global sz_w
        global sz_h
        super().__init__()
        self.title = 'PyQt5 file dialogs - pythonspot.com'
        self.left = 10
        self.top = 10
        self.width = 640*sz_w
        self.height = 480*sz_h

    def initUI(self):
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)
        self.show()


def display_gui(WD_path, WD_name, root_path):
#if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow, WD_path, WD_name, root_path)
    MainWindow.show()
    sys.exit(app.exec_())

