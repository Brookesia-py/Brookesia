def write_conditions():
    import os
    if not os.path.exists('_conditions_input'):  os.mkdir('_conditions_input')
    os.chdir('_conditions_input')

    fichier_data = open('1_reactor.inp', 'w')
    fichier_data.write("""#=============================================
#           Main parameters
#=============================================
main_path         = TEST_1_reactor
mech              = C1_GRI30.cti
verbose           = 4
show_plots        = False
write_ck          = True
tspc              = CH4, CO, CO2
T_check           = True
sp_T              = CO2
Sl_check          = False
sp_Sl             = H
ig_check          = True
sp_ig             = CH3
K_check           = False
sp_K              = H
error_calculation = points
error_coupling    = mean


#=============================================
#           Simulation cases
#=============================================

#======> Case 1
config            = reactor_UV
Ps                = 100000.0
fuel              = CH4
oxidant           = O2
diluent           = N2
diluent_ratio     = N2/O2 3.76
Ts                = 1600.0
phis              = 0.5, 1.0, 1.5
n_pts             = 250.0
delta_npts        = 20.0
t_max_coeff       = 5.0
Scal_ref          = H2O
grad_curv_ratio   = 0.5
tign_nPoints      = 450.0
tign_dt           = 1e-09
tol_ts            = 1e-06, 1e-12



#=============================================
#           Operators
#=============================================


#===========> Op: DRGEP_sp
operator        = DRGEP_sp 
eps             = 0.02
delta_eps       = 0.01
n_points        = 10.0
max_error_sp    = 30, 30, 30
max_error_T     = 30
max_error_ig    = 30
inter_sp_inter  = True
optim           = GA
#====> Genetic Algorithm Optimization
n_gen               = 5
n_indiv             = 5
error_fitness       = mean
Arrh_max_variation  = 5, 5, 5
optim_on_meth       = False
sub_mech_sel        = H2, CO, C1, C2, C3, N
selection_operator  = Roulette
selection_options   = 0.2
Xover_operator      = simple_Xover, multiple_Xover, arith_Xover, heuristic_Xover
Xover_pct           = 10, 20, 20, 20
mut_operator        = uniform_mutation, non_uniform_mutation, boundary_mutation
mut_pct             = 30, 30, 10
mut_opt             = , 3, 
mut_intensity       = 20


#===========> Op: SAR_sp
operator        = SAR_sp 
eps             = 0.02
delta_eps       = 0.01
n_points        = 10.0
max_error_sp    = 30.0, 30.0, 30.0
max_error_T     = 30.0
max_error_ig    = 30.0
inter_sp_inter  = True
ttol_sensi      = 1e-05, 1e-08
optim           = GA
#====> Genetic Algorithm Optimization
n_gen               = 5
n_indiv             = 5
error_fitness       = mean
Arrh_max_variation  = 5, 5, 5
optim_on_meth       = False
sub_mech_sel        = H2, CO, C1, C2, C3, N
selection_operator  = Roulette
selection_options   = 0.2
Xover_operator      = simple_Xover, multiple_Xover, arith_Xover, heuristic_Xover
Xover_pct           = 10, 20, 20, 20
mut_operator        = uniform_mutation, non_uniform_mutation, boundary_mutation
mut_pct             = 30, 30, 10
mut_opt             = , 3, 
mut_intensity       = 20
    """)

    # ================================================================
    # ================================================================
    # ================================================================
    # ================================================================

    fichier_data = open('2_JSR.inp', 'w')
    fichier_data.write("""#=============================================
#           Main parameters
#=============================================
main_path         = TEST_2_JSR
mech              = C1_GRI30.cti
verbose           = 4
show_plots        = False
write_ck          = True
tspc              = OH, CO2, CH4
T_check           = True
sp_T              = CO2
Sl_check          = False
sp_Sl             = H
ig_check          = False
sp_ig             = CH3
K_check           = False
sp_K              = H
error_calculation = points
error_coupling    = mean


#=============================================
#           Simulation cases
#=============================================

#======> Case 1
config            = JSR
Ps                = 100000.0
fuel              = CH4
oxidant           = O2
diluent           = N2
diluent_ratio     = N2/O2 3.76
Ts                = 600.0, 610.0, 620.0, 630.0, 640.0, 650.0, 660.0, 670.0, 680.0, 690.0, 700.0, 710.0, 720.0, 730.0, 740.0, 750.0, 760.0, 770.0, 780.0, 790.0, 800.0, 810.0, 820.0, 830.0, 840.0, 850.0, 860.0, 870.0, 880.0, 890.0, 900.0, 910.0, 920.0, 930.0, 940.0, 950.0, 960.0, 970.0, 980.0, 990.0, 1000.0, 1010.0, 1020.0, 1030.0, 1040.0, 1050.0, 1060.0, 1070.0, 1080.0, 1090.0, 1100.0, 1110.0, 1120.0, 1130.0, 1140.0, 1150.0, 1160.0, 1170.0, 1180.0, 1190.0, 1200.0
phis              = 0.5, 1.0, 1.5
t_max             = 0.2
tol_ts            = 1e-12, 1e-18



#=============================================
#           Operators
#=============================================


#===========> Op: DRGEP_sp
operator        = DRGEP_sp 
eps             = 0.02
delta_eps       = 0.01
n_points        = 10.0
max_error_sp    = 30, 30, 30
max_error_T     = 30
inter_sp_inter  = True

#===========> Op: DRG_r
operator        = DRG_r 
eps             = 0.02
delta_eps       = 0.0
n_points        = 10.0
max_error_sp    = 30.0, 30.0, 30.0
max_error_T     = 30.0
inter_sp_inter  = True
optim           = GA
#====> Particle Swarm Optimization
n_it                = 5
n_indiv             = 5
error_fitness       = mean
Arrh_max_variation  = 5, 5, 5
optim_on_meth       = False
sub_mech_sel        = H2, CO, C1, C2, C3, N
selection_operator  = Roulette
selection_options   = 0.2
Xover_operator      = simple_Xover, multiple_Xover, arith_Xover, heuristic_Xover
Xover_pct           = 10, 20, 20, 20
mut_operator        = uniform_mutation, non_uniform_mutation, boundary_mutation
mut_pct             = 30, 30, 10
mut_opt             = , 3, 
mut_intensity       = 30
    """)

    # ================================================================
    # ================================================================
    # ================================================================
    # ================================================================

    fichier_data = open('3_Freeflame.inp', 'w')
    fichier_data.write("""#=============================================
#           Main parameters
#=============================================
main_path         = TEST_3_Freeflame_pts
mech              = C1_GRI30.cti
verbose           = 4
show_plots        = False
write_ck          = True
tspc              = CH4, CO, CO2
T_check           = True
sp_T              = CO2
Sl_check          = True
sp_Sl             = H
ig_check          = False
sp_ig             = CH3
K_check           = False
sp_K              = H
error_calculation = QoI
error_coupling    = mean


#=============================================
#           Simulation cases
#=============================================

#======> Case 1
config            = free_flame
Ps                = 100000.0
fuel              = CH4
oxidant           = O2
diluent           = N2
diluent_ratio     = N2/O2 3.76
Ts                = 300.0
phis              = 0.5, 1.0, 1.5
xmax              = 0.02
tol_ts            = 1e-05, 1e-08
tol_ss            = 1e-06, 1e-08
transport_model   = Mix
pts_scatter       = [0.0,0.03,0.3,0.5,0.7,1.0]
slope             = 0.05
curve             = 0.05
ratio             = 2.0
prune             = 0.01
restore_flame_folder = ex3_ffl



#=============================================
#           Operators
#=============================================


#===========> Op: DRGEP_sp
operator        = DRGEP_sp 
eps             = 0.02
delta_eps       = 0.01
n_points        = 10.0
max_error_sp    = 30, 30, 30
max_error_T     = 30
max_error_Sl    = 30
inter_sp_inter  = True

#===========> Op: SARGEP_sp
operator        = SARGEP_sp 
eps             = 0.02
delta_eps       = 0.01
n_points        = 10.0
max_error_sp    = 30.0, 30.0, 30.0
max_error_T     = 30.0
max_error_Sl    = 30.0
inter_sp_inter  = True
ttol_sensi      = False, False
optim           = GA
#====> Genetic Algorithm Optimization
n_gen               = 5
n_indiv             = 5
error_fitness       = mean
Arrh_max_variation  = 5, 5, 5
optim_on_meth       = False
sub_mech_sel        = H2, CO, C1, C2, C3, N
selection_operator  = Roulette
selection_options   = 0.2
Xover_operator      = simple_Xover, multiple_Xover, arith_Xover, heuristic_Xover
Xover_pct           = 10, 20, 20, 20
mut_operator        = uniform_mutation, non_uniform_mutation, boundary_mutation
mut_pct             = 30, 30, 10
mut_opt             = , 3, 
mut_intensity       = 20
    """)


    # ================================================================
    # ================================================================
    # ================================================================
    # ================================================================

    fichier_data = open('4_diff.inp', 'w')
    fichier_data.write("""#=============================================
#           Main parameters
#=============================================
main_path         = TEST_4_diff_pts
mech              = C0_Konnov.cti
verbose           = 4
show_plots        = False
write_ck          = True
tspc              = H2, H2O, H
T_check           = True
sp_T              = H2O
Sl_check          = True
sp_Sl             = H
ig_check          = True
sp_ig             = H
K_check           = True
sp_K              = H
error_calculation = QoI
error_coupling    = mean


#=============================================
#           Simulation cases
#=============================================

#======> Case 1
config            = diff_flame
Ps                = 100000.0
fuel_1            = H2
diluent_1         = N2
diluent_ratio_1   = 0.0
Ts_1              = 300.0
phis_1            = 0.5, 1.0, 1.5
mdots_1           = 2.0
oxidant_2         = O2
diluent_2         = N2
diluent_ratio_2   = 79.0
Ts_2              = 300.0
mdots_2           = 2.0
width             = 0.02
tol_ts            = 1e-05, 1e-08
tol_ss            = 1e-06, 1e-08
transport_model   = Mix
pts_scatter       = [0.0,0.2,0.4,0.6,0.8,1.0]
slope             = 0.05
curve             = 0.05
ratio             = 2.0
prune             = 0.01
restore_flame_folder = False



#=============================================
#           Operators
#=============================================


#===========> Op: DRGEP_sp
operator        = DRGEP_sp 
eps             = 0.02
delta_eps       = 0.01
n_points        = 10.0
max_error_sp    = 30, 30, 30
max_error_T     = 30
max_error_ig    = 30
max_error_Sl    = 30
max_error_K     = 30
inter_sp_inter  = True

#===========> Op: DRG_r
operator        = DRG_r 
eps             = 0.02
delta_eps       = 0.01
n_points        = 10.0
max_error_sp    = 30.0, 30.0, 30.0
max_error_T     = 30.0
max_error_ig    = 30.0
max_error_Sl    = 30.0
max_error_K     = 30.0
inter_sp_inter  = True
optim           = GA
#====> Genetic Algorithm Optimization
n_gen               = 5
n_indiv             = 5
error_fitness       = mean
Arrh_max_variation  = 5, 5, 5
optim_on_meth       = False
sub_mech_sel        = H2, CO, N
selection_operator  = Roulette
selection_options   = 0.2
Xover_operator      = simple_Xover, multiple_Xover, arith_Xover, heuristic_Xover
Xover_pct           = 10, 20, 20, 20
mut_operator        = uniform_mutation, non_uniform_mutation, boundary_mutation
mut_pct             = 30, 30, 10
mut_opt             = , 3, 
mut_intensity       = 20
    """)





    # ================================================================
    # ================================================================
    # ================================================================
    # ================================================================

    fichier_data = open('5_pp.inp', 'w')
    fichier_data.write("""#=============================================
#           Main parameters
#=============================================
main_path         = TEST_5_pp_pts
mech              = C1_GRI30.cti
verbose           = 4
show_plots        = False
write_ck          = True
tspc              = CH4, CO, CO2
T_check           = True
sp_T              = CO2
Sl_check          = False
sp_Sl             = H
ig_check          = False
sp_ig             = CH3
K_check           = False
sp_K              = H
error_calculation = points
error_coupling    = mean


#=============================================
#           Simulation cases
#=============================================

#======> Case 1
config            = pp_flame
Ps                = 100000.0
fuel_1            = CH4
oxidant_1         = O2
diluent_1         = N2
diluent_ratio_1   = N2/O2 3.76
Ts_1              = 300.0
phis_1            = 1.0, 1.5
mdots_1           = 1.0
fuel_2            = CH4
oxidant_2         = O2
diluent_2         = N2
diluent_ratio_2   = N2/O2 3.76
Ts_2              = 300.0
phis_2            = 0.5, 0.5
mdots_2           = 1.0
width             = 0.02
tol_ts            = 1e-05, 1e-08
tol_ss            = 1e-06, 1e-08
transport_model   = Mix
pts_scatter       = [0.0,0.03,0.3,0.5,0.7,1.0]
slope             = 0.05
curve             = 0.05
ratio             = 2.0
prune             = 0.01
restore_flame_folder = ex3_ffl



#=============================================
#           Operators
#=============================================


#===========> Op: DRG_sp
operator        = DRG_sp 
eps             = 0.02
delta_eps       = 0.01
n_points        = 10.0
max_error_sp    = 30, 30, 30
max_error_T     = 30
inter_sp_inter  = True

#===========> Op: SAR_r
operator        = SAR_r 
eps             = 0.02
delta_eps       = 0.01
n_points        = 10.0
max_error_sp    = 30.0, 30.0, 30.0
max_error_T     = 30.0
inter_sp_inter  = True
ttol_sensi      = False, False
optim           = GA
#====> Genetic Algorithm Optimization
n_gen               = 5
n_indiv             = 5
error_fitness       = mean
Arrh_max_variation  = 5, 5, 5
optim_on_meth       = True
nb_r2opt            = 30
sub_mech_sel        = H2, CO, C1, C2, C3, N
selection_operator  = Roulette
selection_options   = 0.2
Xover_operator      = simple_Xover, multiple_Xover, arith_Xover, heuristic_Xover
Xover_pct           = 10, 20, 20, 20
mut_operator        = uniform_mutation, non_uniform_mutation, boundary_mutation
mut_pct             = 30, 30, 10
mut_opt             = , 3, 
mut_intensity       = 20     
    """)





    # ================================================================
    # ================================================================
    # ================================================================
    # ================================================================

    fichier_data = open('6_import_Sl.inp', 'w')
    fichier_data.write("""#=============================================
#           Main parameters
#=============================================
main_path         = TEST_6_import_Sl
mech              = C1_GRI30.cti
ext_results_file  = C2H6_Sl_red.csv
conc_units        = Molar_fraction
verbose           = 4
show_plots        = True
write_ck          = True
tspc              = 
T_check           = False
sp_T              = CO2
Sl_check          = True
sp_Sl             = H
ig_check          = False
sp_ig             = CH3
K_check           = False
sp_K              = H
error_calculation = points
error_coupling    = mean


#=============================================
#           Simulation cases
#=============================================



#=============================================
#           Operators
#=============================================

#===========> Op: Optimization without reduction
operator        = NULL 
optim           = GA
#====> Genetic Algorithm Optimization
n_gen              = 5
n_indiv            = 5
error_fitness      = mean
Arrh_max_variation = 10, 5, 5
optim_on_meth      = DRG 
optim_on_meth_pts  = 20
nb_r2opt            = 30
selection_operator  = Roulette
selection_options   = 0.2
Xover_operator      = simple_Xover, multiple_Xover, arith_Xover, heuristic_Xover
Xover_pct           = 10, 20, 20, 20
mut_operator        = uniform_mutation, non_uniform_mutation, boundary_mutation
mut_pct             = 30, 30, 10
mut_opt             = , 3, 
mut_intensity       = 20    
    """)

