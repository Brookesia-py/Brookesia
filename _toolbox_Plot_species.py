'''
se placer à dans le dossier contenant le fichier résultat (typiquement 0_reduction.csv) puis taper :
ipython
import tool_brooks as tb
tb.plot_data('T, CH4, CO, CO2', '0_reduction_results.csv', True)
avec les options :
'T, CH4, CO, CO2'     : nom des espèces à tracer
'0_reduction_results' : nom du fichier résultat
True                  : afficher la légende sur le graph (True : oui / False : non)
'''

#def plot_data_v7(c_path, Var_list = 'T, CH4, CO, CO2, C2H6, C2H4', fileName = ['all.csv'], display_legend = False):
def plot_data_v7(c_path, Var_list = 'CH4, CO, CO2', fileName = '0_reduction_results.csv', display_legend = False):
#def plot_data_v7(c_path, Var_list = 'CO$_2$', fileName = ['all.csv'], display_legend = False):
#def plot_data_v7(c_path, Var_list = 'CO', fileName = ['all.csv'], display_legend = True):


    #%%============================================================================
    ###    Plot options
    #==============================================================================


    # Zoom
    zoom = "n"         # y / n
    x_min = 0.15
    x_max = .27

    Plot_Sl  = False
    Plot_Igt = True #(! ne fonctionne pas correctement si plusieurs pressions)

    #fileName = ["0_reduction_results.csv"] # if fileName not in input options

    plot_raw    = True  # scat for exp
    plot_smooth = False


    txt_size = 12

    fig_width = 4.5
    left_shift = 0.3
    
    plot_style = 'paper' #      'paper' 'paper_latex'  'presentation'  'poster' 


    
    
    
    # pour adapter les polices (si besoin)
    delta_police_legend = -4.5
    delta_police_axes   = -2
    delta_police_ticks  = -2
    delta_police_title  = -2

    

    




    #%%============================================================================
    ###    Library importation
    #==============================================================================
    #System
    import os
    import sys
    import copy

    #Csv files
    import csv

    # Math
    import numpy as np
    from scipy import interpolate

    # Graphs
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    from matplotlib.ticker import (MultipleLocator,FormatStrFormatter,AutoMinorLocator)

    #if 'brookesia.Class_def' not in sys.modules:  
        #del sys.modules["brookesia.Class_def"]
    import brookesia.Class_def as cdef
        
    
    # 3D tools
    #from mpl_toolkits.mplot3d import Axes3D
    #from matplotlib.collections import PolyCollection
    #from matplotlib.colors import colorConverter
    #import matplotlib.pyplot as plt
    #import matplotlib.cm as cm      #color management

    import copy

    os.chdir(c_path)
    if '/' in fileName:
        os.chdir(fileName.split('/')[0])
        fileName = fileName.split('/')[-1]
    try:
        os.mkdir('Plots')
    except: 
        a=2
        
    #print(c_path)

    #%%============================================================================
    ###    Input options (1)
    #==============================================================================
    #sys.argv.append("0_reduction_results.csv")
    #sys.argv.append("references.csv++1")
    #sys.argv.append("CO2")
    #sys.argv.append("CO")
    #sys.argv.append("NO")
    #sys.argv.append("CH")

    Var_name,Var_name_legend = [],[]
    Var_list = Var_list.split(',')
    for var in Var_list:
        var_a = var.replace(' ','')
        var_b = var_a.replace('_','')
        var_c = var_b.replace('$','')
        Var_name_legend.append(var_a)
        Var_name.append(var_c)
    
    delta_polices = [delta_police_legend,delta_police_axes,delta_police_ticks,delta_police_title]


    print(str(sys.argv))

    if len(sys.argv)>1:
        fileName = []
        #if sys.argv[1] == "mult":
            #fileName = []
            #for n in range(len(sys.argv)):
                #if ".csv" in sys.argv[n]:
                    #fileName.append(sys.argv[n])
        #else:
        for n in range(len(sys.argv)):
            if ".csv" in sys.argv[n]:
                fileName.append(sys.argv[n])



    #%%============================================================================
    ###    Get data
    #==============================================================================

    # flags
    case_nb = -1
    step_nb = -1

    # initiating lists
    data               = []
    data_save          = []
    steps_save         = []
    headers_steps      = []
    headers_steps_save = []
    case_titles        = []
    case_data          = []
    case_steps         = []
    Sl                 = []
    Sl_phi             = []
    Igt                = []
    Igt_phi            = []
    Igt_T              = []

    if type(fileName) is str:
        fileName = [fileName]     # obsolete
    for file_n in fileName:
        file_4_ref_only = False
        if '++' in file_n:
            file_n, _ind_st = file_n.split('++')
            print(_ind_st)
            _ind_st = int(_ind_st)
    #        file_n = file_n.replace('++','')
            case_nb=-1
            file_4_ref_only = True
        with open(file_n) as csvfile:
            readcsv = csv.reader(csvfile, delimiter=';')

            read_data=0
            for line in readcsv:
                if line!=[]:
                    if "Case" in line[0]:
                        # if applicable, save previous config data
                        if (case_nb>-1 and not file_4_ref_only) or (case_nb==-1 and file_4_ref_only):
                            data_save.append(data) ; data = [] ;
                            steps_save.append(step_title) ;
                            case_data.append(data_save); data_save = []
                            case_steps.append(steps_save); steps_save = []
                            headers_steps_save.append(headers_steps); headers_steps=[]

                        if file_4_ref_only and case_nb>-1:
                            data_save.append(data)
                            case_data[case_nb].insert(_ind_st,data); data_save = [] ; data = []
                            steps_save.append(step_title)
                            case_steps[case_nb].insert(_ind_st,step_title); steps_save = []
                            headers_steps_save[case_nb].insert(_ind_st,headers_steps[0]); headers_steps=[]

                        case_nb += 1
                        step_nb = -1
                        read_data=0   # if applicable, stop recording data in data list

                    if "reactor" in line[0] or "JSR" in line[0]:
                        c_title = str(case_nb+1)+'_'+line[0]+'_'+line[1]+'_'
                        if '[O2/N2' in line[5] and '3.7' in line[5]:
                            c_title +='air_phi=_'+'%.2f' %float(line[4])+'_Ti='+str(round(float(line[8])))+'K'
                        else :
                            c_title +='O2_dil_phi=_'+'%.2f' %float(line[4])+'_Ti='+str(round(float(line[8])))+'K'
                        case_titles.append(c_title)
                        phi = float(line[4])
                        T_i = float(line[8])
                    elif "flame" in line[0] or "PFR" in line[0]:
                        c_title = str(case_nb+1)+'_'+line[0]+'_'+line[1]+'_'
                        if '[O2/N2' in line[5] and '3.7' in line[5]:
                            c_title +='air_phi=_'+'%.2f' %float(line[4])+'_Ti='+str(round(float(line[8])))+'K'
                        else :
                            c_title +='O2_dil_phi=_'+'%.2f' %float(line[4])+'_Ti='+str(round(float(line[8])))+'K'
                        case_titles.append(c_title)
                        if 'diff' not in line[0]:
                            phi = float(line[4])

    #                try:    Sl_new = float(float(line[1])) ; check_Sl=True
    #                except: check_Sl=False
                    if 'Sl0(cm/s)' in line[0] and Plot_Sl:
                        if step_nb==0:        Sl_phi.append(phi)
                        if len(Sl)==step_nb:  Sl.append([float(line[1])])
                        else:                 Sl[step_nb].append(float(line[1]))

                    if 'Ignition delay time(s)' in line[0] and Plot_Igt:
                        if len(Igt)<step_nb+1: Igt.append([])
                        if phi not in Igt_phi:  Igt_phi.append(phi)
                        idx_phi = Igt_phi.index(phi)
                        if len(Igt_T)<idx_phi+1:
                            Igt_T.append([])
                        if len(Igt[step_nb])<idx_phi+1:
                            Igt[step_nb].append([])
                        if T_i not in Igt_T[idx_phi]: Igt_T[idx_phi].append(T_i)
                        Igt[step_nb][idx_phi].append(float(line[1]))


    #                        if len(Igt)==step_nb+1: Igt.append([float(line[-1])])
    #                        else:
    #                            try:
    #                                Igt[step_nb+1].append(float(line[-1]))



                    if "Step" in line[0]:
                        # if applicable, save previous step data
                        if step_nb>-1:
                            data_save.append(data) ; data = []
                            steps_save.append(step_title)

                        step_title = line[0].split(": ")[1]
                        step_nb +=1
                        read_data=0      # if applicable, stop recording data in data list

                    if line[0]=='':
                        read_data=0

                    if "T(K)" in line or "Tf(K)" in line:   # means the headers line is reached
                        headers_steps.append(line)
                        read_data=1

                    elif read_data == 1:
                        data.append(line) # Construct the data list





    # Save previous config data
    if not file_4_ref_only:
        data_save.append(data)
        case_data.append(data_save)
        steps_save.append(step_title)
        case_steps.append(steps_save)
        headers_steps_save.append(headers_steps)
    else :
        case_data[case_nb].insert(1,data)
        case_steps[case_nb].insert(1,step_title)
        headers_steps_save[case_nb].insert(1,headers_steps[0])


    #%%============================================================================
    ###    Input options (2)
    #==============================================================================


    cases_selected = "all"
    steps_selected = "all"
    DRG1=False ; SA1=False ; LOI1=False

    if len(sys.argv)>1:
        # Multi datafile
        if sys.argv[1] == "mult":
            i=0
            # display cases
            for case_title in case_titles:
                print(str(i) + ":  " + case_title)
                i+=1
            # get the conditions to plot
            case_choices = input("Select the conditions you want to compare (for example: 1 3 5): ")
            cases_selected = [int(case_choices) for case_choices in case_choices.split(" ")]
            steps_choices = input("Select the conditions you want to compare (for example: DRG SA DRG+GA ): ")
            steps_selected = steps_choices.split(" ")
            # If DRG / SA / LOI is to be plotted only for the first condition
            if "DRG1" in steps_selected:
                position_arg = [i+1 for i in range(len(steps_selected)) if steps_selected[i]=="DRG1"]  #find the position of the condition "DRG1"
                steps_selected[position_arg[0]-1]="DRG"
                DRG1=True
            elif "SA1" in steps_selected:
                position_arg = [i+1 for i in range(len(steps_selected)) if steps_selected[i]=="SA1"]  #find the position of the condition "DRG1"
                steps_selected[position_arg[0]-1]="SA"
                SA1=True
            elif "LOI1" in steps_selected:
                position_arg = [i+1 for i in range(len(steps_selected)) if steps_selected[i] =="LOI1"]  #find the position of the condition "DRG1"
                steps_selected[position_arg[0]-1]="LOI"
                LOI1=True



            add_legend = []
            for n_cases in range(len(cases_selected)):
                add_legend.append(input("Additional comments for legend of the condition " + str(cases_selected[n_cases]) + ": "))

            modif_title = input("Do you want to modify the title of the figure (default: last case) (y/n): ")
            if modif_title=="y" or modif_title=="yes":
                new_fig_title = input("New title: ")
    #        while modif_titre == "y":
    #            nb_case = input("Select the number of the case to modify: ")
    #            case_titles[int(nb_case)] = input("New title: ")
    #            modif_titre = input("Do you want to modify the title of a case (y/n): ")

            # get the variables to plot
            Var_name_original = "y"
            for n in range(len(sys.argv)-1):
                n+=1
                if ".csv" not in sys.argv[n] and "mult" not in sys.argv[n] and "modif_title" not in sys.argv[n] and "++" not in  sys.argv[n]:
                    if Var_name_original == "y":
                        Var_name = [] ; Var_name_original = "n"
                    Var_name.append(sys.argv[n])

        # One datafile - get the variables to plot
        elif len(sys.argv)>2 and "noT" in sys.argv:
            Var_name = Var_name[1:]
        elif len(sys.argv)>2:
            Var_name = []
            for n in range(len(sys.argv)-2):
                if ".csv" not in sys.argv[n+2] and "mult" not in sys.argv[n+2] and "modif_title" not in sys.argv[n+2] and "++" not in  sys.argv[n+2]:
                    Var_name.append(sys.argv[n+2])



        # Modif title
        if sys.argv[1] == "modif_title":
            i=1
            # display cases
            for case_title in case_titles:
                print(str(i) + ":  " + case_title)
                i+=1
            modif_title = "y"
            while modif_title == "y" or modif_title == "yes":
                nb_case = input("Select the number of the case to modify: ")
                case_titles[int(nb_case)] = input("New title: ")
                modif_title = input("Do you want to modify the title of a case (y/n): ")





    #%%============================================================================
    ###    Plot graphs
    #==============================================================================
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

    
    os.chdir('Plots')
    
    plt.close("all")

    case_legend=copy.deepcopy(case_steps)
    case_n=False
    if not plot_smooth:
        for case_s in case_legend:
            for step_n in case_s:
                if 'smooth' in step_n:
                    case_n = step_n
        for case_s in case_legend:
            if case_n:  case_s.remove(case_n)
    
    if 'latex' in 'plot_style':
        for cn in range(len(case_legend)):
            for sn in range(len(case_legend[cn])):
                case_legend[cn][sn] = case_legend[cn][sn].replace('_','\_')

    n_case=-1 ; legend_mult=[]
    first_ref_done=0; first_DRG_done=0;first_SA_done=0;first_LOI_done=0;

    if Var_name!=[]:
        for case in range(len(case_data)): #plot every simulation conditions present in the file
            y_min=0; y_max=np.zeros((len(Var_name)))

            if cases_selected=="all":
                fig = plt.figure(case+1)
                fig.set_size_inches(fig_width, len(Var_name)*3)
                fig.subplots_adjust(left = left_shift, bottom = 0.05,
                        right = 0.9, top = 0.95, wspace = 0, hspace = 0.5)
                plt.rc('font', size=txt_size)
            elif case in cases_selected :
                n_case+=1
        #        print(str(n_case))
                fig = plt.figure(1)
                fig.set_size_inches(fig_width, len(Var_name)*3)
                fig.subplots_adjust(left = left_shift, bottom = 0.05,
                            right = 0.9, top = 0.95, wspace = 0, hspace = 0.5)
                plt.rc('font', size=txt_size)
            else:
                continue
            if "mult" not in sys.argv:
                fig_title = case_titles[case].replace("/","")+".png"
                plt,fig = cdef.generate_plot_style(plt,fig, plot_style,delta_polices)


            for steps in range(len(case_steps[case])): #plot every simulation conditions present in the file
                if case_steps[case][steps] in steps_selected or "Reference" in case_steps[case][steps] or steps_selected=="all" :
                    if "Reference" in case_steps[case][steps] and first_ref_done == 1 :
                        continue
                    elif "DRG" in case_steps[case][steps] and first_DRG_done == 1 :
                        continue
                    elif "SA" in case_steps[case][steps] and first_SA_done == 1 :
                        continue
                    elif "LOI" in case_steps[case][steps] and first_LOI_done == 1 :
                        continue
                    else:
                        # ---------------------
                        # Index

                        # ---------------------
                        # Get data
                        data = case_data[case][steps]

                        # Abscissa
                        if "reactor" in case_titles[case]:
                            Idx_abs = headers_steps_save[case][steps].index("Time(s)")
                        elif "JSR" in case_titles[case]:
                            Idx_abs = headers_steps_save[case][steps].index("Ti(K)")
                        elif "flame" in case_titles[case] or "PFR" in case_titles[case]:
                            Idx_abs = headers_steps_save[case][steps].index("Z(m)")
                        abscissa = []
                        for line in range(len(data)):
                            abscissa.append(float(data[line][Idx_abs]))
                        if max(abscissa)<0.1:
                            abscissa = np.array(abscissa)*1000 #abscissa_ms=[i*1000 for i in abscissa]

                        # Ordinate
                        Idx_var = []
                        for var in Var_name:
                            try:
                                if var != 'T' and var != 'T(K)':
                                    Idx_var.append(headers_steps_save[case][steps].index(var))
                                else:
                                    Idx_var.append(1)
                            except:
                                print('Error with var: ' + var)
                                break
                        data2plot = []
                        data2plot_spec = []
                        for var in range(len(Var_name)):
                            for line in range(len(data)):
                                data2plot_spec.append(float(data[line][Idx_var[var]]))

                            data2plot.append(data2plot_spec); data2plot_spec=[]


                        # preparation of the legend for mult cases
                        if "mult" in sys.argv:
                            if first_ref_done==0:
                                legend_mult.append(case_legend[case][steps])
                            elif case_legend[case][steps]=="DRG" and DRG1==True:
                                legend_mult.append(case_legend[case][steps])
                            elif case_legend[case][steps]=="SA" and SA1==True:
                                legend_mult.append(case_legend[case][steps])
                            elif case_legend[case][steps]=="LOI" and LOI1==True:
                                legend_mult.append(case_legend[case][steps])
                            else:
                                legend_mult.append(add_legend[n_case] + case_legend[case][steps])

                        #   Subplot
                #        plt.title(Var_name[var])
                        for var in range(len(Var_name)):
                            axes = plt.subplot(len(Var_name), 1, var+1)

                            # Set the abscissa axis and legend
                            if "reactor" in case_titles[case]:
                                if ("(raw)" in case_steps[case][steps] or "Experimental data" in case_steps[case][steps]) and plot_raw:
                                    axes.scatter(abscissa, data2plot[var], \
                                    marker=scatter_styles[steps], color="black")
                                elif ("(raw)" in case_steps[case][steps] or "Experimental data" in case_steps[case][steps]) and not plot_raw:
                                    a=2
                                elif "smooth" in case_steps[case][steps]:
                                    a=2
                                elif ("(raw)" in str(case_steps[case]) or "Experimental data" in str(case_steps[case])):
                                    axes.plot(abscissa, data2plot[var], \
                                    linestyle=linestyles[steps-1], color=colors_styles[steps], linewidth=2)
                                else:
                                    axes.plot(abscissa, data2plot[var], \
                                    linestyle=linestyles[steps], color=colors_styles[steps], linewidth=2)
                                if max(abscissa)<0.1:
                                    axes.set_xlabel("Time (ms)")
                                else: axes.set_xlabel("Time (s)")

                            if "PFR" in case_titles[case]:
                                if ("(raw)" in case_steps[case][steps] or "Experimental data" in case_steps[case][steps]) and plot_raw:
                                    axes.scatter(abscissa, data2plot[var], \
                                    marker=scatter_styles[steps], color="black")
                                elif ("(raw)" in case_steps[case][steps] or "Experimental data" in case_steps[case][steps]) and not plot_raw:
                                    a=2
                                elif "smooth" in case_steps[case][steps]:
                                    a=2
                                elif ("(raw)" in str(case_steps[case]) or "Experimental data" in str(case_steps[case])):
                                    axes.plot(abscissa, data2plot[var], \
                                    linestyle=linestyles[steps-1], color=colors_styles[steps], linewidth=2)
                                else:
                                    axes.plot(abscissa, data2plot[var], \
                                    linestyle=linestyles[steps], color=colors_styles[steps], linewidth=2)
                                if max(abscissa)<0.1:
                                    axes.set_xlabel("Z (m)")
                                else: axes.set_xlabel("Z (m)")

                            elif "JSR" in case_titles[case]:
                                if ("(raw)" in case_steps[case][steps] or "Experimental data" in case_steps[case][steps]) and plot_raw:
                                    axes.scatter(abscissa, data2plot[var], \
                                    marker=scatter_styles[steps], color="black")
                                elif ("(raw)" in case_steps[case][steps] or "Experimental data" in case_steps[case][steps]) and not plot_raw:
                                    a=2
                                elif "smooth" in case_steps[case][steps]:
                                    a=2
                                elif ("(raw)" in str(case_steps[case]) or "Experimental data" in str(case_steps[case])):
    #                                x=abscissa
    #                                ynew=data2plot[var]
                                    x = np.linspace(abscissa[0], abscissa[-1], num=200, endpoint=True)
                                    f = interpolate.interp1d(abscissa, data2plot[var], kind='quadratic')  # linear slinear quadratic cubic
                                    ynew=f(x)
    #                                f = interpolate.splrep(abscissa, data2plot[var], s=.5)
    #                                ynew = interpolate.splev(x, f, der=0)
                                    axes.plot(x, ynew, \
                                    linestyle=linestyles[steps-1], color=colors_styles[steps], linewidth=2)
                                else:
                                    x=abscissa
                                    ynew=data2plot[var]
    #                                x = np.linspace(abscissa[0], abscissa[-1], num=30, endpoint=True)
    #                                f = interpolate.interp1d(abscissa, data2plot[var], kind='cubic')
    #                                ynew=f(x)
    #                                f = interpolate.splrep(abscissa, data2plot[var], s=.5)
    #                                ynew = interpolate.splev(x, f, der=0)
                                    axes.plot(x, ynew, \
                                    linestyle=linestyles[steps], color=colors_styles[steps], linewidth=2)
                                axes.set_xlabel("T (K)")

                            elif "flame" in case_titles[case]:
                                if ("(raw)" in case_steps[case][steps] or "Experimental data" in case_steps[case][steps]) and plot_raw:
                                    axes.scatter(abscissa, data2plot[var], \
                                    marker=scatter_styles[steps], color="black")
                                elif ("(raw)" in case_steps[case][steps] or "Experimental data" in case_steps[case][steps]) and not plot_raw:
                                    a=2
                                elif "smooth" in case_steps[case][steps]:
                                    a=2
                                elif ("(raw)" in str(case_steps[case]) or "Experimental data" in str(case_steps[case])):
                                    axes.plot(abscissa, data2plot[var], \
                                    linestyle=linestyles[steps-1], color=colors_styles[steps], linewidth=2)
                                else:
                                    axes.plot(abscissa, data2plot[var], \
                                    linestyle=linestyles[steps], color=colors_styles[steps], linewidth=2)
                                axes.set_xlabel("Z (mm)")

                            if "mult" not in sys.argv and display_legend:
                                if ("(raw)" in str(case_steps[case]) or "Experimental data" in str(case_steps[case])):
                                    axes.legend(case_legend[case][-len(case_legend[case])+1:]+[case_legend[case][0]])
                                else:
                                    axes.legend(case_legend[case])

                            #if "Time(s)" in headers_steps_save[case][steps]:
                                #if max(abscissa)<0.1:
                                    #axes.plot(abscissa_ms, data2plot[var])
                                    #if "mult" not in sys.argv:
                                        #axes.legend(case_steps[case])
                                    #axes.set_xlabel("Time(ms)")
                                #else: axes.set_xlabel("Time(s)")
                            #else:
                                #axes.plot(abscissa, data2plot[var])
                                #axes.legend(case_steps[case])
                                #axes.set_xlabel("Z(m)")

                            if zoom == "y":
                                axes.set_xlim(x_min, x_max)
                            # axes.ticklabel_format(axis='y', style='sci', scilimits=(-2,+4))

                            # set the ordinate axis
                            if "T(K)" in Var_name_legend[var] or "T" in Var_name_legend[var]:
                                axes.set_ylabel(Var_name_legend[var])
                            else:
                                axes.set_ylabel("X("+Var_name_legend[var]+")")
                            # y ticks
        #                    start, end = axes.get_ylim()
                            y_max[var]=np.max([y_max[var],max(data2plot[var])*1.1])
                            y_max[var] = float('%.2g' % y_max[var])
                            if "T(K)" in Var_name_legend[var] or "T" in Var_name_legend[var] :
                                axes.yaxis.set_major_formatter(ticker.FormatStrFormatter('%2.f'))
                            else:
                                if y_min!=y_max[var]:
                                    axes.yaxis.set_ticks(np.arange(y_min, y_max[var], (y_max[var]-0)/5))
                                    axes.set_ylim(y_min, y_max[var])
                                axes.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1e'))
                            axes.xaxis.set_minor_locator(AutoMinorLocator())
                            axes.yaxis.set_minor_locator(AutoMinorLocator())


                        # Indicate that the first reference is plotted
                        if "Reference" in case_steps[case][steps] and "mult" in sys.argv:
                            first_ref_done = 1
                        elif case_steps[case][steps]=="DRG" and DRG1==True:
                            first_DRG_done = 1
                        elif case_steps[case][steps]=="SA"  and SA1==True:
                            first_SA_done = 1
                        elif case_steps[case][steps]=="LOI" and LOI1==True:
                            first_LOI_done = 1

                        #last case to use for figure title (mult)
                        last_case_nb = case


            if "mult" not in sys.argv:
                #fig_title = case_titles[case].replace("/","")+".png"
                #plt,fig = cdef.generate_plot_style(plt,fig, plot_style,delta_polices)
                                
                #plt.show()
                plt.savefig(fig_title,dpi=300)



    if Plot_Sl:
        fd = open('Sl_table.csv', 'w')
        # headers
        fd.write('Phi;')
        for steps in range(len(Sl)):
            fd.write(str(case_steps[0][steps])+';')
        fd.write('\n')
        for phi in range(len(Sl_phi)):
            fd.write(str(Sl_phi[phi])+';')
            for steps in range(len(Sl)):
                fd.write(str(Sl[steps][phi])+';')
            fd.write('\n')

    #    plt.rc('font', family='serif') txt_size
        plt.rc('font', size=txt_size)
        fig,axes = plt.subplots()
        for steps in range(len(Sl)):
            if "(raw)" in case_steps[0][steps]:
                axes.scatter(Sl_phi, Sl[steps], \
                marker=scatter_styles[steps], color=colors_styles[steps])
            elif "(smooth" in str(case_steps[0][steps]):
                a=2
            elif "(raw)" in str(case_steps[0]):
                axes.plot(Sl_phi, Sl[steps], \
                linestyle=linestyles[steps-1], color=colors_styles[steps-1], linewidth=2)
            else:
                axes.plot(Sl_phi, Sl[steps], \
                linestyle=linestyles[steps], color=colors_styles[steps], linewidth=2)
            axes.set_xlabel(r'$\Phi$',fontsize=20)
            axes.set_ylabel("Sl (cm/s)",fontsize=14)
            if "(raw)" in str(case_steps[0]) and display_legend:
                axes.legend(case_steps[0][-len(case_steps[0])+2:]+[case_steps[0][0]])
            else:
                axes.legend(case_steps[0])
        axes.xaxis.set_minor_locator(AutoMinorLocator())
        axes.yaxis.set_minor_locator(AutoMinorLocator())
        plt,fig = cdef.generate_plot_style(plt,fig, plot_style,delta_polices)
        #fig.show()
        fig.savefig("Sl.png",dpi=300)

    if Plot_Igt: #(! ne fonctionne pas correctement si plusieurs pressions)
        for _phi in range(len(Igt_phi)):
            fd = open('igt(sec)_table_phi='+str('%.2f' %Igt_phi[_phi])+'.csv', 'w')
            fd.write('T(K);')
            for steps in range(len(Igt)):
                fd.write(str(case_steps[0][steps])+';')
            fd.write('\n')
            for _T in range(len(Igt_T[_phi])):
                fd.write(str(Igt_T[_phi][_T])+';')
                for steps in range(len(Igt)):
                    fd.write(str(Igt[0][_phi][_T])+';')
                fd.write('\n')

    if len(sys.argv)>1:
        # Multi datafile
        if "mult" in sys.argv and display_legend:
            # Add legend
            for var in range(len(Var_name_legend)):
                axes = plt.subplot(len(Var_name_legend), 1, var+1)
                axes.legend(legend_mult)
            # Save figure
            
    plt,fig = cdef.generate_plot_style(plt,fig, plot_style,delta_polices)
    fig.savefig(case_titles[last_case_nb] +".png",dpi=300)
    
    #if 'brookesia.Class_def' in sys.modules:  
        #del sys.modules["brookesia.Class_def"]
    
    os.chdir('..')



