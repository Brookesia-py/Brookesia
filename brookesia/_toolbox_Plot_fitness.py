#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 11:15:39 2019

@author: matynia
"""

def plot_fitness(c_path):

    # Graphs
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    import os
    
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

    os.chdir(c_path)
    try:
        os.mkdir('Plots')
    except:
        a=2            

    worst_list = []
    mean_list  = []
    best_list  = []

    filename = "red_info.txt"

    # lire un fichier 
    fs = open(filename, 'r')

    line = fs.readline()

    while line!='':
        #print(line)
        #  recuperation des donnees
        if 'worst ind' in line:
            worst_list.append(float(line.split('worst ind:')[1].split('best ind:')[0]))
            best_list.append(float(line.split('best ind:')[1].split('mean fitness:')[0]))
            mean_list.append(float(line.split('mean fitness:')[1]))
        line = fs.readline()

        
    fit_evol    = [worst_list,best_list,mean_list]
    list_legend = ['Worst fitness','Best fitness','Mean fitness']
    gen = range(len(fit_evol[0])) 

    os.chdir('Plots')

    # Tracer graph   
    #plt.rc('font', size=14)
    fig,axes = plt.subplots()
    for _curve in range(len(fit_evol)-1):
        axes.plot(gen, fit_evol[_curve+1],\
        linestyle=linestyles[_curve], color=colors_styles[_curve], linewidth=2)
            
        axes.set_xlabel('Iteration number')
        axes.set_ylabel("Fitness")
        
        axes.legend(list_legend[1:])
        fig.show()
        fig.savefig("Fitness.png",dpi=300)

    os.chdir('..')

