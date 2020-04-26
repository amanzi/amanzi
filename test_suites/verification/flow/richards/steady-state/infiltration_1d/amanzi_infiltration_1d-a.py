import os
import sys
import h5py
import numpy as np
import matplotlib
from matplotlib import pyplot as plt

import run_amanzi_standard
from compare_field_results import GetXY_AmanziU_1D
from compare_field_results import GetXY_AmanziS_1D

if __name__ == "__main__":

    # Amanzi: unstructured
    try:
        input_file = "amanzi_infiltration_loam_sand_1d-u.xml"
        path_to_amanzi = "output_2a-u"
        root_amanzi = 'case_2a_plot'
    
        run_amanzi_standard.run_amanzi(input_file, 1, [input_file], path_to_amanzi)

        comp = 'pressure.cell.0'
        x_amanziU, c_amanziU = GetXY_AmanziU_1D(path_to_amanzi,root_amanzi,comp,3)
        unstruct = len(x_amanziU)
    except:
        unstruct = 0


    # Amanzi: structured
    try:
        input_file = "amanzi_infiltration_loam_sand_1d-s.xml"
        path_to_amanzi = "output_2a-s"
        root_amanzi = "plot"

        run_amanzi_standard.run_amanzi(input_file, 1, [input_file], path_to_amanzi)

        comp = "water_Pressure"
        x_amanziS, c_amanziS = GetXY_AmanziS_1D(path_to_amanzi,root_amanzi,comp,2)
        struct = len(x_amanziS)
    except:
        struct = 0
        

    # Amanzi: analytic
    try:
        path_to_amanzi = "golden_output"
        root_amanzi = 'case_2a_plot'

        comp = 'pressure.cell.0'
        x_amanziU_gold, c_amanziU_gold = GetXY_AmanziU_1D(path_to_amanzi,root_amanzi,comp,3)
        unstruct_gold = len(x_amanziU_gold)
    except:
        unstruct_gold = 0


    # subplots
    fig, ax = plt.subplots() 

    if (unstruct>0):
        ax.plot(x_amanziU, c_amanziU,marker='o',color='g',label='AmanziU',linestyle='None')
    if (unstruct_gold>0):
        ax.plot(x_amanziU_gold, c_amanziU_gold,'m-',label='AmanziU-Gold',linewidth=2)
    if (struct>0):
        ax.plot(x_amanziS, c_amanziS,marker='x',label='AmanziS',linestyle='dotted', linewidth=2)

    # axes
    ax.set_xlabel("Distance (m)") #, fontsize=20)
 
    # plot adjustments
    plt.subplots_adjust(left=0.18,bottom=0.12,right=0.98,top=0.9)
    plt.legend(loc='upper left') #, fontsize=13)
    plt.suptitle("Aqueous Pressure",x=0.5,fontsize=20)
    plt.tick_params(axis='both', which='major', labelsize=20)

    # plt.show()
    plt.savefig("infiltration_loam_sand_1d.png",format="png")
    # plt.close()
