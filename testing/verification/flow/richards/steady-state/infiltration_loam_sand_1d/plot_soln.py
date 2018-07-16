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

    try:
        path_to_amanziS = "."
        root_amanziS = "case_2a_plot00001"
        compS = "Aqueous_Pressure"
        x_amanziS, c_amanziS = GetXY_AmanziS_1D(path_to_amanziS,root_amanziS,compS)
        struct = len(x_amanziS)
    except:
        struct = 0
        
    # subplots
    fig, ax = plt.subplots() 
        
    try:
        time = '428'
        comp = 'pressure.cell.0'
        path_to_amanziU = "."
        root_amanziU = 'case_2a_plot'
        x_amanziU, c_amanziU = GetXY_AmanziU_1D(path_to_amanziU,root_amanziU,comp,3)
        unstruct = len(x_amanziU)
     
    except:

        unstruct = 0

    # Do plot
    if (unstruct>0):
        alq = ax.plot(x_amanziU, c_amanziU,'m-',label='AmanziU',linewidth=2)
    if (struct>0):
        sam = ax.plot(x_amanziS, c_amanziS,'g-',label='AmanziS',linewidth=2)

    # Diff
    #sad = ax.plot(x_amanziS, c_amanziS-c_amanziU,'g-',label='AmanziS',linewidth=2) 

    # axes
    ax.set_xlabel("Distance (m)",fontsize=20)
 
    # plot adjustments
    plt.subplots_adjust(left=0.13,bottom=0.1,right=0.91,top=0.9)
    plt.legend(loc='upper left',fontsize=13)
    plt.suptitle("Aqueous Pressure",x=0.5,fontsize=20)
    plt.tick_params(axis='both', which='major', labelsize=20)

    plt.show()
    # plt.savefig(root+"_1d.png",format="png")
    # plt.close()
