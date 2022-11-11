import os
import sys
import h5py
import numpy as np
import matplotlib
from matplotlib import pyplot as plt

import run_amanzi_standard
from compare_field_results import GetXY_AmanziU_Nodes
from compare_field_results import GetXY_AmanziU_Values
from compare_field_results import GetXY_AmanziS_1D


if __name__ == "__main__":

    # mesh
    nx = 164
    ny = 120

    try:
        input_file = "amanzi_unconfined_recharge_1d-u.xml"
        path_to_amanzi = "output-u"
        root_amanzi = 'steady-flow'

        run_amanzi_standard.run_amanzi(input_file, 1, [input_file], path_to_amanzi)

        comp = 'hydraulic_head'
        x_amanziU = GetXY_AmanziU_Nodes(path_to_amanzi, root_amanzi, 0, nx+1, 1, 0)
        c_amanziU = GetXY_AmanziU_Values(path_to_amanzi, root_amanzi, comp, 0, nx*ny, ny)
        unstruct = len(x_amanziU)
    except:
        unstruct = 0

    try:
        input_file = "amanzi_unconfined_recharge_1d-s.xml"
        path_to_amanzi = "output-s"
        root_amanzi = "plot"

        run_amanzi_standard.run_amanzi(input_file, 1, [input_file], path_to_amanzi)

        comp = "Hydraulic_Head"
        x_amanziS, c_amanziS = GetXY_AmanziS_1D(path_to_amanzi, root_amanzi, comp, 1)
        # correct the 1D script (random guess)
        c_amanziS = c_amanziS - 20.9
        struct = len(x_amanziS)
    except:
        struct = 0
        
    try:
        import cmath
        x_analytic = x_amanziU
        ft = 0.3048
        c_analytic = (164*ft) * np.sqrt(1 + x_analytic / (1640*ft) * (2 - x_analytic / (1640*ft)))
        analytic = len(x_analytic)
    except:
        analytic = 0


    # subplots
    fig, ax = plt.subplots() 

    if (unstruct>0):
        # reducing amount of data
        x_amanziU = x_amanziU[2::4]
        c_amanziU = c_amanziU[2::4]
        ax.plot(x_amanziU, c_amanziU,marker='o',color='g',label='AmanziU',linestyle='None')
    if (analytic>0):
        ax.plot(x_analytic, c_analytic,'m-',label='Analytic',linewidth=2)
    if (struct>0):
        x_amanziS = x_amanziS[2::4]
        c_amanziS = c_amanziS[2::4]
        ax.plot(x_amanziS, c_amanziS,marker='x',label='AmanziS',linestyle='dotted', linewidth=2)

    # axes
    ax.set_xlabel("Distance (m)")
 
    # plot adjustments
    plt.subplots_adjust(left=0.18,bottom=0.12,right=0.98,top=0.9)
    plt.legend(loc='upper left') #, fontsize=13)
    plt.suptitle("Hydraulic Head",x=0.5,fontsize=20)
    plt.tick_params(axis='both', which='major', labelsize=20)

    # plt.show()
    plt.savefig("hydraulic_head.png",format="png")
    # plt.close()

