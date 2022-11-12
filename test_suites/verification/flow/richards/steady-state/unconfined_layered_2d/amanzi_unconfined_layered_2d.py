import os
import sys
import h5py
import numpy as np
import matplotlib
from matplotlib import pyplot as plt

import run_amanzi_standard
from compare_field_results import GetXY_AmanziU_Nodes
from compare_field_results import GetXY_AmanziU_Values


if __name__ == "__main__":

    # mesh
    nx = 50 
    ny = 1
    nz = 105
    ft = 0.3048

    try:
        input_file = "amanzi_unconfined_layered_2d-u.xml"
        path_to_amanzi = "output-u"
        root_amanzi = 'steady-flow'

        run_amanzi_standard.run_amanzi(input_file, 1, [input_file], path_to_amanzi)

        comp = 'hydraulic_head'
        x_amanziU = GetXY_AmanziU_Nodes(path_to_amanzi, root_amanzi, 0, nx+1, 1, 0)
        c_amanziU = GetXY_AmanziU_Values(path_to_amanzi, root_amanzi, comp, 55, nx*nz, nz)
        unstruct = len(x_amanziU)
    except:
        unstruct = 0

    try:
        import cmath
        x_analytic = x_amanziU
        c_analytic = ft * np.sqrt(28900 - x_analytic * (12/ft))
        # c_analytic = ft * (170 - x_analytic * (0.04/ft))
        analytic = len(x_analytic)
    except:
        analytic = 0


    # subplots
    fig, ax = plt.subplots() 

    if (unstruct>0):
        # reducing amount of data
        x_amanziU = x_amanziU[1::2]
        c_amanziU = c_amanziU[1::2]
        ax.plot(x_amanziU, c_amanziU,marker='o',color='g',label='AmanziU',linestyle='None')
    if (analytic>0):
        ax.plot(x_analytic, c_analytic,'m-',label='Analytic',linewidth=2)

    # axes
    ax.set_xlabel("Distance (m)")
 
    # plot adjustments
    plt.subplots_adjust(left=0.18,bottom=0.12,right=0.98,top=0.9)
    plt.legend(loc='upper right') #, fontsize=13)
    plt.suptitle("Hydraulic Head",x=0.5,fontsize=20)
    plt.tick_params(axis='both', which='major', labelsize=20)

    # plt.show()
    plt.savefig("hydraulic_head.png",format="png")
    # plt.close()

