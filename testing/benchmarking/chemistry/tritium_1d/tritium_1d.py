# plots calcium concentration along x at last time step 
# benchmark: compares to pflotran simulation results
# author: S.Molins - Sept. 2013

import os
import sys
import h5py
import numpy as np
import matplotlib
#matplotlib.use('Agg')
from matplotlib import pyplot as plt

import run_amanzi_standard
import numpy as np
from compare_field_results import GetXY_AmanziU
from compare_field_results import GetXY_AmanziS
from compare_field_results import GetXY_PFloTran
from compare_field_results import GetXY_CrunchFlow


if __name__ == "__main__":

    try:
        sys.path.append('../../../../tools/amanzi_xml')
    except:
        pass

    try:
        sys.path.append('../../../../MY_TPL_BUILD/ccse/ccse-1.3.4-source/Tools/Py_util')
    except:
        pass
    
    # root name for problem
    root = "tritium"

    # pflotran
    path_to_pflotran = "pflotran"
    root_pflotran = "1d-"+root

    # hardwired for 1d-tritium: time and comp
    time = 'Time:  5.00000E+01 y'
    comp = 'Total_'+root.title()+' [M]'

    x_pflotran, c_pflotran = GetXY_PFloTran(path_to_pflotran,root_pflotran,time,comp)    
    
    # crunchflow GIMRT
    path_to_crunchflow = "crunchflow"

    # hardwired for 1d-tritium-crunch.in: time and comp
    times_CF = 'totcon5.out'
    comp = 0
    ignore = 4

    x_crunchflow, c_crunchflow = GetXY_CrunchFlow(path_to_crunchflow,root,times_CF,comp,ignore)

    CWD = os.getcwd()
    local_path = ""

    # subplots
    fig, ax = plt.subplots() 
        
    # AmanziS + Alqumia + PFlotran chemistry
    try:
        input_file = os.path.join("amanzi-s-1d-"+root+"-alq-pflo.xml")
        path_to_amanzi = "output-s-alq-pflo"
        run_amanzi_standard.run_amanzi(input_file, 1,
                                       ["1d-"+root+"-trim.in",root+".dat",input_file],
                                       path_to_amanzi)
        root_amanzi = "plt00037"
        comp = "Tritium_water_Concentration"
        x_amanziS, c_amanziS = GetXY_AmanziS(path_to_amanzi,root_amanzi,comp)
        struct = len(x_amanziS)
    except:
        struct = 0


    # AmanziS + Alqumia + CrunchFlow chemistry
    try:
        input_file = os.path.join("amanzi-s-1d-"+root+"-alq-crunch.xml")
        path_to_amanzi = "output-s-alq-crunch"
        run_amanzi_standard.run_amanzi(input_file, 1, 
                                       ["1d-"+root+"-crunch.in",root+".dbs","aqueous.dbs",input_file],
                                       path_to_amanzi)
        root_amanzi = "plt00037"
        comp = "Tritium_water_Concentration"
        x_amanziS_crunch, c_amanziS_crunch = GetXY_AmanziS(path_to_amanzi,root_amanzi,comp)
        struct_crunch = len(x_amanziS_crunch)
    except:
        struct_crunch = 0


    # AmanziU + Native chemistry
    try:
        comp = 'total_component_concentration.cell.Tritium conc'
        input_file = os.path.join("amanzi-u-1d-"+root+".xml")
        path_to_amanzi = "output-u"
        run_amanzi_standard.run_amanzi(input_file, 1, [root+".bgd",input_file], path_to_amanzi)

        x_amanziU_native, c_amanziU_native = GetXY_AmanziU(path_to_amanzi,root,comp)
        native = len(x_amanziU_native)
    except:
        native = 0


    # AmanziU + Alqumia + PFloTran chemistry
    try:
        comp = 'total_component_concentration.cell.Tritium conc'
        input_file = os.path.join("amanzi-u-1d-"+root+"-alq-pflo.xml")
        path_to_amanzi = "output-u-alq-pflo"
        run_amanzi_standard.run_amanzi(input_file, 1,
                                       ["1d-"+root+"-trim.in",root+".dat",input_file],
                                       path_to_amanzi)
        x_amanziU, c_amanziU = GetXY_AmanziU(path_to_amanzi,root,comp)
        unstruct = len(x_amanziU)
    except:
        unstruct = 0


    # AmanziU + Alqumia + CrunchFlow chemistry
    try:
        comp = 'total_component_concentration.cell.Tritium conc'
        input_file = os.path.join("amanzi-u-1d-"+root+"-alq-crunch.xml")
        path_to_amanzi = "output-u-alq-crunch"
        run_amanzi_standard.run_amanzi(input_file, 1, 
                                       ["1d-"+root+"-crunch.in",root+".dbs","aqueous"+".dbs",input_file],
                                       path_to_amanzi)
        x_amanziU_crunch, c_amanziU_crunch = GetXY_AmanziU(path_to_amanzi,root,comp)
        unstruct_crunch = len(x_amanziU_crunch)

    except:
        unstruct_crunch = False

    # Do plot
    if (native > 0):
        nat = ax.plot(x_amanziU_native, c_amanziU_native,'rx',label='AmanziU(2nd-Order)+Native',linewidth=2)

    if (unstruct > 0):
        alq = ax.plot(x_amanziU, c_amanziU,'r-',label='AmanziU(2nd-Order)+Alquimia(PFloTran)',linewidth=2)

    if (unstruct_crunch > 0):
        alq_crunch = ax.plot(x_amanziU_crunch, c_amanziU_crunch,'r*',label='AmanziU(2nd-Order)+Alquimia(CrunchFlow)',linewidth=2)

    if (struct > 0):
        sam = ax.plot(x_amanziS, c_amanziS,'g-',label='AmanziS+Alquimia(PFloTran)',linewidth=2) 

    if (struct_crunch > 0):
        sam_crunch = ax.plot(x_amanziS_crunch, c_amanziS_crunch,'g*',label='AmanziS+Alquimia(CrunchFlow)',linewidth=2) 

    pfl = ax.plot(x_pflotran, c_pflotran,'m-',label='PFloTran',linewidth=2)
    crunch = ax.plot(x_crunchflow, c_crunchflow,'m*',label='CrunchFlow(OS3D)',linewidth=2)

    # axes
    ax.set_xlabel("Distance (m)",fontsize=20)
    ax.set_ylabel("Total "+root.title()+" concentration [mol/L]",fontsize=20)

    # ax.set_xlim(43,57)
    # ax.set_ylim(0,.0001)

    # plot adjustments
    plt.subplots_adjust(left=0.20,bottom=0.15,right=0.99,top=0.90)
    plt.legend(loc='upper right',fontsize=10)
    plt.suptitle("Amanzi 1D "+root.title()+" Benchmark at 50 years",x=0.57,fontsize=20)
    plt.tick_params(axis='both', which='major', labelsize=20)

    # Do subplot to zoom in on interface
    a = plt.axes([.65, .25, .3, .35])
    a.set_xlim(43,57)
    a.set_ylim(0,.0001)
    if (native > 0):
        nats = a.plot(x_amanziU_native, c_amanziU_native,'rx',label='Amanzi+Alquimia(PFloTran)',linewidth=2)
    if (unstruct > 0):
        alqs = a.plot(x_amanziU, c_amanziU,'r-',label='Amanzi+Alquimia(PFloTran) - 1st Order',linewidth=2)
    if (unstruct_crunch > 0):
        alqs_crunch = a.plot(x_amanziU_crunch, c_amanziU_crunch,'r*',label='Amanzi+Alquimia(CrunchFlow)',linewidth=2)
    if (struct>0):
        sams = a.plot(x_amanziS, c_amanziS,'g-',label='AmanziS+Alq(PFT)',linewidth=2) 
    if (struct_crunch>0):
        sams_crunch = a.plot(x_amanziS_crunch, c_amanziS_crunch,'g*',label='AmanziS+Alq(CF)',linewidth=2) 

    pfls = a.plot(x_pflotran, c_pflotran,'m-',label='PFloTran',linewidth=2)
    cfs = a.plot(x_crunchflow, c_crunchflow,'m*',label='CrunchFlow',linewidth=2)
    plt.title('Zoom near interface')
    
    # plt.show()
    plt.savefig(root+"_1d.png",format="png")
    # plt.close()

    # finally:
    #     pass
