# plots cation concentration along x at last time step 
# benchmark: compares to pflotran simulation results
# author: S.Molins - Nov. 2013
# modified: E.I.Barker - May 2014
#           - added native chemistry using v2 input spec

import os
import sys
import h5py
import numpy as np
import matplotlib
from matplotlib import pyplot as plt

import run_amanzi_standard
from compare_field_results import GetXY_AmanziU_1D
from compare_field_results import GetXY_AmanziS_1D
from compare_field_results import GetXY_PFloTran_1D
from compare_field_results import GetXY_CrunchFlow_1D

AXES_TICK_SIZE=15
LINE_WIDTH=2

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
    root = "isotherms"
    components = ['A','B','C']
    compcrunch = ['A']

    # times
    # timespfl = ['Time:  0.00000E+00 y', 'Time:  5.00000E+01 y',
    timespfl = ['Time:  5.00000E+01 y',]
    timesama  = ['71']
    timesama2 = ['71']

    # pflotran output
    pflotran_totc_templ = "Total_{0} [M]"
    pflotran_totc = [pflotran_totc_templ.format(x) for x in components]

    pflotran_sorb_templ = "Total_Sorbed_{0} [mol_m^3]"
    pflotran_sorb = [pflotran_sorb_templ.format(x) for x in components]

    # amanzi output
    amanzi_totc_templ = "total_component_concentration.cell.{}" #Component {0} conc"
    amanzi_totc = [amanzi_totc_templ.format(x) for x in components] #range(len(components))]
    amanzi_totc_crunch = [amanzi_totc_templ.format(x) for x in compcrunch] #range(len(components))]

    amanzi_sorb_templ = "total_sorbed.cell.{0}"
    amanzi_sorb = [amanzi_sorb_templ.format(x) for x in range(len(components))]
    amanzi_sorb_crunch = [amanzi_sorb_templ.format(x) for x in range(len(compcrunch))]


# pflotran data --->
    path_to_pflotran = "pflotran"
    root_pflotran = "1d-"+root

    u_pflotran = [[[] for x in range(len(pflotran_totc))] for x in range(len(timespfl))]
    for i, time in enumerate(timespfl):
        for j, comp in enumerate(pflotran_totc):
            x_pflotran, c_pflotran = GetXY_PFloTran_1D(path_to_pflotran,root_pflotran,time,comp)
            u_pflotran[i][j] = c_pflotran
    
    v_pflotran = [[[] for x in range(len(pflotran_sorb))] for x in range(len(timespfl))]
    for i, time in enumerate(timespfl):
        for j, sorb in enumerate(pflotran_sorb):
            x_pflotran, c_pflotran = GetXY_PFloTran_1D(path_to_pflotran,root_pflotran,time,sorb)
            v_pflotran[i][j] = c_pflotran

    CWD = os.getcwd()
    local_path = "" 


# CrunchFlow --->
    path_to_crunchflow = "crunchflow"
    try: 
        times_CF = ['totcon5.out']
        comp = 0
        u_crunchflow = []
        ignore = 4
        for i, time in enumerate(times_CF):
            x_crunchflow, c_crunchflow = GetXY_CrunchFlow_1D(path_to_crunchflow,root,time,comp,ignore)
            u_crunchflow = u_crunchflow + [c_crunchflow]
        crunch = True

    except: 
        crunch = False    


# AmanziU + Native chemistry --->
    try:
        input_file = os.path.join("amanzi-u-1d-"+root+".xml")
        path_to_amanzi = "output-u"
        run_amanzi_standard.run_amanzi(input_file, 1, [root+".bgd",input_file], path_to_amanzi)

        u_native = [[[] for x in range(len(amanzi_totc))] for x in range(len(timesama))]
        for i, time in enumerate(timesama):
            for j, comp in enumerate(amanzi_totc):
                # import pdb; pdb.set_trace()
                x_native, c_native = GetXY_AmanziU_1D(path_to_amanzi,root,comp,1)
                u_native[i][j] = c_native

        v_native = [[[] for x in range(len(amanzi_sorb))] for x in range(len(timesama))]
        for i, time in enumerate(timesama):
            for j, comp in enumerate(amanzi_sorb):
                x_native, c_native = GetXY_AmanziU_1D(path_to_amanzi,root,comp,1)
                v_native[i][j] = c_native

        native = len(x_native)  

    except:
        native = 0 
        pass


# AmanziU + Alquimia + PFloTran --->
    try:  
        input_file = os.path.join("amanzi-u-1d-"+root+"-alq-pflo.xml")
        path_to_amanzi = "output-u-alq-pflo"
        run_amanzi_standard.run_amanzi(input_file, 1, ["1d-"+root+"-trim.in",root+".dat",input_file], path_to_amanzi)

        u_alquimia = [[[] for x in range(len(amanzi_totc))] for x in range(len(timesama))]
        for i, time in enumerate(timesama):
            for j, comp in enumerate(amanzi_totc):
                x_alquimia, c_alquimia = GetXY_AmanziU_1D(path_to_amanzi,root,comp,1)
                u_alquimia[i][j] = c_alquimia
              
        v_alquimia = [[[] for x in range(len(amanzi_sorb))] for x in range(len(timesama))]
        for i, time in enumerate(timesama):
            for j, comp in enumerate(amanzi_sorb):
                x_alquimia, c_alquimia = GetXY_AmanziU_1D(path_to_amanzi,root,comp,1)
                v_alquimia[i][j] = c_alquimia

        alq = True

    except:
        alq = False

# AmanziU + Alquimia + PFloTran WITH WRITER --->
    try:  
        input_file = os.path.join("amanzi-u-1d-"+root+"-alq-pflo-writer.xml")
        path_to_amanzi = "output-u-alq-pflo-writer"
        run_amanzi_standard.run_amanzi(input_file, 1, ["1d-"+root+"-trim.in",root+".dat",input_file], path_to_amanzi)

        u_alquimia_w = [[[] for x in range(len(amanzi_totc))] for x in range(len(timesama))]
        for i, time in enumerate(timesama):
            for j, comp in enumerate(amanzi_totc):
                x_alquimia_w, c_alquimia_w = GetXY_AmanziU_1D(path_to_amanzi,root,comp,1)
                u_alquimia_w[i][j] = c_alquimia_w
              
        v_alquimia_w = [[[] for x in range(len(amanzi_sorb))] for x in range(len(timesama))]
        for i, time in enumerate(timesama):
            for j, comp in enumerate(amanzi_sorb):
                x_alquimia_w, c_alquimia_w = GetXY_AmanziU_1D(path_to_amanzi,root,comp,1)
                v_alquimia_w[i][j] = c_alquimia_w

        alq_writer = True

    except:
        alq_writer = False      

# AmanziU + Alquimia + CrunchFlow chemistry --->
    try:  
        input_file = os.path.join("amanzi-u-1d-"+root+"-alq-crunch.xml")
        path_to_amanzi = "output-u-alq-crunch"
        run_amanzi_standard.run_amanzi(input_file, 1, ["1d-"+root+"-crunch.in",root+".dbs",input_file], path_to_amanzi)

        u_alquimia_crunch = [[[] for x in range(len(amanzi_totc_crunch))] for x in range(len(timesama))]
        for i, time in enumerate(timesama):
            for j, comp in enumerate(amanzi_totc_crunch):
                x_alquimia_crunch, c_alquimia_crunch = GetXY_AmanziU_1D(path_to_amanzi,root,comp,1)
                u_alquimia_crunch[i][j] = c_alquimia_crunch
              
        v_alquimia_crunch = [[[] for x in range(len(amanzi_sorb_crunch))] for x in range(len(timesama))]
        for i, time in enumerate(timesama):
            for j, comp in enumerate(amanzi_sorb_crunch):
                x_alquimia_crunch, c_alquimia_crunch = GetXY_AmanziU_1D(path_to_amanzi,root,comp,1)
                v_alquimia_crunch[i][j] = c_alquimia_crunch

        alqc = True

    except:
        alqc = False


# AmanziS + Alaqumia + PFloTran chemistry --->
    try:
        input_file = os.path.join("amanzi-s-1d-isotherms-alq-pflo.xml")
        path_to_amanzi = "output-s-alq-pflo"
        run_amanzi_standard.run_amanzi(input_file, 1, ["1d-"+root+"-trim.in",root+".dat",input_file], path_to_amanzi)

        root_amanziS = "plt"
        c_amanziS = [ [] for x in range(len(amanzi_totc)) ]
        v_amanziS = [ [] for x in range(len(amanzi_totc)) ]
        for j,comp in enumerate(components):
            compS = "{0}_water_Concentration".format(comp)
            x_amanziS, c_amanziS[j] = GetXY_AmanziS_1D(path_to_amanzi,root_amanziS,compS,1)
            compS = "{0}_Sorbed_Concentration".format(comp)
            x_amanziS, v_amanziS[j] = GetXY_AmanziS_1D(path_to_amanzi,root_amanziS,compS,1)
        struct = len(x_amanziS)
    except:
        struct = 0


# AmanziS + Alqumia + CrunchFlow chemistry --->
    try:
        input_file = os.path.join("amanzi-s-1d-isotherms-alq-crunch.xml")
        path_to_amanzi = "output-s-alq-crunch"
        run_amanzi_standard.run_amanzi(input_file, 1, ["1d-"+root+"-crunch.in",root+".dbs",input_file], path_to_amanzi)

        root_amanziS = "plt"
        compS = "A_water_Concentration"
        x_amanziS_crunch, c_amanziS_crunch = GetXY_AmanziS_1D(path_to_amanzi,root_amanziS,compS,1)
        compS = "A_Sorbed_Concentration"
        x_amanziS_crunch, v_amanziS_crunch = GetXY_AmanziS_1D(path_to_amanzi,root_amanziS,compS,1)
        struct_c = len(x_amanziS_crunch)
    except:
        struct_c = 0


## plotting ---------------------------------
    # subplots
    fig, ax = plt.subplots(4,sharex=True,figsize=(8,12))

    #  colors= ['r'] #,'b','m','g'] # components
    #  styles = ['-','v','o','x'] # codes
    # codes = ['AmanziU (2nd-Ord.)+Alq(PFT)','AmanziU (2nd-Ord.) Native','PFloTran'] + [None,]*9

    # lines on axes

    # for Kd:
    #   ax[0] ---> Aqueous concentrations
    #   ax[1] ---> Sorbed concentrations

    # for Langmuir and Freundlich
    #   ax[2] ---> Aqueous concentrations
    #   ax[3] ---> Sorbed concentrations

    # for i, time in enumerate(times):
    i = 0 # hardwired 50 years -- because the second entry in the list was taken at cycle 71 = 50 years.

#  pflotran
    ax[0].plot(x_pflotran, u_pflotran[i][0],color='m',linestyle='-',linewidth=LINE_WIDTH,label='PFloTran')
    ax[1].plot(x_pflotran, v_pflotran[i][0],color='m',linestyle='-',linewidth=LINE_WIDTH,label='PFloTran')

    ax[2].plot(x_pflotran, u_pflotran[i][1],color='k',linestyle='-',linewidth=LINE_WIDTH,label='Langmuir PFloTran')
    ax[3].plot(x_pflotran, v_pflotran[i][1],color='k',linestyle='-',linewidth=LINE_WIDTH,label='Langmuir PFloTran')

    ax[2].plot(x_pflotran, u_pflotran[i][2],color='c',linestyle='-',linewidth=LINE_WIDTH,label='Freundlich PFloTran')
    ax[3].plot(x_pflotran, v_pflotran[i][2],color='c',linestyle='-',linewidth=LINE_WIDTH,label='Freundlich PFloTran')

# crunchflow
    if crunch:
        ax[0].plot(x_crunchflow, u_crunchflow[0],color='m',linestyle='None',marker='*', label='CrunchFlow OS3D')
        # crunchflow does not output sorbed concentrations
 
# native 
    if native:
        ax[0].plot(x_native, u_native[i][0],color='b',linestyle='None',marker='x',label='AmanziU (2nd-Ord.) Native')
        ax[1].plot(x_native, v_native[i][0],color='b',linestyle='None',marker='x',label='AmanziU (2nd-Ord.) Native')

        ax[2].plot(x_native, u_native[i][1],color='k',linestyle='None',marker='x',label='Langmuir AmanziU (2nd-Ord.) Native')
        ax[3].plot(x_native, v_native[i][1],color='k',linestyle='None',marker='x',label='Langmuir AmanziU (2nd-Ord.) Native')

        ax[2].plot(x_native, u_native[i][2],color='c',linestyle='None',marker='x',label='Freundlich AmanziU (2nd-Ord.) Native')
        ax[3].plot(x_native, v_native[i][2],color='c',linestyle='None',marker='x',label='Freundlich AmanziU (2nd-Ord.) Native')


# unstructured alquimia pflotran
    if alq:
        ax[0].plot(x_alquimia, u_alquimia[i][0],color='r',linestyle='-',linewidth=LINE_WIDTH,label='AmanziU (2nd-Ord.)+Alq(PFT)')
        ax[1].plot(x_alquimia, v_alquimia[i][0],color='r',linestyle='-',linewidth=LINE_WIDTH,label='AmanziU (2nd-Ord.)+Alq(PFT)')

        ax[2].plot(x_alquimia, u_alquimia[i][1],color='k',linestyle='--',linewidth=LINE_WIDTH,label='Langmuir AmanziU (2nd-Ord.)+Alq(PFT)')
        ax[3].plot(x_alquimia, v_alquimia[i][1],color='k',linestyle='--',linewidth=LINE_WIDTH,label='Langmuir AmanziU (2nd-Ord.)+Alq(PFT)')

        ax[2].plot(x_alquimia, u_alquimia[i][2],color='c',linestyle='--',linewidth=LINE_WIDTH,label='Freundlich AmanziU (2nd-Ord.)+Alq(PFT)')
        ax[3].plot(x_alquimia, v_alquimia[i][2],color='c',linestyle='--',linewidth=LINE_WIDTH,label='Freundlich AmanziU (2nd-Ord.)+Alq(PFT)')

# unstructured alquimia pflotran with writer
    if alq_writer:
        ax[0].plot(x_alquimia_w, u_alquimia_w[i][0],color='r',linestyle='--',linewidth=LINE_WIDTH,label='AmanziU (2nd-Ord.)+Alq(PFT)-W')
        ax[1].plot(x_alquimia_w, v_alquimia_w[i][0],color='r',linestyle='--',linewidth=LINE_WIDTH,label='AmanziU (2nd-Ord.)+Alq(PFT)-W')

        ax[2].plot(x_alquimia_w, u_alquimia_w[i][1],color='k',linestyle='-.',linewidth=LINE_WIDTH,label='Langmuir AmanziU (2nd-Ord.)+Alq(PFT)-W')
        ax[3].plot(x_alquimia_w, v_alquimia_w[i][1],color='k',linestyle='-.',linewidth=LINE_WIDTH,label='Langmuir AmanziU (2nd-Ord.)+Alq(PFT)-W')

        ax[2].plot(x_alquimia_w, u_alquimia_w[i][2],color='c',linestyle='-.',linewidth=LINE_WIDTH,label='Freundlich AmanziU (2nd-Ord.)+Alq(PFT)-W')
        ax[3].plot(x_alquimia_w, v_alquimia_w[i][2],color='c',linestyle='-.',linewidth=LINE_WIDTH,label='Freundlich AmanziU (2nd-Ord.)+Alq(PFT)-W')

# unstructured alquimia crunch
    if alqc:
        ax[0].plot(x_alquimia_crunch, u_alquimia_crunch[i][0],color='r',linestyle='None',marker='*',linewidth=LINE_WIDTH,label='AmanziU (2nd-Ord.)+Alq(CF)')
        ax[1].plot(x_alquimia_crunch, v_alquimia_crunch[i][0],color='r',linestyle='None',marker='*',linewidth=LINE_WIDTH,label='AmanziU (2nd-Ord.)+Alq(CF)')

# structured alquimia pflotran
    if (struct>0):
        sam = ax[0].plot(x_amanziS, c_amanziS[0],'g-',label='AmanziS+Alq(PFT)',linewidth=LINE_WIDTH)
        samv = ax[1].plot(x_amanziS, v_amanziS[0],'g-',linewidth=LINE_WIDTH)

        sam1 = ax[2].plot(x_amanziS, c_amanziS[1],'k*',label='Langmuir AmanziS+Alq(PFT)',linewidth=LINE_WIDTH)
        samv1 = ax[3].plot(x_amanziS, v_amanziS[1],'k*',linewidth=LINE_WIDTH)

        sam2 = ax[2].plot(x_amanziS, c_amanziS[2],'c*',label='Freundlich AmanziS+Alq(PFT)',linewidth=LINE_WIDTH)
        samv2 = ax[3].plot(x_amanziS, v_amanziS[2],'c*',linewidth=LINE_WIDTH)

# structured alquimia crunch
    if (struct_c>0):
        samc = ax[0].plot(x_amanziS_crunch, c_amanziS_crunch,'g*',label='AmanziS+Alq(CF)',linewidth=LINE_WIDTH) 
        samcv = ax[1].plot(x_amanziS_crunch, v_amanziS_crunch,'g*',linewidth=LINE_WIDTH) #,markersize=20) 

    # axes
    ax[0].set_title("Kd linear sorption model",fontsize=15)
    ax[1].set_xlabel("Distance (m)",fontsize=15)
    ax[0].set_ylabel("Total A \n Concentration \n [mol/L]",fontsize=15)
    ax[1].set_ylabel("Total A \n Sorbed Concent. \n [mol/m3]",fontsize=15)

    ax[2].set_title("Langmuir and Freundlich sorption models",fontsize=15)
    ax[3].set_xlabel("Distance (m)",fontsize=15)
    ax[2].set_ylabel("Total B, C \n Concentration \n [mol/L]",fontsize=15)
    ax[3].set_ylabel("Total B, C \n Sorbed Concent. \n [mol/m3]",fontsize=15)

    ax[0].legend(loc='upper right',fontsize=10)
    ax[2].legend(loc='upper right',fontsize=10)

    # plot adjustments
    plt.tight_layout() 
    plt.subplots_adjust(left=0.20,bottom=0.08,right=0.95,top=0.90)
    plt.suptitle("Amanzi 1D "+root.title()+" Benchmark at 50 years",x=0.57,fontsize=20)

    for i in range(4):
        ax[i].set_xlim(left=30,right=70)
        ax[i].tick_params(axis='both', which='major', labelsize=AXES_TICK_SIZE)
        ax[i].ticklabel_format(useMathText=True,axis='y',style='sci',scilimits=(-4,-4))

    # pyplot.show()
    plt.savefig(root+"_1d.png",format="png")
    # plt.close()

    # finally:
    #     pass 
