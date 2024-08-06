# plots cation concentration along x at last time step 
# benchmark: compares to pflotran simulation results
# author: S.Molins - Sept. 2013

import os
import sys
import h5py
import numpy as np
import matplotlib
# matplotlib.use('Agg')
from matplotlib import pyplot as plt

import run_amanzi_standard
from compare_field_results import GetXY_AmanziU_1D
from compare_field_results import GetXY_AmanziS_1D
from compare_field_results import GetXY_PFloTran_1D
from compare_field_results import GetXY_CrunchFlow_1D

AXES_TICK_SIZE=15
PFLOTRAN_LINE_COLOR='m'
LINE_WIDTH=4.5
MARKER_SIZE=7

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
    root = "ion-exchange"

    # pflotran
    path_to_pflotran = "pflotran"
    root_pflotran = "1d-"+root

    # hardwired for 1d-calcite: time and comp
    times = ['Time:  5.00000E+01 y',]
    components = ['Na+','Ca++','Mg++','Cl-']

    try:     
        u_pflotran = [[[] for x in range(len(components))] for x in range(len(times))]
        for i, time in enumerate(times):
            for j, comp in enumerate(components):          
                x_pflotran, c_pflotran = GetXY_PFloTran_1D(path_to_pflotran,root_pflotran,time,'Total_'+comp+' [M]')
                u_pflotran[i][j] = c_pflotran
    
        v_pflotran = [[[] for x in range(len(components))] for x in range(len(times))]
        for i, time in enumerate(times):
            for j, comp in enumerate(components):          
                x_pflotran, c_pflotran = GetXY_PFloTran_1D(path_to_pflotran,root_pflotran,time,'Total_Sorbed_'+comp+' [mol_m^3]')
                v_pflotran[i][j] = c_pflotran

        pflotran = True
   
    except:
        pflotran = False


    # CrunchFlow
    path_to_crunchflow = "crunchflow"

    # hardwired for 1d-tritium-crunch.in: time and comp
    times_CF = ['totcon5.out']
    components = [0,1,2,3]
    ignore = 4

    try:
        u_crunchflow = [[[] for x in range(len(components))] for x in range(len(times_CF))]
        for i, time in enumerate(times_CF):
            for j, comp in enumerate(components):          
                x_crunchflow, c_crunchflow = GetXY_CrunchFlow_1D(path_to_crunchflow,root,time,comp,ignore)
                u_crunchflow[i][j] = c_crunchflow

        crunch = True

    except:
        crunch = False
 
    times_CF = ['totexchange5.out']
    components = [0,1,2,3]
    ignore = 4

    try:
      v_crunchflow = [[[] for x in range(len(components))] for x in range(len(times_CF))]
      for i, time in enumerate(times_CF):
         for j, comp in enumerate(components):
            x_crunchflow, c_crunchflow = GetXY_CrunchFlow_1D(path_to_crunchflow,root,time,comp,ignore)
            v_crunchflow[i][j] = c_crunchflow

      crunch = True

    except:
      crunch = False

    # hardwired for 1d-exchange:  last time = '71'
    times = ['71']
    amanzi_components = ['total_component_concentration.Na+', \
                         'total_component_concentration.Ca++', \
                         'total_component_concentration.Mg++', \
                         'total_component_concentration.Cl-']
    amanzi_sorbed     = ['total_sorbed.0', \
                         'total_sorbed.1', \
                         'total_sorbed.2', \
                         'total_sorbed.3']

    amanzi_compS      = ['Na+_water_Concentration', \
                         'Ca++_water_Concentration', \
                         'Mg++_water_Concentration', \
                         'Cl-_water_Concentration']
    amanzi_sorbS      = ['Na+_Sorbed_Concentration', \
                         'Ca++_Sorbed_Concentration', \
                         'Mg++_Sorbed_Concentration', \
                         'Cl-_Sorbed_Concentration']

    # AmanziU + Native chemistry
    try:   
        input_file = os.path.join("amanzi-u-1d-"+root+".xml")
        path_to_amanzi = "output-u"
        run_amanzi_standard.run_amanzi(input_file, 1, [input_file], path_to_amanzi)

        u_amanzi_native = [[[] for x in range(len(amanzi_components))] for x in range(len(times))]
        for i, time in enumerate(times):
            for j, comp in enumerate(amanzi_components):
                x_amanzi_native, c_amanzi_native = GetXY_AmanziU_1D(path_to_amanzi,root,comp,1)
                u_amanzi_native[i][j] = c_amanzi_native

        v_amanzi_native = [[[] for x in range(len(amanzi_sorbed))] for x in range(len(times))]
        for i, time in enumerate(times):
            for j, comp in enumerate(amanzi_sorbed):
                x_amanzi_native, c_amanzi_native = GetXY_AmanziU_1D(path_to_amanzi,root,comp,1)
                v_amanzi_native[i][j] = c_amanzi_native

        native = True

    except:  
        native = False


    # AmanziU + Alquimia + PFloTran chemistry
    try:
        input_file = os.path.join("amanzi-u-1d-"+root+"-alq-pflo.xml")
        path_to_amanzi = "output-u-alq-pflo"
        run_amanzi_standard.run_amanzi(input_file, 1, ["1d-"+root+"-trim.in",root+".dat",input_file], path_to_amanzi)

        u_amanzi_alquimia = [[[] for x in range(len(amanzi_components))] for x in range(len(times))]
        for i, time in enumerate(times):
            for j, comp in enumerate(amanzi_components):
                x_amanzi_alquimia, c_amanzi_alquimia = GetXY_AmanziU_1D(path_to_amanzi,root,comp,1)
                u_amanzi_alquimia[i][j] = c_amanzi_alquimia
              
        v_amanzi_alquimia = [[[] for x in range(len(amanzi_sorbed))] for x in range(len(times))]
        for i, time in enumerate(times):
            for j, comp in enumerate(amanzi_sorbed):
                x_amanzi_alquimia, c_amanzi_alquimia = GetXY_AmanziU_1D(path_to_amanzi,root,comp,1)
                v_amanzi_alquimia[i][j] = c_amanzi_alquimia

        alq = True

    except: 
        alq = False


    # AmanziU + Alquimia + PFloTran chemistry with writer
    try:
        input_file = os.path.join("amanzi-u-1d-"+root+"-alq-pflo-writer.xml")
        path_to_amanzi = "output-u-alq-pflo-writer"
        run_amanzi_standard.run_amanzi(input_file, 1, ["1d-"+root+"-trim.in",root+".dat",input_file], path_to_amanzi)

        u_amanzi_alquimia_w = [[[] for x in range(len(amanzi_components))] for x in range(len(times))]
        for i, time in enumerate(times):
            for j, comp in enumerate(amanzi_components):
                x_amanzi_alquimia_w, c_amanzi_alquimia_w = GetXY_AmanziU_1D(path_to_amanzi,root,comp,1)
                u_amanzi_alquimia_w[i][j] = c_amanzi_alquimia_w
              
        v_amanzi_alquimia_w = [[[] for x in range(len(amanzi_sorbed))] for x in range(len(times))]
        for i, time in enumerate(times):
            for j, comp in enumerate(amanzi_sorbed):
                x_amanzi_alquimia_w, c_amanzi_alquimia_w = GetXY_AmanziU_1D(path_to_amanzi,root,comp,1)
                v_amanzi_alquimia_w[i][j] = c_amanzi_alquimia_w

        alq_writer = True

    except: 
        alq_writer = False
        

    # AmanziU + Alquimia + CrunchFlow chemistry
    try:
        input_file = os.path.join("amanzi-u-1d-"+root+"-alq-crunch.xml")
        path_to_amanzi = "output-u-alq-crunch"
        run_amanzi_standard.run_amanzi(input_file, 1, ["1d-"+root+"-crunch.in",root+".dbs",input_file], path_to_amanzi)

        u_amanzi_alquimia_crunch = [[[] for x in range(len(amanzi_components))] for x in range(len(times))]
        for i, time in enumerate(times):
            for j, comp in enumerate(amanzi_components):
                x_amanzi_alquimia_crunch, c_amanzi_alquimia_crunch = GetXY_AmanziU_1D(path_to_amanzi,root,comp,1)
                u_amanzi_alquimia_crunch[i][j] = c_amanzi_alquimia_crunch
              
        v_amanzi_alquimia_crunch = [[[] for x in range(len(amanzi_sorbed))] for x in range(len(times))]
        for i, time in enumerate(times):
            for j, comp in enumerate(amanzi_sorbed):
                x_amanzi_alquimia_crunch, c_amanzi_alquimia_crunch = GetXY_AmanziU_1D(path_to_amanzi,root,comp,1)
                v_amanzi_alquimia_crunch[i][j] = c_amanzi_alquimia_crunch

        alq_crunch = True

    except: 
        alq_crunch = False


    # AmanziS + Alquimia + PLloTran
    try:
        input_file = os.path.join("amanzi-s-1d-ion-exchange-alq-pflo.xml")
        path_to_amanzi = "output-s-alq-pflo"
        run_amanzi_standard.run_amanzi(input_file, 1, ["1d-"+root+"-trim.in",root+".dat",input_file], path_to_amanzi)

        root_amanziS = "plt"
        # compS = "Na+_Aqueous_Concentration"
        # x_amanziS, c_amanziS = GetXY_AmanziS_1D(path_to_amanzi,root_amanziS,compS,1)

        # import pdb; pdb.set_trace()
        u_amanziS = [[] for x in range(len(amanzi_compS))]
        for j, compS in enumerate(amanzi_compS):
            x_amanziS, c_amanziS = GetXY_AmanziS_1D(path_to_amanzi,root_amanziS,compS,1)
            u_amanziS[j] = c_amanziS

        # compS = "Na+_Sorbed_Concentration"
        # x_amanziS, v_amanziS = GetXY_AmanziS_1D(path_to_amanzi,root_amanziS,compS,1)

        v_amanziS = [[] for x in range(len(amanzi_sorbS))]
        for j, compS in enumerate(amanzi_sorbS):
            x_amanziS, c_amanziS = GetXY_AmanziS_1D(path_to_amanzi,root_amanziS,compS,1)
            v_amanziS[j] = c_amanziS

        struct = len(x_amanziS)

    except:
        struct = 0


    # AmanziS + Alquimia + CrunchFlow
    try:
        input_file = os.path.join("amanzi-s-1d-ion-exchange-alq-crunch.xml")
        path_to_amanzi = "output-s-alq-crunch"
        run_amanzi_standard.run_amanzi(input_file, 1, ["1d-"+root+"-crunch.in",root+".dbs",input_file], path_to_amanzi)

        root_amanziS = "plt"
        # compS = "Na+_Aqueous_Concentration"
        # x_amanziS, c_amanziS = GetXY_AmanziS_1D(path_to_amanzi,root_amanziS,compS,1)

        # import pdb; pdb.set_trace()
        u_amanziS_c = [[] for x in range(len(amanzi_compS))]
        for j, compS in enumerate(amanzi_compS):
            x_amanziS_c, c_amanziS = GetXY_AmanziS_1D(path_to_amanzi,root_amanziS,compS,1)
            u_amanziS_c[j] = c_amanziS

        # compS = "Na+_Sorbed_Concentration"
        # x_amanziS, v_amanziS = GetXY_AmanziS_1D(path_to_amanzi,root_amanziS,compS,1)

        v_amanziS_c = [[] for x in range(len(amanzi_sorbS))]
        for j, compS in enumerate(amanzi_sorbS):
            x_amanziS_c, c_amanziS = GetXY_AmanziS_1D(path_to_amanzi,root_amanziS,compS,1)
            v_amanziS_c[j] = c_amanziS

        struct_c = len(x_amanziS_c)

    except:
        struct_c = 0


    components = ['Na+','Ca++','Mg++','Cl-']

    # define subplots as a subgrid
    fig = plt.figure(figsize=(15,12))

    ax = [] # total concentrations
    bx = [] # sorbed concentrations

    nrows = len(components)
    ncols = 2

    for comp in range(len(components)):
        ax += [plt.subplot(nrows,ncols,2*comp+1)]
        bx += [plt.subplot(nrows,ncols,2*comp+2)]

    if alq:
        i = 0  # hardwired for time '71'
        for j, comp in enumerate(components):
            ax[j].plot(x_amanzi_alquimia, u_amanzi_alquimia[i][j],'r-',label='AmanziU+Alq(PFT)',linewidth=LINE_WIDTH)
            bx[j].plot(x_amanzi_alquimia, v_amanzi_alquimia[i][j],'r-',label='AmanziU+Alq(PFT)',linewidth=LINE_WIDTH)

    if alq_writer:
        i = 0  # hardwired for time '71'
        for j, comp in enumerate(components):
            ax[j].plot(x_amanzi_alquimia_w, u_amanzi_alquimia_w[i][j],'r--',label='AmanziU+Alq(PFT)-W',linewidth=LINE_WIDTH)
            bx[j].plot(x_amanzi_alquimia_w, v_amanzi_alquimia_w[i][j],'r--',label='AmanziU+Alq(PFT)-W',linewidth=LINE_WIDTH)

            
    if alq_crunch:
        i = 0  # hardwired for time '71'
        for j, comp in enumerate(components):
            ax[j].plot(x_amanzi_alquimia_crunch, u_amanzi_alquimia_crunch[i][j],'r*',label='AmanziU+Alq(CF)',markersize=MARKER_SIZE)
            bx[j].plot(x_amanzi_alquimia_crunch, v_amanzi_alquimia_crunch[i][j],'r*',label='AmanziU+Alq(CF)',markersize=MARKER_SIZE)

    if native:
        i = 0  # hardwired for time '71'
        for j, comp in enumerate(components):
            ax[j].plot(x_amanzi_native, u_amanzi_native[i][j],'x',label='AmanziU Native Chem',c=PFLOTRAN_LINE_COLOR)
            bx[j].plot(x_amanzi_native, v_amanzi_native[i][j],'rx',label='AmanziU Native Chem')

    if crunch:
        i = 0  # hardwired for time 50 years
        for j, comp in enumerate(components):
            ax[j].plot(x_crunchflow, u_crunchflow[i][j],'b*',label='CrunchFlow',markersize=MARKER_SIZE)
            bx[j].plot(x_crunchflow, v_crunchflow[i][j],'b*',label='CrunchFlow',markersize=MARKER_SIZE)
    
    if pflotran:
        i = 0  # hardwired for time 50 years
        for j, comp in enumerate(components):
            ax[j].plot(x_pflotran, u_pflotran[i][j],'-',label='PFloTran',c=PFLOTRAN_LINE_COLOR,linewidth=LINE_WIDTH)
            bx[j].plot(x_pflotran, v_pflotran[i][j],'-',label='PFloTran',c=PFLOTRAN_LINE_COLOR,linewidth=LINE_WIDTH)

            ax[j].text(
                0.03,
                0.85,
                comp,
                fontsize=15,
                bbox=dict(facecolor='white', alpha=1.0),
                transform=ax[j].transAxes,
                horizontalalignment='left'
            )

            bx[j].text(
                0.97,
                0.85,
                comp,
                fontsize=15,
                bbox=dict(facecolor='white', alpha=1.0),
                transform=bx[j].transAxes,
                horizontalalignment='right'
            )

    if (struct>0):
        for j in range(len(amanzi_compS)):
            ax[j].plot(x_amanziS, u_amanziS[j],'-',c=PFLOTRAN_LINE_COLOR)#,label='AmanziS+Alq(PFT)')
            bx[j].plot(x_amanziS, v_amanziS[j],'g-',label='AmanziS+Alq(PFT)')

    if (struct_c > 0):
        for j in range(len(amanzi_compS)):
            ax[j].plot(x_amanziS_c, u_amanziS_c[j],'*',c=PFLOTRAN_LINE_COLOR)#,label='AmanziS+Alq(CF)')
            bx[j].plot(x_amanziS_c, v_amanziS_c[j],'g*',label='AmanziS+Alq(CF)')

    # axes
    ax[nrows-1].set_xlabel("Distance (m)",fontsize=15)
    bx[nrows-1].set_xlabel("Distance (m)",fontsize=15)

    # for i,comp in enumerate(components):
    for i,comp in enumerate(components):
        ax[i].set_ylabel("Total Concen. [mol/L]",fontsize=15)
        bx[i].set_ylabel("Total Sorb. [mol/m3]",fontsize=15)
        bx[i].yaxis.tick_right()
        bx[i].yaxis.set_label_position("right")

    # Set axes to be shared + only have ticks on the bottom row
    ax[0].get_shared_x_axes().join(*[ax[i] for i in range(nrows)])
    bx[0].get_shared_x_axes().join(*[bx[i] for i in range(nrows)])

    for i in range(nrows-1):
        ax[i].set_xticklabels([])
        bx[i].set_xticklabels([])

    ax[nrows-1].autoscale()
    bx[nrows-1].autoscale()

    # plot adjustments
    ax[0].legend(fontsize=12,loc='lower right')
    #bx[0].legend(fontsize=12,loc='lower right')

    plt.suptitle("Amanzi 1D "+root.title()+" Benchmark at 50 years",fontsize=20) #,x=0.57,fontsize=20)

    # Set scientific notation ticks and tick size
    from matplotlib.ticker import FormatStrFormatter
    for i in range(nrows):
        ax[i].ticklabel_format(useMathText=True,axis='y',style='sci',scilimits=(-2,-2))
        bx[i].ticklabel_format(useMathText=True,axis='y',style='sci',scilimits=(0,0))
        ax[i].tick_params(axis='both', which='major', labelsize=AXES_TICK_SIZE)
        bx[i].tick_params(axis='both', which='major', labelsize=AXES_TICK_SIZE)
  
    plt.tight_layout() #(pad=0.4, w_pad=0.5, h_pad=1.0)

    plt.subplots_adjust(left=0.10,bottom=0.15,right=0.90,top=0.95)

    # plt.show()
    plt.savefig(root+"_1d.png",format="png")
    # plt.close()

    # finally:
    #     pass 
