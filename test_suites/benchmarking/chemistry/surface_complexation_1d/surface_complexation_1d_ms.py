# plots cation concentration along x at last time step 
# benchmark: compares to pflotran simulation results
# author: S.Molins - Oct. 2013

import os
import sys
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

try:
    sys.path.append('../../../../tools/amanzi_xml')
except:
    pass
    
try:
    sys.path.append('../../../../tools/testing')
except:
    pass

import run_amanzi_standard
import numpy as np
from compare_field_results import GetXY_AmanziU_1D
from compare_field_results import GetXY_AmanziS_1D
from compare_field_results import GetXY_PFloTran_1D
from compare_field_results import GetXY_CrunchFlow_1D


if __name__ == "__main__":

    try:
        sys.path.append('../../../../MY_TPL_BUILD/ccse/ccse-1.3.4-source/Tools/Py_util')
    except:
        pass

    try:
        PF_DIR='/home/smolins/work/PF-Alquimia_verification'
        sys.path.append(PF_DIR)
        from PF_GetXY import GetXY_ParFlow_1D_100
    except:
        print('error: parflow directory not found in ',PF_DIR)
        raise


# root name for problem
    root = "surface-complexation"
    components = ['Na+','NO3-','Zn++'] ## ['H+','Na+','NO3-','Zn++']
    pflotran_totc_templ = "Total_{0} [M]"
    comppflo = [pflotran_totc_templ.format(x) for x in components]
    compcrunch = [1,2,3] 

# times
    timespflo = ['Time:  5.00000E+01 y']
    times = ['71']
    times_CF = ['totcon5.out']
    times_CF_surf = ['totsurface5.out']
    times_CF_pH = ['pH5.out']

# amanzi output (unstructured - native and alquimia)
    amanzi_totc_templ = "total_component_concentration.cell.{0} conc" #Component {0} conc"
    amanzi_totc = [amanzi_totc_templ.format(x) for x in components] #range(len(components))]
    amanzi_totc_crunch = [amanzi_totc_templ.format(x) for x in compcrunch] #range(len(components))]

    amanzi_sorb_templ = "total_sorbed.cell.{0}"
    amanzi_sorb = [amanzi_sorb_templ.format(x+1) for x in range(len(components))]
    amanzi_sorb_crunch = [amanzi_sorb_templ.format(x+1) for x in range(len(compcrunch))]

# amanzi output (structured - alquimia)
    amanzi_totc_templ = "{0}_water_Concentration"
    amanzi_totcS = [amanzi_totc_templ.format(x) for x in components]

    amanzi_sorb_templ = "{0}_Sorbed_Concentration"
    amanzi_sorbS = [amanzi_sorb_templ.format(x) for x in components]

# read pflotran results --->
    path_to_pflotran = "pflotran"
    root_pflotran = "1d-"+root

    u_pflotran = [[[] for x in range(len(comppflo))] for x in range(len(timespflo))]
    for i, time in enumerate(timespflo):
        for j, comp in enumerate(comppflo):          
            x_pflotran, c_pflotran = GetXY_PFloTran_1D(path_to_pflotran,root_pflotran,time,comp)
            u_pflotran[i][j] = c_pflotran
    
    v_pflotran = [[[] for x in range(len(components))] for x in range(len(timespflo))]
    for i, time in enumerate(timespflo):
        for j, comp in enumerate(components):
            x_pflotran, c_pflotran = GetXY_PFloTran_1D(path_to_pflotran,root_pflotran,time,'Total_Sorbed_'+comp+' [mol_m^3]')
            v_pflotran[i][j] = c_pflotran

    pH_pflotran = [[] for x in range(len(timespflo))]
    comp = 'pH'
    for i, time in enumerate(timespflo):
        x_pflotran, c_pflotran = GetXY_PFloTran_1D(path_to_pflotran,root_pflotran,time,comp)
        pH_pflotran[i] = c_pflotran

# read crunchflow results --->
    path_to_crunch = "crunchflow"

    try: 
        u_crunchflow = [ [ [] for x in range(len(compcrunch)) ] for x in range(len(times_CF)) ]
        ignore = 4
        for i, time in enumerate(times_CF):
            for j, comp in enumerate(compcrunch):
                x_crunchflow, c_crunchflow = GetXY_CrunchFlow_1D(path_to_crunch,root,time,comp,ignore)
                u_crunchflow[i][j] = c_crunchflow

        v_crunchflow = [ [ [] for x in range(len(compcrunch)) ] for x in range(len(times_CF_surf)) ]
        ignore = 4
        for i, time in enumerate(times_CF_surf):
            for j, comp in enumerate(compcrunch):
                x_crunchflow, c_crunchflow = GetXY_CrunchFlow_1D(path_to_crunch,root,time,comp,ignore)
                v_crunchflow[i][j] = c_crunchflow

        pH_crunchflow = [ [] for x in range(len(times_CF_pH)) ]
        ignore = 4
        comp = 0
        for i, time in enumerate(times_CF_pH):
            y_crunchflow, c_crunchflow = GetXY_CrunchFlow_1D(path_to_crunch,root,time,comp,ignore)
            pH_crunchflow[i] = c_crunchflow

        crunch = True

    except: 
        crunch = False


# AmanziU  + Native chemistry --->
    try:
        input_file = os.path.join("amanzi-u-1d-"+root+".xml")
        path_to_amanzi = "output-u"
        run_amanzi_standard.run_amanzi(input_file, 1,
                                       [root+".bgd",input_file],
                                       path_to_amanzi)
        
        u_amanzi_native = [[[] for x in range(len(amanzi_totc))] for x in range(len(times))]
        for i, time in enumerate(times):
            for j, comp in enumerate(amanzi_totc):
                x_amanzi_native, c_amanzi_native = GetXY_AmanziU_1D(path_to_amanzi,root,comp,1)
                u_amanzi_native[i][j] = c_amanzi_native

        v_amanzi_native = [[[] for x in range(len(amanzi_sorb))] for x in range(len(times))]
        for i, time in enumerate(times):
            for j, comp in enumerate(amanzi_sorb):
                x_amanzi_native, c_amanzi_native = GetXY_AmanziU_1D(path_to_amanzi,root,comp,1)
                v_amanzi_native[i][j] = c_amanzi_native

        pH_amanzi_native = [ [] for x in range(len(times)) ]
        comp = 'free_ion_species.cell.H+'
        for i, time in enumerate(times):
            x_amanzi_native, c_amanzi_native = GetXY_AmanziU_1D(path_to_amanzi,root,comp,1)
            pH_amanzi_native[i] = -np.log10(c_amanzi_native)

    native = True

    except:
        native = False

# HARDWIRE to exclude native
    native = False
        
# AmanziU +  Alquimia + PFlotran chemistry --->
    try:  
        input_file = os.path.join("amanzi-u-1d-"+root+"-alq-pflo.xml")
        path_to_amanzi = "output-u-alq-pflo"
        run_amanzi_standard.run_amanzi(input_file, 1,
                                       ["1d-"+root+".in",root+".dat",input_file],
                                       path_to_amanzi)

        u_amanzi_alquimia = [[[] for x in range(len(amanzi_totc))] for x in range(len(times))]
        for i, time in enumerate(times):
            for j, comp in enumerate(amanzi_totc):
                x_amanzi_alquimia, c_amanzi_alquimia = GetXY_AmanziU_1D(path_to_amanzi,root,comp,1)
                u_amanzi_alquimia[i][j] = c_amanzi_alquimia

        v_amanzi_alquimia = [[[] for x in range(len(amanzi_sorb))] for x in range(len(times))]
        for i, time in enumerate(times):
            for j, comp in enumerate(amanzi_sorb):
                x_amanzi_alquimia, c_amanzi_alquimia = GetXY_AmanziU_1D(path_to_amanzi,root,comp,1)
                v_amanzi_alquimia[i][j] = c_amanzi_alquimia

        pH_amanzi_alquimia = [ [] for x in range(len(times)) ]
        comp = 'pH.cell.0'
        for i, time in enumerate(times):
             x_amanzi_alquimia, c_amanzi_alquimia = GetXY_AmanziU_1D(path_to_amanzi,root,comp,1)
             pH_amanzi_alquimia[i] = c_amanzi_alquimia ## -np.log10(c_amanzi_native)

        alq = True

    except:
        alq = False


# AmanziU + Alquimia + CrunchFlow chemistry --->
    try:
        input_file = os.path.join("amanzi-u-1d-surface-complexation-alq-crunch.xml")
        path_to_amanzi = "output-u-alq-crunch"
        run_amanzi_standard.run_amanzi(input_file, 1,
                                       ["1d-surface-complexation-crunch.in","surface-complexation.dbs",input_file],
                                       path_to_amanzi)

        u_amanzi_alquimia_crunch = [[[] for x in range(len(amanzi_totc))] for x in range(len(times))]
        for i, time in enumerate(times):
            for j, comp in enumerate(amanzi_totc):
                x_amanzi_alquimia_crunch, c_amanzi_alquimia_crunch = GetXY_AmanziU_1D(path_to_amanzi,root,comp,1)
                u_amanzi_alquimia_crunch[i][j] = c_amanzi_alquimia_crunch
              
        v_amanzi_alquimia_crunch = [[[] for x in range(len(amanzi_sorb))] for x in range(len(times))]
        for i, time in enumerate(times):
            for j, comp in enumerate(amanzi_sorb):
                x_amanzi_alquimia_crunch, c_amanzi_alquimia_crunch = GetXY_AmanziU_1D(path_to_amanzi,root,comp,1)
                v_amanzi_alquimia_crunch[i][j] = c_amanzi_alquimia_crunch

        pH_amanzi_alquimia_crunch = [ [] for x in range(len(times)) ]
        comp = 'pH.cell.0'
        for i, time in enumerate(times):
            x_amanzi_alquimia_crunch, c_amanzi_alquimia_crunch = GetXY_AmanziU_1D(path_to_amanzi,root,comp,1)
            pH_amanzi_alquimia_crunch[i] = c_amanzi_alquimia_crunch ## -np.log10(c_amanzi_native)

        alq_crunch = True

    except:
        alq_crunch = False


# AmanziS + Alquimia + PFlowTran chemistry --->
    try:
        input_file = os.path.join("amanzi-s-1d-surface-complexation-alq-pflo.xml")
        path_to_amanzi = "output-s-alq-pflo"
        run_amanzi_standard.run_amanzi(input_file, 1,
                                       ["1d-"+root+".in",root+".dat",input_file],
                                       path_to_amanzi)
        root_amanziS = "plt"

        u_amanziS = [[[] for x in range(len(amanzi_totcS))] for x in range(len(times))]
        for i, time in enumerate(times):
            for j, comp in enumerate(amanzi_totcS):
                # import pdb; pdb.set_trace()
                x_amanziS, c_amanziS = GetXY_AmanziS_1D(path_to_amanzi,root_amanziS,comp,1)
                u_amanziS[i][j] = c_amanziS
              
        struct = len(x_amanziS)

        # import pdb; pdb.set_trace()

        v_amanziS = [[[] for x in range(len(amanzi_sorbS))] for x in range(len(times))]
        for i, time in enumerate(times):
            for j, comp in enumerate(amanzi_sorbS):
                x_amanziS, c_amanziS = GetXY_AmanziS_1D(path_to_amanzi,root_amanziS,comp,1)
                v_amanziS[i][j] = c_amanziS

        pH_amanziS = [ [] for x in range(len(times)) ]
        comp = 'H+_Free_Ion_Guess'
        for i, time in enumerate(times):
            x_amanziS, c_amanziS = GetXY_AmanziS_1D(path_to_amanzi,root_amanziS,comp,1)
            pH_amanziS[i] = -np.log10(c_amanziS)

    except:
        struct = 0


# AmanziS + Alquimia + CrunchFlow chemistry --->
    try:
        input_file = os.path.join("amanzi-s-1d-surface-complexation-alq-crunch.xml")
        path_to_amanzi = "output-s-alq-crunch"
        # import pdb; pdb.set_trace()
        run_amanzi_standard.run_amanzi(input_file, 1,
                                       ["1d-surface-complexation-crunch.in","surface-complexation.dbs",input_file],
                                       path_to_amanzi)
        root_amanziS = "plt"

        u_amanziS_crunch = [[[] for x in range(len(amanzi_totcS))] for x in range(len(times))]
        for i, time in enumerate(times):
            for j, comp in enumerate(amanzi_totcS):
                # import pdb; pdb.set_trace()
                x_amanziS_crunch, c_amanziS_crunch = GetXY_AmanziS_1D(path_to_amanzi,root_amanziS,comp,1)
                u_amanziS_crunch[i][j] = c_amanziS_crunch
              
        struct_c = len(x_amanziS_crunch)

        v_amanziS_crunch = [[[] for x in range(len(amanzi_sorbS))] for x in range(len(times))]
        for i, time in enumerate(times):
            for j, comp in enumerate(amanzi_sorbS):
                x_amanziS_crunch, c_amanziS_crunch = GetXY_AmanziS_1D(path_to_amanzi,root_amanziS,comp,1)
                v_amanziS_crunch[i][j] = c_amanziS_crunch

        pH_amanziS_crunch = [ [] for x in range(len(times)) ]
        comp = 'H+_Free_Ion_Guess'
        for i, time in enumerate(times):
            x_amanziS_crunch, c_amanziS_crunch = GetXY_AmanziS_1D(path_to_amanzi,root_amanziS,comp,1)
            pH_amanziS_crunch[i] = -np.log10(c_amanziS_crunch)

    except:
        struct_c = 0


    root_ = 'surface_complexation'
    
    # parflow + pflotran
    try:
        path_to_parflow = os.path.join(PF_DIR,root_+'_1d','run_pflotran')

        c_parflow_pflo = [ [] for x in range(len(components)) ]
        v_parflow_pflo = [ [] for x in range(len(components)) ]

        template = "surface_complexation_pf.out.PrimaryMobile.0{}.{}.00005.txt"
        for i,c in enumerate(components):
            compPF = template.format(i+1,c)        
            x_parflow_pflo, c_parflow_pflo[i] = GetXY_ParFlow_1D_100(compPF,path=path_to_parflow)
            parflow_pflo = len(x_parflow_pflo)

        template = "surface_complexation_pf.out.PrimarySorbed.0{}.{}.00005.txt"
        for i,c in enumerate(components):
            compPF = template.format(i+1,c)
            x_parflow_pflo, v_parflow_pflo[i] = GetXY_ParFlow_1D_100(compPF,path=path_to_parflow)

        compPF = "surface_complexation_pf.out.pH.00005.txt"
        x_parflow_pflo, pH_parflow_pflo = GetXY_ParFlow_1D_100(compPF,path=path_to_parflow)


    except:
        parflow_pflo = 0

        
    # parflow + crunch
    try:
        path_to_parflow = os.path.join(PF_DIR,root_+'_1d','run_crunch')

        c_parflow_crunch = [ [] for x in range(len(components)) ]
        v_parflow_crunch = [ [] for x in range(len(components)) ]
        
        template = "surface_complexation_pf.out.PrimaryMobile.0{}.{}.00005.txt"
        for i,c in enumerate(components):
            compPF = template.format(i+1,c)
            x_parflow_crunch, c_parflow_crunch[i] = GetXY_ParFlow_1D_100(compPF,path=path_to_parflow)
            parflow_crunch = len(x_parflow_crunch)

        template = "surface_complexation_pf.out.PrimarySorbed.0{}.{}.00005.txt"
        for i,c in enumerate(components):
            compPF = template.format(i+1,c)
            x_parflow_crunch, v_parflow_crunch[i] = GetXY_ParFlow_1D_100(compPF,path=path_to_parflow)

        compPF = "surface_complexation_pf.out.pH.00005.txt"
        x_parflow_crunch, pH_parflow_crunch = GetXY_ParFlow_1D_100(compPF,path=path_to_parflow)


    except:
        parflow_crunch = 0
        

## ------------------------------------------------------------------------
    # colors= ['r','b','m','g'] # components
    # styles = ['-','--','x'] # codes
    # codes = ['Amanzi+Alquimia(PFloTran)','Amanzi Native Chemistry','PFloTran'] + [None,]*9

    # lines on axes
    # ax[0],b[0] ---> Aqueous concentrations
    # ax[1]      ---> pH
    # ax[2],b[2] ---> Sorbed concentrations

    # first
    # ax[0],b[0] ---> Aqueous concentrations
    # ax[1]      ---> pH

    # define subplots as a subgrid
    fig = plt.figure(figsize=(15,12))

    ax = [] # total concentrations
    bx = [] # sorbed concentrations

    nrows = len(components) + 1
    ncols = 2 

    for comp in range(len(components)):
        ax += [plt.subplot(nrows,ncols,2*comp+1)]
        bx += [plt.subplot(nrows,ncols,2*comp+2)]

    px = plt.subplot(nrows,1,nrows)

    # for i, time in enumerate(times):
    i = 0 # only one time point at 50 years

    # pflotran 
    # if True:
    for j, comp in enumerate(components):
        ax[j].plot(x_pflotran, u_pflotran[i][j],color='m',linestyle='-',linewidth=6,label='PFloTran')
        bx[j].plot(x_pflotran, v_pflotran[i][j],color='m',linestyle='-',linewidth=2)
        ax[j].text(x_pflotran[10],u_pflotran[i][j][10],comp,fontsize=15,bbox=dict(facecolor='white', alpha=1.0))
        bx[j].text(x_pflotran[10],v_pflotran[i][j][10],comp,fontsize=15,bbox=dict(facecolor='white', alpha=1.0))

    px.plot(x_pflotran, pH_pflotran[i],color='m',linestyle='-',linewidth=2,label='PFloTran')


    # crunchflow 
    # if True:
    for j, comp in enumerate(components):
        ax[j].plot(x_crunchflow, u_crunchflow[i][j],'m*',linestyle='None',label='CrunchFlow')
        bx[j].plot(x_crunchflow, v_crunchflow[i][j],'m*',linestyle='None',linewidth=2)

    px.plot(y_crunchflow, pH_crunchflow[i],'m*',linestyle='None',linewidth=2,label='CrunchFlow')


    # amanzi-native
    if native :
        for j, comp in enumerate(components):
            ax[j].plot(x_amanzi_native, u_amanzi_native[i][j],color='b',linestyle='None',marker='x',linewidth=2)
            bx[j].plot(x_amanzi_native, v_amanzi_native[i][j],color='b',linestyle='None',marker='x',linewidth=2,label='AmanziU Native Chemistry')

        px.plot(x_amanzi_native, pH_amanzi_native[i],color='b',linestyle='None',marker='x',linewidth=2,label='AmanziU Native Chemistry')


    # amanzi-unstructured-alquimia-pflotran
    if alq:
        for j, comp in enumerate(components):
            ax[j].plot(x_amanzi_alquimia, u_amanzi_alquimia[i][j],color='r',linestyle='-',linewidth=2)
            bx[j].plot(x_amanzi_alquimia, v_amanzi_alquimia[i][j],color='r',linestyle='-',linewidth=2,label='AmanziU+Alquimia(PFloTran)')

        px.plot(x_amanzi_alquimia, pH_amanzi_alquimia[i],color='r',linestyle='-',linewidth=2,label='AmanziU+Alquimia(PFloTran)')


    # amanzi-unstructured-alquimia-crunchflow
    if alq_crunch:
        for j, comp in enumerate(components):
            ax[j].plot(x_amanzi_alquimia_crunch, u_amanzi_alquimia_crunch[i][j],color='r',linestyle='None',marker='*',linewidth=2)
            bx[j].plot(x_amanzi_alquimia_crunch, v_amanzi_alquimia_crunch[i][j],color='r',linestyle='None',marker='*',linewidth=2,label='AmanziU+Alquimia(CrunchFlow)')

        px.plot(x_amanzi_alquimia_crunch, pH_amanzi_alquimia_crunch[i],color='r',linestyle='None',marker='*',linewidth=2,label='AmanziU+Alquimia(CrunchFlow)')


    # amanzi-structured-alquimia-pflotran
    if struct:
        for j, comp in enumerate(components):
            ax[j].plot(x_amanziS, u_amanziS[i][j],color='g',linestyle='-',linewidth=2,label='AmanziS+Alquimia(PFloTran)')
            bx[j].plot(x_amanziS, v_amanziS[i][j],color='g',linestyle='-',linewidth=2,label='AmanziS+Alquimia(PFloTran)')

        px.plot(x_amanziS, pH_amanziS[i],color='g',linestyle='-',linewidth=2,label='AmanziS+Alquimia(PFloTran)')


    # amanzi-structured-alquimia-crunchflow
    if struct_c:
        for j, comp in enumerate(components):
            ax[j].plot(x_amanziS_crunch, u_amanziS_crunch[i][j],color='g',linestyle='None',marker='*',linewidth=2,label='AmanziS+Alquimia(CrunchFlow)')
            bx[j].plot(x_amanziS_crunch, v_amanziS_crunch[i][j],color='g',linestyle='None',marker='*',linewidth=2,label='AmanziS+Alquimia(CrunchFlow)')

        px.plot(x_amanziS_crunch, pH_amanziS_crunch[i],color='g',linestyle='None',marker='*',linewidth=2,label='AmanziS+Alquimia(CrunchFlow)')

    # parflow + pflotran
    if (parflow_pflo > 0):
        for j in range(len(components)):
            ax[j].plot(x_parflow_pflo, c_parflow_pflo[j],'b-')
            bx[j].plot(x_parflow_pflo, v_parflow_pflo[j],'b-',label='Parflow+Alq(PFT)')

        px.plot(x_parflow_pflo, pH_parflow_pflo,'b-',linewidth=2)

    # parflow + crunch
    if (parflow_crunch > 0):
        for j in range(len(components)):
            ax[j].plot(x_parflow_crunch, c_parflow_crunch[j],'b*')
            bx[j].plot(x_parflow_crunch, v_parflow_crunch[j],'b*',label='Parflow+Alq(CF)')

        px.plot(x_parflow_crunch, pH_parflow_crunch,'b*',linewidth=2)


    # axes
    ax[len(components)-1].set_xlabel("Distance (m)",fontsize=15)
    bx[len(components)-1].set_xlabel("Distance (m)",fontsize=15)

    # for i,comp in enumerate(components):
    i=1
    ax[i].set_ylabel("Total Concentration [mol/L]",fontsize=15)
    bx[i].set_ylabel("Total Sorbed Concent. [mol/m3]",fontsize=15)

    px.set_xlabel("Distance(m)",fontsize=15)
    px.set_ylabel("pH",fontsize=15)

    # for i,comp in enumerate(components):
    #     ax[i].set_ylim(bottom=0)
    #     bx[i].set_ylim(bottom=0)
    px.set_ylim(bottom=4.8)

    # plot adjustments
    ax[0].legend(fontsize=15)
    bx[0].legend(fontsize=15)
    px.legend(fontsize=15,loc='upper right')

    plt.suptitle("Amanzi 1D "+root.title()+" Benchmark at 50 years",fontsize=20) #,x=0.57,fontsize=20)

    plt.tick_params(axis='both', which='major', labelsize=15)
  
    plt.tight_layout() #(pad=0.4, w_pad=0.5, h_pad=1.0)

    plt.subplots_adjust(left=0.10,bottom=0.15,right=0.90,top=0.95)

    # pyplot.show()
    plt.savefig(root+"_1d.png",format="png")
    # plt.close()

    # finally:
    #     pass 
