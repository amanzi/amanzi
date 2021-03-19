# plots calcium concentration along x at last time step 
# benchmark: compares to pflotran simulation results
# author: S.Molins - Sept. 2013

import os
import sys
import h5py
import numpy as np
import matplotlib
# matplotlib.use('Agg')
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
from compare_field_results import GetXY_AmanziU_1D
from compare_field_results import GetXY_AmanziS_1D
from compare_field_results import GetXY_PFloTran_1D
from compare_field_results import GetXY_CrunchFlow_1D


if __name__ == "__main__":

    import os,sys

    try:
        sys.path.append('../../../../MY_TPL_BUILD/ccse/ccse-1.3.4-source/Tools/Py_util')
    except:
        pass
    
    try:
        PF_DIR=os.environ.get('PARFLOW_RESULTS_DIR')
        print('Parflow results directory: ',PF_DIR)
        sys.path.append(PF_DIR)
        from PF_GetXY import GetXY_ParFlow_1D_100
    except:
        print('error: parflow directory not found in ',PF_DIR)
        raise

    # root name for problem
    root = "calcite"

    # pflotran
    path_to_pflotran = "pflotran"
    root_pflotran = "1d-"+root

    # -- hardwired for 1d-calcite: time and comp
    times = ['Time:  5.00000E+01 y']

    comp = 'Total_Ca++ [M]'
    Ca_pflotran = []
    for i, time in enumerate(times):
        x_pflotran, c_pflotran = GetXY_PFloTran_1D(path_to_pflotran,root_pflotran,time,comp)
        Ca_pflotran = Ca_pflotran + [c_pflotran]

    comp = 'pH'
    pH_pflotran = []
    for i, time in enumerate(times):
        x_pflotran, c_pflotran = GetXY_PFloTran_1D(path_to_pflotran,root_pflotran,time,comp)
        pH_pflotran = pH_pflotran + [c_pflotran]

    comp = 'Calcite_VF'
    VF_pflotran = []
    for i, time in enumerate(times):
        x_pflotran, c_pflotran = GetXY_PFloTran_1D(path_to_pflotran,root_pflotran,time,comp)
        VF_pflotran = VF_pflotran + [c_pflotran]

    # -- pflotran: Operator Splitting
    path_to_pflotran = "pflotran/os"

    # hardwired for 1d-calcite: time and comp
    times = ['Time:  5.00000E+01 y']

    comp = 'Total_Ca++ [M]'
    Ca_pflotran_OS = []
    for i, time in enumerate(times):
        x_pflotran_OS, c_pflotran = GetXY_PFloTran_1D(path_to_pflotran,root_pflotran,time,comp)
        Ca_pflotran_OS = Ca_pflotran_OS + [c_pflotran]

    comp = 'pH'
    pH_pflotran_OS = []
    for i, time in enumerate(times):
        x_pflotran_OS, c_pflotran = GetXY_PFloTran_1D(path_to_pflotran,root_pflotran,time,comp)
        pH_pflotran_OS = pH_pflotran_OS + [c_pflotran]

    comp = 'Calcite_VF'
    VF_pflotran_OS = []
    for i, time in enumerate(times):
        x_pflotran_OS, c_pflotran = GetXY_PFloTran_1D(path_to_pflotran,root_pflotran,time,comp)
        VF_pflotran_OS = VF_pflotran_OS + [c_pflotran]

    # crunchflow GIMRT
    path_to_crunchflow = "crunchflow/gimrt"

    # hardwired for calcite_1d_CF.in: time and comp
    times_CF = ['totcon5.out']
    comp = 2
    Ca_crunchflow = []
    ignore = 4
    for i, time in enumerate(times_CF):
        x_crunchflow, c_crunchflow = GetXY_CrunchFlow_1D(path_to_crunchflow,root,time,comp,ignore)
        Ca_crunchflow = Ca_crunchflow + [c_crunchflow]

    # times_CF = ['pH1.out','pH2.out','pH3.out','pH4.out','pH5.out']
    times_CF = ['pH5.out']
    comp = 0
    pH_crunchflow = []
    ignore = 3
    for i, time in enumerate(times_CF):
        x_crunchflow, c_crunchflow = GetXY_CrunchFlow_1D(path_to_crunchflow,root,time,comp,ignore)
        pH_crunchflow = pH_crunchflow + [c_crunchflow]

    # times_CF = ['volume1.out','volume2.out','volume3.out','volume4.out','volume5.out']
    times_CF = ['volume5.out']
    comp = 0
    VF_crunchflow = []
    ignore = 4
    for i, time in enumerate(times_CF):
        x_crunchflow, c_crunchflow = GetXY_CrunchFlow_1D(path_to_crunchflow,root,time,comp,ignore)
        VF_crunchflow = VF_crunchflow + [c_crunchflow/100.0]

    # crunchflow OS3D
    path_to_crunchflow = "crunchflow/os3d"

    # hardwired for calcite_1d_CF.in: time and comp
    # times_CF = ['totcon1.out','totcon2.out','totcon3.out','totcon4.out','totcon5.out']
    times_CF = ['totcon5.out']
    comp = 2
    Ca_crunchOS3D = []
    ignore = 4
    for i, time in enumerate(times_CF):
        x_crunchflow, c_crunchflow = GetXY_CrunchFlow_1D(path_to_crunchflow,root,time,comp,ignore)
        Ca_crunchOS3D = Ca_crunchOS3D + [c_crunchflow]

    # times_CF = ['pH1.out','pH2.out','pH3.out','pH4.out','pH5.out']
    times_CF = ['pH5.out']
    comp = 0
    pH_crunchOS3D = []
    ignore = 3
    for i, time in enumerate(times_CF):
        x_crunchflow, c_crunchflow = GetXY_CrunchFlow_1D(path_to_crunchflow,root,time,comp,ignore)
        pH_crunchOS3D = pH_crunchOS3D + [c_crunchflow]

    # times_CF = ['volume1.out','volume2.out','volume3.out','volume4.out','volume5.out']
    times_CF = ['volume5.out']
    comp = 0
    VF_crunchOS3D = []
    ignore = 4
    for i, time in enumerate(times_CF):
        x_crunchflow, c_crunchflow = GetXY_CrunchFlow_1D(path_to_crunchflow,root,time,comp,ignore)
        VF_crunchOS3D = VF_crunchOS3D + [c_crunchflow/100.0]

    CWD = os.getcwd()
    local_path = "" 


    # subplots ----------------------------------------------------------
    fig, ax = plt.subplots(3,figsize=(8,10))

    # Amanzi U + Native chemistry
    try:
        input_file = os.path.join("amanzi-u-1d-calcite.xml")
        path_to_amanzi = "output-u"
        run_amanzi_standard.run_amanzi(input_file, 1, ["calcite.bgd",input_file], path_to_amanzi)
        
        comp = 'total_component_concentration.cell.Ca++ conc'
        Ca_amanzi_native = []
        for i, time in enumerate(times):
            x_amanzi_native, c_amanzi_native = GetXY_AmanziU_1D(path_to_amanzi,root,comp,1)
            Ca_amanzi_native = Ca_amanzi_native +[c_amanzi_native]

        comp = 'free_ion_species.cell.H+'
        pH_amanzi_native = []
        for i, time in enumerate(times):
            x_amanzi_native, c_amanzi_native = GetXY_AmanziU_1D(path_to_amanzi,root,comp,1)
            pH_amanzi_native = pH_amanzi_native +[-np.log10(c_amanzi_native)]

        comp = 'mineral_volume_fractions.cell.Calcite vol frac'
        VF_amanzi_native = []
        for i, time in enumerate(times):
            x_amanzi_native, c_amanzi_native = GetXY_AmanziU_1D(path_to_amanzi,root,comp,1)
            VF_amanzi_native = VF_amanzi_native +[c_amanzi_native]

        native = True

    except:
        native = False
        pass    

    #HARDWIRE to exclude native chemistry
    native=False

    # Amanzi U + Alquimia + PFloTran chemistry
    try:
        input_file = os.path.join("amanzi-u-1d-calcite-alq-pflo.xml")
        path_to_amanzi = "output-u-alq-pflo"
        run_amanzi_standard.run_amanzi(input_file, 1, ["1d-calcite-trim.in","calcite.dat",input_file], path_to_amanzi)

        comp = 'total_component_concentration.cell.Ca++ conc'
        Ca_amanzi_alquimia = []
        for i, time in enumerate(times):
            x_amanzi_alquimia, c_amanzi_alquimia = GetXY_AmanziU_1D(path_to_amanzi,root,comp,1)
            Ca_amanzi_alquimia = Ca_amanzi_alquimia +[c_amanzi_alquimia]

        comp = 'free_ion_species.cell.H+'
        pH_amanzi_alquimia = []
        for i, time in enumerate(times):
            x_amanzi_alquimia, c_amanzi_alquimia = GetXY_AmanziU_1D(path_to_amanzi,root,comp,1)
            pH_amanzi_alquimia = pH_amanzi_alquimia +[-np.log10(c_amanzi_alquimia)]

        comp = 'mineral_volume_fractions.cell.Calcite vol frac'
        VF_amanzi_alquimia = []
        for i, time in enumerate(times):
            x_amanzi_alquimia, c_amanzi_alquimia = GetXY_AmanziU_1D(path_to_amanzi,root,comp,1)
            VF_amanzi_alquimia = VF_amanzi_alquimia +[c_amanzi_alquimia]

        alq = True

    except:
        alq = False


    # Amanzi U + Alquimia + CruchFlow chemistry
    try:
        input_file = os.path.join("amanzi-u-1d-calcite-alq-crunch.xml")
        path_to_amanzi = "output-u-alq-crunch"
        run_amanzi_standard.run_amanzi(input_file, 1, ["1d-calcite-crunch.in","calcite.dbs",input_file], path_to_amanzi)

        comp = 'total_component_concentration.cell.Ca++ conc'
        Ca_amanzi_alquimia_crunch = []
        for i, time in enumerate(times):
            x_amanzi_alquimia_crunch, c_amanzi_alquimia_crunch = GetXY_AmanziU_1D(path_to_amanzi,root,comp,1)
            Ca_amanzi_alquimia_crunch = Ca_amanzi_alquimia_crunch +[c_amanzi_alquimia_crunch]

        comp = 'free_ion_species.cell.H+'
        pH_amanzi_alquimia_crunch = []
        for i, time in enumerate(times):
           x_amanzi_alquimia_crunch, c_amanzi_alquimia_crunch = GetXY_AmanziU_1D(path_to_amanzi,root,comp,1)
           pH_amanzi_alquimia_crunch = pH_amanzi_alquimia_crunch +[-np.log10(c_amanzi_alquimia_crunch)]

        comp = 'mineral_volume_fractions.cell.Calcite vol frac'
        VF_amanzi_alquimia_crunch = []
        for i, time in enumerate(times):
           x_amanzi_alquimia_crunch, c_amanzi_alquimia_crunch = GetXY_AmanziU_1D(path_to_amanzi,root,comp,1)
           VF_amanzi_alquimia_crunch = VF_amanzi_alquimia_crunch +[c_amanzi_alquimia_crunch]

        alq_crunch = True

    except:
        alq_crunch = False

    
    # Amanzi S + Alquimia + PFloTran chemistry
    try:
        # import pdb; pdb.set_trace()
        input_file = os.path.join("amanzi-s-1d-calcite-alq-pflo.xml")
        path_to_amanzi = "output-s-alq-pflo"
        run_amanzi_standard.run_amanzi(input_file, 1, ["1d-calcite-trim.in","calcite.dat",input_file], path_to_amanzi)
        
        root_amanziS = "plt"
        compS = "Ca++_water_Concentration"
        x_amanziS, c_amanziS = GetXY_AmanziS_1D(path_to_amanzi,root_amanziS,compS,1)
        struct = len(x_amanziS)
        compS = "H+_Free_Ion_Guess"
        x_amanziS, pH_amanziS = GetXY_AmanziS_1D(path_to_amanzi,root_amanziS,compS,1)
        pH_amanziS =  -np.log10(pH_amanziS)
        compS = "Calcite_Volume_Fraction"
        x_amanziS, VF_amanziS = GetXY_AmanziS_1D(path_to_amanzi,root_amanziS,compS,1)
    except:
        struct = 0


    # Amanzi S + Alquimia + CrunchFlow chemistry
    try:
        # import pdb; pdb.set_trace()
        input_file = os.path.join("amanzi-s-1d-calcite-alq-crunch.xml")
        path_to_amanzi = "output-s-alq-crunch"
        run_amanzi_standard.run_amanzi(input_file, 1, ["1d-calcite-crunch.in","calcite.dbs",input_file], path_to_amanzi)

        root_amanziS = "plt"
        compS = "Ca++_water_Concentration"
        x_amanziS_crunch, c_amanziS_crunch = GetXY_AmanziS_1D(path_to_amanzi,root_amanziS,compS,1)
        struct_c = len(x_amanziS_crunch)
        compS = "H+_Free_Ion_Guess"
        x_amanziS_crunch, pH_amanziS_crunch = GetXY_AmanziS_1D(path_to_amanzi,root_amanziS,compS,1)
        pH_amanziS_crunch = -np.log10(pH_amanziS_crunch)
        compS = "Calcite_Volume_Fraction"
        x_amanziS_crunch, VF_amanziS_crunch = GetXY_AmanziS_1D(path_to_amanzi,root_amanziS,compS,1)
    except:
        struct_c = 0

#    import pdb; pdb.set_trace()
        
    # parflow + pflotran
    try:
        path_to_parflow = os.path.join(PF_DIR,root+'_1d','pflotran')
        
        compPF = "calcite_pf.out.PrimaryMobile.02.Ca++.00005.txt"
        x_parflow_pflo, c_parflow_pflo = GetXY_ParFlow_1D_100(compPF,path=path_to_parflow)
        parflow_pflo = len(x_parflow_pflo)

        compPF = "calcite_pf.out.pH.00005.txt"
        x_parflow_pflo, pH_parflow_pflo = GetXY_ParFlow_1D_100(compPF,path=path_to_parflow)

        compPF = "calcite_pf.out.MineralVolfx.00.Calcite.00005.txt"
        x_parflow_pflo, VF_parflow_pflo = GetXY_ParFlow_1D_100(compPF,path=path_to_parflow)

    except:
        parflow_pflo = 0

    # parflow + crunch
    try:
        path_to_parflow = os.path.join(PF_DIR,root+'_1d','crunch')
        
        compPF = "calcite_pf.out.PrimaryMobile.02.Ca++.00005.txt"
        x_parflow_crunch, c_parflow_crunch = GetXY_ParFlow_1D_100(compPF,path=path_to_parflow)
        parflow_crunch = len(x_parflow_crunch)

        compPF = "calcite_pf.out.pH.00005.txt"
        x_parflow_crunch, pH_parflow_crunch = GetXY_ParFlow_1D_100(compPF,path=path_to_parflow)

        compPF = "calcite_pf.out.MineralVolfx.00.Calcite.00005.txt"
        x_parflow_crunch, VF_parflow_crunch = GetXY_ParFlow_1D_100(compPF,path=path_to_parflow)

    except:
        parflow_crunch = 0


    # subplots
    # fig, ax = plt.subplots(3,sharex=True,figsize=(8,10))
      
    for i, time in enumerate(times):
        # lines on axes
        if alq: 
            ax[0].plot(x_amanzi_alquimia, Ca_amanzi_alquimia[i],'r-',linewidth=2)
        if alq_crunch:
            ax[0].plot(x_amanzi_alquimia_crunch, Ca_amanzi_alquimia_crunch[i],'r*',linewidth=2,markersize=12)
        if native:
            ax[0].plot(x_amanzi_native, Ca_amanzi_native[i],'rx')

        ax[0].plot(x_pflotran_OS, Ca_pflotran_OS[i],'m-',linewidth=2)
        ax[0].plot(x_crunchflow, Ca_crunchOS3D[i],'m*')

        if alq:
            ax[1].plot(x_amanzi_alquimia, pH_amanzi_alquimia[i],'r-',linewidth=2)#,label='AmanziU(2nd-O)+Alq(PFLOTRAN)')
        if alq_crunch:
            ax[1].plot(x_amanzi_alquimia_crunch, pH_amanzi_alquimia_crunch[i],'r*')#,linewidth=2,markersize=12,label='AmanziU(2nd-O)+Alq(CrunchFlow)')
        if native:
            ax[1].plot(x_amanzi_native, pH_amanzi_native[i],'rx',label='AmanziU(2nd-O) Native Chem.')

        ax[1].plot(x_pflotran_OS, pH_pflotran_OS[i],'m-',linewidth=2)
        ax[1].plot(x_crunchflow, pH_crunchOS3D[i],'m*')

        if i==0:
            if alq:
                ax[2].plot(x_amanzi_alquimia, VF_amanzi_alquimia[i],'r-',linewidth=2,label='AmanziU+Alq(PFLOTRAN)')
            if alq_crunch:
                ax[2].plot(x_amanzi_alquimia_crunch, VF_amanzi_alquimia_crunch[i],'r*',linewidth=2,markersize=12,label='AmanziU+Alq(CrunchFlow)')
            if native:
                ax[2].plot(x_amanzi_native, VF_amanzi_native[i],'rx',label='AmanziU(2nd-O)+Native')

            ax[2].plot(x_pflotran_OS, VF_pflotran_OS[i],'m-',label='PFLOTRAN',linewidth=2)
            ax[2].plot(x_crunchflow, VF_crunchOS3D[i],'m*',label='CrunchFlow OS3D',markersize=7)
        else:
            if alq:
                ax[2].plot(x_amanzi_alquimia, VF_amanzi_alquimia[i],'r-',linewidth=2)
            if alq_crunch:
                ax[2].plot(x_amanzi_alquimia_crunch, VF_amanzi_alquimia_crunch[i],'r*',linewidth=2)
            if native:
                ax[2].plot(x_amanzi_native, VF_amanzi_native[i],'rx')

            ax[2].plot(x_pflotran_OS, VF_pflotran_OS[i],'m-',linewidth=2)

    #import pdb; pdb.set_trace()
    if (struct>0):
        sam = ax[0].plot(x_amanziS, c_amanziS,'g-',linewidth=2)#label='AmanziS+Alq(PFLOTRAN)',linewidth=2)
        sampH = ax[1].plot(x_amanziS, pH_amanziS,'g-',linewidth=2)
        samVF = ax[2].plot(x_amanziS, VF_amanziS,'g-',linewidth=2,label='AmanziS+Alq(PFLOTRAN)')

    if (struct_c>0):
        samc = ax[0].plot(x_amanziS_crunch, c_amanziS_crunch,'g*',linewidth=2)#,label='AmanziS+Alq(CrunchFlow)',linewidth=2)     
        samcpH = ax[1].plot(x_amanziS_crunch, pH_amanziS_crunch,'g*',linewidth=2)     
        samcVF = ax[2].plot(x_amanziS_crunch, VF_amanziS_crunch,'g*',linewidth=2,label='AmanziS+Alq(CrunchFlow)')

    # parflow + pflotran
    if (parflow_pflo>0):
        pfpfC  = ax[0].plot(x_parflow_pflo, c_parflow_pflo,'b-',linewidth=2)#label='Parflow+Alq(PFLOTRAN)',linewidth=2)
        pfpfpH = ax[1].plot(x_parflow_pflo, pH_parflow_pflo,'b-',linewidth=2)
        pfpfVF = ax[2].plot(x_parflow_pflo, VF_parflow_pflo,'b-',linewidth=2,label='Parflow+Alq(PFLOTRAN)')
        
    # parflow + crunch    
    if (parflow_crunch>0):
        pfcfC  = ax[0].plot(x_parflow_crunch, c_parflow_crunch,'b*',linewidth=2)#label='Parflow+Alq(CrunchFlow)',linewidth=2)
        pfcfpH = ax[1].plot(x_parflow_crunch, pH_parflow_crunch,'b*',linewidth=2)
        pfcfVF = ax[2].plot(x_parflow_crunch, VF_parflow_crunch,'b*',linewidth=2,label='Parflow+Alq(CrunchFlow)')

    # set x lim
    # ax[0].set_xlim((18,32))
  
    # axes
    ax[0].set_xlabel("Distance (m)",fontsize=15)
    ax[1].set_xlabel("Distance (m)",fontsize=15)
    ax[2].set_xlabel("Distance (m)",fontsize=15)
    ax[0].set_ylabel("Total Ca conc. [mol/L]",fontsize=15)
    ax[1].set_ylabel("pH",fontsize=15)
    ax[2].set_ylabel("Calcite volume fraction",fontsize=15)

    # plot adjustments
#    plt.subplots_adjust(left=0.15,bottom=0.15,right=0.99,top=0.90)
    #import pdb; pdb.set_trace()
#    if (struct>0 or struct_c>0):
#        ax[0].legend(loc='lower right',fontsize=14)
#
#    if (alq>0 or alq_crunch>0):
#        ax[1].legend(loc='lower right',fontsize=14) 

    ax[2].legend(loc='lower right',fontsize=13)
#    plt.suptitle("Amanzi 1D Calcite Benchmark",x=0.57,fontsize=20)
    ax[0].tick_params(axis='x', which='major', labelsize=15)
    ax[1].tick_params(axis='x', which='major', labelsize=15)
    ax[2].tick_params(axis='x', which='major', labelsize=15)
    ax[0].tick_params(axis='y', which='major', labelsize=15)
    ax[1].tick_params(axis='y', which='major', labelsize=15)
    ax[2].tick_params(axis='y', which='major', labelsize=15)

    ax[0].ticklabel_format(axis='y',style='scientific',scilimits=(0,0))
    ax[2].ticklabel_format(axis='y',style='scientific',scilimits=(0,0))
    
#    ax[2].set_ylim((0,1e-5))
    plt.tight_layout(pad=0.25)
    
    # pyplot.show()
    plt.savefig(local_path+"calcite_1d.png",format="png")
    # plt.close()

    # set x lim
    # ax[0].set_xlim((30,70))
    # plt.savefig(local_path+"calcite_1d_2.png",format="png")

    # finally:
    #     pass 
