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
from compare_field_results import GetXY_AmanziU_1D
from compare_field_results import GetXY_AmanziS_1D
from compare_field_results import GetXY_PFloTran_1D
from compare_field_results import GetXY_CrunchFlow_1D

if __name__ == "__main__":

    import os, sys
    try:
        sys.path.append('../../../../tools/amanzi_xml')
    except:
        pass
    import run_amanzi_standard
    import numpy as np

    try:
        sys.path.append('../../../../MY_TPL_BUILD/ccse/ccse-1.3.4-source/Tools/Py_util')
    except:
        pass

    
    # root name for problem
    root = "1d-calcite"

    # pflotran
    path_to_pflotran = "pflotran"
    # path_to_pflotran = "/home/scratch/smolins/amanzi/examples/phase2/chemistry/1d-calcite/pflotran-run/"

     # hardwired for 1d-calcite: time and comp
    # times = ['Time:  1.00000E+01 y','Time:  2.00000E+01 y','Time:  3.00000E+01 y','Time:  4.00000E+01 y','Time:  5.00000E+01 y']
    times = ['Time:  5.00000E+01 y']

    comp = 'Total_Ca++ [M]'
    Ca_pflotran = []
    for i, time in enumerate(times):
       x_pflotran, c_pflotran = GetXY_PFloTran_1D(path_to_pflotran,root,time,comp)
       Ca_pflotran = Ca_pflotran + [c_pflotran]

    comp = 'pH'
    pH_pflotran = []
    for i, time in enumerate(times):
       x_pflotran, c_pflotran = GetXY_PFloTran_1D(path_to_pflotran,root,time,comp)
       pH_pflotran = pH_pflotran + [c_pflotran]

    comp = 'Calcite_VF'
    VF_pflotran = []
    for i, time in enumerate(times):
       x_pflotran, c_pflotran = GetXY_PFloTran_1D(path_to_pflotran,root,time,comp)
       VF_pflotran = VF_pflotran + [c_pflotran]

    # pflotran Operator Splitting
    path_to_pflotran = "pflotran"

     # hardwired for 1d-calcite: time and comp
    times = ['Time:  1.00000E+01 y','Time:  2.00000E+01 y','Time:  3.00000E+01 y','Time:  4.00000E+01 y','Time:  5.00000E+01 y']
    times = ['Time:  5.00000E+01 y']

    comp = 'Total_Ca++ [M]'
    Ca_pflotran_OS = []
    for i, time in enumerate(times):
       x_pflotran_OS, c_pflotran = GetXY_PFloTran_1D(path_to_pflotran,root,time,comp)
       Ca_pflotran_OS = Ca_pflotran_OS + [c_pflotran]

    comp = 'pH'
    pH_pflotran_OS = []
    for i, time in enumerate(times):
       x_pflotran_OS, c_pflotran = GetXY_PFloTran_1D(path_to_pflotran,root,time,comp)
       pH_pflotran_OS = pH_pflotran_OS + [c_pflotran]

    comp = 'Calcite_VF'
    VF_pflotran_OS = []
    for i, time in enumerate(times):
       x_pflotran_OS, c_pflotran = GetXY_PFloTran_1D(path_to_pflotran,root,time,comp)
       VF_pflotran_OS = VF_pflotran_OS + [c_pflotran]

    # crunchflow OS3D
    path_to_crunchflow = "crunchflow"

     # hardwired for calcite_1d_CF.in: time and comp
#    times_CF = ['totcon1.out','totcon2.out','totcon3.out','totcon4.out','totcon5.out']
    times_CF = ['totcon5.out']
    comp = 2
    Ca_crunchOS3D = []
    ignore = 4
    for i, time in enumerate(times_CF):
       x_crunchflow, c_crunchflow = GetXY_CrunchFlow_1D(path_to_crunchflow,root,time,comp,ignore)
       Ca_crunchOS3D = Ca_crunchOS3D + [c_crunchflow]

#    times_CF = ['pH1.out','pH2.out','pH3.out','pH4.out','pH5.out']
    times_CF = ['pH5.out']
    comp = 0
    pH_crunchOS3D = []
    ignore = 3
    for i, time in enumerate(times_CF):
       x_crunchflow, c_crunchflow = GetXY_CrunchFlow_1D(path_to_crunchflow,root,time,comp,ignore)
       pH_crunchOS3D = pH_crunchOS3D + [c_crunchflow]

#    times_CF = ['volume1.out','volume2.out','volume3.out','volume4.out','volume5.out']
    times_CF = ['volume5.out']
    comp = 0
    VF_crunchOS3D = []
    ignore = 4
    for i, time in enumerate(times_CF):
       x_crunchflow, c_crunchflow = GetXY_CrunchFlow_1D(path_to_crunchflow,root,time,comp,ignore)
       VF_crunchOS3D = VF_crunchOS3D + [c_crunchflow]

    CWD = os.getcwd()
    local_path = "" 

    # subplots
    fig, ax = plt.subplots(3,sharex=True,figsize=(8,10))
    
    try:
        # Amanzi native chemistry
        input_file = os.path.join("amanzi-u-1d-calcite.xml")
        path_to_amanzi = "output-u"
        run_amanzi_standard.run_amanzi(input_file, 1,
                                       ["calcite.bgd",input_file], path_to_amanzi)
        
        comp = 'total_component_concentration.Ca++ conc'
        Ca_amanzi_native = []
        for i, time in enumerate(times):
           x_amanzi_native, c_amanzi_native = GetXY_AmanziU_1D(path_to_amanzi,root,comp,1)
           Ca_amanzi_native = Ca_amanzi_native +[c_amanzi_native]

        comp = 'free_ion_species.H+'
        pH_amanzi_native = []
        for i, time in enumerate(times):
           x_amanzi_native, c_amanzi_native = GetXY_AmanziU_1D(path_to_amanzi,root,comp,1)
           pH_amanzi_native = pH_amanzi_native +[-np.log10(c_amanzi_native)]

        comp = 'mineral_volume_fractions.Calcite vol frac'
        VF_amanzi_native = []
        for i, time in enumerate(times):
           x_amanzi_native, c_amanzi_native = GetXY_AmanziU_1D(path_to_amanzi,root,comp,1)
           VF_amanzi_native = VF_amanzi_native +[c_amanzi_native]

        native = True

    except:
        native = False

        pass    

    try:
        # Amanzi-Alquimia
        input_file = os.path.join("amanzi-u-1d-calcite-het-alq-pflo.xml")
        path_to_amanzi = "output-u-alq-pflo"
        run_amanzi_standard.run_amanzi(input_file, 1, 
                                       ["1d-calcite.in","calcite.dat",input_file],
                                       path_to_amanzi)

        comp = 'total_component_concentration.Ca++ conc'
        Ca_amanzi_alquimia = []
        for i, time in enumerate(times):
           x_amanzi_alquimia, c_amanzi_alquimia = GetXY_AmanziU_1D(path_to_amanzi,root,comp,1)
           Ca_amanzi_alquimia = Ca_amanzi_alquimia +[c_amanzi_alquimia]

        comp = 'free_ion_species.H+'
        pH_amanzi_alquimia = []
        for i, time in enumerate(times):
           x_amanzi_alquimia, c_amanzi_alquimia = GetXY_AmanziU_1D(path_to_amanzi,root,comp,1)
           pH_amanzi_alquimia = pH_amanzi_alquimia +[-np.log10(c_amanzi_alquimia)]

        comp = 'mineral_volume_fractions.Calcite vol frac'
        VF_amanzi_alquimia = []
        for i, time in enumerate(times):
           x_amanzi_alquimia, c_amanzi_alquimia = GetXY_AmanziU_1D(path_to_amanzi,root,comp,1)
           VF_amanzi_alquimia = VF_amanzi_alquimia +[c_amanzi_alquimia]

        alq = True

    except:
        
        alq = False

    try:
        # Amanzi-Alquimia-Crunch
        input_file = os.path.join("amanzi-u-1d-calcite-het-alq-crunch.xml")
        path_to_amanzi = "output-u-alq-crunch"
        run_amanzi_standard.run_amanzi(input_file, 1, 
                                       ["1d-calcite-crunch.in","calcite.dbs",input_file],
                                       path_to_amanzi)

        comp = 'total_component_concentration.Ca++ conc'
        Ca_amanzi_alquimia_crunch = []
        for i, time in enumerate(times):
           x_amanzi_alquimia_crunch, c_amanzi_alquimia_crunch = GetXY_AmanziU_1D(path_to_amanzi,root,comp,1)
           Ca_amanzi_alquimia_crunch = Ca_amanzi_alquimia_crunch +[c_amanzi_alquimia_crunch]

        comp = 'free_ion_species.H+'
        pH_amanzi_alquimia_crunch = []
        for i, time in enumerate(times):
           x_amanzi_alquimia_crunch, c_amanzi_alquimia_crunch = GetXY_AmanziU_1D(path_to_amanzi,root,comp,1)
           pH_amanzi_alquimia_crunch = pH_amanzi_alquimia_crunch +[-np.log10(c_amanzi_alquimia_crunch)]

        comp = 'mineral_volume_fractions.Calcite vol frac'
        VF_amanzi_alquimia_crunch = []
        for i, time in enumerate(times):
           x_amanzi_alquimia_crunch, c_amanzi_alquimia_crunch = GetXY_AmanziU_1D(path_to_amanzi,root,comp,1)
           VF_amanzi_alquimia_crunch = VF_amanzi_alquimia_crunch +[c_amanzi_alquimia_crunch]

        alq_crunch = True

    except:
        
        alq_crunch = False

    # amanziS data
    
    # +pflotran
    try:
        # import pdb; pdb.set_trace()
        input_file = os.path.join("amanzi-s-1d-calcite-het-alq.xml")
        path_to_amanzi = "output-s-alq-pflo"
        run_amanzi_standard.run_amanzi(input_file, 1,
                                       ["1d-calcite.in","calcite.dat",input_file],
                                       path_to_amanzi)

        root_amanziS = "plt"
        compS = "Ca++_Aqueous_Concentration"
        x_amanziS, c_amanziS = GetXY_AmanziS_1D(path_to_amanzi,root_amanziS,compS,1)
        struct = len(x_amanziS)
        compS = "H+_Free_Ion_Guess"
        x_amanziS, pH_amanziS = GetXY_AmanziS_1D(path_to_amanzi,root_amanziS,compS,1)
        pH_amanziS =  -np.log10(pH_amanziS)
        compS = "Calcite_Volume_Fraction"
        x_amanziS, VF_amanziS = GetXY_AmanziS_1D(path_to_amanzi,root_amanziS,compS,1)
    except:
        struct = 0

    # +crunchflow
    try:
        # import pdb; pdb.set_trace()
        input_file = os.path.join("amanzi-s-1d-calcite-het-alq-crunch.xml")
        path_to_amanzi = "output-s-alq-crunch"
        run_amanzi_standard.run_amanzi(input_file, 1,
                                       ["1d-calcite-crunch.in","calcite.dbs",input_file],
                                       path_to_amanzi)
        root_amanziS = "plt"
        compS = "Ca++_Aqueous_Concentration"
        x_amanziS_crunch, c_amanziS_crunch = GetXY_AmanziS_1D(path_to_amanzi,root_amanziS,compS,1)
        struct_c = len(x_amanziS_crunch)
        compS = "H+_Free_Ion_Guess"
        x_amanziS_crunch, pH_amanziS_crunch = GetXY_AmanziS_1D(path_to_amanzi,root_amanziS,compS,1)
        pH_amanziS_crunch = -np.log10(pH_amanziS_crunch)
        compS = "Calcite_Volume_Fraction"
        x_amanziS_crunch, VF_amanziS_crunch = GetXY_AmanziS_1D(path_to_amanzi,root_amanziS,compS,1)
    except:
        struct_c = 0

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
##          ax[0].plot(x_pflotran, Ca_pflotran[i],'bx',linewidth=2)
          ax[0].plot(x_pflotran_OS, Ca_pflotran_OS[i],'m-',linewidth=2)
#          if i>1:
#                  ax[0].plot(x_crunchflow, Ca_crunchflow[i-1],'y.')
#                  ax[0].plot(x_crunchflow, Ca_crunchOS3D[i-1],'g.')
##          ax[0].plot(x_crunchflow, Ca_crunchflow[i],'y.')
          ax[0].plot(x_crunchflow, Ca_crunchOS3D[i],'m*')


          if alq:
                 ax[1].plot(x_amanzi_alquimia, pH_amanzi_alquimia[i],'r-',linewidth=2,label='AmanziU(2nd-O)+Alq(PFT)')
          if alq_crunch:
                 ax[1].plot(x_amanzi_alquimia_crunch, pH_amanzi_alquimia_crunch[i],'r*',linewidth=2,markersize=12,label='AmanziU(2nd-O)+Alq(CF)')
          if native:
                 ax[1].plot(x_amanzi_native, pH_amanzi_native[i],'rx',label='AmanziU(2nd-O) Native Chem.')
##          ax[1].plot(x_pflotran, pH_pflotran[i],'bx',linewidth=2)
          ax[1].plot(x_pflotran_OS, pH_pflotran_OS[i],'m-',linewidth=2)
#          if i>0:
#                  ax[1].plot(x_crunchflow, pH_crunchflow[i-1],'y.')
#                  ax[1].plot(x_crunchflow, pH_crunchOS3D[i-1],'g.')
##          ax[1].plot(x_crunchflow, pH_crunchflow[i],'y.')
          ax[1].plot(x_crunchflow, pH_crunchOS3D[i],'m*')

          if i==0:
            if alq:
                   ax[2].plot(x_amanzi_alquimia, VF_amanzi_alquimia[i],'r-',linewidth=2)
            if alq_crunch:
                   ax[2].plot(x_amanzi_alquimia_crunch, VF_amanzi_alquimia_crunch[i],'r*',linewidth=2,markersize=12)
            if native:
                   ax[2].plot(x_amanzi_native, VF_amanzi_native[i],'rx')
##            ax[2].plot(x_pflotran, VF_pflotran[i],'bx',label='PFloTran',linewidth=2)
            ax[2].plot(x_pflotran_OS, VF_pflotran_OS[i],'m-',label='PFloTran OS',linewidth=2)
##            ax[2].plot(x_crunchflow, VF_crunchflow[i],'y.',label='CrunchFlow GIMRT')
            ax[2].plot(x_crunchflow, VF_crunchOS3D[i],'m*',label='CrunchFlow OS3D',markersize=7)
          else:
            if alq:
                   ax[2].plot(x_amanzi_alquimia, VF_amanzi_alquimia[i],'r-',linewidth=2)
            if alq_crunch:
                   ax[2].plot(x_amanzi_alquimia_crunch, VF_amanzi_alquimia_crunch[i],'r*',linewidth=2,markersize=12)
            if native:
                  ax[2].plot(x_amanzi_native, VF_amanzi_native[i],'rx')
##            ax[2].plot(x_pflotran, VF_pflotran[i],'bx',linewidth=2)
            ax[2].plot(x_pflotran_OS, VF_pflotran_OS[i],'m-',linewidth=2)
#            if i==1:
#                  ax[2].plot(x_crunchflow, VF_crunchflow[i-1],'y.',label='CrunchFlow GIMRT')
#                  ax[2].plot(x_crunchflow, VF_crunchOS3D[i-1],'g.',label='CrunchFlow OS3D')
#            else: 
#                  ax[2].plot(x_crunchflow, VF_crunchflow[i-1],'y.')
#                  ax[2].plot(x_crunchflow, VF_crunchOS3D[i-1],'g.')

    #import pdb; pdb.set_trace()
    if (struct>0):
        sam = ax[0].plot(x_amanziS, c_amanziS,'g-',label='AmanziS+Alq(PFT)',linewidth=2)     
        sampH = ax[1].plot(x_amanziS, pH_amanziS,'g-',linewidth=2)
        samVF = ax[2].plot(x_amanziS, VF_amanziS,'g-',linewidth=2)

    if (struct_c>0):
        samc = ax[0].plot(x_amanziS_crunch, c_amanziS_crunch,'g*',label='AmanziS+Alq(CF)',linewidth=2)     
        samcpH = ax[1].plot(x_amanziS_crunch, pH_amanziS_crunch,'g*',linewidth=2)     
        samcVF = ax[2].plot(x_amanziS_crunch, VF_amanziS_crunch,'g*',linewidth=2)     

    # set x lim
    ax[0].set_xlim((30,65))
  
    # axes
    ax[2].set_xlabel("Distance (m)",fontsize=20)
    ax[0].set_ylabel("Total Ca concentration [mol/L]",fontsize=15)
    ax[1].set_ylabel("pH",fontsize=20)
    ax[2].set_ylabel("Calcite volume fraction",fontsize=15)

    # plot adjustments
    plt.subplots_adjust(left=0.15,bottom=0.15,right=0.99,top=0.90)
    #import pdb; pdb.set_trace()
    if (struct>0 or struct_c>0):
        ax[0].legend(loc='lower right',fontsize=13)

    if (alq>0 or alq_crunch>0):
        ax[1].legend(loc='lower right')#,fontsize=13) 

    ax[2].legend(loc='lower right')#,fontsize=13)
    plt.suptitle("Amanzi 1D Calcite Benchmark",x=0.57,fontsize=20)
    plt.tick_params(axis='x', which='major', labelsize=20)

    #pyplot.show()
    plt.savefig(local_path+"calcite_1d.png",format="png")
    #plt.close()

    # set x lim
#    ax[0].set_xlim((30,70))
#    plt.savefig(local_path+"calcite_1d_2.png",format="png")

    #finally:
    #    pass
