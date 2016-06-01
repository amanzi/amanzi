# plots calcium concentration along x at last time step 
# benchmark: compares to pflotran simulation results
# author: S.Molins - Sept. 2013

import os
import sys
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt


# ----------- AMANZI + ALQUIMIA -----------------------------------------------------------------

def GetXY_Amanzi(path,root,comp):

    # open amanzi concentration and mesh files
    dataname = os.path.join(path,root+"_data.h5")
    amanzi_file = h5py.File(dataname,'r')
    meshname = os.path.join(path,root+"_mesh.h5")
    amanzi_mesh = h5py.File(meshname,'r')

    # extract cell coordinates
    y = np.array(amanzi_mesh['0']['Mesh']["Nodes"][0:len(amanzi_mesh['0']['Mesh']["Nodes"])/4,0])
    # y = np.array(amanzi_mesh['Mesh']["Nodes"][0:len(amanzi_mesh['Mesh']["Nodes"])/4,0]) # old style

    # center of cell
    x_amanzi_alquimia  = np.diff(y)/2+y[0:-1]

    # extract concentration array
    time = max(amanzi_file[comp].keys())
    c_amanzi_alquimia = np.array(amanzi_file[comp][time])
    amanzi_file.close()
    amanzi_mesh.close()
    
    return (x_amanzi_alquimia, c_amanzi_alquimia)

def GetXY_AmanziS(path,root,time,comp):
    try:
        import fsnapshot
        fsnok = True
    except:
        fsnok = False

    #import pdb; pdb.set_trace()

    plotfile = os.path.join(path,root)
    if os.path.isdir(plotfile) & fsnok:
        (nx, ny, nz) = fsnapshot.fplotfile_get_size(plotfile)
        x = np.zeros( (nx), dtype=np.float64)
        y = np.zeros( (nx), dtype=np.float64)
        (y, x, npts, err) = fsnapshot.fplotfile_get_data_1d(plotfile, comp, y, x)
    else:
        x = np.zeros( (0), dtype=np.float64)
        y = np.zeros( (0), dtype=np.float64)
    
    return (x, y)

# ----------- PFLOTRAN STANDALONE ------------------------------------------------------------

def GetXY_PFloTran(path,root,time,comp):

    # read pflotran data
    filename = os.path.join(path,"1d-"+root+".h5")
    pfdata = h5py.File(filename,'r')

    # extract coordinates
    y = np.array(pfdata['Coordinates']['X [m]'])
    x_pflotran = np.diff(y)/2+y[0:-1]

    # extract concentrations
    c_pflotran = np.array(pfdata[time][comp])
    c_pflotran = c_pflotran.flatten()
    pfdata.close()

    return (x_pflotran, c_pflotran)

# ------------- CRUNCHFLOW ------------------------------------------------------------------
def GetXY_CrunchFlow(path,root,cf_file,comp,ignore):

    # read CrunchFlow data
    filename = os.path.join(path,cf_file)
    f = open(filename,'r')
    lines = f.readlines()
    f.close()

    # ignore couple of lines
    for i in range(ignore):
      lines.pop(0)

    # extract data x0, x1, ..., xN-1 per line, keep only two columns
    xv=[]
    yv=[] 
    for line in lines:
      xv = xv + [float(line.split()[0])]
      yv = yv + [float(line.split()[comp+1])]
    
    xv = np.array(xv)
    yv = np.array(yv)

    return (xv, yv)

# Main -------------------------------------------------------------------------------------
if __name__ == "__main__":

    import os
    import run_amanzi_standard
    import numpy as np

    # root name for problem
    root = "calcite"

    # pflotran
    path_to_pflotran = "pflotran"
    # path_to_pflotran = "/home/scratch/smolins/amanzi/examples/phase2/chemistry/1d-calcite/pflotran-run/"

     # hardwired for 1d-calcite: time and comp
    # times = ['Time:  1.00000E+01 y','Time:  2.00000E+01 y','Time:  3.00000E+01 y','Time:  4.00000E+01 y','Time:  5.00000E+01 y']
    times = ['Time:  5.00000E+01 y']

    comp = 'Total_Ca++ [M]'
    Ca_pflotran = []
    for i, time in enumerate(times):
       x_pflotran, c_pflotran = GetXY_PFloTran(path_to_pflotran,root,time,comp)
       Ca_pflotran = Ca_pflotran + [c_pflotran]

    comp = 'pH'
    pH_pflotran = []
    for i, time in enumerate(times):
       x_pflotran, c_pflotran = GetXY_PFloTran(path_to_pflotran,root,time,comp)
       pH_pflotran = pH_pflotran + [c_pflotran]

    comp = 'Calcite_VF'
    VF_pflotran = []
    for i, time in enumerate(times):
       x_pflotran, c_pflotran = GetXY_PFloTran(path_to_pflotran,root,time,comp)
       VF_pflotran = VF_pflotran + [c_pflotran]

    # pflotran Operator Splitting
    path_to_pflotran = "pflotran/os"

     # hardwired for 1d-calcite: time and comp
    times = ['Time:  1.00000E+01 y','Time:  2.00000E+01 y','Time:  3.00000E+01 y','Time:  4.00000E+01 y','Time:  5.00000E+01 y']
    times = ['Time:  5.00000E+01 y']

    comp = 'Total_Ca++ [M]'
    Ca_pflotran_OS = []
    for i, time in enumerate(times):
       x_pflotran_OS, c_pflotran = GetXY_PFloTran(path_to_pflotran,root,time,comp)
       Ca_pflotran_OS = Ca_pflotran_OS + [c_pflotran]

    comp = 'pH'
    pH_pflotran_OS = []
    for i, time in enumerate(times):
       x_pflotran_OS, c_pflotran = GetXY_PFloTran(path_to_pflotran,root,time,comp)
       pH_pflotran_OS = pH_pflotran_OS + [c_pflotran]

    comp = 'Calcite_VF'
    VF_pflotran_OS = []
    for i, time in enumerate(times):
       x_pflotran_OS, c_pflotran = GetXY_PFloTran(path_to_pflotran,root,time,comp)
       VF_pflotran_OS = VF_pflotran_OS + [c_pflotran]

    # crunchflow GIMRT
    path_to_crunchflow = "crunchflow/gimrt"

     # hardwired for calcite_1d_CF.in: time and comp
#    times_CF = ['totcon1.out','totcon2.out','totcon3.out','totcon4.out','totcon5.out']
    times_CF = ['totcon5.out']
    comp = 2
    Ca_crunchflow = []
    ignore = 4
    for i, time in enumerate(times_CF):
       x_crunchflow, c_crunchflow = GetXY_CrunchFlow(path_to_crunchflow,root,time,comp,ignore)
       Ca_crunchflow = Ca_crunchflow + [c_crunchflow]

#    times_CF = ['pH1.out','pH2.out','pH3.out','pH4.out','pH5.out']
    times_CF = ['pH5.out']
    comp = 0
    pH_crunchflow = []
    ignore = 3
    for i, time in enumerate(times_CF):
       x_crunchflow, c_crunchflow = GetXY_CrunchFlow(path_to_crunchflow,root,time,comp,ignore)
       pH_crunchflow = pH_crunchflow + [c_crunchflow]

#    times_CF = ['volume1.out','volume2.out','volume3.out','volume4.out','volume5.out']
    times_CF = ['volume5.out']
    comp = 0
    VF_crunchflow = []
    ignore = 4
    for i, time in enumerate(times_CF):
       x_crunchflow, c_crunchflow = GetXY_CrunchFlow(path_to_crunchflow,root,time,comp,ignore)
       VF_crunchflow = VF_crunchflow + [c_crunchflow/100.0]

    # crunchflow OS3D
    path_to_crunchflow = "crunchflow/os3d"

     # hardwired for calcite_1d_CF.in: time and comp
#    times_CF = ['totcon1.out','totcon2.out','totcon3.out','totcon4.out','totcon5.out']
    times_CF = ['totcon5.out']
    comp = 2
    Ca_crunchOS3D = []
    ignore = 4
    for i, time in enumerate(times_CF):
       x_crunchflow, c_crunchflow = GetXY_CrunchFlow(path_to_crunchflow,root,time,comp,ignore)
       Ca_crunchOS3D = Ca_crunchOS3D + [c_crunchflow]

#    times_CF = ['pH1.out','pH2.out','pH3.out','pH4.out','pH5.out']
    times_CF = ['pH5.out']
    comp = 0
    pH_crunchOS3D = []
    ignore = 3
    for i, time in enumerate(times_CF):
       x_crunchflow, c_crunchflow = GetXY_CrunchFlow(path_to_crunchflow,root,time,comp,ignore)
       pH_crunchOS3D = pH_crunchOS3D + [c_crunchflow]

#    times_CF = ['volume1.out','volume2.out','volume3.out','volume4.out','volume5.out']
    times_CF = ['volume5.out']
    comp = 0
    VF_crunchOS3D = []
    ignore = 4
    for i, time in enumerate(times_CF):
       x_crunchflow, c_crunchflow = GetXY_CrunchFlow(path_to_crunchflow,root,time,comp,ignore)
       VF_crunchOS3D = VF_crunchOS3D + [c_crunchflow/100.0]

    CWD = os.getcwd()
    local_path = "" 

    # subplots
    fig, ax = plt.subplots(3,sharex=True,figsize=(8,10))
    
    try:
        # Amanzi native chemistry
        input_filename = os.path.join("amanzi-u-1d-calcite.xml")
        path_to_amanzi = "amanzi-native-output"
        run_amanzi_standard.run_amanzi(input_filename, 1, ["calcite.bgd"], path_to_amanzi)
        
        comp = 'total_component_concentration.cell.Ca++ conc'
        Ca_amanzi_native = []
        for i, time in enumerate(times):
           x_amanzi_native, c_amanzi_native = GetXY_Amanzi(path_to_amanzi,root,comp)
           Ca_amanzi_native = Ca_amanzi_native +[c_amanzi_native]

        comp = 'free_ion_species.cell.H+'
        pH_amanzi_native = []
        for i, time in enumerate(times):
           x_amanzi_native, c_amanzi_native = GetXY_Amanzi(path_to_amanzi,root,comp)
           pH_amanzi_native = pH_amanzi_native +[-np.log10(c_amanzi_native)]

        comp = 'mineral_volume_fractions.cell.Calcite vol frac'
        VF_amanzi_native = []
        for i, time in enumerate(times):
           x_amanzi_native, c_amanzi_native = GetXY_Amanzi(path_to_amanzi,root,comp)
           VF_amanzi_native = VF_amanzi_native +[c_amanzi_native]

        native = True

    except:
        native = False

        pass    

    try:
        # Amanzi-Alquimia
        input_filename = os.path.join("amanzi-u-1d-calcite-alq.xml")
        path_to_amanzi = "amanzi-alquimia-output"
        run_amanzi_standard.run_amanzi(input_filename, 1, ["1d-calcite-trim.in","calcite.dat"], path_to_amanzi)

        comp = 'total_component_concentration.cell.Ca++ conc'
        Ca_amanzi_alquimia = []
        for i, time in enumerate(times):
           x_amanzi_alquimia, c_amanzi_alquimia = GetXY_Amanzi(path_to_amanzi,root,comp)
           Ca_amanzi_alquimia = Ca_amanzi_alquimia +[c_amanzi_alquimia]

        comp = 'free_ion_species.cell.H+'
        pH_amanzi_alquimia = []
        for i, time in enumerate(times):
           x_amanzi_alquimia, c_amanzi_alquimia = GetXY_Amanzi(path_to_amanzi,root,comp)
           pH_amanzi_alquimia = pH_amanzi_alquimia +[-np.log10(c_amanzi_alquimia)]

        comp = 'mineral_volume_fractions.cell.Calcite vol frac'
        VF_amanzi_alquimia = []
        for i, time in enumerate(times):
           x_amanzi_alquimia, c_amanzi_alquimia = GetXY_Amanzi(path_to_amanzi,root,comp)
           VF_amanzi_alquimia = VF_amanzi_alquimia +[c_amanzi_alquimia]

        alq = True

    except:
        
        alq = False

    try:
        # Amanzi-Alquimia-Crunch
        input_filename = os.path.join("amanzi-u-1d-calcite-alq-crunch.xml")
        path_to_amanzi = "amanzi-alquimia-crunch-output"
        run_amanzi_standard.run_amanzi(input_filename, 1, ["1d-calcite-crunch.in","calcite.dbs"], path_to_amanzi)

        comp = 'total_component_concentration.cell.Ca++ conc'
        Ca_amanzi_alquimia_crunch = []
        for i, time in enumerate(times):
           x_amanzi_alquimia_crunch, c_amanzi_alquimia_crunch = GetXY_Amanzi(path_to_amanzi,root,comp)
           Ca_amanzi_alquimia_crunch = Ca_amanzi_alquimia_crunch +[c_amanzi_alquimia_crunch]

        comp = 'free_ion_species.cell.H+'
        pH_amanzi_alquimia_crunch = []
        for i, time in enumerate(times):
           x_amanzi_alquimia_crunch, c_amanzi_alquimia_crunch = GetXY_Amanzi(path_to_amanzi,root,comp)
           pH_amanzi_alquimia_crunch = pH_amanzi_alquimia_crunch +[-np.log10(c_amanzi_alquimia_crunch)]

        comp = 'mineral_volume_fractions.cell.Calcite vol frac'
        VF_amanzi_alquimia_crunch = []
        for i, time in enumerate(times):
           x_amanzi_alquimia_crunch, c_amanzi_alquimia_crunch = GetXY_Amanzi(path_to_amanzi,root,comp)
           VF_amanzi_alquimia_crunch = VF_amanzi_alquimia_crunch +[c_amanzi_alquimia_crunch]

        alq_crunch = True

    except:
        
        alq_crunch = False

    # amanziS data
    
    # +pflotran
    try:
        # import pdb; pdb.set_trace()
        input_filename = os.path.join("amanzi-s-1d-calcite-alq.xml")
        path_to_amanziS = "struct_amanzi-output-pflo"
        run_amanzi_standard.run_amanzi(input_filename, 1, [], path_to_amanziS)
        root_amanziS = "plt00501"
        compS = "Ca++_Aqueous_Concentration"
        x_amanziS, c_amanziS = GetXY_AmanziS(path_to_amanziS,root_amanziS,time,compS)
        struct = len(x_amanziS)
        compS = "H+_Free_Ion_Guess"
        x_amanziS, pH_amanziS = GetXY_AmanziS(path_to_amanziS,root_amanziS,time,compS)
        pH_amanziS =  -np.log10(pH_amanziS)
        compS = "Calcite_Volume_Fraction"
        x_amanziS, VF_amanziS = GetXY_AmanziS(path_to_amanziS,root_amanziS,time,compS)
    except:
        struct = 0

    # +crunchflow
    try:
        # import pdb; pdb.set_trace()
        input_filename = os.path.join("amanzi-s-1d-calcite-alq-crunch.xml")
        path_to_amanziS = "struct_amanzi-output-crunch"
        run_amanzi_chem.run_amanzi_chem(input_filename,run_path=path_to_amanziS,chemfiles=None)
        root_amanziS = "plt00501"
        compS = "Ca++_Aqueous_Concentration"
        x_amanziS_crunch, c_amanziS_crunch = GetXY_AmanziS(path_to_amanziS,root_amanziS,time,compS)
        struct_c = len(x_amanziS_crunch)
        compS = "H+_Free_Ion_Guess"
        x_amanziS_crunch, pH_amanziS_crunch = GetXY_AmanziS(path_to_amanziS,root_amanziS,time,compS)
        pH_amanziS_crunch = -np.log10(pH_amanziS_crunch)
        compS = "Calcite_Volume_Fraction"
        x_amanziS_crunch, VF_amanziS_crunch = GetXY_AmanziS(path_to_amanziS,root_amanziS,time,compS)
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
    ax[0].set_xlim((18,32))
  
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
