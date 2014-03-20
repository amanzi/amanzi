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

def GetXY_Amanzi(path,root,time,comp):

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
    c_amanzi_alquimia = np.array(amanzi_file[comp][time])
    amanzi_file.close()
    amanzi_mesh.close()
    
    return (x_amanzi_alquimia, c_amanzi_alquimia)

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
    import run_amanzi_chem
    import numpy as np

    # root name for problem
    root = "calcite"

    # pflotran
    path_to_pflotran = "pflotran"
    # path_to_pflotran = "/home/scratch/smolins/amanzi/examples/phase2/chemistry/1d-calcite/pflotran-run/"

     # hardwired for 1d-calcite: time and comp
    times = ['Time:  1.00000E+01 y','Time:  2.00000E+01 y','Time:  3.00000E+01 y','Time:  4.00000E+01 y','Time:  5.00000E+01 y']

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
    times_CF = ['totcon1.out','totcon2.out','totcon3.out','totcon4.out','totcon5.out']
    comp = 2
    Ca_crunchflow = []
    ignore = 4
    for i, time in enumerate(times_CF):
       x_crunchflow, c_crunchflow = GetXY_CrunchFlow(path_to_crunchflow,root,time,comp,ignore)
       Ca_crunchflow = Ca_crunchflow + [c_crunchflow]

    times_CF = ['pH1.out','pH2.out','pH3.out','pH4.out','pH5.out']
    comp = 0
    pH_crunchflow = []
    ignore = 3
    for i, time in enumerate(times_CF):
       x_crunchflow, c_crunchflow = GetXY_CrunchFlow(path_to_crunchflow,root,time,comp,ignore)
       pH_crunchflow = pH_crunchflow + [c_crunchflow]

    times_CF = ['volume1.out','volume2.out','volume3.out','volume4.out','volume5.out']
    comp = 0
    VF_crunchflow = []
    ignore = 4
    for i, time in enumerate(times_CF):
       x_crunchflow, c_crunchflow = GetXY_CrunchFlow(path_to_crunchflow,root,time,comp,ignore)
       VF_crunchflow = VF_crunchflow + [c_crunchflow/100.0]

    # crunchflow OS3D
    path_to_crunchflow = "crunchflow/os3d"

     # hardwired for calcite_1d_CF.in: time and comp
    times_CF = ['totcon1.out','totcon2.out','totcon3.out','totcon4.out','totcon5.out']
    comp = 2
    Ca_crunchOS3D = []
    ignore = 4
    for i, time in enumerate(times_CF):
       x_crunchflow, c_crunchflow = GetXY_CrunchFlow(path_to_crunchflow,root,time,comp,ignore)
       Ca_crunchOS3D = Ca_crunchOS3D + [c_crunchflow]

    times_CF = ['pH1.out','pH2.out','pH3.out','pH4.out','pH5.out']
    comp = 0
    pH_crunchOS3D = []
    ignore = 3
    for i, time in enumerate(times_CF):
       x_crunchflow, c_crunchflow = GetXY_CrunchFlow(path_to_crunchflow,root,time,comp,ignore)
       pH_crunchOS3D = pH_crunchOS3D + [c_crunchflow]

    times_CF = ['volume1.out','volume2.out','volume3.out','volume4.out','volume5.out']
    comp = 0
    VF_crunchOS3D = []
    ignore = 4
    for i, time in enumerate(times_CF):
       x_crunchflow, c_crunchflow = GetXY_CrunchFlow(path_to_crunchflow,root,time,comp,ignore)
       VF_crunchOS3D = VF_crunchOS3D + [c_crunchflow/100.0]

    CWD = os.getcwd()
    local_path = "" 

    # local_path="/home/scratch/smolins/amanzi-fresh/demos/phase2/chemistry/1d-calcite/"
    # path_to_amanzi="/home/scratch/smolins/amanzi-alquimia/examples/phase2/chemistry/1d-calcite/"

    # subplots
    fig, ax = plt.subplots(3,sharex=True,figsize=(8,10))
    
    try:
        # hardwired for 1d-calcite: Ca = component 2, last time = '71'
        times = ['31','41','51','61','71']

        # Amanzi native chemistry
        input_filename = os.path.join("amanzi-u-1d-calcite.xml")
        path_to_amanzi = "amanzi-native-output"
        run_amanzi_chem.run_amanzi_chem("../"+input_filename,run_path=path_to_amanzi,chemfiles=["calcite.bgd"])
        
        comp = 'total_component_concentration.cell.Ca++ conc'
        Ca_amanzi_native = []
        for i, time in enumerate(times):
           x_amanzi_native, c_amanzi_native = GetXY_Amanzi(path_to_amanzi,root,time,comp)
           Ca_amanzi_native = Ca_amanzi_native +[c_amanzi_native]

        comp = 'free_ion_species.cell.0'
        pH_amanzi_native = []
        for i, time in enumerate(times):
           x_amanzi_native, c_amanzi_native = GetXY_Amanzi(path_to_amanzi,root,time,comp)
           pH_amanzi_native = pH_amanzi_native +[-np.log10(c_amanzi_native)]

        comp = 'mineral_volume_fractions.cell.Calcite vol frac'
        VF_amanzi_native = []
        for i, time in enumerate(times):
           x_amanzi_native, c_amanzi_native = GetXY_Amanzi(path_to_amanzi,root,time,comp)
           VF_amanzi_native = VF_amanzi_native +[c_amanzi_native]
    except:

        pass    

    try:
        # Amanzi-Alquimia
        input_filename = os.path.join("amanzi-u-1d-calcite-alq.xml")
        path_to_amanzi = "amanzi-alquimia-output"
        run_amanzi_chem.run_amanzi_chem("../"+input_filename,run_path=path_to_amanzi,chemfiles=["1d-calcite.in","calcite.dat"])

        comp = 'total_component_concentration.cell.Ca++ conc'
        Ca_amanzi_alquimia = []
        for i, time in enumerate(times):
           x_amanzi_alquimia, c_amanzi_alquimia = GetXY_Amanzi(path_to_amanzi,root,time,comp)
           Ca_amanzi_alquimia = Ca_amanzi_alquimia +[c_amanzi_alquimia]

        comp = 'free_ion_species.cell.0'
        pH_amanzi_alquimia = []
        for i, time in enumerate(times):
           x_amanzi_alquimia, c_amanzi_alquimia = GetXY_Amanzi(path_to_amanzi,root,time,comp)
           pH_amanzi_alquimia = pH_amanzi_alquimia +[-np.log10(c_amanzi_alquimia)]

        comp = 'mineral_volume_fractions.cell.Calcite vol frac'
        VF_amanzi_alquimia = []
        for i, time in enumerate(times):
           x_amanzi_alquimia, c_amanzi_alquimia = GetXY_Amanzi(path_to_amanzi,root,time,comp)
           VF_amanzi_alquimia = VF_amanzi_alquimia +[c_amanzi_alquimia]

        alq = True

    except:
        
        alq = False

        # subplots
       # fig, ax = plt.subplots(3,sharex=True,figsize=(8,10))
      
    for i, time in enumerate(times):

          # lines on axes
          if alq: 
                 ax[0].plot(x_amanzi_alquimia, Ca_amanzi_alquimia[i],'r-',linewidth=2)
          ax[0].plot(x_amanzi_native, Ca_amanzi_native[i],'r--')
          ax[0].plot(x_pflotran, Ca_pflotran[i],'bx',linewidth=2)
          ax[0].plot(x_pflotran_OS, Ca_pflotran_OS[i],'b*',linewidth=2)
          if i>1:
                  ax[0].plot(x_crunchflow, Ca_crunchflow[i-1],'y.')
                  ax[0].plot(x_crunchflow, Ca_crunchOS3D[i-1],'g.')

          if alq:
                 ax[1].plot(x_amanzi_alquimia, pH_amanzi_alquimia[i],'r-',linewidth=2)
          ax[1].plot(x_amanzi_native, pH_amanzi_native[i],'r--')
          ax[1].plot(x_pflotran, pH_pflotran[i],'bx',linewidth=2)
          ax[1].plot(x_pflotran_OS, pH_pflotran_OS[i],'b*',linewidth=2)
          if i>0:
                  ax[1].plot(x_crunchflow, pH_crunchflow[i-1],'y.')
                  ax[1].plot(x_crunchflow, pH_crunchOS3D[i-1],'g.')

          if i==0:
            if alq:
                   ax[2].plot(x_amanzi_alquimia, VF_amanzi_alquimia[i],'r-',label='Amanzi+Alquimia(PFloTran)',linewidth=2)
            ax[2].plot(x_amanzi_native, VF_amanzi_native[i],'r--',label='Amanzi Native Chemistry')
            ax[2].plot(x_pflotran, VF_pflotran[i],'bx',label='PFloTran',linewidth=2)
            ax[2].plot(x_pflotran_OS, VF_pflotran_OS[i],'b*',label='PFloTran OS',linewidth=2)
          else:
            if alq:
                   ax[2].plot(x_amanzi_alquimia, VF_amanzi_alquimia[i],'r-',linewidth=2)
            ax[2].plot(x_amanzi_native, VF_amanzi_native[i],'r--')
            ax[2].plot(x_pflotran, VF_pflotran[i],'bx',linewidth=2)
            ax[2].plot(x_pflotran_OS, VF_pflotran_OS[i],'b*',linewidth=2)
            if i==1:
                  ax[2].plot(x_crunchflow, VF_crunchflow[i-1],'y.',label='CrunchFlow GIMRT')
                  ax[2].plot(x_crunchflow, VF_crunchOS3D[i-1],'g.',label='CrunchFlow OS3D')
            else: 
                  ax[2].plot(x_crunchflow, VF_crunchflow[i-1],'y.')
                  ax[2].plot(x_crunchflow, VF_crunchOS3D[i-1],'g.')
          
    # axes
    ax[2].set_xlabel("Distance (m)",fontsize=20)
    ax[0].set_ylabel("Total Ca concentration [mol/L]",fontsize=15)
    ax[1].set_ylabel("pH",fontsize=20)
    ax[2].set_ylabel("Calcite volume fraction",fontsize=15)

    # plot adjustments
    plt.subplots_adjust(left=0.15,bottom=0.15,right=0.99,top=0.90)
    ax[2].legend(loc='lower right',fontsize=13)
    plt.suptitle("Amanzi 1D Calcite Benchmark",x=0.57,fontsize=20)
    plt.tick_params(axis='x', which='major', labelsize=20)

    #pyplot.show()
    plt.savefig(local_path+"calcite_1d.png",format="png")
    #plt.close()

    #finally:
    #    pass 
