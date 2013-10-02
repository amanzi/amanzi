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
    times = ['Time:  0.00000E+00 y','Time:  1.00000E+01 y','Time:  2.00000E+01 y','Time:  3.00000E+01 y','Time:  4.00000E+01 y','Time:  5.00000E+01 y']

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

    CWD = os.getcwd()
    local_path = "" 

    # local_path="/home/scratch/smolins/amanzi-fresh/demos/phase2/chemistry/1d-calcite/"
    # path_to_amanzi="/home/scratch/smolins/amanzi-alquimia/examples/phase2/chemistry/1d-calcite/"
        
    try:
        # hardwired for 1d-calcite: Ca = component 2, last time = '71'
        times = ['0','31','41','51','61','71']

        # Amanzi native chemistry
        input_filename = os.path.join("amanzi-u-1d-calcite.xml")
        path_to_amanzi = "amanzi-native-output"
        run_amanzi_chem.run_amanzi_chem("../"+input_filename,run_path=path_to_amanzi,chemfiles=["calcite.bgd"])
        
        comp = 'total_component_concentration.cell.Component 2 conc'
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

        # Amanzi-Alquimia
        input_filename = os.path.join("amanzi-u-1d-calcite-alq.xml")
        path_to_amanzi = "amanzi-alquimia-output"
        run_amanzi_chem.run_amanzi_chem("../"+input_filename,run_path=path_to_amanzi,chemfiles=["1d-calcite.in","calcite.dat"])

        comp = 'total_component_concentration.cell.Component 2 conc'
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

        # subplots
        fig, ax = plt.subplots()

        alq = [0.0] * len(times)
        ama = [0.0] * len(times)
        pfl = [0.0] * len(times)
        for i, time in enumerate(times):

          # lines on axes
          alq[i] = ax.plot(x_amanzi_alquimia, Ca_amanzi_alquimia[i],'r-',label='Amanzi+Alquimia(PFloTran)',linewidth=2)
          ama[i] = ax.plot(x_amanzi_native, Ca_amanzi_native[i],'ro',label='Amanzi Native Chemistry')
          pfl[i] = ax.plot(x_pflotran, Ca_pflotran[i],'b-',label='PFloTran',linewidth=2)

        # axes
        ax.set_xlabel("Distance (m)",fontsize=20)
        ax.set_ylabel("Total Ca concentration [mol/L]",fontsize=20)

        # plot adjustments
        plt.subplots_adjust(left=0.20,bottom=0.15,right=0.95,top=0.90)
#        plt.legend(loc='upper left',fontsize=13)
        plt.suptitle("Amanzi 1D Calcite Benchmark",x=0.57,fontsize=20)
        plt.tick_params(axis='both', which='major', labelsize=20)

        #pyplot.show()
        plt.savefig(local_path+"calcite_1d.png",format="png")
        plt.close()

    finally:
        pass 
