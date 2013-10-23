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
    c_amanzi_alquimia = np.array(amanzi_file[comp][time]).flatten()
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
    c_pflotran = np.array(pfdata[time][comp]).flatten()
#    c_pflotran = c_pflotran.flatten()
    pfdata.close()

    return (x_pflotran, c_pflotran)

if __name__ == "__main__":

    import os
    import run_amanzi_chem
    import numpy as np

    # root name for problem
    root = "surface-complexation"
    components = ['H+','Na+','NO3-','Zn++']

    # pflotran
    path_to_pflotran = "pflotran"

     # hardwired for 1d-calcite: time and comp
    times = ['Time:  0.00000E+00 y', 'Time:  5.00000E+01 y',]
    comppflo = ['pH','Total_Na+ [M]','Total_NO3- [M]','Total_Zn++ [M]']
    
    u_pflotran = [[[] for x in range(len(comppflo))] for x in range(len(times))]
    for i, time in enumerate(times):
       for j, comp in enumerate(comppflo):          
          x_pflotran, c_pflotran = GetXY_PFloTran(path_to_pflotran,root,time,comp)
          u_pflotran[i][j] = c_pflotran
    
    v_pflotran = [[[] for x in range(len(components))] for x in range(len(times))]
    for i, time in enumerate(times):
       for j, comp in enumerate(components):
          x_pflotran, c_pflotran = GetXY_PFloTran(path_to_pflotran,root,time,'Total_Sorbed_'+comp+' [mol_m^3]')
          v_pflotran[i][j] = c_pflotran

    CWD = os.getcwd()
    local_path = "" 

    # subplots
    fig, ax = plt.subplots(3,sharex=True,figsize=(15,12))
    bx =[None,]*3
    bx[0] = ax[0].twinx()
    bx[2] = ax[2].twinx()
    
    try:
        # hardwired for 1d-calcite: Tritium = component 0, last time = '71'
        times = ['0','71']
        amanzi_components = ['free_ion_species.cell.0', \
                             'total_component_concentration.cell.Component 1 conc', \
                             'total_component_concentration.cell.Component 2 conc', \
                             'total_component_concentration.cell.Component 3 conc']
        amanzi_sorbed     = ['total_sorbed.cell.0', \
                             'total_sorbed.cell.1', \
                             'total_sorbed.cell.2', \
                             'total_sorbed.cell.3']
 
        # Amanzi native chemistry
        input_filename = os.path.join("amanzi-u-1d-"+root+".xml")
        path_to_amanzi = "amanzi-native-output"
        run_amanzi_chem.run_amanzi_chem("../"+input_filename,run_path=path_to_amanzi,chemfiles=[root+".bgd"])

        u_amanzi_native = [[[] for x in range(len(amanzi_components))] for x in range(len(times))]
        for i, time in enumerate(times):
           for j, comp in enumerate(amanzi_components):
              x_amanzi_native, c_amanzi_native = GetXY_Amanzi(path_to_amanzi,root,time,comp)
              if j ==0:
                       u_amanzi_native[i][j] = -np.log10(c_amanzi_native)
              else:
                       u_amanzi_native[i][j] = c_amanzi_native

        v_amanzi_native = [[[] for x in range(len(amanzi_sorbed))] for x in range(len(times))]
        for i, time in enumerate(times):
           for j, comp in enumerate(amanzi_sorbed):
              x_amanzi_native, c_amanzi_native = GetXY_Amanzi(path_to_amanzi,root,time,comp)
              v_amanzi_native[i][j] = c_amanzi_native

    except:
        
        pass

    try:  
        # Amanzi-Alquimia
        input_filename = os.path.join("amanzi-u-1d-"+root+"-alq.xml")
        path_to_amanzi = "amanzi-alquimia-output"
        run_amanzi_chem.run_amanzi_chem("../"+input_filename,run_path=path_to_amanzi,chemfiles=["1d-"+root+".in",root+".dat"])

        u_amanzi_alquimia = [[[] for x in range(len(amanzi_components))] for x in range(len(times))]
        for i, time in enumerate(times):
           for j, comp in enumerate(amanzi_components):
              x_amanzi_alquimia, c_amanzi_alquimia = GetXY_Amanzi(path_to_amanzi,root,time,comp)
              if j == 0:
                       u_amanzi_alquimia[i][j] = -np.log10(c_amanzi_alquimia)
              else:         
                       u_amanzi_alquimia[i][j] = c_amanzi_alquimia
              
        v_amanzi_alquimia = [[[] for x in range(len(amanzi_sorbed))] for x in range(len(times))]
        for i, time in enumerate(times):
           for j, comp in enumerate(amanzi_sorbed):
              x_amanzi_alquimia, c_amanzi_alquimia = GetXY_Amanzi(path_to_amanzi,root,time,comp)
              v_amanzi_alquimia[i][j] = c_amanzi_alquimia

        alq = True

    except:

        alq = False


    colors= ['r','b','m','g'] # components
    styles = ['-','--','x'] # codes
    codes = ['Amanzi+Alquimia(PFloTran)','Amanzi Native Chemistry','PFloTran'] + [None,]*9

    # lines on axes
    # ax[0],b[0] ---> Aqueous concentrations
    # ax[1]      ---> pH
    # ax[2],b[2] ---> Sorbed concentrations

    # first
    # ax[0],b[0] ---> Aqueous concentrations
    # ax[1]      ---> pH
    # for i, time in enumerate(times):
    i = 0
    for j, comp in enumerate(components):
        if j == 0:
            if alq:
                   ax[1].plot(x_amanzi_alquimia, u_amanzi_alquimia[i][j],color=colors[j],linestyle=styles[0],linewidth=2)
            ax[1].plot(x_amanzi_native, u_amanzi_native[i][j],color=colors[j],linestyle=styles[1],linewidth=2,label='pH')
            ax[1].plot(x_pflotran, u_pflotran[i][j],color=colors[j],linestyle='None',marker=styles[2],linewidth=2)
        elif j == 3:
            if alq:
                   bx[0].plot(x_amanzi_alquimia, u_amanzi_alquimia[i][j],color=colors[j],linestyle=styles[0],linewidth=2)
            bx[0].plot(x_amanzi_native, u_amanzi_native[i][j],color=colors[j],linestyle=styles[1],linewidth=2,label=comp)
            bx[0].plot(x_pflotran, u_pflotran[i][j],color=colors[j],linestyle='None',marker=styles[2],linewidth=2)
        else:
            if alq:
                   ax[0].plot(x_amanzi_alquimia, u_amanzi_alquimia[i][j],color=colors[j],linestyle=styles[0],linewidth=2)
            ax[0].plot(x_amanzi_native, u_amanzi_native[i][j],color=colors[j],linestyle=styles[1],linewidth=2,label=comp)
            ax[0].plot(x_pflotran, u_pflotran[i][j],color=colors[j],linestyle='None',marker=styles[2],linewidth=2)

    # second
    # ax[2],b[2] ---> Sorbed concentrations
    # for i, time in enumerate(times):
    i = 0
    for j, comp in enumerate(components):

        if j == 3:
            if alq:
                   bx[2].plot(x_amanzi_alquimia, v_amanzi_alquimia[i][j],color=colors[j],linestyle=styles[0],linewidth=2,label=codes[j*len(styles)])
            bx[2].plot(x_amanzi_native, v_amanzi_native[i][j],color=colors[j],linestyle=styles[1],linewidth=2,label=codes[j*len(styles)+1])
            bx[2].plot(x_pflotran, v_pflotran[i][j],color=colors[j],linestyle='None',marker=styles[2],linewidth=2,label=codes[j*len(styles)+2])
        else:
            if alq:
                   ax[2].plot(x_amanzi_alquimia, v_amanzi_alquimia[i][j],color=colors[j],linestyle=styles[0],linewidth=2,label=codes[j*len(styles)])
            ax[2].plot(x_amanzi_native, v_amanzi_native[i][j],color=colors[j],linestyle=styles[1],linewidth=2,label=codes[j*len(styles)+1])
            ax[2].plot(x_pflotran, v_pflotran[i][j],color=colors[j],linestyle='None',marker=styles[2],linewidth=2,label=codes[j*len(styles)+2])

    # axes
    ax[2].set_xlabel("Distance (m)",fontsize=15)
    ax[0].set_ylabel("Total Concentration [mol/L]",fontsize=15)
    bx[0].set_ylabel("Total Concentration [mol/L]",fontsize=15,color=colors[3])
    ax[1].set_ylabel("pH",fontsize=15)
    ax[2].set_ylabel("Total Sorbed Concent. [mol/m3]",fontsize=15)
    bx[2].set_ylabel("Total Sorbed Concent. [mol/m3]",fontsize=15,color=colors[3])

    # plot adjustments
    plt.subplots_adjust(left=0.10,bottom=0.15,right=0.90,top=0.90)
    ax[0].legend(loc='center right',fontsize=15)
    bx[0].legend(loc='center left',fontsize=15)
    ax[1].legend(loc='center right',fontsize=15)
    ax[2].legend(loc='center right',fontsize=15)
    ax[2].set_ylim(bottom=-0.01)
    bx[2].set_ylim(bottom=-2.0e-6)
    #bx[2].legend(loc='center',fontsize=15)
    plt.suptitle("Amanzi 1D "+root.title()+" Benchmark at 50 years",x=0.57,fontsize=20)
    plt.tick_params(axis='both', which='major', labelsize=15)

    #pyplot.show()
    plt.savefig(root+"_1d.png",format="png")
    #plt.close()

    #finally:
    #    pass 
