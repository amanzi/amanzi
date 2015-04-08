# plots cation concentration along x at last time step 
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
    c_amanzi_alquimia = np.array(amanzi_file[comp][time]).flatten()
    amanzi_file.close()
    amanzi_mesh.close()
    
    return (x_amanzi_alquimia, c_amanzi_alquimia)

def GetXY_AmanziS(path,root,time,comp):
    try:
        import fsnapshot
        fsnok = True
    except:
        fsnok = False

   # import pdb; pdb.set_trace()

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
    c_pflotran = np.array(pfdata[time][comp]).flatten()
#    c_pflotran = c_pflotran.flatten()
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
#-------------------

if __name__ == "__main__":

    import os
    import run_amanzi_chem
    import numpy as np

    # root name for problem
    root = "ion-exchange"

    # pflotran
    path_to_pflotran = "pflotran"

     # hardwired for 1d-calcite: time and comp
    times = ['Time:  0.00000E+00 y', 'Time:  5.00000E+01 y',]
    components = ['Na+','Ca++','Mg++','Cl-']

    # comp = 'Total_'+root.title()+' [M]'
    
    u_pflotran = [[[] for x in range(len(components))] for x in range(len(times))]
    for i, time in enumerate(times):
       for j, comp in enumerate(components):          
          x_pflotran, c_pflotran = GetXY_PFloTran(path_to_pflotran,root,time,'Total_'+comp+' [M]')
          u_pflotran[i][j] = c_pflotran
    
    v_pflotran = [[[] for x in range(len(components))] for x in range(len(times))]
    for i, time in enumerate(times):
       for j, comp in enumerate(components):          
          x_pflotran, c_pflotran = GetXY_PFloTran(path_to_pflotran,root,time,'Total_Sorbed_'+comp+' [mol_m^3]')
          v_pflotran[i][j] = c_pflotran

    # crunchflow
    path_to_crunchflow = "crunchflow"

     # hardwired for 1d-tritium-crunch.in: time and comp
    times_CF = ['totcon1.out','totcon5.out']
    components = [0,1,2,3]
    ignore = 4

    u_crunchflow = [[[] for x in range(len(components))] for x in range(len(times_CF))]
    for i, time in enumerate(times_CF):
       for j, comp in enumerate(components):          
          x_crunchflow, c_crunchflow = GetXY_CrunchFlow(path_to_crunchflow,root,time,comp,ignore)
          u_crunchflow[i][j] = c_crunchflow

    times_CF = ['totexchange1.out','totexchange5.out']
    components = [0,1,2,3]
    ignore = 4
    
    v_crunchflow = [[[] for x in range(len(components))] for x in range(len(times_CF))]
    for i, time in enumerate(times_CF):
       for j, comp in enumerate(components):
          x_crunchflow, c_crunchflow = GetXY_CrunchFlow(path_to_crunchflow,root,time,comp,ignore)
          v_crunchflow[i][j] = c_crunchflow

    CWD = os.getcwd()
    local_path = "" 
    
    # subplots
    fig, ax = plt.subplots(2,sharex=True,figsize=(8,8))
    
    
    # hardwired for 1d-exchange:  last time = '71'
    times = ['1','72']
    amanzi_components = ['total_component_concentration.cell.Na+ conc', \
                         'total_component_concentration.cell.Ca++ conc', \
                         'total_component_concentration.cell.Mg++ conc', \
                         'total_component_concentration.cell.Cl- conc']
    amanzi_sorbed     = ['total_sorbed.cell.0', \
                         'total_sorbed.cell.1', \
                         'total_sorbed.cell.2', \
                         'total_sorbed.cell.3']

    amanzi_compS      = ['Na+_Aqueous_Concentration', \
                         'Ca++_Aqueous_Concentration', \
                         'Mg++_Aqueous_Concentration', \
                         'Cl-_Aqueous_Concentration']
    amanzi_sorbS      = ['Na+_Sorbed_Concentration', \
                         'Ca++_Sorbed_Concentration', \
                         'Mg++_Sorbed_Concentration', \
                         'Cl-_Sorbed_Concentration']

    try:   
        # Amanzi native chemistry
        input_filename = os.path.join("amanzi-u-1d-"+root+".xml")
        path_to_amanzi = "amanzi-native-output"
        run_amanzi_chem.run_amanzi_chem("../"+input_filename,run_path=path_to_amanzi,chemfiles=[root+".bgd"])

        u_amanzi_native = [[[] for x in range(len(amanzi_components))] for x in range(len(times))]
        for i, time in enumerate(times):
           for j, comp in enumerate(amanzi_components):
              x_amanzi_native, c_amanzi_native = GetXY_Amanzi(path_to_amanzi,root,time,comp)
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
        run_amanzi_chem.run_amanzi_chem("../"+input_filename,run_path=path_to_amanzi,chemfiles=["1d-"+root+"-trim.in",root+".dat"])

        u_amanzi_alquimia = [[[] for x in range(len(amanzi_components))] for x in range(len(times))]
        for i, time in enumerate(times):
           for j, comp in enumerate(amanzi_components):
              x_amanzi_alquimia, c_amanzi_alquimia = GetXY_Amanzi(path_to_amanzi,root,time,comp)
              u_amanzi_alquimia[i][j] = c_amanzi_alquimia
              
        v_amanzi_alquimia = [[[] for x in range(len(amanzi_sorbed))] for x in range(len(times))]
        for i, time in enumerate(times):
           for j, comp in enumerate(amanzi_sorbed):
              x_amanzi_alquimia, c_amanzi_alquimia = GetXY_Amanzi(path_to_amanzi,root,time,comp)
              v_amanzi_alquimia[i][j] = c_amanzi_alquimia

        alq = True

    except: 

        alq = False
  
    try:

        # Amanzi-Alquimia-Crunch
        input_filename = os.path.join("amanzi-u-1d-"+root+"-alq-crunch.xml")
        path_to_amanzi = "amanzi-alquimia-crunch-output"
        run_amanzi_chem.run_amanzi_chem("../"+input_filename,run_path=path_to_amanzi,chemfiles=["1d-"+root+"-crunch.in",root+".dbs"])

        u_amanzi_alquimia_crunch = [[[] for x in range(len(amanzi_components))] for x in range(len(times))]
        for i, time in enumerate(times):
           for j, comp in enumerate(amanzi_components):
              x_amanzi_alquimia_crunch, c_amanzi_alquimia_crunch = GetXY_Amanzi(path_to_amanzi,root,time,comp)
              u_amanzi_alquimia_crunch[i][j] = c_amanzi_alquimia_crunch
              
        v_amanzi_alquimia_crunch = [[[] for x in range(len(amanzi_sorbed))] for x in range(len(times))]
        for i, time in enumerate(times):
           for j, comp in enumerate(amanzi_sorbed):
              x_amanzi_alquimia_crunch, c_amanzi_alquimia_crunch = GetXY_Amanzi(path_to_amanzi,root,time,comp)
              v_amanzi_alquimia_crunch[i][j] = c_amanzi_alquimia_crunch

        alq_crunch = True

    except: 

        alq_crunch = False

        # subplots
        #fig, ax = plt.subplots(2,sharex=True,figsize=(15,8))

    # amanziS data
    
    # +pflotran
    try:
        input_filename = os.path.join("amanzi-s-1d-ion-exchange-alq.xml")
        path_to_amanziS = "struct_amanzi-output-pflo"
        run_amanzi_chem.run_amanzi_chem(input_filename,run_path=path_to_amanziS,chemfiles=None)
        root_amanziS = "plt00051"
        #compS = "Na+_Aqueous_Concentration"
        #x_amanziS, c_amanziS = GetXY_AmanziS(path_to_amanziS,root_amanziS,time,compS)

        #import pdb; pdb.set_trace()
        u_amanziS = [[] for x in range(len(amanzi_compS))]
        for j, compS in enumerate(amanzi_compS):
           x_amanziS, c_amanziS = GetXY_AmanziS(path_to_amanziS,root_amanziS,time,compS)
           u_amanziS[j] = c_amanziS

        #compS = "Na+_Sorbed_Concentration"
        #x_amanziS, v_amanziS = GetXY_AmanziS(path_to_amanziS,root_amanziS,time,compS)

        v_amanziS = [[] for x in range(len(amanzi_sorbS))]
        for j, compS in enumerate(amanzi_sorbS):
           x_amanziS, c_amanziS = GetXY_AmanziS(path_to_amanziS,root_amanziS,time,compS)
           v_amanziS[j] = c_amanziS

        struct = len(x_amanziS)
    except:
        struct = 0

    # +crunch
    try:
        input_filename = os.path.join("amanzi-s-1d-ion-exchange-alq-crunch.xml")
        path_to_amanziS = "struct_amanzi-output-crunch"
        run_amanzi_chem.run_amanzi_chem(input_filename,run_path=path_to_amanziS,chemfiles=None)
        root_amanziS = "plt00051"
        #compS = "Na+_Aqueous_Concentration"
        #x_amanziS, c_amanziS = GetXY_AmanziS(path_to_amanziS,root_amanziS,time,compS)

        #import pdb; pdb.set_trace()
        u_amanziS_c = [[] for x in range(len(amanzi_compS))]
        for j, compS in enumerate(amanzi_compS):
           x_amanziS_c, c_amanziS = GetXY_AmanziS(path_to_amanziS,root_amanziS,time,compS)
           u_amanziS_c[j] = c_amanziS

        #compS = "Na+_Sorbed_Concentration"
        #x_amanziS, v_amanziS = GetXY_AmanziS(path_to_amanziS,root_amanziS,time,compS)

        v_amanziS_c = [[] for x in range(len(amanzi_sorbS))]
        for j, compS in enumerate(amanzi_sorbS):
           x_amanziS_c, c_amanziS = GetXY_AmanziS(path_to_amanziS,root_amanziS,time,compS)
           v_amanziS_c[j] = c_amanziS

        struct_c = len(x_amanziS_c)
    except:
        struct_c = 0


    components = ['Na+','Ca++','Mg++','Cl-']
    colors= ['r','b','m','g'] # components
    styles = ['-','*','+','x'] # codes
#    codes = ['Amanzi+Alquimia(PFloTran)','Amanzi+Alquimia(CrunchFlow)','Amanzi Native Chemistry','PFloTran'] + [None,]*12
    codes = ['Amanzi(2nd-Ord)+Alquimia(PFloTran)','Amanzi(2nd-Ord)+Alquimia(CrunchFlow)','CrunchFlow(OS3D)','PFloTran'] + [None,]*12

    # lines on axes
    # for i, time in enumerate(times):
    i = 1
    for j, comp in enumerate(components):
        if alq:
               ax[0].plot(x_amanzi_alquimia, u_amanzi_alquimia[i][j],color=colors[j],linestyle=styles[0],linewidth=2)
        if alq_crunch:
               ax[0].plot(x_amanzi_alquimia_crunch, u_amanzi_alquimia_crunch[i][j],color=colors[j],linestyle='None',marker=styles[1],linewidth=2)
        # ax[0].plot(x_amanzi_native, u_amanzi_native[i][j],color=colors[j],marker=styles[2],linestyle='None',linewidth=2,label=comp)
        ax[0].plot(x_crunchflow, u_crunchflow[i][j],color=colors[j],linestyle='None',marker=styles[2],linewidth=2,label=comp)
        ax[0].plot(x_pflotran, u_pflotran[i][j],color=colors[j],linestyle='None',marker=styles[3],linewidth=2)

    # for i, time in enumerate(times):
    i = 1
    for j, comp in enumerate(components):
        if alq:
               ax[1].plot(x_amanzi_alquimia, v_amanzi_alquimia[i][j],color=colors[j],linestyle=styles[0],linewidth=2,label=codes[j*len(styles)])
        if alq_crunch:
               ax[1].plot(x_amanzi_alquimia_crunch, v_amanzi_alquimia_crunch[i][j],color=colors[j],linestyle='None',marker=styles[1],linewidth=2,label=codes[j*len(styles)+1])
      #  ax[1].plot(x_amanzi_native, v_amanzi_native[i][j],color=colors[j],marker=styles[2],linestyle='None',linewidth=2,label=codes[j*len(styles)+2])
      #  import pdb; pdb.set_trace()
        ax[1].plot(x_crunchflow, v_crunchflow[i][j],color=colors[j],linestyle='None',marker=styles[2],linewidth=2,label=codes[j*len(styles)+2])
        ax[1].plot(x_pflotran, v_pflotran[i][j],color=colors[j],linestyle='None',marker=styles[3],linewidth=2,label=codes[j*len(styles)+3])

    if (struct>0):
        
        for j in range(len(amanzi_compS)):
            ax[0].plot(x_amanziS, u_amanziS[j],color=colors[j],linestyle='--',linewidth=2)      
        
        for j in range(len(amanzi_sorbS)-1):
            ax[1].plot(x_amanziS, v_amanziS[j],color=colors[j],linestyle='--',linewidth=2) 

        ax[1].plot(x_amanziS, v_amanziS[len(amanzi_sorbS)-1],color=colors[len(amanzi_sorbS)-1],linestyle='--',label='AmanziS+Alquimia(PFloTran)',linewidth=2)

    if (struct_c > 0):
        
        for j in range(len(amanzi_compS)):
            ax[0].plot(x_amanziS_c, u_amanziS_c[j],color=colors[j],linestyle='None',marker='o',markerfacecolor='None',linewidth=2)      
        
        for j in range(len(amanzi_sorbS)-1):
            ax[1].plot(x_amanziS_c, v_amanziS_c[j],color=colors[j],linestyle='None',marker='o',markerfacecolor='None',linewidth=2) 

        ax[1].plot(x_amanziS_c, v_amanziS_c[len(amanzi_sorbS)-1],color=colors[len(amanzi_sorbS)-1],linestyle='None',marker='o',label='AmanziS+Alquimia(Crunch)',linewidth=2,markerfacecolor='None')

    # axes
    ax[1].set_xlabel("Distance (m)",fontsize=15)
    ax[0].set_ylabel("Total Concentration [mol/L]",fontsize=15)
    ax[1].set_ylabel("Total Sorbed Concent. [mol/m3]",fontsize=15)

    # ax[0].set_xlim(47,70)
    # ax[0].set_ylim(0.06,0.10)

    # ax[1].set_xlim(47,70)
    # ax[1].set_ylim(80,170)

    # plot adjustments
    plt.subplots_adjust(left=0.15,bottom=0.15,right=0.99,top=0.90)
    ax[0].legend(loc='upper left',fontsize=15)
    ax[1].legend(loc='lower right',fontsize=15)
    plt.suptitle("Amanzi 1D "+root.title()+" Benchmark at 50 years",x=0.57,fontsize=20)
    plt.tick_params(axis='both', which='major', labelsize=15)

    #pyplot.show()
    plt.savefig(root+"_1d.png",format="png")
    #plt.close()

    # finally:
    #    pass 
