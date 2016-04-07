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
#   c_pflotran = c_pflotran.flatten()
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
    import run_amanzi_standard
    import numpy as np

    # root name for problem
    root = "ion-exchange"

    # pflotran
    path_to_pflotran = "pflotran"

     # hardwired for 1d-calcite: time and comp
    times = ['Time:  5.00000E+01 y',]
    components = ['Na+','Ca++','Mg++','Cl-']

    # comp = 'Total_'+root.title()+' [M]'

    try:     
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

      pflotran = True
   
    except:
      pflotran = False

    # crunchflow
    path_to_crunchflow = "crunchflow"

     # hardwired for 1d-tritium-crunch.in: time and comp
    times_CF = ['totcon5.out']
    components = [0,1,2,3]
    ignore = 4

    try:
      u_crunchflow = [[[] for x in range(len(components))] for x in range(len(times_CF))]
      for i, time in enumerate(times_CF):
         for j, comp in enumerate(components):          
            x_crunchflow, c_crunchflow = GetXY_CrunchFlow(path_to_crunchflow,root,time,comp,ignore)
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
            x_crunchflow, c_crunchflow = GetXY_CrunchFlow(path_to_crunchflow,root,time,comp,ignore)
            v_crunchflow[i][j] = c_crunchflow

      crunch = True

    except:
      crunch = False

    # hardwired for 1d-exchange:  last time = '71'
    times = ['71']
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
        run_amanzi_standard.run_amanzi(input_filename, 1, [root+".bgd"], path_to_amanzi)

        u_amanzi_native = [[[] for x in range(len(amanzi_components))] for x in range(len(times))]
        for i, time in enumerate(times):
           for j, comp in enumerate(amanzi_components):
              x_amanzi_native, c_amanzi_native = GetXY_Amanzi(path_to_amanzi,root,comp)
              u_amanzi_native[i][j] = c_amanzi_native

        v_amanzi_native = [[[] for x in range(len(amanzi_sorbed))] for x in range(len(times))]
        for i, time in enumerate(times):
           for j, comp in enumerate(amanzi_sorbed):
              x_amanzi_native, c_amanzi_native = GetXY_Amanzi(path_to_amanzi,root,comp)
              v_amanzi_native[i][j] = c_amanzi_native

        native = True

    except:  
        native = False

    try:
        # Amanzi-Alquimia
        input_filename = os.path.join("amanzi-u-1d-"+root+"-alq.xml")
        path_to_amanzi = "amanzi-alquimia-output"
        run_amanzi_standard.run_amanzi(input_filename, 1, ["1d-"+root+"-trim.in",root+".dat"], path_to_amanzi)

        u_amanzi_alquimia = [[[] for x in range(len(amanzi_components))] for x in range(len(times))]
        for i, time in enumerate(times):
           for j, comp in enumerate(amanzi_components):
              x_amanzi_alquimia, c_amanzi_alquimia = GetXY_Amanzi(path_to_amanzi,root,comp)
              u_amanzi_alquimia[i][j] = c_amanzi_alquimia
              
        v_amanzi_alquimia = [[[] for x in range(len(amanzi_sorbed))] for x in range(len(times))]
        for i, time in enumerate(times):
           for j, comp in enumerate(amanzi_sorbed):
              x_amanzi_alquimia, c_amanzi_alquimia = GetXY_Amanzi(path_to_amanzi,root,comp)
              v_amanzi_alquimia[i][j] = c_amanzi_alquimia

        alq = True

    except: 
        alq = False
  
    try:

        # Amanzi-Alquimia-Crunch
        input_filename = os.path.join("amanzi-u-1d-"+root+"-alq-crunch.xml")
        path_to_amanzi = "amanzi-alquimia-crunch-output"
        run_amanzi_standard.run_amanzi(input_filename, 1, ["1d-"+root+"-crunch.in",root+".dbs"], path_to_amanzi)

        u_amanzi_alquimia_crunch = [[[] for x in range(len(amanzi_components))] for x in range(len(times))]
        for i, time in enumerate(times):
           for j, comp in enumerate(amanzi_components):
              x_amanzi_alquimia_crunch, c_amanzi_alquimia_crunch = GetXY_Amanzi(path_to_amanzi,root,comp)
              u_amanzi_alquimia_crunch[i][j] = c_amanzi_alquimia_crunch
              
        v_amanzi_alquimia_crunch = [[[] for x in range(len(amanzi_sorbed))] for x in range(len(times))]
        for i, time in enumerate(times):
           for j, comp in enumerate(amanzi_sorbed):
              x_amanzi_alquimia_crunch, c_amanzi_alquimia_crunch = GetXY_Amanzi(path_to_amanzi,root,comp)
              v_amanzi_alquimia_crunch[i][j] = c_amanzi_alquimia_crunch

        alq_crunch = True

    except: 
        alq_crunch = False


    # amanziS data
    
    # +pflotran
    try:
        input_filename = os.path.join("amanzi-s-1d-ion-exchange-alq.xml")
        path_to_amanziS = "struct_amanzi-output-pflo"
        run_amanzi_standard.run_amanzi(input_filename, 1, [], path_to_amanziS)
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
        run_amanzi_standard.run_amanzi(input_filename, 1, [], path_to_amanziS)
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

           ax[j].plot(x_amanzi_alquimia, u_amanzi_alquimia[i][j],'r-',label='AmanziU+Alq(PFT)')
           bx[j].plot(x_amanzi_alquimia, v_amanzi_alquimia[i][j],'r-')

    if alq_crunch:
      i = 0  # hardwired for time '71'

      for j, comp in enumerate(components):

           ax[j].plot(x_amanzi_alquimia_crunch, u_amanzi_alquimia_crunch[i][j],'r*',label='AmanziU+Alq(CF)')
           bx[j].plot(x_amanzi_alquimia_crunch, v_amanzi_alquimia_crunch[i][j],'r*')

    if native:
      i = 0  # hardwired for time '71'

      for j, comp in enumerate(components):
           ax[j].plot(x_amanzi_native, u_amanzi_native[i][j],'rx')
           bx[j].plot(x_amanzi_native, v_amanzi_native[i][j],'rx',label='AmanziU Native Chem')

    if crunch:
      i = 0  # hardwired for time 50 years

      for j, comp in enumerate(components):
           ax[j].plot(x_crunchflow, u_crunchflow[i][j],'m*',label='CrunchFlow')
           bx[j].plot(x_crunchflow, v_crunchflow[i][j],'m*')
    
    if pflotran:
      i = 0  # hardwired for time 50 years

      for j, comp in enumerate(components):
           ax[j].plot(x_pflotran, u_pflotran[i][j],'m-',label='PFloTran')
           bx[j].plot(x_pflotran, v_pflotran[i][j],'m-')
           ax[j].text(x_pflotran[10],u_pflotran[i][j][10],comp,fontsize=15,bbox=dict(facecolor='white', alpha=1.0))
           bx[j].text(x_pflotran[10],v_pflotran[i][j][10],comp,fontsize=15,bbox=dict(facecolor='white', alpha=1.0))

    if (struct>0):

        for j in range(len(amanzi_compS)):
            ax[j].plot(x_amanziS, u_amanziS[j],'g-')
            bx[j].plot(x_amanziS, v_amanziS[j],'g-',label='AmanziS+Alq(PFT)')

    if (struct_c > 0):
        
        for j in range(len(amanzi_compS)):
            ax[j].plot(x_amanziS_c, u_amanziS_c[j],'g*')      
            bx[j].plot(x_amanziS_c, v_amanziS_c[j],'g*',label='AmanziS+Alq(CF)')

    # axes
    ax[len(components)-1].set_xlabel("Distance (m)",fontsize=15)
    bx[len(components)-1].set_xlabel("Distance (m)",fontsize=15)

##    for i,comp in enumerate(components):
    for i,comp in enumerate(components):
      ax[i].set_ylabel("Total Concen. [mol/L]",fontsize=15)
      bx[i].set_ylabel("Total Sorb. [mol/m3]",fontsize=15)

    # plot adjustments
    ax[0].legend(fontsize=12,loc='lower right')
    bx[0].legend(fontsize=12,loc='lower right')

    plt.suptitle("Amanzi 1D "+root.title()+" Benchmark at 50 years",fontsize=20) #,x=0.57,fontsize=20)

    plt.tick_params(axis='both', which='major', labelsize=15)
  
    plt.tight_layout() #(pad=0.4, w_pad=0.5, h_pad=1.0)

    plt.subplots_adjust(left=0.10,bottom=0.15,right=0.90,top=0.95)

    #pyplot.show()
    plt.savefig(root+"_1d.png",format="png")
    #plt.close()

    #finally:
    #    pass 
