# plots tracer concentration along x at last time step 
# benchmark: compares to pflotran simulation results
# author: S.Molins - Sept. 2013

import os
import sys
import h5py
import numpy as np
import matplotlib
# matplotlib.use('Agg')
from matplotlib import pyplot as plt


# ----------- AMANZI + ALQUIMIA -----------------------------------------------------------------

def GetXY_AmanziU(path,root,comp):
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

def GetXY_AmanziS(path,root,comp):
    try:
        import fsnapshot
        fsnok = True
    except:
        fsnok = False

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

if __name__ == "__main__":

    import os
    import run_amanzi_standard
    import numpy as np

    # root name for problem
    root = "tracer"

    # pflotran
    path_to_pflotran = "pflotran"

     # hardwired for 1d-calcite: time and comp
    time = 'Time:  5.00000E+01 y'
    comp = 'Total_'+root.title()+' [M]'

    x_pflotran, c_pflotran = GetXY_PFloTran(path_to_pflotran,root,time,comp)    
    
    # CrunchFlow: hardwired for calcite_1d_CF.in: time and comp
    times_CF = 'totcon5.out'
    comp = 0
    ignore = 4

    # crunchflow GIMRT
    path_to_crunchflow = "crunchflow/gimrt"
    x_crunchflow, c_crunchflow = GetXY_CrunchFlow(path_to_crunchflow,root,times_CF,comp,ignore)

    # crunchflow OS3D
    path_to_crunchflow = "crunchflow/os3d"
    x_crunchOS3D, c_crunchOS3D = GetXY_CrunchFlow(path_to_crunchflow,root,times_CF,comp,ignore)
    
    CWD = os.getcwd()
    local_path = "" 

    # amanziU - native

    try:
      comp = 'total_component_concentration.cell.tracer conc'
      input_filename = os.path.join("amanzi-u-1d-"+root+".xml")
      path_to_amanzi = "amanzi-u-native-output"
      run_amanzi_standard.run_amanzi(input_filename, 1, [], path_to_amanzi)
      x_amanzi_native, c_amanzi_native = GetXY_AmanziU(path_to_amanzi,root,comp)
      native = len(x_amanzi_native)
  
    except:
      native = 0

    # amanziU - alquimia pflotran

    try:
        comp = 'total_component_concentration.cell.tracer conc'
        input_filename = os.path.join("amanzi-u-1d-"+root+"-alq-pflo.xml")
        path_to_amanzi = "amanzi-u-alq-pflo-output"
        run_amanzi_standard.run_amanzi(input_filename, 1, ["1d-"+root+".in",root+".dat"], path_to_amanzi)
        x_amanzi_alquimia, c_amanzi_alquimia = GetXY_AmanziU(path_to_amanzi,root,comp)
        alq = len(x_amanzi_alquimia)

    except:
        alq = 0

    # amanziS
    
    try:
        input_filename = os.path.join("amanzi-s-1d-"+root+"-alq-pflo.xml")
        path_to_amanziS = "amanzi-s-alq-pflo-output"
        run_amanzi_standard.run_amanzi(input_filename, 1, ["1d-"+root+".in",root+".dat"], path_to_amanziS)
        root_amanziS = "plt00051"
        compS = "Tracer_water_Concentration"
        x_amanziS, c_amanziS = GetXY_AmanziS(path_to_amanziS,root_amanziS,compS)
        struct = len(x_amanziS)
    except:
        struct = 0

# plotting --------------------------------------------------------

# subplots
    fig, ax = plt.subplots()

# pflotran
    ax.plot(x_pflotran, c_pflotran,'m-',label='PFloTran',linewidth=2)

# crunchflow
    ax.plot(x_crunchflow, c_crunchflow,'m--',label='CrunchFlow GIMRT',linewidth=2)
    ax.plot(x_crunchOS3D, c_crunchOS3D,'m*',label='CrunchFlow OS3D',linewidth=2) 

# unstruct amanzi native
    if native>0:
        ax.plot(x_amanzi_native, c_amanzi_native,'rx',label='AmanziU+Native',linewidth=2)

# unstruct amanzi alquimia + pflotran
    if alq>0:
        ax.plot(x_amanzi_alquimia, c_amanzi_alquimia,'r-',label='AmanziU+Alq(PFT)',linewidth=2)

# struct amanzi alquimia + pflotran
    if (struct>0):
        sam = ax.plot(x_amanziS, c_amanziS,'g-',label='AmanziS+Alq(PFT)',linewidth=2)     

# figure look
    # axes
    ax.set_xlabel("Distance (m)",fontsize=20)
    ax.set_ylabel("Total "+root.title()+" concentration [mol/L]",fontsize=20)

    # plot adjustments
    plt.subplots_adjust(left=0.20,bottom=0.15,right=0.95,top=0.90)
    plt.legend(loc='upper right',fontsize=13)
    plt.suptitle("Amanzi 1D "+root.title()+" Benchmark at 50 years",x=0.57,fontsize=20)
    plt.tick_params(axis='both', which='major', labelsize=20)

    plt.savefig(root+"_1d.png",format="png")
#   plt.show()
#   plt.close()

