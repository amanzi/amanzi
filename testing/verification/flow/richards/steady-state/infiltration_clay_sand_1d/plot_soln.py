import os
import sys
import h5py
import numpy as np
import matplotlib
from matplotlib import pyplot as plt


def GetXY_AmanziU(path,root,time,comp):

    # open amanzi concentration and mesh files
    dataname = os.path.join(path,root+"_data.h5")
    amanzi_file = h5py.File(dataname,'r')
    meshname = os.path.join(path,root+"_mesh.h5")
    amanzi_mesh = h5py.File(meshname,'r')

    # extract cell coordinates (z = comp:2)
    y=np.array(amanzi_mesh['0']['Mesh']['Nodes'][0:len(amanzi_mesh['0']['Mesh']['Nodes']),2])

    # center of cell
    yy = np.array([y[4*i] for i in range(len(y)/4)])
    x_amanziU = yy[0:-1]+np.diff(yy)/2

    # extract concentration array
    c_amanziU = np.array(amanzi_file[comp][time])
    c_amanziU = c_amanziU.reshape(len(c_amanziU))

    amanzi_file.close()
    amanzi_mesh.close()
    
    return (x_amanziU, c_amanziU)

def GetXY_AmanziS(path,root,comp):
    try:
        import fsnapshot
        fsnok = True
    except:
        fsnok = False

    plotfile = os.path.join(path,root)

    if os.path.isdir(plotfile) & fsnok:
        (nx, ny, nz) = fsnapshot.fplotfile_get_size(plotfile)
        v = np.zeros( (nx,ny), dtype=np.float64)
        (v, err) = fsnapshot.fplotfile_get_data_2d(plotfile, comp, v)

        (xmin, xmax, ymin, ymax, zmin, zmax) = fsnapshot.fplotfile_get_limits(plotfile)
        dy = (ymax - ymin)/ny
        y = ymin + dy*0.5 + np.arange( (ny), dtype=np.float64 )*dy
        v = v[0]
        
    else:
        y = np.zeros( (0), dtype=np.float64)
        v = np.zeros( (0), dtype=np.float64)
    
    return (y, v)

if __name__ == "__main__":

    try:
        path_to_amanziS = "."
        root_amanziS = "case_2c_plot00001"
        compS = "Aqueous_Pressure"
        x_amanziS, c_amanziS = GetXY_AmanziS(path_to_amanziS,root_amanziS,compS)
        struct = len(x_amanziS)
    except:
        struct = 0
        
    # subplots
    fig, ax = plt.subplots() 
        
    try:
        #time = '385'
        time = '404'
        comp = 'pressure.cell.0'
        path_to_amanziU = "."
        root_amanziU = 'case_2c_plot'
        x_amanziU, c_amanziU = GetXY_AmanziU(path_to_amanziU,root_amanziU,time,comp)
        unstruct = len(x_amanziU)
    except:
        unstruct = 0

    try:
        time = '407'
        comp = 'pressure.cell.0'
        path_to_amanziU = "golden_output"
        root_amanziU = 'case_2c_plot'
        x_amanziU_gold, c_amanziU_gold = GetXY_AmanziU(path_to_amanziU,root_amanziU,time,comp)
        unstruct_gold = len(x_amanziU_gold)
    except:
        unstruct_gold = 0
    print "EIB>> unstruct      = ",unstruct
    print "EIB>> unstruct_gold = ",unstruct_gold
    print "EIB>> struct        = ",struct

    # Do plot
    if (unstruct>0):
        alq = ax.plot(x_amanziU, c_amanziU,'m-',label='AmanziU',linewidth=2)
    if (unstruct_gold>0):
        alqG = ax.plot(x_amanziU_gold, c_amanziU_gold,'m-',label='AmanziU-Gold',linewidth=2)
    if (struct>0):
        sam = ax.plot(x_amanziS, c_amanziS,'g-',label='AmanziS',linewidth=2)

    # Diff
    #sad = ax.plot(x_amanziS, c_amanziS-c_amanziU,'g-',label='AmanziS',linewidth=2) 

    # axes
    ax.set_xlabel("Distance (m)") #,fontsize=20)
 
    # plot adjustments
    plt.subplots_adjust(left=0.13,bottom=0.1,right=0.91,top=0.9)
    #plt.legend(loc='upper left',fontsize=13)
    plt.legend(loc='upper left') #,fontsize=13)
    plt.suptitle("Aqueous Pressure",x=0.5,fontsize=20)
    plt.tick_params(axis='both', which='major', labelsize=20)

    plt.show()
    # plt.savefig(root+"_1d.png",format="png")
    # plt.close()
