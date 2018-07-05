import os
import sys
import h5py
import numpy as np


def GetXY_AmanziU(path,root,comp):

    # open amanzi concentration and mesh files
    dataname = os.path.join(path,root+"_data.h5")
    amanzi_file = h5py.File(dataname,'r')
    meshname = os.path.join(path,root+"_mesh.h5")
    amanzi_mesh = h5py.File(meshname,'r')

    # extract cell coordinates (z = comp:2)
    y = np.array(amanzi_mesh['0']['Mesh']['Nodes'][0:len(amanzi_mesh['0']['Mesh']['Nodes']),2])

    # center of cell
    yy = np.array([y[4*i] for i in range(len(y)/4)])
    x_amanziU = yy[0:-1]+np.diff(yy)/2


    # determine 'time'
    times = amanzi_file[comp].keys()
    time = times[len(times)-1]

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

    path_to_golden = "golden_output"
    if len(sys.argv) > 1:
	path_to_golden = sys.argv[1]
        path_to_golden += "/golden_output"

    try:
        path_to_amanziS = "."
        root_amanziS = "plot00053"
        compS = "Aqueous_Pressure"
        x_amanziS, c_amanziS = GetXY_AmanziS(path_to_amanziS,root_amanziS,compS)
        struct = len(x_amanziS)
    except:
        struct = 0
        
    try:
        comp = 'hydraulic_head.cell.0'
        path_to_amanziU = "."
        root_amanziU = 'plot'
        #x_amanziU, c_amanziU = GetXY_AmanziU(path_to_amanziU,root_amanziU,time,comp)
        x_amanziU, c_amanziU = GetXY_AmanziU(path_to_amanziU,root_amanziU,comp)
        unstruct = len(x_amanziU)
    except:
        unstruct = 0

    try:
        comp = 'hydraulic_head.cell.0'
        path_to_amanziU = path_to_golden
        root_amanziU = 'plot'
        x_amanziU_gold, c_amanziU_gold = GetXY_AmanziU(path_to_amanziU,root_amanziU,comp)
        unstruct_gold = len(x_amanziU_gold)
    except:
        unstruct_gold = 0

    # Diff
    msg = ""
    if (unstruct and unstruct_gold):
        diff = c_amanziU_gold - c_amanziU
        error = np.linalg.norm(diff)

    # Report
        tol = 1e-8
        if error < tol:
	    msg = msg + "Comparison Passed"
            msg = msg + "\n  error = " + str(error)
	    print msg
        else:
	    msg = msg + "Comparison Failed"
            msg = msg + "\n  error = " + str(error)
	    sys.exit(msg)
    else:
	msg = msg + "Comparison Failed"
	msg = msg + "\n  tests results or golden_output missing"
        sys.exit(msg)

