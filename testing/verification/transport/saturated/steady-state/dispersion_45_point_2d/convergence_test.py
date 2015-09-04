import os
import sys
import h5py
import numpy as np


def GetOBS_AmanziU(path,obs_name,comp):

    # open amanzi observations file
    dataname = os.path.join(path,obs_name)
    datafile = open(dataname,'r')

    # parse file to extra location and value
    # throw away 2 header lines
    datafile.readline()
    datafile.readline()
    # read the rest
    obs = []
    for line in datafile.readlines():
	vals = line.split(",")
	obs.append(vals)
    # get individual observations
    ob_names = []
    for ob in obs:
        if (ob[1] not in ob_names):
	    ob_names.append(ob[1])
    obs_coords = {}
    obs_values = {}
    for name in ob_names:
	tmp_coords = []
	tmp_values = []
	for ob in obs:
	    if (ob[1] == name):
		tmp_coords.append(float(ob[len(ob)-2]))
		tmp_values.append(float(ob[len(ob)-1]))
	obs_coords[name] = np.array(tmp_coords)
	obs_values[name] = np.array(tmp_values)

    datafile.close()

    return (obs_coords, obs_values)

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
        root_amanziS = "plot00001"
        compS = "Aqueous_Pressure"
        x_amanziS, c_amanziS = GetXY_AmanziS(path_to_amanziS,root_amanziS,compS)
        struct = len(x_amanziS)
    except:
        struct = 0
        
    try:
        comp = 'drawdown'
        path_to_amanziU = "."
        root_amanziU = 'obs5_2_1_r5.out'
        #x_amanziU, c_amanziU = GetOBS_AmanziU(path_to_amanziU,obs_name,comp)
        x_amanziU, c_amanziU = GetOBS_AmanziU(path_to_amanziU,root_amanziU,comp)
        unstruct = len(x_amanziU.keys()[0])
    except:
        unstruct = 0

    try:
        comp = 'drawdown'
        path_to_amanziU = path_to_golden
        root_amanziU = 'obs5_2_1_r5.out'
        x_amanziU_gold, c_amanziU_gold = GetOBS_AmanziU(path_to_amanziU,root_amanziU,comp)
        unstruct_gold = len(x_amanziU_gold.keys()[0])
    except:
        unstruct_gold = 0

    # Diff
    msg = ""
    if (unstruct and unstruct_gold):
	err = {}
	for ob in c_amanziU_gold.keys():
            diff = c_amanziU_gold[ob] - c_amanziU[ob]
            error = np.linalg.norm(diff)
	    err[ob] = error

    # Report
        tol = 1e-8
	total_error = 0
	for ob in err:
	    if (err[ob] > tol):
		total_error += 1

        if total_error > 0 :
	    msg = msg + "Comparison Failed"
	    for ob in err:
	        msg = msg +  "\n  error("+str(ob)+") = " + str(err[ob])
            sys.exit(msg)
        else:
	    msg = msg + "Comparison Passed"
	    for ob in err:
	        msg = msg +  "\n  error("+str(ob)+") = " + str(err[ob])
	    print msg
    else:
	msg = msg + "Comparison Failed"
	msg = msg + "\n  tests results or golden_output missing"
        sys.exit(msg)
