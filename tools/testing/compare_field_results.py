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

def find_output_data_path(input_xml):
    import xml.etree.ElementTree as ET
    tree = ET.parse(input_xml)
    root = tree.getroot()
    if root.tag != 'amanzi_input':
        raise RuntimeError, 'The given XML file is not an Amanzi input.'
    output = root.find('output')
    observations = output.find('observations')
    obs_file = observations.find('filename').text
    return obs_file

def compute_error_norms(x_output, c_output, x_ref, c_ref):
    errors = {}
    for region in x_ref.keys():
        diff = c_ref[region] - c_output[region]
        errors[region] = np.linalg.norm(diff)
    return errors

if __name__ == "__main__":

    # Check argument number.
    if len(sys.argv) != 5:
        print('%s: usage:'%sys.argv[0])
        print('%s observation input.xml reference.obs tolerance'%sys.argv[0])
        exit(0)

    # Get arguments.
    observation = sys.argv[1]
    input_xml = sys.argv[2]
    ref_file = sys.argv[3]
    tolerance = float(sys.argv[4])

    # Check arguments.
    if not os.path.exists(input_xml):
        print('%s: input file %s does not exist.'%(sys.argv[0], input_xml))
        exit(0)
    if not os.path.exists(ref_file):
        print('%s: reference file %s does not exist.'%(sys.argv[0], ref_file))
        exit(0)
    if tolerance <= 0.0:
        print('%s: tolerance must be positive.'%sys.argv[0])
        exit(0)

    # Root through the input file and find out where the output lives.
    obs_file = find_output_observation_path(input_xml)
    if not os.path.exists(obs_file):
        print('%s: output file %s does not exist.'%(sys.argv[0], obs_file))

    # Now fetch the data from the output and reference files.
    x_output, c_output = get_observation_data(obs_file, observation)
    x_ref, c_ref = get_observation_data(ref_file, observation)

    # Compute the error norms for each of the observations in their respective regions.
    errors = compute_error_norms(x_output, c_output, x_ref, c_ref)
    msg = ""

    # Report.
    total_errors = 0
    error_lines = []
    for region, error in errors.items():
        if (error > tolerance):
            total_errors += 1
            error_lines.append('  error %g > tolerance %g in region %s'%(error, tolerance, region))

    if total_errors == 0:
        sys.exit('Comparison Passed')
    else:
        sys.exit('\n'.join(['Comparison Failed for observation %s:'%observation] + error_lines))
