import os, os.path, sys
import h5py
import numpy as np
import glob

## Extract 1d array along dimension dim = 1,2,3 
def GetXY_AmanziU_1D(path,root,comp,dim):

    # open amanzi concentration and mesh files
    dataname = os.path.join(path,root+"_data.h5")
    amanzi_file = h5py.File(dataname,'r')
    meshname = os.path.join(path,root+"_mesh.h5")
    amanzi_mesh = h5py.File(meshname,'r')

    # extract cell coordinates
    if (dim == 1):
        y = np.array(amanzi_mesh['0']['Mesh']["Nodes"][0:len(amanzi_mesh['0']['Mesh']["Nodes"])//4,0])
    if (dim == 3):
        y = np.array(amanzi_mesh['0']['Mesh']["Nodes"][0::4,2])

    # shift array by dy/2
    x_amanziU = np.diff(y)/2 + y[0:-1]

    # extract concentration array
    alltimes = [int(n) for n in amanzi_file[comp].keys()]
    time = list(amanzi_file[comp].keys())[alltimes.index(max(alltimes))]

    c_amanziU = np.array(amanzi_file[comp][time]).flatten()
    amanzi_file.close()
    amanzi_mesh.close()

    return (x_amanziU, c_amanziU)


def GetXY_AmanziS_1D(path,root,comp,dim):
    try:
        import fsnapshot
        fsnok = True
    except:
        fsnok = False

    # plotfile = os.path.join(path,root)
    plotfile = max(glob.glob(path+"/"+root+'*'))

    if os.path.isdir(plotfile) & fsnok:
        (nx, ny, nz) = fsnapshot.fplotfile_get_size(plotfile)
        if (dim == 1):
            nn = nx
        if (dim == 2):
            nn = ny
        x = np.zeros( (nn), dtype=np.float64)
        y = np.zeros( (nn), dtype=np.float64)
        (y, x, npts, err) = fsnapshot.fplotfile_get_data_1d(plotfile, comp, y, x, dim)
    else:
        y = np.zeros( (0), dtype=np.float64)
        v = np.zeros( (0), dtype=np.float64)
    
    return (x, y)


def GetXY_PFloTran_1D(path,root,time,comp):

    # read pflotran data
    filename = os.path.join(path,root+".h5")
    pfdata = h5py.File(filename,'r')

    # extract coordinates
    y = np.array(pfdata['Coordinates']['X [m]'])
    x_pflotran = np.diff(y)/2+y[0:-1]

    # extract concentrations
    c_pflotran = np.array(pfdata[time][comp]).flatten()
    pfdata.close()

    return (x_pflotran, c_pflotran)


def GetXY_CrunchFlow_1D(path,root,cf_file,comp,ignore):

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


## Cut logically structured array in logical direction d.
## Node that mesh and data may be cut in diffrerent directions.
##
## Unstructured data are placed in a long array. If this array
## has an underlying structure, we can slice it using step. 
## The slice range is controlled by start and stop.
def GetXY_AmanziU_Nodes(path, root, start, stop, step, d):

    # open mesh file
    name = os.path.join(path,root+"_mesh.h5")
    amanzi_file = h5py.File(name,'r')

    # extract cell d-coordinates
    tmp = np.array(amanzi_file['0']['Mesh']["Nodes"][start:stop,d])
    xp = tmp[::step]
    xc = np.diff(xp) + xp[0:-1] 
    amanzi_file.close()

    return xc


def GetXY_AmanziU_Values(path, root, comp, start, stop, step):

    # open data file
    name = os.path.join(path,root+"_data.h5")
    amanzi_file = h5py.File(name,'r')

    # extract data array with maximum time stamp
    alltimes = [int(n) for n in amanzi_file[comp].keys()]
    time = list(amanzi_file[comp].keys())[alltimes.index(max(alltimes))]

    tmp = np.array(amanzi_file[comp][time][start:stop]).flatten()
    values = tmp[::step]
    amanzi_file.close()

    return values


def parse_input_xml(input_xml):
    """Parses the Amanzi input XML file and returns a tuple (type, output_file), 
where type is 'structured' or 'unstructured', and output_file is the full path 
of the most recent output file to be used in the comparison."""
    import xml.etree.ElementTree as ET
    tree = ET.parse(input_xml)
    root = tree.getroot()
    if root.tag != 'amanzi_input':
        raise RuntimeError('The given XML file is not a valid Amanzi input.')
    if 'type' not in root.attrib.keys():
        raise RuntimeError('Could not find a type (structured/unstructured) in the given Amanzi input.')

    # Find the simulation type.
    file_type = root.attrib['type'].lower()
    if file_type not in ['structured', 'unstructured']:
        raise RuntimeError('Invalid simulation type in given Amanzi input: %s' % file_type)

    # Now find the output prefix.
    output = root.find('output')
    vis = output.find('vis')
    vis_prefix = vis.find('base_filename').text

    # Search on the disk for the most recent vis file.
    class Walker:
        pass
    walker = Walker()
    walker.prefix = vis_prefix
    walker.vis_files = []

    def append_vis_file(walker, dir_name, files):
        walker.vis_files.extend([os.path.join(dir_name, f) for f in files if (walker.prefix in f and '_mesh' not in f)])

    os.path.walk(os.path.dirname(input_xml), append_vis_file, walker)
    walker.vis_files.sort()

    return file_type, walker.vis_files[-1]

def compute_error_norm(x_output, c_output, x_ref, c_ref):
    diff = c_ref - c_output
    return np.linalg.norm(diff)

def compare_structured_data(output_file, ref_file, field):
    out_dir = os.path.dirname(output_file)
    out_prefix = os.path.basename(output_file).replace('_data.h5', '')

    ref_dir = os.path.dirname(ref_file)
    ref_prefix = os.path.basename(ref_file).replace('_data.h5', '')

    # Now fetch the data from the output and reference files.
    x_output, c_output = GetXY_AmanziS(out_dir, out_prefix, field)
    x_ref, c_ref = GetXY_AmanziS(ref_dir, ref_prefix, field)

    # Compute the error norm for the field.
    error = compute_error_norm(x_output, c_output, x_ref, c_ref)

    # Report.
    if (error > tolerance):
        sys.exit('error %g > tolerance %g for field %s'%(error, tolerance, field))
    else:
        sys.exit('Comparison Passed')

def compare_unstructured_data(output_file, ref_file, field):
    out_dir = os.path.dirname(output_file)
    out_prefix = os.path.basename(output_file).replace('_data.h5', '')

    ref_dir = os.path.dirname(ref_file)
    ref_prefix = os.path.basename(ref_file).replace('_data.h5', '')

    # Now fetch the data from the output and reference files.
    x_output, c_output = GetXY_AmanziU(out_dir, out_prefix, field)
    x_ref, c_ref = GetXY_AmanziU(ref_dir, ref_prefix, field)

    # Compute the error norm for the field.
    error = compute_error_norm(x_output, c_output, x_ref, c_ref)

    # Report.
    if (error > tolerance):
        sys.exit('error %g > tolerance %g for field %s'%(error, tolerance, field))
    else:
        sys.exit('Comparison Passed')

if __name__ == "__main__":

    # Check argument number.
    if len(sys.argv) != 5:
        print('%s: usage:'%sys.argv[0])
        print('%s field input.xml reference.h5 tolerance'%sys.argv[0])
        exit(0)

    # Get arguments.
    field = sys.argv[1]
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
    sim_type, output_file = parse_input_xml(input_xml)
    if not os.path.exists(output_file):
        print('%s: output file %s does not exist.'%(sys.argv[0], output_file))

    if sim_type == 'structured':
        compare_structured_data(output_file, ref_file, field)
    else:
        compare_unstructured_data(output_file, ref_file, field)

