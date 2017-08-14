# rename: convert_modflow_mesh_to_exodus.py
import argparse, sys

# PyLaGriT must be installed for this script to work
# For more information, visit: https://lanl.github.io/LaGriT/gettingstarted.html
from pylagrit import PyLaGriT

# User-defined variables
# Edit for your use-case
PARAMETER_FILE="PARAMS" # Input file holding parameters
MESH_NAME="mo_modflow" # Internal LaGriT mesh name

def load_params(infile):
    '''
    Loads parameters from file into a parameter dictionary
    File must contain parameters in the format:
    PARAMETER=VALUE
    Comments can be added with # this character
    
    :param infile: Parameter file to load
    :type infile: string
    
    :returns dict: Parameter dictionary
    '''
    
    # Create empty variable/value dict
    dict = {}
    
    # Load parameters file
    try:
        f = open(infile,'r')
    except:
        print("ERROR: Could not open parameters file "+infile)
        print("Verify the filename is correct and the filename is located relative")
        print("to the directory this script was executed")
        exit(1)
    
    for line in f:
        # If not a comment line...
        if (line[0] != "#") and (line.strip() != ''):
            if "#" in line:
                nline = line.split("#")[0]
                
                try:
                    key = nline.split("=")[0].strip()
                    val = nline.split("=")[1].strip()
                except:
                    print("ERROR IN LINE: " + line.strip("\n"))
                    print("Parameter line must be in the format: PARAM=VALUE")
                    exit(1)
            else:
                try:
                    key = line.split("=")[0].strip()
                    val = line.split("=")[1].strip()
                except:
                    print("ERROR IN LINE: " + line.strip("\n"))
                    print("Parameter line must be in the format: PARAM=VALUE")
                    exit(1)
            
            # Populate the dictionary
            dict[key.lower()] = val
    
    # Close file and return dictionary
    f.close()
    return dict

def verify_params(dict):
    '''
    Verifies that the necessary parameters were loaded.
    If not, the user is prompted to enter them.
    Defaults are available if the user does not know.
    
    :param dict: Parameter dictionary
    :type dict: dictionary
    '''
    
    # List of parameter names to search for
    params = ['nrows','ncols','input_bnds','export_name','dx','dy','height','minx','miny']
    
    # Default values for select parameters
    default = {'dx':100,'dy':100,'minx':0,'miny':0,'export_name':'modflow_out.inp','input_bnds':'init_bnds.inf'}
    
    # Search dictionary and prompt for missing values
    for key in params:
        if (not key in dict) or (key in dict and dict[key]==''):
            prompt = raw_input("Enter the value for {}: ".format(key))
            
            if prompt != '':
                dict[key] = prompt
            else:
                if key in default:
                    print("Using default value for {} as {}".format(key,default[key]))
                    dict[key] = default[key]
                else:
                    print("ERROR: Unable to continue without {} value".format(key))
                    exit(1)
                    

# Main function
def main(argv=None):
    
    if argv == None:
        argv = sys.argv
    
    parser = argparse.ArgumentParser(description = "Convert Modflow to Amanzi")
    parser.add_argument("-p", "--pfile", help="name of parameter file", default=PARAMETER_FILE)
    parser.add_argument("-v", "--view", help="view in Paraview", action="store_true", default = False)
    args = parser.parse_args()
    
    PFILE = args.pfile
    VIEW_IN_PARAVIEW = args.view
    
    # Initialize LaGriT
    l = PyLaGriT()
    
    # Load parameters file & check for missing values
    params = load_params(PFILE)
    verify_params(params)
    
    # Assign parameters to variables
    input_bnds = params['input_bnds']
    export_name = params['export_name']
    nrows = int(params['nrows'])
    ncols = int(params['ncols'])
    llcorner = [float(params['minx']),float(params['miny'])]
    D_XY=[float(params['dx']),float(params['dy'])]
    
    # Edit later - initialize faceset variables
    ibnd = {"active": 1, "noflow": 0, "head": -1, "edge": -2, "halite": -3}
    imat = {"active": 1, "noflow": 2, "head":  3, "edge":  4, "halite":  3}
    
    '''
    ibnd_active = 1
    ibnd_noflow = 0
    ibnd_head   = -1
    ibnd_edge   = -2
    ibnd_halite = -3
    
    imat_active = 1
    imat_noflow = 2
    imat_head   = 3
    imat_halite = 3
    imat_edge   = 4
    '''
    
    fs_bottom = 1
    fs_top    = 2
    fs_east   = 3
    fs_north  = 4
    fs_west   = 5
    fs_south  = 6
    fs_halite = 7
    fs_head   = 7
    fs_noflow = 8
    
    rmmat_edge = 4
    maxmat = 1
    
    # Generate hexmesh with cell-centered materials and optional elevation
    hexmesh = l.read_modflow(input_bnds, nrows, ncols, DXY=D_XY, name=MESH_NAME)
    
    # Eltset stuff
    
    # Cycle through all IMAT/IBND values and create & set eltsets
    i = 0
    for ibnd_key in ibnd:
        for imat_key in imat:
            if ibnd_key == imat_key:
                i += 1
                hexmesh.eltset_attribute("mod_bnds", ibnd[key], boolstr='eq', name="e{}".format(i))
                hexmesh.setatt("itetclr",imat[key],stride=["eltset","get","e{}".format(i)])
    
    # Remove edge cells
    hexmesh.eltset_attribute("itetclr",rmmat_edge,boolstr='eq',name="eremove") # Find edge cells
    hexmesh.rmpoint_eltset("eremove") # Remove edge cells
    hexmesh.rmpoint_compress(filter_bool=False, resetpts_itp=True) # compress and reset
    
    # Temp facesets for internal side interfaces
    hexmesh.resetpts_itp()
    
    l.sendline("extract/surfmesh/1 0 0/mo_sides/{}/-all-".format(hexmesh.name))
    
    
    hexmesh.dump_avs2(export_name)
    
    if (VIEW_IN_PARAVIEW == True):
        hexmesh.paraview()

if __name__ == "__main__":
    main()