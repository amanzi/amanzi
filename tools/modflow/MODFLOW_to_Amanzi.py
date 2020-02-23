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
    except Exception:
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
                except Exception:
                    print("ERROR IN LINE: " + line.strip("\n"))
                    print("Parameter line must be in the format: PARAM=VALUE")
                    exit(1)
            else:
                try:
                    key = line.split("=")[0].strip()
                    val = line.split("=")[1].strip()
                except Exception:
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
    params = ['nrows','ncols','input_bnds','export_name','dx','dy','height','minx','miny',\
                'exo_file']
    
    # Default values for select parameters
    default = {'dx':100,'dy':100,'minx':0,'miny':0,'export_name':'modflow_out.inp','input_bnds':'init_bnds.inf'}
    
    # Search dictionary and prompt for missing values
    for key in params:
        if (key not in dict) or (key in dict and dict[key]==''):
            prompt = input("Enter the value for {}: ".format(key))
            
            if prompt != '':
                dict[key] = prompt
            else:
                if key in default:
                    print("Using default value for {} as {}".format(key,default[key]))
                    dict[key] = default[key]
                else:
                    print("ERROR: Unable to continue without {} value".format(key))
                    exit(1)
                    
def write_fs_file(l):
    commands = ["cmo / copy / mo_tmp / mo_surf",\
    "cmo / select / mo_tmp",\
    "eltset / e_keep / id_side / eq / SS_ID",\
    "eltset / e_delete / not / e_keep",\
    "rmpoint / element / eltset get e_delete",\
    "rmpoint / compress",\
    "cmo / DELATT / mo_tmp / id_side",\
    "dump / avs2 / FILENAME / mo_tmp / 0 0 0 2",\
    "cmo / delete / mo_tmp"]
    
    for command in commands:
        l.sendline(command)

# Main function
def main(argv=None):
    
    if argv is None:
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
    
    # IBND and IMAT materials properties
    ibnd = {"active": int(params['ibnd_active']), "noflow": int(params['ibnd_noflow']),\
            "head": int(params['ibnd_head']), "edge": int(params['ibnd_edge']), "halite": int(params['ibnd_halite'])}
    imat = {"active": int(params['imat_active']), "noflow": int(params['imat_noflow']),\
            "head": int(params['imat_head']), "edge": int(params['imat_edge']), "halite": int(params['imat_halite'])}

    # Cycle through these keys instead of dictionary to preserve sequence
    mat_keys = ['active','noflow','head','edge','halite']

    fs_bottom = int(params['fs_bottom'])
    fs_top = int(params['fs_top'])
    fs_east = int(params['fs_east'])
    fs_north = int(params['fs_north'])
    fs_west = int(params['fs_west'])
    fs_south = int(params['fs_south'])
    fs_halite = int(params['fs_halite'])
    fs_head = int(params['fs_head'])
    fs_noflow = int(params['fs_noflow'])
    
    rmmat_edge = int(params['rmmat_edge'])
    maxmat = int(params['maxmat'])
    
    EXO_FILE=params['exo_file']
    
    # Generate hexmesh with cell-centered materials and optional elevation
    hexmesh = l.read_modflow(input_bnds, nrows, ncols, DXY=D_XY, name=MESH_NAME)
    
    # Eltset stuff
    for i, key in enumerate(mat_keys):
        hexmesh.eltset_attribute("mod_bnds", ibnd[key], boolstr='eq', name="e{}".format(i+1))
        hexmesh.setatt("itetclr",imat[key],stride=["eltset","get","e{}".format(i+1)])
    
    l.sendline('cmo/select/{}'.format(hexmesh.name)) # tmp ---------------------------------------------------
    
    # Remove edge cells
    hexmesh.eltset_attribute("itetclr",rmmat_edge,boolstr='eq',name="eremove") # Find edge cells
    hexmesh.rmpoint_eltset("eremove") # Remove edge cells
    hexmesh.rmpoint_compress(filter_bool=False, resetpts_itp=True) # compress and reset
    
    # Temp facesets for internal side interfaces
    mo_sides = l.extract_surfmesh(name="mo_sides",cmo_in=hexmesh,stride=[1,0,0],append='-all-',reorder=False,resetpts_itp=False)
    l.sendline("cmo/select/mo_sides")
    
    # Create psets from attributes 1 & 2
    ptop = mo_sides.pset_attribute("pts_topbot",2,comparison='eq',stride=(1,0,0),name="ptop")
    pbot = mo_sides.pset_attribute("pts_topbot",1,comparison='eq',stride=(1,0,0),name="pbot")
    
    # Generate eltsets from psets
    etop = ptop.eltset(membership='exclusive',name='etop')
    ebot = pbot.eltset(membership='exclusive',name='ebot')
    
    # Create a single eltset from the above two & remove
    e_delete = mo_sides.eltset_bool([etop,ebot],boolstr='union',name='e_delete')
    mo_sides.rmpoint_eltset("e_delete",compress=True,resetpts_itp=False)
    
    # NOFLOW faces (west corner)
    mo_tmp = l.copy(mo_sides, name="mo_tmp")
    edel = mo_tmp.eltset_attribute("itetclr0",imat["noflow"],boolstr='ne',name='edel')
    mo_tmp.rmpoint_eltset("edel",compress=True,resetpts_itp=True)
    mo_fs1 = l.copy(mo_tmp,name="mo_fs1")
    mo_tmp.delete()
    
    # Halite faces (eastern curve)
    mo_tmp = l.copy(mo_sides,name="mo_tmp")
    l.sendline("cmo/select/mo_tmp")
    edel = mo_tmp.eltset_attribute("itetclr0",imat["halite"],boolstr='ne',name='edel')
    mo_tmp.rmpoint_eltset("edel",compress=True,resetpts_itp=True)
    mo_fs2 = l.copy(mo_tmp,name="mo_fs2")
    mo_tmp.delete()
    
    mo_sides.delete()
    
    l.sendline("cmo/select/{}".format(hexmesh.name))
    edel = hexmesh.eltset_attribute("itetclr",maxmat,boolstr='gt',name='edel')
    hexmesh.rmpoint_eltset("edel",compress=True,resetpts_itp=True)
    
    # Final facesets for cropped mesh
    mo_surf = l.extract_surfmesh(name="mo_surf",cmo_in=hexmesh,stride=[1,0,0],append='external',reorder=False,resetpts_itp=False)
    mo_surf.addatt("id_normal",vtype='vint',rank='scalar',length='nelements',interpolate='',persistence='',ioflag='',value='')
    mo_surf.addatt("id_tmp",vtype='vint',rank='scalar',length='nelements',interpolate='',persistence='',ioflag='',value='')
    l.sendline("cmo/select/mo_surf")
    mo_surf.settets(method='normal')
    mo_surf.copyatt('itetclr',attname_sink='id_normal',mo_src=mo_surf)
    
    # Tag top and bottom face sets
    ptop = mo_surf.pset_attribute("pts_topbot",2,comparison='eq',stride=(1,0,0),name="ptop")
    pbot = mo_surf.pset_attribute("pts_topbot",1,comparison='eq',stride=(1,0,0),name="pbot")
    etop = ptop.eltset(membership='exclusive',name='etop')
    ebot = pbot.eltset(membership='exclusive',name='ebot')
    mo_surf.setatt("itetclr",2,stride=["eltset","get","etop"])
    mo_surf.setatt("itetclr",1,stride=["eltset","get","ebot"])
    
    mo_fs1.setatt("itetclr",1)
    mo_surf.setatt("id_tmp",0)
    mo_surf.interpolate('map','id_tmp',mo_fs1,'itetclr',stride=[1,0,0])
    etmp = mo_surf.eltset_attribute('id_tmp',1,boolstr='eq',name='etmp')
    mo_surf.setatt('itetclr',fs_noflow,stride=['eltset','get','etmp'])
    etmp.delete()

    # Interface between material 1 and halite and head cells
    mo_fs2.setatt('itetclr',1)
    mo_surf.setatt('id_tmp',0)
    mo_surf.interpolate('map','id_tmp',mo_fs2,'itetclr',stride=[1,0,0])
    etmp = mo_surf.eltset_attribute('id_tmp',1,boolstr='eq',name='etmp')
    mo_surf.setatt('itetclr',fs_halite,stride=['eltset','get','etmp'])
    etmp.delete()

    mo_surf.addatt('id_side',vtype='VINT',rank='scalar',length='nelements',interpolate='',persistence='',ioflag='',value='')
    mo_surf.copyatt('itetclr',attname_sink='id_side',mo_src=mo_surf)
    
    # Check for continuous connected boundary of 1
    mo_chk = l.copy(mo_surf,name='mo_chk')
    l.boundary_components(style='element',reset=False)
    mo_chk.delete()
    
    # Make sure to remove all attributes except idelem1 and idface1
    bad_atts = ['id_normal','itetclr0','idnode0','idelem0','facecol','itetclr1','idface0','id_tmp','ncon50','nconbnd','icontab']
    for att in bad_atts:
        mo_surf.delatt(att)
    
    l.sendline("cmo/select/mo_surf")
    filenames = ['bottom','top','east','north','west','south','head','noflow']
    ss_ids = [1,2,3,4,5,6,7,8]
    
    # WRITE ALL FACESET FILES
    for i in range(0,len(ss_ids)):
        l.sendline('define/FILENAME/fs_{}.faceset'.format(filenames[i]))
        l.sendline('define/SS_ID/{}'.format(ss_ids[i]))
        write_fs_file(l)
    
    
    #hexmesh.dump_exo(EXO_FILE,facesets=['fs_bottom.faceset','fs_top.faceset','fs_east.faceset'])
    l.sendline('dump / exo / {} / {} / / / facesets &\n fs_bottom.faceset fs_top.faceset fs_east.faceset &\n fs_north.faceset fs_west.faceset fs_south.faceset &\n fs_head.faceset fs_noflow.faceset'.format(EXO_FILE, hexmesh.name))
    
    hexmesh.dump_avs2(export_name)
    
    if VIEW_IN_PARAVIEW:
        hexmesh.paraview()

if __name__ == "__main__":
    main()
