# plots calcium concentration along x at last time step 
# benchmark: compares to pflotran simulation results
# author: S.Molins - Sept. 2013

import os
import sys
import h5py
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import cm as cm

# break output arrays in chunks
def get_chunks(arr, chunk_size = 10):
    chunks  = [ arr[start:start+chunk_size] for start in range(0, len(arr), chunk_size)]
    return chunks

# ----------- AMANZI + ALQUIMIA -----------------------------------------------------------------

def GetXYZ_Amanzi(path,root,comp,nodesx=121,nodesy=101):

    # open amanzi concentration and mesh files
    dataname = os.path.join(path,root+"_data.h5")
    amanzi_file = h5py.File(dataname,'r')

    meshname = os.path.join(path,root+"_mesh.h5")
    amanzi_mesh = h5py.File(meshname,'r')

    # nodal x, y
    x = [r[0] for r in amanzi_mesh['0']['Mesh']['Nodes']]
    y = [r[1] for r in amanzi_mesh['0']['Mesh']['Nodes']]

    # element x
    x_ = get_chunks(x,chunk_size=nodesx)
    x_amanzi = np.array( [np.diff(xcoord[0:nodesx])/2+xcoord[0:nodesx-1] for xcoord in x_[0:-1]] )

    # element y
    y_ = get_chunks(y,chunk_size=nodesx)
    y_amanzi = np.diff(y_,axis=0)/2 + y_[0:-1]
    y_amanzi = np.array( [ycoord[0:-1] for ycoord in y_amanzi] )

    # element data for x, y -- not sure why this thing is transposed
    time = max(amanzi_file[comp].keys())
    v = [v[0] for v in amanzi_file[comp][time]]
    z_amanzi = np.array( get_chunks(v,chunk_size=nodesy-1) )
    z_amanzi = z_amanzi.T

    amanzi_file.close()
    amanzi_mesh.close()
    
    return (x_amanzi, y_amanzi, z_amanzi)

# Main -------------------------------------------------------------------------------------
if __name__ == "__main__":

    import run_amanzi_standard

    # root name for problem
    root = "non_grid_aligned"
    nodesx = 121
    nodesy = 101

    local_path = "" 

    # subplots
    plt.subplots(1,figsize=(11,8))
    
    # Amanzi + Native chemistry
    try:
        input_file = os.path.join("non_grid_aligned-u.xml")
        path_to_amanzi = "output-u"
        run_amanzi_standard.run_amanzi(input_file, 1, [input_file], path_to_amanzi)
        
        comp = 'mineral_volume_fractions.cell.Calcite'
        x_native, y_native, z_native = GetXYZ_Amanzi(path_to_amanzi,root,comp,nodesx=nodesx,nodesy=nodesy)

    except Exception:
        pass    
      
    extent = [0.0, 0.60, 0.0, 0.50]
    
    # plot adjustments
    #plt.subplots_adjust(left=0.0,bottom=0.15,right=0.99,top=0.90)
    plt.suptitle("Amanzi 2D Non-grid Aligned Flow and Transport",fontsize=20)
    plt.tick_params(axis='both', which='major', labelsize=20)

    plt.xlabel("X (m)",fontsize=20)
    plt.ylabel("Y (m)",fontsize=20)

    plt.imshow(z_native, vmin=z_native.min(), vmax=z_native.max(), origin='lower',extent=extent,cmap=cm.bwr)
    cbar = plt.colorbar(shrink=.8, format='%.2e')

    cbar.ax.tick_params(axis='both', which='both', labelsize=20)
    cbar.ax.set_ylabel('Calcite volume\nfraction [-]',fontsize=20,rotation=0)
    cbar.ax.yaxis.set_label_coords(4., 1.15)

    # plt.show()
    plt.savefig(local_path+"non_grid_aligned_2d.png",format="png")

    
