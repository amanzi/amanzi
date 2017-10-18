"""Plots 2D data on quadrilaterals using matplotlib.
"""

import sys,os
import numpy as np
import matplotlib.collections
from matplotlib import pyplot as plt
import h5py
import mesh
import colors

def fullname(varname):
    fullname = varname
    if not '.cell.' in fullname:
        fullname = fullname+'.cell.0'
    return fullname


def transect_data(varnames, keys='all', directory=".", filename="visdump_data.h5",
                mesh_filename="visdump_mesh.h5", coord_order=None, deformable=False):
    """Pulls simulation output into structured 2D arrays for transect-based, (i,j) indexing.

    Input:
      varnames       | A list of variable names to pull, e.g.
                     |  ['saturation_liquid', 'saturation_ice']
      keys           | Indices of timesteps to pull.  Either an int (i.e. 0, -1, etc) 
                     |  for the kth timestep, or a list of ints, or 'all'.
      directory      | Directory of the run.  Defaults to '.'
      filename       | Filename of the run.  Defaults to 'visdump_data.h5'
      mesh_filename  | Filename of the mesh.  Defaults to 'visdump_mesh.h5'
      coord_order    | Order of the transect coordinates.  Defaults to ['x','z'].  The 
                     |  mesh is sorted in this order.
      deformable     | Is the mesh deforming?
   
    Output:
      Output is an array of shape:
      ( len(varnames+2), len(keys), n_cells_coord_order[0], n_cells_coord_order[1] )
    
      data[0,0,:,:] is the coord_order[0] centroid
      data[1,0,:,:] is the coord_order[1] centroid
      data[i+2,k,:,:] is the ith varname data at the kth requested timestep, sorted in 
                      the same way as the centroids.

      Note that the data is re-ordered in INCREASING coordinate, i.e. bottom to top in z.

    Example usage:  
      Calculate and plot the thaw depth at step 5.

      // Pull saturation ice -- TD is where sat ice = 0."
      data = transect_data(['saturation_ice', 5)

      // x coordinate for plotting
      x = data[0,0,:,0]
    
      // for each column, find highest z where sat_ice > 0.
      td_i = np.array([np.where(data[2,0,i,:] > 0.)[0][-1] for i in range(data.shape[2])])

      // now that we have an index into the highest cell with ice, determine td as the 
      // mean of the highest cell with ice and the one above that.  Note this assumes
      // all columns have some thawing.
      td_z = np.array( [  (dat[1,0,i,td_i[i]] + dat[1,0,i,td_i[i+1]]) / 2. 
                              for i in range(len(td_i)) ] )

      plt.plot(x, td_z)

    """
    if coord_order is None:
        coord_order = ['x','z']

    if type(varnames) is str:
        varnames = [varnames,]

    # get centroids
    xyz = mesh.meshElemCentroids(mesh_filename, directory)

    # round to the nearest 0.1
    xyz = 0.1 * np.round(10*xyz)

    # get ordering of centroids
    dtype = [(coord_order[0], float), (coord_order[1], float)]
    num_order = []
    for i in coord_order:
        if i == 'x':
            num_order.append(0)
        elif i == 'y':
            num_order.append(1)
        elif i == 'z':
            num_order.append(2)

    xyz_sort_order = np.array([tuple([xyz[i,x] for x in num_order]) for i in range(len(xyz))], dtype=dtype)
    xyz_sorting = xyz_sort_order.argsort(order=coord_order)

    with h5py.File(os.path.join(directory,filename),'r') as dat:
        keys_avail = dat[fullname(varnames[0])].keys()
        keys_avail.sort(lambda a,b: int.__cmp__(int(a),int(b)))

        if keys == 'all':
            keys = keys_avail
        elif type(keys) is str:
            keys = [keys,]
        elif type(keys) is int:
            keys = [keys_avail[keys],]
        elif type(keys) is slice:
            keys = keys_avail[keys]
        elif type(keys) is list:
            if all(type(k) is int for k in keys):
                keys = [keys_avail[k] for k in keys]
            elif all(type(k) is str for k in keys):
                pass
            else:
                raise RuntimeError("Keys requested cannot be processed -- should be 'all', int, or str key, or list of ints or strs.")
                

        # get data
        vals = np.zeros((len(varnames)+2, len(keys), len(xyz)), 'd')

        for i,key in enumerate(keys):
            if deformable:
                xyz = mesh.meshElemCentroids(mesh_filename, directory)
            vals[0,i,:] = xyz[xyz_sorting,num_order[0]]
            vals[1,i,:] = xyz[xyz_sorting,num_order[1]]
            for j,varname in enumerate(varnames):
                vals[j+2,i,:] = dat[fullname(varname)][key][:,0][xyz_sorting]

    # reshape the data
    # determine nx
    nx = len(set(vals[0,0,:]))
    nz = vals.shape[2] / nx
    if (nx * nz != vals.shape[2]):
        raise RuntimeError("Assumption about first coordinate being cleanly binnable is falling apart -- ask Ethan to rethink this algorithm!")
    shp = vals.shape
    return vals.reshape(shp[0], shp[1], nx, nz)


def plot(dataset, ax, cax=None, vmin=None, vmax=None, cmap="jet",
         label=None,
    mesh_filename="visdump_mesh.h5", directory="."):
    """Draws a dataset on an ax."""
    if vmin is None:
        vmin = dataset.min()
    if vmax is None:
        vmax = dataset.max()

    # get the mesh and collapse to 2D
    etype, coords, conn = mesh.meshElemXYZ(filename=mesh_filename, directory=directory)
    if etype is not 'HEX':
        raise RuntimeError("Only works for Hexs")

    coords2 = np.array([[coords[i][0::2] for i in c if coords[i][1] == 0.0] for c in conn])
    assert coords2.shape[2] == 2
    assert coords2.shape[1] == 4

    polygons = matplotlib.collections.PolyCollection(coords2, edgecolor='k', cmap=cmap)
    polygons.set_array(dataset)
    polygons.set_clim(vmin,vmax)
    ax.add_collection(polygons)

    xmin = min(c[0] for c in coords.itervalues())
    xmax = max(c[0] for c in coords.itervalues())
    zmin = min(c[2] for c in coords.itervalues())
    zmax = max(c[2] for c in coords.itervalues())

    ax.set_xlim(xmin,xmax)
    ax.set_ylim(zmin,zmax)

    if cax is not None:
        cb = plt.colorbar(polygons, cax=cax)
        if label is not None:
            cb.set_label(label)
        
    return ((xmin,xmax),(zmin,zmax))
    
