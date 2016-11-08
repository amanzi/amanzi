"""Plots 2D data on quadrilaterals using matplotlib.
"""

import sys,os
import numpy as np
import matplotlib.collections
from matplotlib import pyplot as plt
import h5py
import mesh
import colors

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
    
