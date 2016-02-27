#!/usr/bin/env python

"""Loads, sorts, and writes column data from a 1D (potentially disordered) simulation into a H5 file.

Usage: column_data.py

Loads ATS simulation data in the current directory and writes a sorted
file for pressure and temperature named column_data.h5 which can then
be read by other simulations looking to interpolate a 1D initial
condition.
"""

import sys,os
import numpy as np
import h5py
import mesh

def meshX(coord=2, filename="visdump_mesh.h5", directory="."):
    return mesh.meshElemCentroids(filename,directory)[:,coord]


def fullname(varname):
    fullname = varname
    if not '.cell.' in fullname:
        fullname = fullname+'.cell.0'
    return fullname

def column_data(varnames, keys='all', directory=".", filename="visdump_data.h5",
                mesh_filename="visdump_mesh.h5", coord=2):
    """Returns data of shape ( len(varnames+1), len(keys), n_cells )"""
    z = meshX(coord, mesh_filename, directory)
    if type(varnames) is str:
        varnames = [varnames,]

    with h5py.File(os.path.join(directory,filename),'r') as dat:
        keys_avail = dat[fullname(varnames[0])].keys()
        keys_avail.sort(lambda a,b: int.__cmp__(int(a),int(b)))

        if keys == 'all':
            keys = keys_avail
        elif keys == '-1' or keys == -1:
            keys = [keys_avail[-1]]
        elif type(keys) is str:
            keys = [keys]
        elif type(keys) is int:
            keys = [str(keys)]
        elif type(keys) is slice:
            keys = keys_avail[keys]

        vals = np.zeros((len(varnames)+1, len(keys), len(z)), 'd')
        for i,key in enumerate(keys):
            vals[0,i,:] = z
            for j,varname in enumerate(varnames):
                vals[j+1,i,:] = dat[fullname(varname)][key][:,0]

    # sort in z coordinate
    return vals[:,:,vals[0,0,:].argsort()]



if __name__ == "__main__":
    if sys.argv[-1] == "column_data.py":
        z0 = 0.0
    else:
        z0 = float(sys.argv[-1])
    
    with h5py.File("column_data.h5", 'w') as fout:
        dat = column_data(["pressure", "temperature"], keys=-1)
        #dat = column_data(["pressure",], keys=-1)
        z_depth = z0 - dat[0,0,:] # correction to get surface set, depth coordinate
        fout.create_dataset("z", data=np.flipud(z_depth))
        fout.create_dataset("pressure", data=np.flipud(dat[1,0,:]))
        fout.create_dataset("temperature", data=np.flipud(dat[2,0,:]))
