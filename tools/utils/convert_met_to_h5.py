#!/usr/bin/env python

"""Converts text-based Met data to H5 files for use with ATS.

Usage: convert_met_to_h5.py MY_MET_DATA.txt

creates MY_MET_DATA.h5
"""

import numpy as np
import h5py


def convert(filename):


    # grab headers
    with open(filename,'r') as fid:
        header = fid.readline().strip().split()

    # grab data
    dtxt = np.loadtxt(filename, skiprows=2)

    # make new data
    dnew = np.zeros((dtxt.shape[0], dtxt.shape[1] - 3),'d')
    dnew[:,1:] = dtxt[:,4:]
    dnew[:,0] = dtxt[:,0] * 24 * 3600

    # make new headers
    headers_new = ['time']
    headers_new.extend(header[4:])

    for i, header in enumerate(headers_new):
        # convert to Kelvin
        if header == "Ta":
            dnew[:,i] = dnew[:,i] + 273.15

        # convert preciptiation into rates
        if header == "Ps" or header == "Pr":
            dnew[:,i] = dnew[:,i] / (24 * 3600)

    # turn into hdf5
    h5filename = filename[:-4]+".h5"
    with h5py.File(h5filename,'w') as dat:
        for i in range(len(headers_new)):
            dat[headers_new[i]] = dnew[:,i]

if __name__ == "__main__":
    import sys
    filename = sys.argv[-1]
    convert(filename)
