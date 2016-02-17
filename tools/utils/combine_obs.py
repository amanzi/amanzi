#!/usr/bin/env python
# -------------------------------------------------------------
# file: combine_xmf.py
#
# Combines xmf files from restarted runs into one VisIt readable file.
#

import sys,os
import numpy as np


def loadFile(filename, directory):
    return np.loadtxt(os.path.join(directory, filename))

def combine(filename, directory_list):
    dat = [loadFile(filename, directory) for directory in directory_list]
    for i in range(len(dat)-1):
        first_included = np.where(dat[i][:,0] >= dat[i+1][0,0])[0]
        if len(first_included) > 0:
            dat[i] = dat[i][0:first_included[0],:]

    return np.concatenate(dat)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Combine *.VisIt.xmf files into a single VisIt readable file.")
    parser.add_argument("DIRECTORIES", nargs="+", type=str,
                        help="List of directories to combine.")

    parser.add_argument("--filename", "-f", type=str,
                        help="File name to be combined")

    args = parser.parse_args()
    if len(args.DIRECTORIES) < 1:
        parser.print_help()
        raise RuntimeError("Specify nonzero length list of directories.")

    # create the directory objects
    data = combine(args.filename, args.DIRECTORIES)

    # get the header
    with open(os.path.join(args.DIRECTORIES[0],args.filename), 'r') as fid:
        line = fid.readline()
        header = []
        while line.strip().startswith("#"):
            header.append(line)
            line = fid.readline()
    headerlines = "\n".join(header)

    # write the file
    np.savetxt(args.filename, data, header=headerlines)


