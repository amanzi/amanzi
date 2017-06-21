#!/usr/bin/env python

"""Loads and plots timestep data for a given run.

Usage: parse_logfile.py out.log
"""

import numpy as np

def print_headers():
    print "cycle, time, dt, iteration count, wallclock avg (s)"

def parse_file(fid, wallclock=False):
    """Reads a file, and returns a list of good and bad timesteps.

    Each are a list of 3-tuples: (step number, time of step, dt)
    """
    
    data = []
    faildata = []

    # find number of "cycles"
    total_cycles = 0
    for line in fid:
        if "Cycle =" in line:
            sline = line.split()
            cyc = int(sline[4][:-1])
            time = float(sline[8][:-1])
            dt = float(sline[12])

            if len(data) > 0 and data[-1][0] == cyc:
                faildata.append(data.pop())
            data.append([cyc,time,dt])

    return np.array(data), np.array(faildata)


def float_list_type(mystring):
    colors = []
    for f in mystring.strip("(").strip(")").strip("[").strip("]").split(","):
        try:
            colors.append(float(f))
        except:
            colors.append(f)
    return colors

if __name__ == "__main__":
    import sys
    from matplotlib import pyplot as plt

    import argparse
    parser = argparse.ArgumentParser(description="Plot timestep histories from an ATS run logfile.")
    parser.add_argument("LOG_FILES", nargs="+", type=str,
                        help="List of logfiles to parse.")
    parser.add_argument("--colors", "-c", type=float_list_type,
                        default=None,
                        help="List of color indices to use, of the form: --colors=[0,0.1,1], where the doubles are in the range (0,1) and are mapped to a color using the colormap.")
    parser.add_argument("--colormap", "-m", type=str,
                        default="jet",
                        help="Colormap used to pick a color.")
    args = parser.parse_args()

    import colors
    cm = colors.cm_mapper(0,1,args.colormap)
        
    fnames = args.LOG_FILES

    for i,fname in enumerate(fnames):
        with open(fname,'r') as fid:
            data = parse_file(fid)
            plt.subplot(121)
            if args.colors is None:
                if len(fnames) > 1:
                    c = cm(float(i)/(len(fnames)-1))
                else:
                    c = 'b'
            else:
                if type(args.colors[i]) is float:
                    c = cm(args.colors[i])
                else:
                    c = args.colors[i]
            plt.semilogy(data[0][:,1], data[0][:,2], '-x', color=c, label=fname)
            if data[1].shape != (0,):
                plt.semilogy(data[1][:,1], data[1][:,2], 'x', color=c)
            plt.xlabel("time [days]")
            plt.ylabel("dt [days]")
            plt.legend(loc='lower left')
            plt.subplot(122)
            plt.semilogy(data[0][:,0], data[0][:,2], '-x', color=c)
            if data[1].shape != (0,):
                plt.semilogy(data[1][:,0], data[1][:,2], 'x', color=c)
            plt.xlabel("cycles [-]")
            plt.ylabel("dt [days]")
    plt.show()


