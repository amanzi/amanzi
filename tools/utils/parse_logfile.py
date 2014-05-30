#!/usr/bin/env python

"""Loads and plots timestep data for a given run.

Usage: parse_logfile.py out.log
"""

import numpy as np

def print_headers():
    print "cycle, time, dt, iteration count, wallclock avg (s)"

def parse_file(fid, wallclock=False):
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


if __name__ == "__main__":
    import sys
    from matplotlib import pyplot as plt
    fname = sys.argv[-1]
    with open(fname,'r') as fid:
        data = parse_file(fid)
    plt.subplot(121)
    plt.semilogy(data[0][:,1], data[0][:,2], 'b-x')
    plt.semilogy(data[1][:,1], data[1][:,2], 'bx')
    plt.xlabel("time [days]")
    plt.ylabel("dt [days]")
    plt.subplot(122)
    plt.semilogy(data[0][:,0], data[0][:,2], 'b-x')
    plt.semilogy(data[1][:,0], data[1][:,2], 'bx')
    plt.xlabel("cycles [-]")
    plt.ylabel("dt [days]")
    plt.show()


