#!/usr/bin/env python

"""Loads and plots timestep data for a given run.

Usage: parse_logfile.py out.log

NOTE: This is deprecated -- please use plot_timestep_history.py instead.
"""
from __future__ import print_function
from __future__ import division

import numpy as np

def print_headers():
    print("cycle, time, dt, iteration count, wallclock avg (s)")

def split_header(line):
    """Cleans a line and splits off the VerboseObject header and remainder of the line"""
    line_comp = line.split("|")
    if len(line_comp) < 2:
        return "", line.strip()
    else:
        header = line_comp[0].strip()
        remainder = "|".join(line_comp[1:]).strip()
        return header, remainder

    
    
def parse_logfile2(fid):
    """Detailed reader"""
    steps = []
    inside_cycle = False
    
    for line in fid:
        header, line = split_header(line)
        if line.startswith("Cycle ="):
            inside_cycle = True
            sline = line.split()
            cyc = int(sline[2][:-1])
            time = float(sline[6][:-1])
            dt = float(sline[10])
            steps.append(dict(cycle=cyc, time=time, dt=dt, failed=False, nonlinear_solves=list()))

            nonlinear_itr = -1
            backtracking_total_count = 0
            backtracking_num_itrs = 0
            max_backtrack_count = 0
            error_hist = []

        elif inside_cycle:
            if line.startswith("Solve succeeded"):
                solver = dict(nonlinear_iterations=nonlinear_itr,
                              backtracking_total_count=backtracking_total_count,
                              backtracking_iterations=backtracking_num_itrs,
                              backtracking_max_backtracks=max_backtrack_count,
                              error_history=error_hist,
                              failed=False)
                steps[-1]['nonlinear_solves'].append(solver)

            elif line.startswith("Solve failed") or "Solution iterate is not admissible, FAIL" in line:
                solver = dict(nonlinear_iterations=nonlinear_itr,
                              backtracking_total_count=backtracking_total_count,
                              backtracking_iterations=backtracking_num_itrs,
                              backtracking_max_backtracks=max_backtrack_count,
                              error_history=error_hist,
                              failed=True)
                steps[-1]['nonlinear_solves'].append(solver)
                steps[-1]['failed'] = True

            elif "error(res)" in line and "L2" not in line:
                sline = line.split(":")
                nonlinear_itr = int(sline[0])
                error_hist.append(float(sline[-1].split()[-1]))

                if "backtrack" in sline[1]:
                    backtracking_total_count += 1
                    max_backtrack_count = max(max_backtrack_count, int(sline[1].split()[1]))
                    error_hist[-1] = float(sline[-1].split()[-1])
                if sline[1].strip() == "backtrack 1":
                    backtracking_num_itrs += 1
    return steps
                    
            
    
    
def parse_logfile(fid):
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
    """Convert string-form list of doubles into list of doubles."""
    colors = []
    for f in mystring.strip("(").strip(")").strip("[").strip("]").split(","):
        try:
            colors.append(float(f))
        except:
            colors.append(f)
    return colors


def get_axs():
    """Gets a figure and list of axes for plotting."""
    return plt.subplots(1,2)

def decorate_axs(axs):
    """Adds legends, labels, limits."""
    axs[0].set_xlabel("time [days]")
    axs[0].set_ylabel("dt [days]")
    axs[0].legend(loc='lower left')
    axs[1].set_xlabel("cycles [-]")
    axs[1].set_ylabel("dt [days]")
    return

def plot(data, axs, color, label, symbol='x'):
    """Plot the data."""
    axs[0].semilogy(data[0][:,1], data[0][:,2], '-'+symbol, color=color, label=label)
    if data[1].shape != (0,):
        axs[0].semilogy(data[1][:,1], data[1][:,2], symbol, color=color)
    axs[1].semilogy(data[0][:,0], data[0][:,2], '-'+symbol, color=color)
    if data[1].shape != (0,):
        axs[1].semilogy(data[1][:,0], data[1][:,2], symbol, color=color)

def write_to_file(data, fnamebase):
    """Writes the data to a file for future reading"""
    fname = fnamebase+".npz"
    np.savez(fname, good_timesteps=data[0], bad_timesteps=data[1])

def read_from_file(fname):
    """Reads a .npz file"""
    read = np.load(fname)
    return [read["good_timesteps"], read["bad_timesteps"]]
            
if __name__ == "__main__":
    import sys,os
    from matplotlib import pyplot as plt

    import argparse
    parser = argparse.ArgumentParser(description="Plot timestep histories from an ATS run logfile.  Store the results in a .npz file for future faster reading.")
    parser.add_argument("LOG_FILES", nargs="+", type=str,
                        help="List of logfiles to parse.")
    parser.add_argument("--colors", "-c", type=float_list_type,
                        default=None,
                        help="List of color indices to use, of the form: --colors=[0,0.1,1], where the doubles are in the range (0,1) and are mapped to a color using the colormap.")
    parser.add_argument("--colormap", "-m", type=str,
                        default="jet",
                        help="Colormap used to pick a color.")
    parser.add_argument("--overwrite", "-o", action="store_true",
                        help="Do not use any existing .npz file -- instead reload from the logfile and overwrite the .npz file.")
    args = parser.parse_args()

    import colors
    cm = colors.cm_mapper(0,1,args.colormap)
    fig, axs = get_axs()
        
    fnames = args.LOG_FILES
    for i,fname in enumerate(fnames):
        if fname.endswith(".npz"):
            data = read_from_file(fname)
        elif os.path.isfile(fname+".npz") and not args.overwrite:
            data = read_from_file(fname+".npz")
        else:
            with open(fname,'r') as fid:
                data = parse_logfile(fid)
            write_to_file(data, fname)
                
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
            
        plot(data, axs, c, fname)

    decorate_axs(axs)
    plt.show()
    sys.exit(0)


