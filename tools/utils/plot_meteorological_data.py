#!/usr/bin/env python
"""Plots meteorological data, using expected default names, from an ATS h5 input file.

Usage: plot_meteorological_data.py INPUT.h5
"""

import h5py
import numpy as np
from matplotlib import pyplot as plt

import colors


names1 = ["air temperature [K]",
         "relative humidity [-]",
         "incoming shortwave radiation [W m^-2]",
         "incoming longwave radiation [W m^-2]",
         "wind speed [m s^-1]"]

precip1 = ["precipitation rain [m s^-1]",
          "precipitation snow [m SWE s^-1]"]


names2 = ["Ta", "RH", "Qswin", "Qlwin", "Us"]
precip2 = ["Pr", "Ps"]

def plot_met(fname, axs, color='b', end_time_in_years=np.inf, style='-'):
    axs = axs.ravel()

    with h5py.File(fname, 'r') as fid:
        try:
            time = fid['time [s]'][:]/86400.0/365.0
        except KeyError:
            try:
                time = fid['time'][:]/86400.0/365.0
            except KeyError:
                raise KeyError('Missing time entry "time [s]"')
            else:
                names = names2
                precip = precip2
        else:
            names = names1
            precip = precip1

        time = np.squeeze(time)            

        if end_time_in_years > time[-1]:
            end = len(time)
        else:
            end = np.where(time > end_time_in_years)[0][0]
        
        for i,n in enumerate(names):
            if n in fid.keys():
                axs[i].plot(time[0:end], np.squeeze(fid[n][0:end]), style, color=color)
            axs[i].set_title(n)

        for p, s2 in zip(precip, ['-','--']):
            if p in fid.keys():
                axs[i+1].plot(time[0:end], np.squeeze(fid[p][0:end]), s2, color=color)
        axs[i+1].set_title("precip (solid=rain, dash=snow) [m SWE s^-1]")


def get_axs():
    return plt.subplots(3,2, figsize=(10,10))
        
if __name__ == "__main__":
    import sys

    import argparse
    parser = argparse.ArgumentParser(description="Plot met data from an ATS input h5 file using default names.")
    parser.add_argument("INPUT_FILES", nargs="+", type=str,
                        help="List of logfiles to parse.")
    parser.add_argument("--colors", "-c", type=colors.float_list_type,
                        default=None,
                        help="List of color indices to use, of the form: --colors=[0,0.1,1], where the doubles are in the range (0,1) and are mapped to a color using the colormap.")
    parser.add_argument("--colormap", "-m", type=str,
                        default="jet",
                        help="Colormap used to pick a color.")
    args = parser.parse_args()

    import colors
    cm = colors.cm_mapper(0,1,args.colormap)
    fig, axs = get_axs()

    fnames = args.INPUT_FILES
    for i,fname in enumerate(fnames):
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
            
        plot_met(fname, axs, c)

    plt.show()
        
