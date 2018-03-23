#!/usr/bin/env python
"""
Plot met data from an ATS input h5 file using default names.

This is currently only useful on a single, 1D column.
"""

import h5py
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.cm
import parse_ats
    

def plot_surface((keys,times,dat), ax, color, varname, style='-', v86=False, negate=False, label=None):
    if v86:
        prefix = ""
    else:
        prefix = "surface-"
    sd = parse_ats.getSurfaceData(keys, dat, prefix+varname)
    if negate:
        sd = -sd

    if sd[0] == 0.:
        sd[0] = np.nan

    ax.plot(times, sd, style, color=color, label=label)
    ax.set_ylabel(varname)
    ax.set_xlabel("time [yr]")


def plot_surface_balance(ktd, ax, color, v86=False, label=None):
    """Plot data provided as a (keys,times,dat) tuple on a set of axes."""
    ax = ax.ravel()
    
    plot_surface(ktd,ax[0], color, 'incoming_shortwave_radiation', v86=v86, label=label)
    try:
        plot_surface(ktd, ax[1], color, 'incoming_longwave_radiation', v86=v86)
    except KeyErrror:
        print "Cannot find longwave radiation (old runs may not have this).  Skipping..."

    plot_surface(ktd, ax[2], color, 'qE_lw_out', v86=v86)
    plot_surface(ktd, ax[3], color, 'qE_latent_heat', v86=v86)
    plot_surface(ktd, ax[4], color, 'qE_sensible_heat', v86=v86)
    plot_surface(ktd, ax[5], color, 'conducted_energy_source')
        
    plot_surface(ktd, ax[6][0], color, 'snow_temperature')
    plot_surface(ktd, ax[6][1], colors.darken(color), 'snow_depth', style='--', v86=v86)

    plot_surface(ktd, ax[7][0], color, 'temperature')
    plot_surface(ktd, ax[7][1], colors.darken(color), 'ponded_depth', style='--', v86=v86)

    plot_surface(ktd, ax[8], color, 'albedo', v86=v86)

def get_axs():
    f,axs = plt.subplots(3,3, figsize=(14.4,8.3))
    axs = axs.ravel()
    axs[6] = (axs[6], axs[6].twinx())
    axs[7] = (axs[7], axs[7].twinx())
    return f,axs
        
if __name__ == "__main__":
    import sys

    import argparse
    import colors
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("INPUT_DIRECTORIES", nargs="+", type=str,
                        help="List of directories to plot.")
    parser.add_argument("--colors", "-c", type=colors.float_list_type,
                        default=None,
                        help="List of color indices to use, of the form: --colors=[0,0.1,1], where the doubles are in the range (0,1) and are mapped to a color using the colormap.")
    parser.add_argument("--colormap", "-m", type=str,
                        default="jet",
                        help="Colormap used to pick a color.")
    args = parser.parse_args()

    cm = colors.cm_mapper(0,1,args.colormap)
    fig, axs = get_axs()

    fnames = args.INPUT_DIRECTORIES
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

        ktd = parse_ats.readATS(fname, "visdump_surface_data.h5")
        plot_surface_balance(ktd, axs, c, label=fname)
        ktd[2].close()

    plt.tight_layout()
    axs[0].legend(bbox_to_anchor=(0., 1., 3.6, .05), loc=3,
               ncol=len(fnames), mode="expand", borderaxespad=0.)
    plt.show()
