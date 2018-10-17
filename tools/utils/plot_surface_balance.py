#!/usr/bin/env python
"""
Plot met data from an ATS input h5 file using default names.

This is currently only useful on a single, 1D column.
"""

import os,sys
import h5py
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.cm
import parse_ats
import itertools

import colors

def plot_surface((keys,times,dat), ax, color, varname, style='-', version='dev-new',
                 negate=False, label=None, hackfactor=None, divide_by_cell_volume=False):
    if version == '0.86':
        if varname == 'pressure' or varname == 'temperature':
            prefix = "surface-"
        else:
            prefix = ""
    elif version == 'dev':
        prefix = "surface-"
    elif version == 'dev-new':
        if varname.startswith('snow'):
            varname = varname[5:]
            prefix = 'snow-'
        else:
            prefix = 'surface-'
    else:
        raise RuntimeError("Unrecognized version '%s'"%version)
    sd = parse_ats.getSurfaceData(keys, dat, prefix+varname)

    if divide_by_cell_volume:
        cv = parse_ats.getSurfaceData(keys, dat, 'surface-cell_volume')
        sd = sd / cv
    
    if negate:
        sd = -sd

    if sd[0] == 0.:
        sd[0] = np.nan

    if hackfactor is not None:
        sd = sd * hackfactor

    ax.plot(times, sd, style, color=color, label=label)
    ax.set_ylabel(varname)
    ax.set_xlabel("time [yr]")


def plot_surface_balance(ktd, ax, color, marker=None, version='dev-new', label=None, hackfactor=None):
    """Plot data provided as a (keys,times,dat) tuple on a set of axes."""
    if marker is not None:
        style = '-'+marker
        style2 = '--'+marker
    else:
        style = '-'
        style2 = '--'

    ax = ax.ravel()

    if ktd[2][1] is None:
        ktds = (ktd[0], ktd[1], ktd[2][0])
    else:
        ktds = (ktd[0], ktd[1], ktd[2][1])
    ktd = (ktd[0], ktd[1], ktd[2][0])

    plot_surface(ktd,ax[0], color, 'incoming_shortwave_radiation', style=style, version=version, label=label)
    try:
        plot_surface(ktd, ax[1], color, 'incoming_longwave_radiation', style=style, version=version)
    except KeyError:
        print "Cannot find longwave radiation (old runs may not have this).  Skipping..."

    if version == '0.86':
        plot_surface(ktd, ax[2], color, 'qE_lw_out', style=style, version=version, negate=True)
    else:
        plot_surface(ktd, ax[2], color, 'qE_lw_out', style=style, version=version)
    plot_surface(ktd, ax[3], color, 'qE_latent_heat', style=style, version=version)
    plot_surface(ktd, ax[4], color, 'qE_sensible_heat', style=style, version=version)

    if version == '0.86':
        plot_surface(ktd, ax[5], color, 'total_energy_source', style=style, hackfactor=hackfactor, divide_by_cell_volume=True)
    else:
        plot_surface(ktd, ax[5], color, 'total_energy_source', style=style, hackfactor=hackfactor)
        
    plot_surface(ktds, ax[6][0], color, 'snow_temperature', style=style, version=version)
    plot_surface(ktds, ax[6][1], colors.darken(color), 'snow_depth', style=style2, version=version)

    plot_surface(ktd, ax[7][0], color, 'temperature', style=style, version=version)
    plot_surface(ktd, ax[7][1], colors.darken(color), 'ponded_depth', style=style2, version=version)

    plot_surface(ktd, ax[8], color, 'unfrozen_fraction', style=style, version=version)


def get_version(dirname):
    if os.path.isfile(os.path.join(dirname, 'visdump_snow_data.h5')):
        version = 'dev-new'
    else:
        d = h5py.File(os.path.join(dirname, 'visdump_surface_data.h5'), 'r')
        if 'surface-precipitation_rain.cell.0' in d.keys():
            version = 'dev'
        else:
            version = '0.86'
        d.close()
    return version

    
def get_axs():
    f,axs = plt.subplots(3,3, figsize=(14.4,8.3))
    axs = axs.ravel()
    axs[6] = (axs[6], axs[6].twinx())
    axs[7] = (axs[7], axs[7].twinx())
    return f,axs
        
if __name__ == "__main__":
    import sys
    import argparse
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
    markers = reversed(['x', '+', 'o', '^', 'v', 's'][0:len(fnames)])
    for i,(fname,marker) in enumerate(zip(fnames, itertools.cycle(markers))):
        version = get_version(fname)
        print("Directory %s got version: %s"%(fname,version))
        hackfactor = None
        if fname == "run-transient4":
            hackfactor = 2.0
            
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
        if version == 'dev-new':
            ktds = parse_ats.readATS(fname, 'visdump_snow_data.h5')[2]
        else:
            ktds = None
        plot_surface_balance((ktd[0], ktd[1], (ktd[2],ktds)), axs, c, label=fname, hackfactor=hackfactor, version=version, marker=marker)
        if ktds is not None:
            ktds.close()
        ktd[2].close()

    plt.tight_layout()
    axs[0].legend(bbox_to_anchor=(0., 1., 3.6, .05), loc=3,
               ncol=len(fnames), mode="expand", borderaxespad=0.)
    plt.show()
