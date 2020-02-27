"""Plotting functions for met data"""

import h5py
import numpy as np
from matplotlib import pyplot as plt

names = ["air temperature [K]",
         "relative humidity [-]",
         "incoming shortwave radiation [W m^-2]",
         "incoming longwave radiation [W m^-2]",
         "wind speed [m s^-1]"]

precip = ["precipitation rain [m s^-1]",
          "precipitation snow [m SWE s^-1]"]

params = {'backend':'ps',
        'axes.labelsize':16,
        'text.fontsize':24,
        'text.usetex':False,
        'font.family':'serif',
        'font.size':24,
        'font.serif':['Times','Palatino','serif'],
        'ps.usedistiller':'xpdf',
        'lines.linewidth':2,
        'figure.subplot.left':0.07,
        'figure.subplot.right':0.97,
        'figure.subplot.bottom':0.07,
        'figure.subplot.top':0.97}

        



def plot_met(fname, end_time_in_years=np.inf, figsize=(10,10), style='b'):
    fax = plt.subplots(3,2,figsize=figsize)
    axs = fax[1].ravel()

    plt.rcParams.update(params)
    with h5py.File(fname, 'r') as fid:
        time = fid['time [s]'][:]/86400.0/365.0

        if end_time_in_years > time[-1]:
            end = len(time)
        else:
            end = np.where(time > end_time_in_years)[0][0]
        
        for i,n in enumerate(names):
            if n in fid.keys():
                axs[i].plot(time[0:end], fid[n][0:end], style)
            
            axs[i].set_title(n, fontsize=24)

        # for p, s2 in zip(precip, ['r','b']):
        #     if p in fid.keys():
        axs[i+1].plot(time[0:end], fid[precip[0]][0:end], 'r')
        axs[i+1].plot(time[0:end], fid[precip[1]][0:end], 'b')
        axs[i+1].set_title("precip (r=rain, b=snow) [m SWE s^-1]", fontsize=24)
        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
        
if __name__ == "__main__":
    import sys    
    plot_met(sys.argv[-1], figsize=(18,12))
    plt.savefig("met_data.png")
    plt.show()
        
