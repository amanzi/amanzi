from matplotlib import pyplot as plt
import numpy as np

pc = np.linspace(0,7, 1000)
pc = 10**pc

def plot(wrm, ax=None, color='b', label=None):
    if ax is None:
        fig,ax = plt.subplots(1,1,squeeze=True)

    s = np.array([wrm.saturation(apc) for apc in pc])
    ax.semilogy(s, pc, color=color, label=label)
    ax.set_xlabel("saturation [-]")
    ax.set_ylabel("capillary pressure [Pa]")
    return ax
