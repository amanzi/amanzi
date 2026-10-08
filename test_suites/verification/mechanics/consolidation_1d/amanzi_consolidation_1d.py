import os
import sys
import h5py
import numpy as np
import matplotlib
from matplotlib import pyplot as plt

import run_amanzi_standard

if __name__ == "__main__":

    # Amanzi: unstructured
    try:
        input_file = "amanzi_consolidation_1d.xml"
        path_to_amanzi = "output"
        root_amanzi = 'plot'
    
        run_amanzi_standard.run_amanzi(input_file, 1, [input_file], path_to_amanzi)

        dataname = os.path.join(path_to_amanzi,root_amanzi+"_data.h5")
        amanzi_file = h5py.File(dataname,'r')
        meshname = os.path.join(path_to_amanzi,root_amanzi+"_mesh.h5")
        amanzi_mesh = h5py.File(meshname,'r')

        x = np.array(amanzi_mesh['0']['Mesh']["Nodes"][0:len(amanzi_mesh['0']['Mesh']["Nodes"]),0])
        x_amanziU = np.diff(x)/2 + x[0:-1]
        x_amanziU = x_amanziU[:40]

        comp = 'pressure'
        alltimes = [n for n in amanzi_file[comp].keys()]
        remove = []
        for i in alltimes:
            try:
                a = int(i)
            except:
                remove += [i]

        for i in remove:
            alltimes.remove(i)
        alltimes = [int(n) for n in alltimes]   
        time = list(amanzi_file[comp].keys())[alltimes.index(max(alltimes))]

        c_amanziU = np.array(amanzi_file[comp][time]).flatten()
        amanzi_file.close()
        amanzi_mesh.close()
 
        c_amanziU = c_amanziU[::6]
        unstruct = len(x_amanziU)
    except:
        unstruct = 0


    # Amanzi: analytic
    g = 10.0
    p0 = 1e+7
    rho = 1e+3
    mu = 1e-3
    k = 6e-18
    phi_cp = 1e-11

    E = 3e+10
    nu = 0.0

    t = 36000.0
    L = 10.0

    Ss = rho * g * (phi_cp + 1.0 / E)
    cv = rho * g / mu * k / Ss  # effective consolidation coefficient

    x_analytic = np.arange(0.0, L, 0.1)
    c_analytic = [0.0] * len(x_analytic)

    for m in range(0,100):
        k = 2 * m + 1
        factor = np.pi / 2 / L
        for n in range(len(x_analytic)):
            xc = x_analytic[n]
            c_analytic[n] = c_analytic[n] + (1.0 / k) * np.sin(k * factor * xc) * np.exp(-k * k * factor * factor * cv * t)

    for n in range(len(x_analytic)):
        c_analytic[n] = c_analytic[n] * (4 * p0 / np.pi)


    # plots
    fig = plt.figure()
    axes = fig.add_axes([.1,.1,.8,.8])
       
    if (unstruct>0):
        axes.plot(x_amanziU,c_amanziU,marker='o',color='g',label='AmanziU',linestyle='None')
    axes.plot(x_analytic,c_analytic,color='b',label='analytic')

    # axes
    axes.set_xlabel("Distance (m)") #, fontsize=20)
    axes.set_ylabel("Excess Pressure (Pa)") #, fontsize=20)
 
    # plot adjustments
    plt.legend(loc='upper left') #, fontsize=13)

    #plt.show()
    plt.savefig("consolidation_1d.png",format="png")
    # plt.close()
