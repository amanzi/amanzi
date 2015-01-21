# plots cation concentration along x at last time step 
# benchmark: compares to pflotran simulation results
# author: S.Molins - Nov. 2013
# modified: E.I.Barker - May 2014
#           - added native chemistry using v2 input spec

import os
import sys
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt


# ----------- AMANZI + ALQUIMIA -----------------------------------------------------------------

def GetXY_Amanzi(path,root,time,comp):

    # open amanzi concentration and mesh files
    dataname = os.path.join(path,root+"_data.h5")
    amanzi_file = h5py.File(dataname,'r')
    meshname = os.path.join(path,root+"_mesh.h5")
    amanzi_mesh = h5py.File(meshname,'r')

    # extract cell coordinates
    y = np.array(amanzi_mesh['0']['Mesh']["Nodes"][0:len(amanzi_mesh['0']['Mesh']["Nodes"])/4,0])
    # y = np.array(amanzi_mesh['Mesh']["Nodes"][0:len(amanzi_mesh['Mesh']["Nodes"])/4,0]) # old style
    # center of cell
    x_amanzi_alquimia  = np.diff(y)/2+y[0:-1]

    # extract concentration array
    c_amanzi_alquimia = np.array(amanzi_file[comp][time]).flatten()
    amanzi_file.close()
    amanzi_mesh.close()

    return (x_amanzi_alquimia, c_amanzi_alquimia)

# ----------- PFLOTRAN STANDALONE ------------------------------------------------------------

def GetXY_PFloTran(path,root,time,comp):

    # read pflotran data
    filename = os.path.join(path,"1d-"+root+".h5")
    pfdata = h5py.File(filename,'r')

    # extract coordinates
    y = np.array(pfdata['Coordinates']['X [m]'])
    x_pflotran = np.diff(y)/2+y[0:-1]

    # extract concentrations
    c_pflotran = np.array(pfdata[time][comp]).flatten()
#    c_pflotran = c_pflotran.flatten()
    pfdata.close()

    return (x_pflotran, c_pflotran)

# ------------- CRUNCHFLOW ------------------------------------------------------------------
def GetXY_CrunchFlow(path,root,cf_file,comp,ignore):

    # read CrunchFlow data
    filename = os.path.join(path,cf_file)
    f = open(filename,'r')
    lines = f.readlines()
    f.close()

    # ignore couple of lines
    for i in range(ignore):
      lines.pop(0)

    # extract data x0, x1, ..., xN-1 per line, keep only two columns
    xv=[]
    yv=[] 
    for line in lines:
      xv = xv + [float(line.split()[0])]
      yv = yv + [float(line.split()[comp+1])]
    
    xv = np.array(xv)
    yv = np.array(yv)

    return (xv, yv)

if __name__ == "__main__":

    import os
    import run_amanzi_chem
    import numpy as np

    # root name for problem
    root = "isotherms"
    components = ['A']#,'B','C']
    compcrunch = ['A']

    # times
    timespfl = ['Time:  0.00000E+00 y', 'Time:  5.00000E+01 y',]
    timesama  = ['0','71']
    timesama2 = ['0','71']

    # pflotran output
    pflotran_totc_templ = "Total_{0} [M]"
    pflotran_totc = [pflotran_totc_templ.format(x) for x in components]

    pflotran_sorb_templ = "Total_Sorbed_{0} [mol_m^3]"
    pflotran_sorb = [pflotran_sorb_templ.format(x) for x in components]

    # amanzi output
    amanzi_totc_templ = "total_component_concentration.cell.{} conc" #Component {0} conc"
    amanzi_totc = [amanzi_totc_templ.format(x) for x in components] #range(len(components))]
    amanzi_totc_crunch = [amanzi_totc_templ.format(x) for x in compcrunch] #range(len(components))]

    amanzi_sorb_templ = "total_sorbed.cell.{0}"
    amanzi_sorb = [amanzi_sorb_templ.format(x) for x in range(len(components))]
    amanzi_sorb_crunch = [amanzi_sorb_templ.format(x) for x in range(len(compcrunch))]

    # pflotran data
    path_to_pflotran = "pflotran"
    
    u_pflotran = [[[] for x in range(len(pflotran_totc))] for x in range(len(timespfl))]
    for i, time in enumerate(timespfl):
       for j, comp in enumerate(pflotran_totc):
          x_pflotran, c_pflotran = GetXY_PFloTran(path_to_pflotran,root,time,comp)
          u_pflotran[i][j] = c_pflotran
    
    v_pflotran = [[[] for x in range(len(pflotran_sorb))] for x in range(len(timespfl))]
    for i, time in enumerate(timespfl):
       for j, sorb in enumerate(pflotran_sorb):
          x_pflotran, c_pflotran = GetXY_PFloTran(path_to_pflotran,root,time,sorb)
          v_pflotran[i][j] = c_pflotran

    CWD = os.getcwd()
    local_path = "" 

    path_to_crunchflow = "crunchflow"

     # hardwired for 1d-isotherms-crunch.in: time and comp
    times_CF = ['totcon5.out']
    comp = 0
    u_crunchflow = []
    ignore = 4
    for i, time in enumerate(times_CF):
       x_crunchflow, c_crunchflow = GetXY_CrunchFlow(path_to_crunchflow,root,time,comp,ignore)
       u_crunchflow = u_crunchflow + [c_crunchflow]
    crunch = True

#   crunchflow does not explicitly calculate sorbed concs in Kd approach   

    # subplots
    fig, ax = plt.subplots(2,sharex=True,figsize=(8,8))
    
    try:

        # Amanzi native chemistry
        input_filename = os.path.join("amanzi-u-1d-"+root+".xml")
        path_to_amanzi = "amanzi-native-output"
        run_amanzi_chem.run_amanzi_chem("../"+input_filename,run_path=path_to_amanzi,chemfiles=[root+".bgd"])

        u_amanzi_native = [[[] for x in range(len(amanzi_totc))] for x in range(len(timesama))]
        for i, time in enumerate(timesama):
           for j, comp in enumerate(amanzi_totc):
              x_amanzi_native, c_amanzi_native = GetXY_Amanzi(path_to_amanzi,root,time,comp)
              u_amanzi_native[i][j] = c_amanzi_native

        v_amanzi_native = [[[] for x in range(len(amanzi_sorb))] for x in range(len(timesama))]
        for i, time in enumerate(timesama):
           for j, comp in enumerate(amanzi_sorb):
              x_amanzi_native, c_amanzi_native = GetXY_Amanzi(path_to_amanzi,root,time,comp)
              v_amanzi_native[i][j] = c_amanzi_native

    except:
        
        pass

    try:  
        # Amanzi-Alquimia
        input_filename = os.path.join("amanzi-u-1d-"+root+"-alq.xml")
        path_to_amanzi = "amanzi-alquimia-output"
        run_amanzi_chem.run_amanzi_chem("../"+input_filename,run_path=path_to_amanzi,chemfiles=["1d-"+root+"-trim.in",root+".dat"])

        u_amanzi_alquimia = [[[] for x in range(len(amanzi_totc))] for x in range(len(timesama))]
        for i, time in enumerate(timesama):
           for j, comp in enumerate(amanzi_totc):
              x_amanzi_alquimia, c_amanzi_alquimia = GetXY_Amanzi(path_to_amanzi,root,time,comp)
              u_amanzi_alquimia[i][j] = c_amanzi_alquimia
              
        v_amanzi_alquimia = [[[] for x in range(len(amanzi_sorb))] for x in range(len(timesama))]
        for i, time in enumerate(timesama):
           for j, comp in enumerate(amanzi_sorb):
              x_amanzi_alquimia, c_amanzi_alquimia = GetXY_Amanzi(path_to_amanzi,root,time,comp)
              v_amanzi_alquimia[i][j] = c_amanzi_alquimia

        alq = True

    except:

        alq = False

    try:  
        # Amanzi-ISV2
        input_filename = os.path.join("amanzi-u-1d-"+root+"-isv2.xml")
        path_to_amanzi = "amanzi-native-isv2-output"
        run_amanzi_chem.run_amanzi_chem("../"+input_filename,run_path=path_to_amanzi,isV2=True)

        u_amanzi_native_isv2 = [[[] for x in range(len(amanzi_totc))] for x in range(len(timesama2))]
        for i, time in enumerate(timesama2):
           for j, comp in enumerate(amanzi_totc):
              x_amanzi_native_isv2, c_amanzi_native_isv2 = GetXY_Amanzi(path_to_amanzi,root,time,comp)
              u_amanzi_native_isv2[i][j] = c_amanzi_native_isv2
              
        v_amanzi_native_isv2 = [[[] for x in range(len(amanzi_sorb))] for x in range(len(timesama2))]
        for i, time in enumerate(timesama2):
           for j, comp in enumerate(amanzi_sorb):
              x_amanzi_native_isv2, c_amanzi_native_isv2 = GetXY_Amanzi(path_to_amanzi,root,time,comp)
              v_amanzi_native_isv2[i][j] = c_amanzi_native_isv2

        isv2 = True

    except:

        isv2 = False

    try:  
        # Amanzi-Alquimia-Crunch
        input_filename = os.path.join("amanzi-u-1d-"+root+"-alq-crunch.xml")
        path_to_amanzi = "amanzi-alquimia-crunch-output"
        run_amanzi_chem.run_amanzi_chem("../"+input_filename,run_path=path_to_amanzi,chemfiles=["1d-"+root+"-crunch.in",root+".dbs"])

        u_amanzi_alquimia_crunch = [[[] for x in range(len(amanzi_totc_crunch))] for x in range(len(timesama))]
        for i, time in enumerate(timesama):
           for j, comp in enumerate(amanzi_totc_crunch):
              x_amanzi_alquimia_crunch, c_amanzi_alquimia_crunch = GetXY_Amanzi(path_to_amanzi,root,time,comp)
              u_amanzi_alquimia_crunch[i][j] = c_amanzi_alquimia_crunch
              
        v_amanzi_alquimia_crunch = [[[] for x in range(len(amanzi_sorb_crunch))] for x in range(len(timesama))]
        for i, time in enumerate(timesama):
           for j, comp in enumerate(amanzi_sorb_crunch):
              x_amanzi_alquimia_crunch, c_amanzi_alquimia_crunch = GetXY_Amanzi(path_to_amanzi,root,time,comp)
              v_amanzi_alquimia_crunch[i][j] = c_amanzi_alquimia_crunch

        alqc = True

    except:

        alqc = False



    colors= ['r'] #,'b','m','g'] # components
    styles = ['-','v','o','x'] # codes
    codes = ['AmanziU (2nd-Ord.)+Alq(PFT)','AmanziU (2nd-Ord.) Native(v2)','Amanzi (2nd-Ord.) Native(v1.2)','PFloTran'] + [None,]*9

    # lines on axes
    # ax[0] ---> Aqueous concentrations
    # ax[1] ---> Sorbed concentrations

    # for i, time in enumerate(times):
    i = 1 # hardwired 50 years
    for j, comp in enumerate(components):
            if alq:
                   ax[0].plot(x_amanzi_alquimia, u_amanzi_alquimia[i][j],color=colors[j],linestyle=styles[0],linewidth=2)
            if isv2:
                   ax[0].plot(x_amanzi_native_isv2, u_amanzi_native_isv2[i][j],color=colors[j],marker=styles[1],markeredgecolor=colors[j],linestyle='None')
            ax[0].plot(x_amanzi_native, u_amanzi_native[i][j],markeredgecolor=colors[j],marker=styles[2],linewidth=2,label=comp,mfc='None',markeredgewidth=2,linestyle='None')
            ax[0].plot(x_pflotran, u_pflotran[i][j],color=colors[j],linestyle='None',marker=styles[3],linewidth=2)

    #import pdb; pdb.set_trace()
    if crunch:
             ax[0].plot(x_crunchflow, u_crunchflow[0],color=colors[0],linestyle='--',marker='|',linewidth=5, label='CrunchFlow')

    if alqc:
            comp=='A'
            ax[0].plot(x_amanzi_alquimia_crunch, u_amanzi_alquimia_crunch[i][0],color=colors[0],linestyle='None',marker='*',linewidth=2)

    # for i, time in enumerate(times):
    for j, comp in enumerate(components):
            if alq:
                   ax[1].plot(x_amanzi_alquimia, v_amanzi_alquimia[i][j],color=colors[j],linestyle=styles[0],linewidth=2,label=codes[j*len(styles)])
            if isv2:
                   ax[1].plot(x_amanzi_native_isv2, v_amanzi_native_isv2[i][j],color=colors[j],linestyle='None',marker=styles[1],markeredgecolor=colors[j],label=codes[j*len(styles)+1])
            ax[1].plot(x_amanzi_native, v_amanzi_native[i][j],markeredgecolor=colors[j],marker=styles[2],linewidth=2,label=codes[j*len(styles)+2],mfc='None',markeredgewidth=2,linestyle='None')
            ax[1].plot(x_pflotran, v_pflotran[i][j],color=colors[j],linestyle='None',marker=styles[3],linewidth=2,label=codes[j*len(styles)+3])

    if alqc:
            comp=='A'
            ax[1].plot(x_amanzi_alquimia_crunch, v_amanzi_alquimia_crunch[i][0],color=colors[0],linestyle='None',marker='*',linewidth=2,label='AmanziU (2nd-Ord.)+Alq(CF)')

    # axes
    ax[1].set_xlabel("Distance (m)",fontsize=15)
    ax[0].set_ylabel("Total Concentration [mol/L]",fontsize=15)
    ax[1].set_ylabel("Total Sorbed Concent. [mol/m3]",fontsize=15)

    # plot adjustments
    plt.subplots_adjust(left=0.15,bottom=0.15,right=0.90,top=0.90)
    ax[0].legend(loc='center right',fontsize=15)
    ax[1].legend(loc='upper right',fontsize=10)
    ax[0].set_xlim(left=30,right=70)
    ax[1].set_xlim(left=30,right=70)
    #bx[2].legend(loc='center',fontsize=15)
    plt.suptitle("Amanzi 1D "+root.title()+" Benchmark at 50 years",x=0.57,fontsize=20)
    plt.tick_params(axis='both', which='major', labelsize=15)

    #pyplot.show()
    plt.savefig(root+"_1d.png",format="png")
    #plt.close()

    #finally:
    #    pass 
