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
from matplotlib import pyplot as plt


# ----------- AMANZI + ALQUIMIA -----------------------------------------------------------------

def GetXY_Amanzi(path,root,comp):

    # open amanzi concentration and mesh files
    dataname = os.path.join(path,root+"_data.h5")
    amanzi_file = h5py.File(dataname,'r')
    meshname = os.path.join(path,root+"_mesh.h5")
    amanzi_mesh = h5py.File(meshname,'r')

    # extract cell coordinates
    y = np.array(amanzi_mesh['0']['Mesh']["Nodes"][0:len(amanzi_mesh['0']['Mesh']["Nodes"])/4,0])
    # y = np.array(amanzi_mesh['Mesh']["Nodes"][0:len(amanzi_mesh['Mesh']["Nodes"])/4,0]) # old style
    # center of cell
    x_alquimia = np.diff(y)/2+y[0:-1]

    # extract concentration array
    time = max(amanzi_file[comp].keys())
    c_alquimia = np.array(amanzi_file[comp][time]).flatten()
    amanzi_file.close()
    amanzi_mesh.close()

    return (x_alquimia, c_alquimia)

# Amanzi Struct

def GetXY_AmanziS(path,root,time,comp):
    try:
        import fsnapshot
        fsnok = True
    except:
        fsnok = False

    #import pdb; pdb.set_trace()

    plotfile = os.path.join(path,root)
    if os.path.isdir(plotfile) & fsnok:
        (nx, ny, nz) = fsnapshot.fplotfile_get_size(plotfile)
        x = np.zeros( (nx), dtype=np.float64)
        y = np.zeros( (nx), dtype=np.float64)
        (y, x, npts, err) = fsnapshot.fplotfile_get_data_1d(plotfile, comp, y, x)
    else:
        x = np.zeros( (0), dtype=np.float64)
        y = np.zeros( (0), dtype=np.float64)
    
    return (x, y)


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
    # c_pflotran = c_pflotran.flatten()
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
    import run_amanzi_standard
    import numpy as np

    # root name for problem
    root = "isotherms"
    components = ['A','B','C']
    compcrunch = ['A']

    # times
#    timespfl = ['Time:  0.00000E+00 y', 'Time:  5.00000E+01 y',
    timespfl = ['Time:  5.00000E+01 y',]
    timesama  = ['71']
    timesama2 = ['71']

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


# pflotran data --->
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


# CrunchFlow --->
    path_to_crunchflow = "crunchflow"
    try: 
        times_CF = ['totcon5.out']
        comp = 0
        u_crunchflow = []
        ignore = 4
        for i, time in enumerate(times_CF):
           x_crunchflow, c_crunchflow = GetXY_CrunchFlow(path_to_crunchflow,root,time,comp,ignore)
           u_crunchflow = u_crunchflow + [c_crunchflow]
        crunch = True

    except: 
        crunch = False    


# Amanzi native chemistry --->
    try:
        input_filename = os.path.join("amanzi-u-1d-"+root+".xml")
        path_to_amanzi = "amanzi-native-output"
        run_amanzi_standard.run_amanzi(input_filename, 1, [root+".bgd"], path_to_amanzi)

        u_native = [[[] for x in range(len(amanzi_totc))] for x in range(len(timesama))]
        for i, time in enumerate(timesama):
           for j, comp in enumerate(amanzi_totc):
#              import pdb; pdb.set_trace()
              x_native, c_native = GetXY_Amanzi(path_to_amanzi,root,comp)
              u_native[i][j] = c_native

        v_native = [[[] for x in range(len(amanzi_sorb))] for x in range(len(timesama))]
        for i, time in enumerate(timesama):
           for j, comp in enumerate(amanzi_sorb):
              x_native, c_native = GetXY_Amanzi(path_to_amanzi,root,comp)
              v_native[i][j] = c_native

        native = len(x_native)  

    except:
        native = 0 
        pass


# Amanzi-Alquimia-PFloTran --->
    try:  
        input_filename = os.path.join("amanzi-u-1d-"+root+"-alq.xml")
        path_to_amanzi = "amanzi-alquimia-output"
        run_amanzi_standard.run_amanzi(input_filename, 1, ["1d-"+root+"-trim.in",root+".dat"], path_to_amanzi)

        u_alquimia = [[[] for x in range(len(amanzi_totc))] for x in range(len(timesama))]
        for i, time in enumerate(timesama):
           for j, comp in enumerate(amanzi_totc):
              x_alquimia, c_alquimia = GetXY_Amanzi(path_to_amanzi,root,comp)
              u_alquimia[i][j] = c_alquimia
              
        v_alquimia = [[[] for x in range(len(amanzi_sorb))] for x in range(len(timesama))]
        for i, time in enumerate(timesama):
           for j, comp in enumerate(amanzi_sorb):
              x_alquimia, c_alquimia = GetXY_Amanzi(path_to_amanzi,root,comp)
              v_alquimia[i][j] = c_alquimia

        alq = True

    except:
        alq = False


# Amanzi-Alquimia-Crunch --->
    try:  
        input_filename = os.path.join("amanzi-u-1d-"+root+"-alq-crunch.xml")
        path_to_amanzi = "amanzi-alquimia-crunch-output"
        run_amanzi_standard.run_amanzi(input_filename, 1, ["1d-"+root+"-crunch.in",root+".dbs"], path_to_amanzi)

        u_alquimia_crunch = [[[] for x in range(len(amanzi_totc_crunch))] for x in range(len(timesama))]
        for i, time in enumerate(timesama):
           for j, comp in enumerate(amanzi_totc_crunch):
              x_alquimia_crunch, c_alquimia_crunch = GetXY_Amanzi(path_to_amanzi,root,comp)
              u_alquimia_crunch[i][j] = c_alquimia_crunch
              
        v_alquimia_crunch = [[[] for x in range(len(amanzi_sorb_crunch))] for x in range(len(timesama))]
        for i, time in enumerate(timesama):
           for j, comp in enumerate(amanzi_sorb_crunch):
              x_alquimia_crunch, c_alquimia_crunch = GetXY_Amanzi(path_to_amanzi,root,comp)
              v_alquimia_crunch[i][j] = c_alquimia_crunch

        alqc = True

    except:
        alqc = False


# Amanzi-structured --->

    # +pflotran
    try:
        input_filename = os.path.join("amanzi-s-1d-isotherms-alq-pflo.xml")
        path_to_amanziS = "struct_amanzi-output-pflo"
        run_amanzi_standard.run_amanzi(input_filename, 1, [], path_to_amanziS)
        root_amanziS = "plt00051"
        c_amanziS = [ [] for x in range(len(amanzi_totc)) ]
        v_amanziS = [ [] for x in range(len(amanzi_totc)) ]
        for j,comp in enumerate(components):
           compS = "{0}_Aqueous_Concentration".format(comp)
           x_amanziS, c_amanziS[j] = GetXY_AmanziS(path_to_amanziS,root_amanziS,time,compS)
           compS = "{0}_Sorbed_Concentration".format(comp)
           x_amanziS, v_amanziS[j] = GetXY_AmanziS(path_to_amanziS,root_amanziS,time,compS)
        struct = len(x_amanziS)
    except:
        struct = 0

    # +crunchflow
    try:
        input_filename = os.path.join("amanzi-s-1d-isotherms-alq-crunch.xml")
        path_to_amanziS = "struct_amanzi-output-crunch"
        run_amanzi_standard.run_amanzi(input_filename, 1, [], path_to_amanziS)
        root_amanziS = "plt00051"
        compS = "A_Aqueous_Concentration"
        x_amanziS_crunch, c_amanziS_crunch = GetXY_AmanziS(path_to_amanziS,root_amanziS,time,compS)
        compS = "A_Sorbed_Concentration"
        x_amanziS_crunch, v_amanziS_crunch = GetXY_AmanziS(path_to_amanziS,root_amanziS,time,compS)        
        struct_c = len(x_amanziS_crunch)
    except:
        struct_c = 0

## plotting ---------------------------------

    # subplots
    fig, ax = plt.subplots(4,sharex=True,figsize=(8,12))

#    colors= ['r'] #,'b','m','g'] # components
#    styles = ['-','v','o','x'] # codes
#    codes = ['AmanziU (2nd-Ord.)+Alq(PFT)','AmanziU (2nd-Ord.) Native','PFloTran'] + [None,]*9

    # lines on axes

    # for Kd:
    #   ax[0] ---> Aqueous concentrations
    #   ax[1] ---> Sorbed concentrations

    # for Langmuir and Freundlich
    #   ax[2] ---> Aqueous concentrations
    #   ax[3] ---> Sorbed concentrations

    # for i, time in enumerate(times):
    i = 0 # hardwired 50 years -- because the second entry in the list was taken at cycle 71 = 50 years.

#  pflotran
    ax[0].plot(x_pflotran, u_pflotran[i][0],color='m',linestyle='-',linewidth=2,label='PFloTran')    
    ax[1].plot(x_pflotran, v_pflotran[i][0],color='m',linestyle='-',linewidth=2)

    ax[2].plot(x_pflotran, u_pflotran[i][1],color='k',linestyle='-',linewidth=2,label='Langmuir PFloTran ')    
    ax[3].plot(x_pflotran, v_pflotran[i][1],color='k',linestyle='-',linewidth=2)

    ax[2].plot(x_pflotran, u_pflotran[i][2],color='c',linestyle='-',linewidth=2,label='Freundlich PFloTran')    
    ax[3].plot(x_pflotran, v_pflotran[i][2],color='c',linestyle='-',linewidth=2)

# crunchflow
    if crunch:
        ax[0].plot(x_crunchflow, u_crunchflow[0],color='m',linestyle='None',marker='*', label='CrunchFlow OS3D')
        # crunchflow does not output sorbed concentrations
 
# native 
    if native:
       ax[0].plot(x_native, u_native[i][0],color='b',linestyle='None',marker='x')
       ax[1].plot(x_native, v_native[i][0],color='b',linestyle='None',marker='x',label='AmanziU (2nd-Ord.) Native')

       ax[2].plot(x_native, u_native[i][1],color='k',linestyle='None',marker='x')
       ax[3].plot(x_native, v_native[i][1],color='k',linestyle='None',marker='x',label='Langmuir AmanziU (2nd-Ord.) Native')

       ax[2].plot(x_native, u_native[i][2],color='c',linestyle='None',marker='x')
       ax[3].plot(x_native, v_native[i][2],color='c',linestyle='None',marker='x',label='Freundlich AmanziU (2nd-Ord.) Native')


# unstructured alquimia pflotran
    if alq:
       ax[0].plot(x_alquimia, u_alquimia[i][0],color='r',linestyle='-',linewidth=2)
       ax[1].plot(x_alquimia, v_alquimia[i][0],color='r',linestyle='-',linewidth=2,label='AmanziU (2nd-Ord.)+Alq(PFT)')

       ax[2].plot(x_alquimia, u_alquimia[i][1],color='k',linestyle='--',linewidth=2)
       ax[3].plot(x_alquimia, v_alquimia[i][1],color='k',linestyle='--',linewidth=2,label='Langmuir AmanziU (2nd-Ord.)+Alq(PFT)')

       ax[2].plot(x_alquimia, u_alquimia[i][2],color='c',linestyle='--',linewidth=2)
       ax[3].plot(x_alquimia, v_alquimia[i][2],color='c',linestyle='--',linewidth=2,label='Freundlich AmanziU (2nd-Ord.)+Alq(PFT)')

# unstructured alquimia crunch
    if alqc:
       ax[0].plot(x_alquimia_crunch, u_alquimia_crunch[i][0],color='r',linestyle='None',marker='*',linewidth=2)
       ax[1].plot(x_alquimia_crunch, v_alquimia_crunch[i][0],color='r',linestyle='None',marker='*',linewidth=2,label='AmanziU (2nd-Ord.)+Alq(CF)')

# structured alquimia pflotran
    if (struct>0):
        sam = ax[0].plot(x_amanziS, c_amanziS[0],'g-',label='AmanziS+Alq(PFT)',linewidth=2)
        samv = ax[1].plot(x_amanziS, v_amanziS[0],'g-',linewidth=2)

        sam1 = ax[2].plot(x_amanziS, c_amanziS[1],'k*',label='Langmuir AmanziS+Alq(PFT)',linewidth=2)
        samv1 = ax[3].plot(x_amanziS, v_amanziS[1],'k*',linewidth=2)

        sam2 = ax[2].plot(x_amanziS, c_amanziS[2],'c*',label='Freundlich AmanziS+Alq(PFT)',linewidth=2)
        samv2 = ax[3].plot(x_amanziS, v_amanziS[2],'c*',linewidth=2)

# structured alquimia crunch
    if (struct_c>0):
        samc = ax[0].plot(x_amanziS_crunch, c_amanziS_crunch,'g*',label='AmanziS+Alq(CF)',linewidth=2) 
        samcv = ax[1].plot(x_amanziS_crunch, v_amanziS_crunch,'g*',linewidth=2) #,markersize=20) 

    # for i, time in enumerate(times):
    ##for j, comp in enumerate(components):
    ##        if alq:
    ##               ax[1].plot(x_alquimia, v_alquimia[i][j],color=colors[j],linestyle=styles[0],linewidth=2,label=codes[j*len(styles)])
     ##       if isv2:
     ##              
     ##       ax[1].plot(x_native, v_native[i][j],markeredgecolor=colors[j],marker=styles[2],linewidth=2,label=codes[j*len(styles)+2],mfc='None',markeredgewidth=2,linestyle='None')
     ##       ax[1].plot(x_pflotran, v_pflotran[i][j],color=colors[j],linestyle='None',marker=styles[3],linewidth=2,label=codes[j*len(styles)+3])

    # axes
    ax[0].set_title("Kd linear sorption model",fontsize=15)
    ax[1].set_xlabel("Distance (m)",fontsize=15)
    ax[0].set_ylabel("Total A \n Concentration \n [mol/L]",fontsize=15)
    ax[1].set_ylabel("Total A \n Sorbed Concent. \n [mol/m3]",fontsize=15)

    ax[2].set_title("Langmuir and Freundlich sorption models",fontsize=15)
    ax[3].set_xlabel("Distance (m)",fontsize=15)
    ax[2].set_ylabel("Total B, C \n Concentration \n [mol/L]",fontsize=15)
    ax[3].set_ylabel("Total B, C \n Sorbed Concent. \n [mol/m3]",fontsize=15)

    ax[0].legend(loc='upper right',fontsize=10)
    ax[1].legend(loc='upper right',fontsize=10)
    ax[0].set_xlim(left=30,right=70)
    ax[1].set_xlim(left=30,right=70)

    ax[2].legend(loc='upper right',fontsize=10)
    ax[3].legend(loc='upper right',fontsize=8)
    ax[2].set_xlim(left=30,right=70)
    ax[3].set_xlim(left=30,right=70)

    # plot adjustments
    plt.tight_layout() 
    plt.subplots_adjust(left=0.20,bottom=0.15,right=0.95,top=0.90)
    plt.suptitle("Amanzi 1D "+root.title()+" Benchmark at 50 years",x=0.57,fontsize=20)
    plt.tick_params(axis='both', which='major', labelsize=15)

    # pyplot.show()
    plt.savefig(root+"_1d.png",format="png")
    # plt.close()

    # finally:
    #     pass 
