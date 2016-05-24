# plots cation concentration along x at last time step 
# benchmark: compares to pflotran simulation results
# author: S.Molins - Oct. 2013

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
    root = "surface-complexation"
    components = ['Na+','NO3-','Zn++'] ## ['H+','Na+','NO3-','Zn++']
    pflotran_totc_templ = "Total_{0} [M]"
    comppflo = [pflotran_totc_templ.format(x) for x in components]
    compcrunch = [1,2,3] 

# times
    timespflo = ['Time:  5.00000E+01 y']          # ['Time:  0.00000E+00 y', 'Time:  5.00000E+01 y',]
    times = ['71'] #  ['1','72']
    times_CF = ['totcon5.out']
    times_CF_surf = ['totsurface5.out']
    times_CF_pH = ['pH5.out']

# amanzi output (unstructured - native and alquimia)
    amanzi_totc_templ = "total_component_concentration.cell.{0} conc" #Component {0} conc"
    amanzi_totc = [amanzi_totc_templ.format(x) for x in components] #range(len(components))]
    amanzi_totc_crunch = [amanzi_totc_templ.format(x) for x in compcrunch] #range(len(components))]

    amanzi_sorb_templ = "total_sorbed.cell.{0}"
    amanzi_sorb = [amanzi_sorb_templ.format(x+1) for x in range(len(components))]
    amanzi_sorb_crunch = [amanzi_sorb_templ.format(x+1) for x in range(len(compcrunch))]

# amanzi output (structured - alquimia)
    amanzi_totc_templ = "{0}_Aqueous_Concentration"
    amanzi_totcS = [amanzi_totc_templ.format(x) for x in components]

    amanzi_sorb_templ = "{0}_Sorbed_Concentration"
    amanzi_sorbS = [amanzi_sorb_templ.format(x) for x in components]

# read pflotran results --->
    path_to_pflotran = "pflotran"

    u_pflotran = [[[] for x in range(len(comppflo))] for x in range(len(timespflo))]
    for i, time in enumerate(timespflo):
       for j, comp in enumerate(comppflo):          
          x_pflotran, c_pflotran = GetXY_PFloTran(path_to_pflotran,root,time,comp)
          u_pflotran[i][j] = c_pflotran
    
    v_pflotran = [[[] for x in range(len(components))] for x in range(len(timespflo))]
    for i, time in enumerate(timespflo):
       for j, comp in enumerate(components):
          x_pflotran, c_pflotran = GetXY_PFloTran(path_to_pflotran,root,time,'Total_Sorbed_'+comp+' [mol_m^3]')
          v_pflotran[i][j] = c_pflotran

    pH_pflotran = [[] for x in range(len(timespflo))]
    comp = 'pH'
    for i, time in enumerate(timespflo):
          x_pflotran, c_pflotran = GetXY_PFloTran(path_to_pflotran,root,time,comp)
          pH_pflotran[i] = c_pflotran

# read crunchflow results --->
    path_to_crunch = "crunchflow"

    try: 
        u_crunchflow = [ [ [] for x in range(len(compcrunch)) ] for x in range(len(times_CF)) ]
        ignore = 4
        for i, time in enumerate(times_CF):
           for j, comp in enumerate(compcrunch):
              x_crunchflow, c_crunchflow = GetXY_CrunchFlow(path_to_crunch,root,time,comp,ignore)
              u_crunchflow[i][j] = c_crunchflow

        v_crunchflow = [ [ [] for x in range(len(compcrunch)) ] for x in range(len(times_CF_surf)) ]
        ignore = 4
        for i, time in enumerate(times_CF_surf):
           for j, comp in enumerate(compcrunch):
              x_crunchflow, c_crunchflow = GetXY_CrunchFlow(path_to_crunch,root,time,comp,ignore)
              v_crunchflow[i][j] = c_crunchflow

        pH_crunchflow = [ [] for x in range(len(times_CF_pH)) ]
        ignore = 4
        comp = 0
        for i, time in enumerate(times_CF_pH):
            y_crunchflow, c_crunchflow = GetXY_CrunchFlow(path_to_crunch,root,time,comp,ignore)
            pH_crunchflow[i] = c_crunchflow

        crunch = True

    except: 
        crunch = False

# Amanzi native chemistry --->
    try:

        input_filename = os.path.join("amanzi-u-1d-"+root+".xml")
        path_to_amanzi = "amanzi-native-output"
        #run_amanzi_chem.run_amanzi_chem("../"+input_filename,run_path=path_to_amanzi,chemfiles=[root+".bgd"])
        
        u_amanzi_native = [[[] for x in range(len(amanzi_totc))] for x in range(len(times))]
        for i, time in enumerate(times):
           for j, comp in enumerate(amanzi_totc):
              x_amanzi_native, c_amanzi_native = GetXY_Amanzi(path_to_amanzi,root,time,comp)
              u_amanzi_native[i][j] = c_amanzi_native

        v_amanzi_native = [[[] for x in range(len(amanzi_sorb))] for x in range(len(times))]
        for i, time in enumerate(times):
           for j, comp in enumerate(amanzi_sorb):
              x_amanzi_native, c_amanzi_native = GetXY_Amanzi(path_to_amanzi,root,time,comp)
              v_amanzi_native[i][j] = c_amanzi_native

        pH_amanzi_native = [ [] for x in range(len(times)) ]
        comp = 'free_ion_species.cell.H+'
        for i, time in enumerate(times):
              x_amanzi_native, c_amanzi_native = GetXY_Amanzi(path_to_amanzi,root,time,comp)
              pH_amanzi_native[i] = -np.log10(c_amanzi_native)

	native = True

    except:
      
        native = False

# Amanzi-Alquimia-PFlotran --->
    try:  

        input_filename = os.path.join("amanzi-u-1d-"+root+"-alq.xml")
        path_to_amanzi = "amanzi-alquimia-output"
        run_amanzi_chem.run_amanzi_chem("../"+input_filename,run_path=path_to_amanzi,chemfiles=["1d-"+root+".in",root+".dat"])

        u_amanzi_alquimia = [[[] for x in range(len(amanzi_totc))] for x in range(len(times))]
        for i, time in enumerate(times):
            for j, comp in enumerate(amanzi_totc):
                x_amanzi_alquimia, c_amanzi_alquimia = GetXY_Amanzi(path_to_amanzi,root,time,comp)
                u_amanzi_alquimia[i][j] = c_amanzi_alquimia

        v_amanzi_alquimia = [[[] for x in range(len(amanzi_sorb))] for x in range(len(times))]
        for i, time in enumerate(times):
           for j, comp in enumerate(amanzi_sorb):
              x_amanzi_alquimia, c_amanzi_alquimia = GetXY_Amanzi(path_to_amanzi,root,time,comp)
              v_amanzi_alquimia[i][j] = c_amanzi_alquimia

        pH_amanzi_alquimia = [ [] for x in range(len(times)) ]
        comp = 'pH.cell.0'
        for i, time in enumerate(times):
              x_amanzi_alquimia, c_amanzi_alquimia = GetXY_Amanzi(path_to_amanzi,root,time,comp)
              pH_amanzi_alquimia[i] = c_amanzi_alquimia ## -np.log10(c_amanzi_native)

        alq = True

    except:

        alq = False

# Amanzi-Alquimia-Crunch --->
    try:
        input_filename = os.path.join("amanzi-u-1d-surface-complexation-alq-crunch.xml")
        path_to_amanzi = "amanzi-alquimia-crunch-output"
        run_amanzi_chem.run_amanzi_chem("../"+input_filename,run_path=path_to_amanzi,chemfiles=["1d-surface-complexation-crunch.in","surface-complexation.dbs"])

        u_amanzi_alquimia_crunch = [[[] for x in range(len(amanzi_totc))] for x in range(len(times))]
        for i, time in enumerate(times):
           for j, comp in enumerate(amanzi_totc):
                x_amanzi_alquimia_crunch, c_amanzi_alquimia_crunch = GetXY_Amanzi(path_to_amanzi,root,time,comp)
                u_amanzi_alquimia_crunch[i][j] = c_amanzi_alquimia_crunch
              
        v_amanzi_alquimia_crunch = [[[] for x in range(len(amanzi_sorb))] for x in range(len(times))]
        for i, time in enumerate(times):
           for j, comp in enumerate(amanzi_sorb):
              x_amanzi_alquimia_crunch, c_amanzi_alquimia_crunch = GetXY_Amanzi(path_to_amanzi,root,time,comp)
              v_amanzi_alquimia_crunch[i][j] = c_amanzi_alquimia_crunch

        pH_amanzi_alquimia_crunch = [ [] for x in range(len(times)) ]
        comp = 'pH.cell.0'
        for i, time in enumerate(times):
              x_amanzi_alquimia_crunch, c_amanzi_alquimia_crunch = GetXY_Amanzi(path_to_amanzi,root,time,comp)
              pH_amanzi_alquimia_crunch[i] = c_amanzi_alquimia_crunch ## -np.log10(c_amanzi_native)

        alq_crunch = True

    except:
        
        alq_crunch = False


    # amanziS data  ------>
    
    # +pflotran
    try:
        input_filename = os.path.join("amanzi-s-1d-surface-complexation-alq.xml")
        path_to_amanziS = "struct_amanzi-output-pflo"
        run_amanzi_chem.run_amanzi_chem(input_filename,run_path=path_to_amanziS,chemfiles=None)
        root_amanziS = "plt00501"

        u_amanziS = [[[] for x in range(len(amanzi_totcS))] for x in range(len(times))]
        for i, time in enumerate(times):
           for j, comp in enumerate(amanzi_totcS):
                x_amanziS, c_amanziS = GetXY_AmanziS(path_to_amanziS,root_amanziS,time,comp)
                u_amanziS[i][j] = c_amanziS
              
        struct = len(x_amanziS)

#        import pdb; pdb.set_trace()

        v_amanziS = [[[] for x in range(len(amanzi_sorbS))] for x in range(len(times))]
        for i, time in enumerate(times):
           for j, comp in enumerate(amanzi_sorbS):
              x_amanziS, c_amanziS = GetXY_AmanziS(path_to_amanziS,root_amanziS,time,comp)
              v_amanziS[i][j] = c_amanziS

        pH_amanziS = [ [] for x in range(len(times)) ]
        comp = 'H+_Free_Ion_Guess'
        for i, time in enumerate(times):
              x_amanziS, c_amanziS = GetXY_AmanziS(path_to_amanziS,root_amanziS,time,comp)
              pH_amanziS[i] = -np.log10(c_amanziS)

    except:
        struct = 0

    # +crunchflow
    try:
        input_filename = os.path.join("amanzi-s-1d-surface-complexation-alq-crunch.xml")
        path_to_amanziS = "struct_amanzi-output-crunch"
#        import pdb; pdb.set_trace()
        run_amanzi_chem.run_amanzi_chem(input_filename,run_path=path_to_amanziS,chemfiles=None)
        root_amanziS = "plt00501"

        u_amanziS_crunch = [[[] for x in range(len(amanzi_totcS))] for x in range(len(times))]
        for i, time in enumerate(times):
           for j, comp in enumerate(amanzi_totcS):
                x_amanziS_crunch, c_amanziS_crunch = GetXY_AmanziS(path_to_amanziS,root_amanziS,time,comp)
                u_amanziS_crunch[i][j] = c_amanziS_crunch
              
        struct_c = len(x_amanziS_crunch)

        v_amanziS_crunch = [[[] for x in range(len(amanzi_sorbS))] for x in range(len(times))]
        for i, time in enumerate(times):
           for j, comp in enumerate(amanzi_sorbS):
              x_amanziS_crunch, c_amanziS_crunch = GetXY_AmanziS(path_to_amanziS,root_amanziS,time,comp)
              v_amanziS_crunch[i][j] = c_amanziS_crunch

        pH_amanziS_crunch = [ [] for x in range(len(times)) ]
        comp = 'H+_Free_Ion_Guess'
        for i, time in enumerate(times):
              x_amanziS_crunch, c_amanziS_crunch = GetXY_AmanziS(path_to_amanziS,root_amanziS,time,comp)
              pH_amanziS_crunch[i] = -np.log10(c_amanziS_crunch)

    except:
        struct_c = 0

## ------------------------------------------------------------------------

    ## colors= ['r','b','m','g'] # components
    ## styles = ['-','--','x'] # codes
    ## codes = ['Amanzi+Alquimia(PFloTran)','Amanzi Native Chemistry','PFloTran'] + [None,]*9

    # lines on axes
    # ax[0],b[0] ---> Aqueous concentrations
    # ax[1]      ---> pH
    # ax[2],b[2] ---> Sorbed concentrations

    # first
    # ax[0],b[0] ---> Aqueous concentrations
    # ax[1]      ---> pH

# define subplots as a subgrid
    fig = plt.figure(figsize=(15,12))

    ax = [] # total concentrations
    bx = [] # sorbed concentrations

    nrows = len(components) + 1
    ncols = 2 

    for comp in range(len(components)):

        ax += [plt.subplot(nrows,ncols,2*comp+1)]
        bx += [plt.subplot(nrows,ncols,2*comp+2)]

    px = plt.subplot(nrows,1,nrows)

    # for i, time in enumerate(times):
    i = 0 # only one time point at 50 years

# pflotran 
#    if True:
    for j, comp in enumerate(components):

           ax[j].plot(x_pflotran, u_pflotran[i][j],color='m',linestyle='-',linewidth=6,label='PFloTran')
           bx[j].plot(x_pflotran, v_pflotran[i][j],color='m',linestyle='-',linewidth=2)
           ax[j].text(x_pflotran[10],u_pflotran[i][j][10],comp,fontsize=15,bbox=dict(facecolor='white', alpha=1.0))
           bx[j].text(x_pflotran[10],v_pflotran[i][j][10],comp,fontsize=15,bbox=dict(facecolor='white', alpha=1.0))

    px.plot(x_pflotran, pH_pflotran[i],color='m',linestyle='-',linewidth=2,label='PFloTran')

# crunchflow 
#    if True:
    for j, comp in enumerate(components):

           ax[j].plot(x_crunchflow, u_crunchflow[i][j],'m*',linestyle='None',label='CrunchFlow')
           bx[j].plot(x_crunchflow, v_crunchflow[i][j],'m*',linestyle='None',linewidth=2)

    px.plot(y_crunchflow, pH_crunchflow[i],'m*',linestyle='None',linewidth=2,label='CrunchFlow')

# amanzi-native
    if native :
        for j, comp in enumerate(components):

            ax[j].plot(x_amanzi_native, u_amanzi_native[i][j],color='b',linestyle='None',marker='x',linewidth=2)
            bx[j].plot(x_amanzi_native, v_amanzi_native[i][j],color='b',linestyle='None',marker='x',linewidth=2,label='AmanziU Native Chemistry')

        px.plot(x_amanzi_native, pH_amanzi_native[i],color='b',linestyle='None',marker='x',linewidth=2,label='AmanziU Native Chemistry')

# amanzi-unstructured-alquimia-pflotran
    if alq:

        for j, comp in enumerate(components):

            ax[j].plot(x_amanzi_alquimia, u_amanzi_alquimia[i][j],color='r',linestyle='-',linewidth=2)
            bx[j].plot(x_amanzi_alquimia, v_amanzi_alquimia[i][j],color='r',linestyle='-',linewidth=2,label='AmanziU+Alquimia(PFloTran)')

        px.plot(x_amanzi_alquimia, pH_amanzi_alquimia[i],color='r',linestyle='-',linewidth=2,label='AmanziU+Alquimia(PFloTran)')

# amanzi-unstructured-alquimia-crunchflow
    if alq_crunch:

        for j, comp in enumerate(components):

            ax[j].plot(x_amanzi_alquimia_crunch, u_amanzi_alquimia_crunch[i][j],color='r',linestyle='None',marker='*',linewidth=2)
            bx[j].plot(x_amanzi_alquimia_crunch, v_amanzi_alquimia_crunch[i][j],color='r',linestyle='None',marker='*',linewidth=2,label='AmanziU+Alquimia(CrunchFlow)')

        px.plot(x_amanzi_alquimia_crunch, pH_amanzi_alquimia_crunch[i],color='r',linestyle='None',marker='*',linewidth=2,label='AmanziU+Alquimia(CrunchFlow)')

# amanzi-structured-alquimia-pflotran
    if struct:

        for j, comp in enumerate(components):

            ax[j].plot(x_amanziS, u_amanziS[i][j],color='g',linestyle='-',linewidth=2)
            bx[j].plot(x_amanziS, v_amanziS[i][j],color='g',linestyle='-',linewidth=2,label='AmanziS+Alquimia(PFloTran)')

        px.plot(x_amanziS, pH_amanziS[i],color='g',linestyle='-',linewidth=2,label='AmanziS+Alquimia(PFloTran)')

# amanzi-structured-alquimia-crunchflow
    if struct_c:

        for j, comp in enumerate(components):

            ax[j].plot(x_amanziS_crunch, u_amanziS_crunch[i][j],color='g',linestyle='None',marker='*',linewidth=2)
            bx[j].plot(x_amanziS_crunch, v_amanziS_crunch[i][j],color='g',linestyle='None',marker='*',linewidth=2,label='AmanziS+Alquimia(CrunchFlow)')

        px.plot(x_amanziS_crunch, pH_amanziS_crunch[i],color='g',linestyle='None',marker='*',linewidth=2,label='AmanziS+Alquimia(CrunchFlow)')

    # axes
    ax[len(components)-1].set_xlabel("Distance (m)",fontsize=15)
    bx[len(components)-1].set_xlabel("Distance (m)",fontsize=15)

##    for i,comp in enumerate(components):
    i=1
    ax[i].set_ylabel("Total Concentration [mol/L]",fontsize=15)
    bx[i].set_ylabel("Total Sorbed Concent. [mol/m3]",fontsize=15)

    px.set_xlabel("Distance(m)",fontsize=15)
    px.set_ylabel("pH",fontsize=15)

#    for i,comp in enumerate(components):
#        ax[i].set_ylim(bottom=0)
#        bx[i].set_ylim(bottom=0)
    px.set_ylim(bottom=4.8)

    # plot adjustments
    ax[0].legend(fontsize=15)
    bx[0].legend(fontsize=15)
    px.legend(fontsize=15,loc='upper right')

    plt.suptitle("Amanzi 1D "+root.title()+" Benchmark at 50 years",fontsize=20) #,x=0.57,fontsize=20)

    plt.tick_params(axis='both', which='major', labelsize=15)
  
    plt.tight_layout() #(pad=0.4, w_pad=0.5, h_pad=1.0)

    plt.subplots_adjust(left=0.10,bottom=0.15,right=0.90,top=0.95)

    #pyplot.show()
    plt.savefig(root+"_1d.png",format="png")
    #plt.close()

    #finally:
    #    pass 
