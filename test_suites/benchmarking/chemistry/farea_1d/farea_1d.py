# plots cation concentration along x at last time step 
# benchmark: compares to pflotran simulation results
# author: S.Molins - Oct. 2013

import os
import sys
import h5py
import numpy as np
import matplotlib
from matplotlib import pyplot as plt

import run_amanzi_standard
from compare_field_results import GetXY_AmanziU_1D
from compare_field_results import GetXY_PFloTran_1D


if __name__ == "__main__":

    try:
        sys.path.append('../../../../tools/amanzi_xml')
    except:
        pass

    try:
        sys.path.append('../../../../MY_TPL_BUILD/ccse/ccse-1.3.4-source/Tools/Py_util')
    except:
        pass

    # root name for problem
    root = "farea-full"

    # components and minerals
    components = ['H+', 'Al+++', 'Ca++', 'Cl-', 'Fe+++', 'CO2(aq)', 'K+', 'Mg++', 'Na+', 'SiO2(aq)', 'SO4--', 'Tritium', 'NO3-', 'UO2++']
    minerals =['Quartz', 'Goethite', 'Kaolinite', 'Schoepite', 'Gibbsite', 'Jurbanite', 'Basaluminite', 'Opal']

    # amanzi output
    amanzi_totc_templ = "total_component_concentration.cell.%s"
    amanzi_totc = [amanzi_totc_templ%comp for comp in components]

    amanzi_sorb_templ = "total_sorbed.cell.{0}"
    amanzi_sorb = [amanzi_sorb_templ.format(x) for x in range(len(components))]

    amanzi_vf_templ = "mineral_volume_fractions.cell.{0}"
    amanzi_vf = [amanzi_vf_templ.format(x) for x in minerals]

    # pflotran output
    pflotran_totc_templ = "Total_{0} [M]"
    pflotran_totc = [pflotran_totc_templ.format(x) for x in components]

    pflotran_sorb_templ = "Total_Sorbed_{0} [mol_m^3]"
    pflotran_sorb = [pflotran_sorb_templ.format(x) for x in components]

    pflotran_vf_templ = "{0}_VF"
    pflotran_vf = [pflotran_vf_templ.format(x) for x in minerals]

    # hardwired time / add or remove here
    timespflo = ['Time:  5.00000E+01 y']
    timesama  = ['71']
    
    # hardwired selected components / add or remove here
    search = ['Ca++','Mg++','K+','Cl-'] # ['UO2++'] # ['Ca++','Mg++','K+','Cl-'] # ['Ca++', 'Mg++','SiO2(aq)']
    index = [components.index(comp) for comp in search]

    # hardwired selected minerals / add or remove here
    searchm = ['Kaolinite'] # 'Goethite', 'Kaolinite', 'Schoepite', 'Gibbsite'
    indexm = [minerals.index(min) for min in searchm] 

    # pflotran selected output
    totcpflo = [pflotran_totc[i] for i in index]
    sorbpflo = [pflotran_sorb[i] for i in index]
    vfpflo   = [pflotran_vf[i] for i in indexm]

    # amanzi selected output
    totcama = [amanzi_totc[i] for i in index]
    sorbama = [amanzi_sorb[i] for i in index]
    vfama   = [amanzi_vf[i] for i in indexm]

    # start with pflotran results
    path_to_pflotran = "pflotran"
    root_pflotran = "ascem-2012-1d-"+root
    time = timespflo[0]

    # tot concentration
    u_pflotran = [[] for x in range(len(totcpflo))]
    for j, comp in enumerate(totcpflo):          
        x_pflotran, c_pflotran = GetXY_PFloTran_1D(path_to_pflotran,root_pflotran,time,comp)
        u_pflotran[j] = c_pflotran

    # sorbed concentration    
    v_pflotran = [[] for x in range(len(sorbpflo))]
    for j, sorb in enumerate(sorbpflo):
        x_pflotran, c_pflotran = GetXY_PFloTran_1D(path_to_pflotran,root_pflotran,time,sorb)
        v_pflotran[j] = c_pflotran

    # mineral volume fraction
    w_pflotran = [[] for x in range(len(vfpflo))]
    for j, vf in enumerate(vfpflo):
        x_pflotran, c_pflotran = GetXY_PFloTran_1D(path_to_pflotran,root_pflotran,time,vf)
        w_pflotran[j] = c_pflotran

    CWD = os.getcwd()
    local_path = "" 
    

    # Amanzi + Native chemistry
    try:
        input_file = os.path.join("amanzi-u-1d-"+root+".xml")
        path_to_amanzi = "output-u"
        run_amanzi_standard.run_amanzi(input_file, 1, [root+".bgd",input_file], path_to_amanzi)

        time = timesama[0]

        # tot conc
        u_amanzi_native = [[] for x in range(len(totcama))]
        for j, comp in enumerate(totcama):
            x_amanzi_native, c_amanzi_native = GetXY_AmanziU_1D(path_to_amanzi,'plot',comp,1)
            u_amanzi_native[j] = c_amanzi_native

        # sorb conc
        v_amanzi_native = [[] for x in range(len(sorbama))]
        for j, sorb in enumerate(sorbama):
            x_amanzi_native, c_amanzi_native = GetXY_AmanziU_1D(path_to_amanzi,'plot',sorb,1)
            v_amanzi_native[j] = c_amanzi_native

        # mineral volume fraction
        w_amanzi_native = [[] for x in range(len(vfama))]
        for j, vf in enumerate(vfama):
            x_amanzi_native, c_amanzi_native = GetXY_AmanziU_1D(path_to_amanzi,'plot',vf,1)
            w_amanzi_native[j] = c_amanzi_native

        native = True

    except:
        native = False


    # Amanzi + Alquimia + PFloTran chemistry
    try:  
        input_file = os.path.join("amanzi-u-1d-"+root+"-alq-pflo.xml")
        path_to_amanzi = "output-u-alq-pflo"
        run_amanzi_standard.run_amanzi(input_file, 1, 
                                       ["1d-"+root+"-trim.in", root+".dat",input_file],
                                       path_to_amanzi)

        # import pdb; pdb.set_trace()
        time = timesama[0]

        # tot concentration
        u_amanzi_alquimia = [[] for x in range(len(totcama))]
        for j, comp in enumerate(totcama):
            x_amanzi_alquimia, c_amanzi_alquimia = GetXY_AmanziU_1D(path_to_amanzi,"plot",comp,1)
            u_amanzi_alquimia[j] = c_amanzi_alquimia  

        # sorbed concentration
        v_amanzi_alquimia = [[] for x in range(len(sorbama))]
        for j, sorb in enumerate(sorbama):
            x_amanzi_alquimia, c_amanzi_alquimia = GetXY_AmanziU_1D(path_to_amanzi,"plot",sorb,1)
            v_amanzi_alquimia[j] = c_amanzi_alquimia

        # mineral volume fraction
        w_amanzi_alquimia = [[] for x in range(len(vfama))]
        for j, vf in enumerate(vfama):
            x_amanzi_alquimia, c_amanzi_alquimia = GetXY_AmanziU_1D(path_to_amanzi,"plot",vf,1)
            w_amanzi_alquimia[j] = c_amanzi_alquimia

        alq = True

    except:
        alq = False

    # Amanzi + Alquimia + PFloTran chemistry with writer
    try:  
        input_file = os.path.join("amanzi-u-1d-"+root+"-alq-pflo-writer.xml")
        path_to_amanzi = "output-u-alq-pflo-writer"
        run_amanzi_standard.run_amanzi(input_file, 1, 
                                       ["1d-"+root+"-trim.in", root+".dat",input_file],
                                       path_to_amanzi)

        # import pdb; pdb.set_trace()
        time = timesama[0]

        # tot concentration
        u_amanzi_alquimia_w = [[] for x in range(len(totcama))]
        for j, comp in enumerate(totcama):
            x_amanzi_alquimia_w, c_amanzi_alquimia_w = GetXY_AmanziU_1D(path_to_amanzi,"plot",comp,1)
            u_amanzi_alquimia_w[j] = c_amanzi_alquimia_w  

        # sorbed concentration
        v_amanzi_alquimia_w = [[] for x in range(len(sorbama))]
        for j, sorb in enumerate(sorbama):
            x_amanzi_alquimia_w, c_amanzi_alquimia_w = GetXY_AmanziU_1D(path_to_amanzi,"plot",sorb,1)
            v_amanzi_alquimia_w[j] = c_amanzi_alquimia_w

        # mineral volume fraction
        w_amanzi_alquimia_w = [[] for x in range(len(vfama))]
        for j, vf in enumerate(vfama):
            x_amanzi_alquimia_w, c_amanzi_alquimia_w = GetXY_AmanziU_1D(path_to_amanzi,"plot",vf,1)
            w_amanzi_alquimia_w[j] = c_amanzi_alquimia_w

        alq_writer = True

    except:
        alq_writer = False
        
    # initialize subplots
    fig, ax = plt.subplots(3,sharex=True,figsize=(8,15))
    # bx =[None,]*3
    # bx[0] = ax[0].twinx()
    # bx[2] = ax[2].twinx()

    colors= ['r','b','m','g'] # components
    colors2= ['c','k','g','y'] # components
    styles = ['-','o','x','--'] # codes
    codes = ['Amanzi+Alquimia(PFloTran)','Amanzi Native Chemistry','PFloTran'] + [None,]*9

    # lines on axes
    # ax[0],b[0] ---> Aqueous concentrations
    # ax[1],b[1] ---> Sorbed concentrations

    for j, comp in enumerate(search):
        if alq:
            ax[0].plot(x_amanzi_alquimia, u_amanzi_alquimia[j],color=colors[j],linestyle=styles[0],linewidth=2)
        if alq_writer:
            ax[0].plot(x_amanzi_alquimia_w, u_amanzi_alquimia_w[j],color=colors[j],linestyle=styles[3],linewidth=2)
        if native:
            ax[0].plot(x_amanzi_native, u_amanzi_native[j],color=colors[j],marker=styles[1],linestyle='None',linewidth=2,label=comp)

        ax[0].plot(x_pflotran, u_pflotran[j],color=colors[j],linestyle='None',marker=styles[2],linewidth=2)
        
        if alq:
            label='Amanzi+Alquimia(PFloTran)' if j==0 else None
            ax[1].plot(x_amanzi_alquimia, v_amanzi_alquimia[j],color=colors[j],linestyle=styles[0],linewidth=2,label=label)
        if alq_writer:
            label='Amanzi+Alquimia(PFloTran)-W' if j==0 else None
            ax[1].plot(x_amanzi_alquimia_w, v_amanzi_alquimia_w[j],color=colors[j],linestyle=styles[3],linewidth=2,label=label)
        if native:
            label='Amanzi Native Chemistry' if j==0 else None
            ax[1].plot(x_amanzi_native, v_amanzi_native[j],color=colors[j],marker=styles[1],linestyle='None',linewidth=2,label=label)
        label='PFloTran' if j==0 else None
        ax[1].plot(x_pflotran, v_pflotran[j],color=colors[j],linestyle='None',marker=styles[2],linewidth=2,label=label)

    # ax[2],b[2] ---> Mineral Volume Fractions

    for j, vf in enumerate(searchm):
        if alq:
            ax[2].plot(x_amanzi_alquimia, w_amanzi_alquimia[j],color=colors2[j],linestyle=styles[0],linewidth=2)
        if alq_writer:
            ax[2].plot(x_amanzi_alquimia_w, w_amanzi_alquimia_w[j],color=colors2[j],linestyle=styles[3],linewidth=2)
        if native:
            ax[2].plot(x_amanzi_native, w_amanzi_native[j],color=colors2[j],marker=styles[1],linestyle='None',linewidth=2,label=vf) 
            ax[2].plot(x_pflotran, w_pflotran[j],color=colors2[j],linestyle='None',marker=styles[2],linewidth=2)

    # axes
    ax[2].set_xlabel("Distance (m)",fontsize=15)

    ax[0].set_ylabel("Total Concentration [mol/L]",fontsize=15)
    ax[1].set_ylabel("Total Sorbed Concent. [mol/m3]",fontsize=15)
    ax[2].set_ylabel("Mineral Volume Fraction [m3/m3]",fontsize=15)

    # plot adjustments
    plt.subplots_adjust(left=0.10,bottom=0.05,right=0.90,top=0.90)
    ax[0].legend(loc='center right',fontsize=15)
    ax[1].legend(loc='center right',fontsize=15)
    ax[2].legend(loc='center right',fontsize=15)

    # ax[2].set_ylim(bottom=-0.01)
    # bx[2].set_ylim(bottom=-2.0e-6)
 
    plt.suptitle("Amanzi 1D "+root.title()+" Benchmark at 50 years",x=0.57,fontsize=20)
    plt.tick_params(axis='both', which='major', labelsize=15)

    # pyplot.show()
    plt.savefig(root+"_1d.png",format="png")
    # plt.close()

    # finally:
    #     pass 
