import calcite_1d
import os
import sys
import h5py
import numpy as np

times = ['Time:  5.00000E+01 y']
root = "calcite"


def read_pflotran(path_to_pflotran,root,time):

    comp = 'Total_Ca++ [M]'
    Ca_pflotran = []
    for i, time in enumerate(times):
        x_pflotran, c_pflotran = calcite_1d.GetXY_PFloTran(path_to_pflotran,root,time,comp)
        Ca_pflotran = Ca_pflotran + [c_pflotran]

    comp = 'pH'
    pH_pflotran = []
    for i, time in enumerate(times):
        x_pflotran, c_pflotran = calcite_1d.GetXY_PFloTran(path_to_pflotran,root,time,comp)
        pH_pflotran = pH_pflotran + [c_pflotran]

    comp = 'Calcite_VF'
    VF_pflotran = []
    for i, time in enumerate(times):
        x_pflotran, c_pflotran =calcite_1d.GetXY_PFloTran(path_to_pflotran,root,time,comp)
        VF_pflotran = VF_pflotran + [c_pflotran]

    return (x_pflotran, Ca_pflotran, pH_pflotran, VF_pflotran)

    

def test_calcite_1d_AmanziUNative():
    input_file = os.path.join("amanzi-u-1d-calcite.xml")
    path_to_amanzi = "output-u"
    # complete_native, x_amanzi, Ca_amanzi_native, pH_amanzi_native, VF_amanzi_native = calcite_1d.RunAmanziUTest(input_file, path_to_amanzi, times)
    result = calcite_1d.RunAmanziUTest(root, input_file, path_to_amanzi, times)

    complete_native = result[0]
    x_amanzi = result[1]
    Ca_amanzi_native = result[2]
    pH_amanzi_native = result[3]
    VF_amanzi_native = result[4]

    assert complete_native==True
    
    path_to_pflotran = "pflotran"    
    x_pflotran, Ca_pflotran, pH_pflotran, VF_pflotran = read_pflotran(path_to_pflotran, root, times)
       
    err_ca = 0.
    err_ca_max = 0.
    err_ca_min = 1e+32
    for i, time in enumerate(times):       
        err_ca = err_ca + np.linalg.norm(Ca_pflotran[i] - Ca_amanzi_native[i].T)
        err_ca_max = max(err_ca_max, np.linalg.norm(Ca_pflotran[i] - Ca_amanzi_native[i].T, 1))
        #err_ca_min = min(err_ca_min, np.linalg.norm(Ca_pflotran[i] - Ca_amanzi_native[i].T, -1))
        sol_norm = np.linalg.norm(Ca_pflotran[i])
        sol_norm_max = np.linalg.norm(Ca_pflotran[i], 1)

    print err_ca/sol_norm,  err_ca_max/sol_norm_max
        
    assert err_ca/sol_norm <= 8e-2
    assert err_ca_max/sol_norm_max <= 4e-3

    err_ph = 0.
    err_ph_max = 0.
    err_ph_min = 1e+32   
    for i, time in enumerate(times):       
        err_ph = err_ph + np.linalg.norm(pH_pflotran[i] - pH_amanzi_native[i].T)
        err_ph_max = max(err_ca_max, np.linalg.norm(pH_pflotran[i] - pH_amanzi_native[i].T, 1))
        sol_norm = np.linalg.norm(pH_pflotran[i])
        sol_norm_max = np.linalg.norm(pH_pflotran[i], 1)
        #err_ph_min = min(err_ph_min, np.linalg.norm(pH_pflotran[i] - pH_amanzi_native[i].T, -1))

    # print pH_pflotran[i]
    # print pH_amanzi_native[i].T
    print err_ph/sol_norm, err_ph_max/sol_norm_max
        
    assert err_ph/sol_norm <= 2e-2
    assert err_ph_max/sol_norm_max <= 2e-3

    err_VF = 0.
    err_VF_max = 0.
    err_VF_min = 1e+32
    for i, time in enumerate(times):       
        err_VF = err_VF + np.linalg.norm(VF_pflotran[i] - VF_amanzi_native[i].T)
        err_VF_max = max(err_ca_max, np.linalg.norm(VF_pflotran[i] - VF_amanzi_native[i].T, 1))
        sol_norm = np.linalg.norm(pH_pflotran[i])
        sol_norm_max = np.linalg.norm(pH_pflotran[i], 1)
        #err_VF_min = min(err_VF_min, np.linalg.norm(VF_pflotran[i] - VF_amanzi_native[i].T, -1))

    # print VF_pflotran[i]
    # print VF_amanzi_native[i].T
    print err_VF/sol_norm, err_VF_max/sol_norm_max
        
    assert err_VF/sol_norm <= 1e-8
    assert err_VF_max/sol_norm_max <= 3e-7    


# def test_calcite_1d():
    
#     calcite_1d.RunAllTests()

    
