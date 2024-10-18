import os
import sys
import h5py
import numpy as np

import run_amanzi_standard
from compare_field_results import GetXY_AmanziU_1D
from compare_field_results import GetXY_AmanziS_1D


if __name__ == "__main__":

    path_to_golden = "golden_output"
    if len(sys.argv) > 1:
        path_to_golden = sys.argv[1]
        path_to_golden += "/golden_output"

    try:
        path_to_amanziS = "."
        root_amanziS = "steady-flow-isv2"
        compS = "Aqueous_Pressure"
        x_amanziS, c_amanziS = GetXY_AmanziS_1D(path_to_amanziS,root_amanziS,compS)
        struct = len(x_amanziS)
    except:
        struct = 0
        
    try:
        comp = 'hydraulic_head.cell.0'
        path_to_amanziU = "."
        root_amanziU = 'steady-flow-isv2'
        x_amanziU, c_amanziU = GetXY_AmanziU_1D(path_to_amanziU,root_amanziU,comp,3)
        unstruct = len(x_amanziU)
    except:
        unstruct = 0

    try:
        comp = 'hydraulic_head.cell.0'
        path_to_amanziU = path_to_golden
        root_amanziU = 'steady-flow-isv2'
        x_amanziU_gold, c_amanziU_gold = GetXY_AmanziU_1D(path_to_amanziU,root_amanziU,comp,3)
        unstruct_gold = len(x_amanziU_gold)
    except:
        unstruct_gold = 0

    # Diff
    msg = ""
    if (unstruct and unstruct_gold):
        diff = c_amanziU_gold - c_amanziU
        error = np.linalg.norm(diff)

    # Report
        tol = 1e-8
        if error < tol:
            msg = msg + "Comparison Passed"
            msg = msg + "\n  error = " + str(error)
            print(msg)
        else:
            msg = msg + "Comparison Failed"
            msg = msg + "\n  error = " + str(error)
            sys.exit(msg)
    else:
        msg = msg + "Comparison Failed"
        msg = msg + "\n  tests results or golden_output missing"
        sys.exit(msg)

