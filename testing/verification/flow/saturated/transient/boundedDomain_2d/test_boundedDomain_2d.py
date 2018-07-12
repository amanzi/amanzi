import amanzi_boundedDomain_2d

import os
import sys
import h5py
import numpy as np
import pytest


def test_boundedDomain_2d():

    try:
        sys.path.append('../../../../../../tools/amanzi_xml')
    except:
        pass

    import run_amanzi_standard
    input_file =os.path.join("amanzi_boundedDomain_2d.xml")
    run_dir = "amanzi-output"
    golden_dir = "golden_output"
    
    print "run_dir", run_dir
    
    CWD = os.getcwd()
    try:
        print "CWD", CWD
        print "input_file", input_file
        run_amanzi_standard.run_amanzi(input_file, 16, [input_file], run_dir)
        obs_data = amanzi_boundedDomain_2d.load_amanzi_obs(run_dir + "/observations.out")
        gold_data = amanzi_boundedDomain_2d.load_amanzi_obs(golden_dir + "/observations.out")

        err_obs = 0.

        assert err_obs < 1e-5
        
    finally:
        os.chdir(CWD)
