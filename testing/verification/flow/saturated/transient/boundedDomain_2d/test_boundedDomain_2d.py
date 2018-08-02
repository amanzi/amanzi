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

    try:
        sys.path.append('../../../../../../tools/testing')
    except:
        pass

    import compare_observation_results as comp_obs
    
    CWD = os.getcwd()
    try:
        print "CWD", CWD
        print "input_file", input_file
        #run_amanzi_standard.run_amanzi(input_file, 16, [input_file], run_dir)
        obs_name = "drawdown"
        obs_name1 = "Drawdown"        
        run_obs_coord, run_obs_val = comp_obs.get_observation_data(run_dir + "/observations.gold.out", obs_name1)
        gold_obs_coord, gold_obs_val = comp_obs.get_observation_data(golden_dir + "/observations.out", obs_name1)

        errors = comp_obs.compute_error_norms(run_obs_coord, run_obs_val, gold_obs_coord, gold_obs_val)

        for region, error in errors.items():        
        
            assert error < 1e-5
        
    finally:
        os.chdir(CWD)
