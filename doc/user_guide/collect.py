#!/usr/bin/env python

import os, sys, utils, shutil, subprocess
import optparse
import fnmatch

#  Create dictionary that describes:
#
#  - layout of directories
#  - index files to be created
#  - subdirectories (tests/tutorials) to be copied
#

#
# support routines
#
def find(pattern, path):
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return result

#
# Install
#
install={}
install['index'] = {
  'index_title' : 'Building Amanzi',
  'index_file' : 'doc/user_guide/install/index.rst',
  'index_list' : ['bootstrap','tpls','amanzi']
}
    
install['tpls'] = {
  'from_file' : 'config/SuperBuild/INSTALL',
  'dest_file' : 'doc/user_guide/install/building_tpls.rst',
  'index_entry' : 'building_tpls.rst'
}

install['amanzi'] = {
  'from_file' : 'doc/build_instructions/build_instructions.rst',
  'dest_file' : 'doc/user_guide/install/build_instructions.rst',
  'index_entry' : 'build_instructions.rst'
}

install['bootstrap'] = {
  'from_file' : 'doc/build_instructions/building_bootstrap.rst',
  'dest_file' : 'doc/user_guide/install/building_bootstrap.rst',
  'index_entry' : 'building_bootstrap.rst'
}

 
#
# Tutorials
#
tutorial = {}
tutorial['index'] = {
  'index_title' : 'Tutorial',
  'index_file' : 'doc/user_guide/tutorial/index.rst',
  'index_list' : ['flow_1','flow_2']
}
    
tutorial['flow_1'] = {
  'from_dir' : 'doc/tutorial/saturated',
  'dest_dir' : 'doc/user_guide/tutorial/steady_confined',
  'index_entry' : 'steady_confined/Tutorial_steady_linear.rst'
}
    
tutorial['flow_2'] = {
  'from_dir' : 'doc/tutorial/richards',
  'dest_dir' : 'doc/user_guide/tutorial/transient_infiltration',
  'index_entry' : 'transient_infiltration/Tutorial_dvz.rst'
}

#
# Verificaiton Tests
#
verification={}
verification['index'] = {
  'index_title' : 'Verification Testing',
  'index_file' : 'doc/user_guide/verification/index.rst',
  'index_list' : [ 'confined_flow',
                   'unconfined_flow',
                   'infiltration',
                   'transport'
                 ],
}

verification['confined_flow'] = {
  'index_entry' : 'confined_flow/index.rst',
  'index' : {'index_title' : 'Confined Flow Tests',
             'index_file' : 'doc/user_guide/verification/confined_flow/index.rst',
             'index_list' : [ 'linear_head_head_1d',
                              'linear_flux_head_1d',
                              'linear_materials_serial_1d',
                              'linear_materials_parallel_1d',
                              'theis_isotropic_1d',
                              'butler_pod_2d',
                              'butler_strip_2d',
                              'hantush_anisotropic_2d',
                              'boundedDomain_2d', 
                            ],
            },
  'linear_head_head_1d' : {
    'from_dir' : 'test_suites/verification/flow/saturated/steady-state/linear_head_head_1d',
    'dest_dir' : 'doc/user_guide/verification/confined_flow/linear_head_head_1d',
    'index_entry' : 'linear_head_head_1d/amanzi_linear_head_head_1d.rst'
  },
  'linear_flux_head_1d' : {
    'from_dir' : 'test_suites/verification/flow/saturated/steady-state/linear_flux_head_1d',
    'dest_dir' : 'doc/user_guide/verification/confined_flow/linear_flux_head_1d',
    'index_entry' : 'linear_flux_head_1d/amanzi_linear_flux_head_1d.rst'
  },
  'linear_materials_serial_1d' : {
    'from_dir' : 'test_suites/verification/flow/saturated/steady-state/linear_materials_serial_1d',
    'dest_dir' : 'doc/user_guide/verification/confined_flow/linear_materials_serial_1d',
    'index_entry' : 'linear_materials_serial_1d/amanzi_linear_materials_serial_1d.rst'
  },
  'linear_materials_parallel_1d' : {
    'from_dir' : 'test_suites/verification/flow/saturated/steady-state/linear_materials_parallel_1d',
    'dest_dir' : 'doc/user_guide/verification/confined_flow/linear_materials_parallel_1d',
    'index_entry' : 'linear_materials_parallel_1d/amanzi_linear_materials_parallel_1d.rst'
  },
  'theis_isotropic_1d' : {
    'from_dir' : 'test_suites/verification/flow/saturated/transient/theis_isotropic_1d',
    'dest_dir' : 'doc/user_guide/verification/confined_flow/theis_isotropic_1d',
    'index_entry' : 'theis_isotropic_1d/amanzi_theis_isotropic_1d.rst',
  },
  'butler_pod_2d' : {
    'from_dir' : 'test_suites/verification/flow/saturated/transient/butler_pod_2d',
    'dest_dir' : 'doc/user_guide/verification/confined_flow/butler_pod_2d',
    'index_entry' : 'butler_pod_2d/amanzi_butler_pod_2d.rst',
  },
  'butler_strip_2d' : {
    'from_dir' : 'test_suites/verification/flow/saturated/transient/butler_strip_2d',
    'dest_dir' : 'doc/user_guide/verification/confined_flow/butler_strip_2d',
    'index_entry' : 'butler_strip_2d/amanzi_butler_strip_2d.rst',
  },
  'hantush_anisotropic_2d' : {
    'from_dir' : 'test_suites/verification/flow/saturated/transient/hantush_anisotropic_2d',
    'dest_dir' : 'doc/user_guide/verification/confined_flow/hantush_anisotropic_2d',
    'index_entry' : 'hantush_anisotropic_2d/amanzi_hantush_anisotropic_2d.rst',
  },
  'boundedDomain_2d' : {
     'from_dir' : 'test_suites/verification/flow/saturated/transient/boundedDomain_2d',
     'dest_dir' : 'doc/user_guide/verification/confined_flow/boundedDomain_2d',
     'index_entry' : 'boundedDomain_2d/amanzi_boundedDomain_2d.rst',
  },
}

verification['unconfined_flow'] = {
  'index_entry': 'unconfined_flow/index.rst',
  'index' : {'index_title' : 'Unconfined Flow Tests',
             'index_file'  : 'doc/user_guide/verification/unconfined_flow/index.rst',
             'index_list'  : [ 'unconfined_no_recharge_1d', 
                               'unconfined_seepage_1d',
                               'unconfined_layered_2d',
                               'unconfined_recharge_1d',
                             ], 
            },
  'unconfined_no_recharge_1d': {
    'from_dir' : 'test_suites/verification/flow/richards/steady-state/unconfined_no_recharge_1d',
    'dest_dir' : 'doc/user_guide/verification/unconfined_flow/unconfined_no_recharge_1d',
    'index_entry' : 'unconfined_no_recharge_1d/amanzi_unconfined_no_recharge_1d.rst',
  },
  'unconfined_seepage_1d': {
    'from_dir' : 'test_suites/verification/flow/richards/steady-state/unconfined_seepage_1d',
    'dest_dir' : 'doc/user_guide/verification/unconfined_flow/unconfined_seepage_1d',
    'index_entry' : 'unconfined_seepage_1d/amanzi_unconfined_seepage_1d.rst',
  },
  'unconfined_layered_2d': {
    'from_dir' : 'test_suites/verification/flow/richards/steady-state/unconfined_layered_2d',
    'dest_dir' : 'doc/user_guide/verification/unconfined_flow/unconfined_layered_2d',
    'index_entry' : 'unconfined_layered_2d/amanzi_unconfined_layered_2d.rst'
  },
  'unconfined_recharge_1d': {
    'from_dir' : 'test_suites/verification/flow/richards/steady-state/unconfined_recharge_1d',
    'dest_dir' : 'doc/user_guide/verification/unconfined_flow/unconfined_recharge_1d',
    'index_entry' : 'unconfined_recharge_1d/amanzi_unconfined_recharge_1d.rst',
  },
}

verification['infiltration'] = {
  'index_entry': 'infiltration/index.rst',
  'index' : {'index_title' : 'Infiltration Flow Tests',
             'index_file'  : 'doc/user_guide/verification/infiltration/index.rst',
             'index_list'  : [ 'infiltration_1d', 
                             ], 
            },
  'infiltration_1d': {
    'from_dir' : 'test_suites/verification/flow/richards/steady-state/infiltration_1d',
    'dest_dir' : 'doc/user_guide/verification/infiltration/infiltration_1d',
    'index_entry' : 'infiltration_1d/amanzi_infiltration_1d.rst',
  },
}

verification['transport'] = {
  'index_entry': 'transport/index.rst',
  'index' : {'index_title' : 'Transport of Solutes: Advection, Dispersion and Diffusion',
             'index_file' : 'doc/user_guide/verification/transport/index.rst',
             'index_list' : [ 'dispersion_aligned_point_2d',
                              'dispersion_45_point_2d',
                              'dual_porosity_1d'
                            ],
            },
  'dispersion_aligned_point_2d' : {
    'from_dir' : 'test_suites/verification/transport/saturated/steady-state/dispersion_aligned_point_2d',
    'dest_dir' : 'doc/user_guide/verification/transport/dispersion_aligned_point_2d',
    'index_entry' : 'dispersion_aligned_point_2d/amanzi_dispersion_aligned_point_2d.rst',
  },
  'dispersion_45_point_2d' : {
    'from_dir' : 'test_suites/verification/transport/saturated/steady-state/dispersion_45_point_2d',
    'dest_dir' : 'doc/user_guide/verification/transport/dispersion_45_point_2d',
    'index_entry' : 'dispersion_45_point_2d/amanzi_dispersion_45_point_2d.rst',
  },
  'dual_porosity_1d' : {
    'from_dir' : 'test_suites/verification/transport/saturated/transient/dual_porosity_1d',
    'dest_dir' : 'doc/user_guide/verification/transport/dual_porosity_1d',
    'index_entry' : 'dual_porosity_1d/amanzi_dual_porosity_1d.rst',
  },
}

#
#  Benchmarks
#
benchmark = {}
benchmark['index'] = {
  'index_title' : 'Benchmark Testing',
  'index_file' : 'doc/user_guide/benchmarking/index.rst',
  'index_list' : [ 'chemistry',
                   'transport'
                 ],
}

benchmark['chemistry'] = {
  'index_entry' : 'chemistry/index.rst',
  'index' : { 'index_title' : 'Chemistry',
              'index_label' : 'sec-benchmarks-chemistry',
              'index_file' : 'doc/user_guide/benchmarking/chemistry/index.rst',
              'index_list' : [ 'tracer_1d',
                               'tritium_1d',
                               'calcite_1d',
                               'isotherms_1d',
                               'ion_exchange_1d',
                               'surface_complexation_1d',
                               'farea_1d',
                             ]
            },
  'tracer_1d' : {
    'from_dir' : 'test_suites/benchmarking/chemistry/tracer_1d',
    'dest_dir' : 'doc/user_guide/benchmarking/chemistry/tracer_1d',
    'index_entry' : 'tracer_1d/amanzi_u-1d-tracer.rst'
  },
  'tritium_1d': {
    'from_dir' : 'test_suites/benchmarking/chemistry/tritium_1d',
    'dest_dir' : 'doc/user_guide/benchmarking/chemistry/tritium_1d',
    'index_entry' : 'tritium_1d/amanzi_u-1d-tritium.rst'
  },
  'calcite_1d' : {
    'from_dir' : 'test_suites/benchmarking/chemistry/calcite_1d',
    'dest_dir' : 'doc/user_guide/benchmarking/chemistry/calcite_1d',
    'index_entry' : 'calcite_1d/amanzi_u-1d-calcite.rst'
  },
  'isotherms_1d' : {
    'from_dir' : 'test_suites/benchmarking/chemistry/isotherms_1d',
    'dest_dir' : 'doc/user_guide/benchmarking/chemistry/isotherms_1d',
    'index_entry' : 'isotherms_1d/amanzi_u-1d-isotherms.rst'
  },
  'ion_exchange_1d' : {
    'from_dir' : 'test_suites/benchmarking/chemistry/ion_exchange_1d',
    'dest_dir' : 'doc/user_guide/benchmarking/chemistry/ion_exchange_1d',
    'index_entry' : 'ion_exchange_1d/amanzi_u-1d-ion-exchange.rst'
  },
  'surface_complexation_1d': {
    'from_dir' : 'test_suites/benchmarking/chemistry/surface_complexation_1d',
    'dest_dir' : 'doc/user_guide/benchmarking/chemistry/surface_complexation_1d',
    'index_entry' : 'surface_complexation_1d/amanzi_u-1d-surface-complexation.rst'
  },
  'farea_1d': {
    'from_dir' : 'test_suites/benchmarking/chemistry/farea_1d',
    'dest_dir' : 'doc/user_guide/benchmarking/chemistry/farea_1d',
    'index_entry' : 'farea_1d/amanzi_u-1d-farea.rst'
  },
}

benchmark['transport'] = {
  'index_entry' : 'transport/index.rst',
  'index' : { 'index_title' : 'Transport',
              'index_file' : 'doc/user_guide/benchmarking/transport/index.rst',
              'index_list' : [ 'single_fracture',
                               'non_grid_aligned',
                             ]
            },
  'single_fracture' : {
    'from_dir' : 'test_suites/benchmarking/coupled_flow_transport/single_fracture',
    'dest_dir' : 'doc/user_guide/benchmarking/transport/single_fracture',
    'index_entry' : 'single_fracture/amanzi_single_fracture.rst'
  },
  'non_grid_aligned' : {
    'from_dir' : 'test_suites/benchmarking/chemistry/non_grid_aligned',
    'dest_dir' : 'doc/user_guide/benchmarking/transport/non_grid_aligned',
    'index_entry' : 'non_grid_aligned/non_grid_aligned.rst'
  },
}

#
#  Personal tests
#
mycase = {}
mycase['index'] = {
  'index_title' : 'My Prototype Documentation',
  'index_file' : 'doc/user_guide/mycase/index.rst',
  'index_list' : [ 'my_test', 
                 ],
}

mycase['my_test'] = {
  'from_dir' : 'test_suites/verification/flow/saturated/steady-state/my_test',
  'dest_dir' : 'doc/user_guide/mycase/my_test',
  'index_entry' : 'linear_head_head/amanzi_my_test.rst'
}

# =======================================================================================

#
#  Create parser and options
#
p = optparse.OptionParser()
p.add_option('--full-guide', help='Build the full User Guide', default=False, dest='full_guide', action='store_true')
p.add_option('--mycase', help='Build the "mycase" test', default=False, dest='mycase', action='store_true')
p.add_option('--install', default=False, dest='install', action='store_true')
p.add_option('--tutorial', default=False, dest='tutorial', action='store_true')
p.add_option('--verification', default=False, dest='verification', action='store_true')
p.add_option('--parallel', default=False, dest='parallel', action='store_true')
p.add_option('--benchmarking', default=False, dest='benchmarking', action='store_true')
p.add_option('--run-tests', default=False, dest='run_tests', action='store_true')

(opts,args) = p.parse_args()

#
#  Create dictionary for sections
#
sections = {}

#
#  Table of Contents (Top Level)
#
toc_user_guide = {
  'index_list' : [ 'background', 'quickstart', 'capabilities', 'input' ],
  'background'   : { 'index_entry' : 'background/index.rst' },
  'quickstart'   : { 'index_entry' : 'background/getting_started.rst'},
  'capabilities' : { 'index_entry' : 'capabilities/index.rst' },
  'input'        : { 'index_entry' : 'input/index.rst'  },
}

if ( opts.install or opts.full_guide ):
    toc_user_guide['index_list'].append('install')
    toc_user_guide['install'] = { 'index_entry' : 'install/index.rst' }
    sections['install'] = install

if ( opts.tutorial or opts.full_guide ):
    toc_user_guide['index_list'].append('tutorial')
    toc_user_guide['tutorial'] = { 'index_entry' : 'tutorial/index.rst' }
    sections['tutorial'] = tutorial

if ( opts.verification or opts.full_guide ):
    toc_user_guide['index_list'].append('verification')
    toc_user_guide['verification'] = { 'index_entry' : 'verification/index.rst' }
    sections['verification'] = verification

if ( opts.benchmarking or opts.full_guide ):
    toc_user_guide['index_list'].append('benchmarking')
    toc_user_guide['benchmarking'] = {'index_entry' : 'benchmarking/index.rst'}
    sections['benchmarking'] = benchmark

if (opts.mycase):
    toc_user_guide['index_list'].append('mycase')
    toc_user_guide['mycase'] = {'index_entry' : 'mycase/index.rst'}
    sections['mycase'] = mycase
    
# =======================================================================================

# Set the logfile
logfile=sys.stdout

# Set Amanzi source directory
amanzi_home=utils.AmanziHome(logfile)

# Set level counter
level=1

# Copy top-level base index file
shutil.copyfile('index.in','index.rst')

# Create index files
utils.RecurseIndex(amanzi_home,sections,level,logfile)

# Copy content 
utils.RecurseCopy(amanzi_home,sections,level,logfile)

# Fix top-level index 
index_entries = [ ]
for e in toc_user_guide['index_list']:
    index_entries.append(toc_user_guide[e]['index_entry'])
utils.IndexInsert('index.rst',index_entries)

# Fix paths on plot directives
utils.WalkRstFiles(amanzi_home,sections,logfile)

# =======================================================================================

# Run the tests
print("\nRunning verification and benchmarking tests...")

mpi_exec = os.getenv('AMANZI_MPI_EXEC', 'mpirun')
mpi_np = os.getenv('AMANZI_MPI_NP', '1')
print(">>> Enviromental variable or default 'AMANZI_INSTALL_DIR' = " 
      + os.getenv('AMANZI_INSTALL_DIR', 'missing'))
print(">>> Enviromental variable or default 'AMANZI_MPI_EXEC' = " + mpi_exec)
print(">>> Enviromental variable or default 'AMANZI_MPI_NP' = " + mpi_np)
print

suffices = {"", "-a", "-b", "-c"}

if ( opts.verification or opts.full_guide and opts.run_tests):

    for name in verification['index']['index_list']:
        for key in verification[name]['index']['index_list']:
            cwd = os.getcwd()
            run_directory= os.path.join(cwd, 'verification/' + name + '/' + key)
            print("\nTEST: " + run_directory)
            os.chdir(run_directory)

            for s in suffices:
                filename = "amanzi_" + key + s + ".py"
                if os.path.isfile(filename):
                    stdout_file = open("collect.log", "w")
                    print("   script: " + filename)
                    subprocess.call(["python", filename], stdout=stdout_file, stderr=subprocess.STDOUT)
                    
                    for el in find("stdout.out", run_directory):
                        el2 = el[el.find("verification"):]
                        print("   file: " + el2)
                        if ("Amanzi::SIMULATION_SUCCESSFUL" in open(el).read()):
                            print("   result: SIMULATION_SUCCESSFUL")
                        else:
                            print("   ERROR: " + el)
            os.chdir(cwd)

if ( opts.benchmarking or opts.full_guide and opts.run_tests):
            
    for name in benchmark['index']['index_list']:
        for key in benchmark[name]['index']['index_list']:
            cwd = os.getcwd()
            run_directory= os.path.join(cwd, 'benchmarking/' + name + '/' + key)
            print("\nTEST: " + run_directory)
            os.chdir(run_directory)

            for s in suffices:
                filename = key + s + ".py"
                if os.path.isfile(filename):
                    stdout_file = open("collect.log", "w")
                    print("   script: " + filename)
                    subprocess.call(["python", filename], stdout=stdout_file, stderr=subprocess.STDOUT)
                    
                    for el in find("stdout.out", run_directory):
                        el2 = el[el.find("benchmarking"):]
                        print("   file: " + el2)
                        if ("Amanzi::SIMULATION_SUCCESSFUL" in open(el).read()):
                            print("   result: SIMULATION_SUCCESSFUL")
                        else:
                            print("   ERROR:" + open(el).readline())
            os.chdir(cwd)
