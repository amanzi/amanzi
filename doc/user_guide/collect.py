#!/usr/bin/env python

import os, sys, utils
import optparse

#  Create dictionary that describes:
#
#  - layout of directories
#  - index files to be created
#  - subdirectories (tests/tutorials) to be copied
#

#
# Install
#
install={}
install['index']={'index_title' : 'Building Amanzi',
                   'index_file' : 'doc/user_guide/install/index.rst',
                   'index_list' : ['bootstrap','tpls','amanzi','quickstart'],
                   }
    
   
install['tpls']={'from_file' : 'config/SuperBuild/INSTALL',
                 'dest_file' : 'doc/user_guide/install/building_tpls.rst',
                 'index_entry' : 'building_tpls.rst',
                 }

install['amanzi']={'from_dir' : 'doc/build_instructions',
                   'dest_dir' : 'doc/user_guide/install',
                   'index_entry' : 'build_instructions.rst',
                  }
install['bootstrap']={'index_entry' : 'building_bootstrap.rst', }
install['quickstart']={'index_entry' : 'cmake_quickstart.rst', }

    
#
# Tutorials
#
tutorial={}
tutorial['index']={'index_title' : 'Tutorial',
                   'index_file' : 'doc/user_guide/tutorial/index.rst',
                   'index_list' : ['flow_1','flow_2'],
                   }
    
tutorial['flow_1']={'from_dir' : 'doc/tutorial/saturated',
                    'dest_dir' : 'doc/user_guide/tutorial/steady_confined',
                    'index_entry' : 'steady_confined/Tutorial_steady_linear.rst',
                    }
    
tutorial['flow_2']={'from_dir' : 'doc/tutorial/richards',
                    'dest_dir' : 'doc/user_guide/tutorial/transient_infiltration',
                    'index_entry' : 'transient_infiltration/Tutorial_dvz.rst',
                    }

#
# Verificaiton Tests
#
verification={}
verification['index']={'index_title':'Verification Testing',
                       'index_file':'doc/user_guide/verification/index.rst',
                       'index_list':['confined_flow','unconfined_flow','transport'],
}

verification['confined_flow']={'index_entry' : 'confined_flow/index.rst',
                               'index' : 
                               {'index_title' : 'Confined Flow Tests',
                                'index_file' : 'doc/user_guide/verification/confined_flow/index.rst',
                                'index_list' : ['linear_head_head', 'linear_flux_head',
                                                'linear_materials_serial','linear_materials_parallel',
                                                'theis_isotropic', 'hantush_anisotropic',
                                                ],
                                },
                               'linear_head_head' :
                                   {'from_dir' : 'testing/verification/flow/saturated/steady-state/linear_head_head_1d',
                                    'dest_dir' : 'doc/user_guide/verification/confined_flow/linear_head_head',
                                    'index_entry' : 'linear_head_head/amanzi_linear_head_head_1d.rst'
                                    },
                               'linear_flux_head' :
                                   {'from_dir' : 'testing/verification/flow/saturated/steady-state/linear_flux_head_1d',
                                    'dest_dir' : 'doc/user_guide/verification/confined_flow/linear_flux_head',
                                    'index_entry' : 'linear_flux_head/amanzi_linear_flux_head_1d.rst'
                                    },
                               'linear_materials_serial' :
                                   {'from_dir' : 'testing/verification/flow/saturated/steady-state/linear_materials_serial_1d',
                                    'dest_dir' : 'doc/user_guide/verification/confined_flow/linear_materials_serial',
                                    'index_entry' : 'linear_materials_serial/amanzi_linear_materials_serial_1d.rst'
                                    },
                               'linear_materials_parallel' :
                                   {'from_dir' : 'testing/verification/flow/saturated/steady-state/linear_materials_parallel_1d',
                                    'dest_dir' : 'doc/user_guide/verification/confined_flow/linear_materials_parallel',
                                    'index_entry' : 'linear_materials_parallel/amanzi_linear_materials_parallel_1d.rst'
                                    },
                               'theis_isotropic' :
                                   {'from_dir' : 'testing/verification/flow/saturated/transient/theis_isotropic_1d',
                                    'dest_dir' : 'doc/user_guide/verification/confined_flow/theis_isotropic',
                                    'index_entry' : 'theis_isotropic/amanzi_theis_isotropic_1d.rst',
                                    },
                               'hantush_anisotropic' :
                                   {'from_dir' : 'testing/verification/flow/saturated/transient/hantush_anisotropic_2d',
                                    'dest_dir' : 'doc/user_guide/verification/confined_flow/hantush_anisotropic_2d',
                                    'index_entry' : 'hantush_anisotropic_2d/amanzi_hantush_anisotropic_2d.rst',
                                    },
                               }


verification['unconfined_flow']={'index_entry': 'unconfined_flow/index.rst',
                                 'index' : 
                                 {'index_title' : 'Unconfined Flow Tests',
                                  'index_file' : 'doc/user_guide/verification/unconfined_flow/index.rst',
                                  'index_list' : ['unconfined_no_recharge','unconfined_layered','unconfined_seepage'], 
                                  },
                                 'unconfined_no_recharge':
                                     {'from_dir' : 'testing/verification/flow/richards/steady-state/unconfined_no_recharge_1d',
                                      'dest_dir' : 'doc/user_guide/verification/unconfined_flow/unconfined_no_recharge',
                                      'index_entry' : 'unconfined_no_recharge/amanzi_unconfined_no_recharge_1d.rst'
                                      },
                                 'unconfined_recharge':
                                     {'from_dir' : 'testing/verification/flow/richards/steady-state/unconfined_recharge_1d',
                                      'dest_dir' : 'doc/user_guide/verification/unconfined_flow/unconfined_recharge',
                                      },
                                 'unconfined_layered':
                                     {'from_dir' : 'testing/verification/flow/richards/steady-state/unconfined_layered_2d',
                                      'dest_dir' : 'doc/user_guide/verification/unconfined_flow/unconfined_layered',
                                      'index_entry' : 'unconfined_layered/amanzi_unconfined_layered_2d.rst'
                                      },
                                 'unconfined_seepage':
                                     {'from_dir' : 'testing/verification/flow/richards/steady-state/unconfined_seepage_1d',
                                      'dest_dir' : 'doc/user_guide/verification/unconfined_flow/unconfined_seepage',
                                      'index_entry' : 'unconfined_seepage/amanzi_unconfined_seepage_1d.rst',
                                      },
                                 }

verification['transport']={'index_entry': 'transport/index.rst',
                           'index' : 
                           {'index_title' : 'Transport of Solutes: Advection, Dispersion and Diffusion',
                            'index_file' : 'doc/user_guide/verification/transport/index.rst',
                            'index_list' : ['dispersion_aligned_point_2d'],
                            },
                           'dispersion_aligned_point_2d' :
                               {'from_dir' : 'testing/verification/transport/saturated/steady-state/dispersion_aligned_point_2d',
                                'dest_dir' : 'doc/user_guide/verification/transport/dispersion_aligned_point_2d',
                                'index_entry' : 'dispersion_aligned_point_2d/amanzi_dispersion_aligned_point_2d.rst',
                                },
                           }

#
#  Benchmarks
#
benchmark={}
benchmark['index']={'index_title' : 'Benchmark Testing',
                    'index_file' : 'doc/user_guide/benchmarking/index.rst',
                    'index_list' : ['chemistry'],
                }

benchmark['chemistry']={'index_entry' : 'chemistry/index.rst',
                        'index' : 
                        {'index_title' : 'Chemistry',
                         'index_file' : 'doc/user_guide/benchmarking/chemistry/index.rst',
                         'index_list' : ['tracer', 'tritium', 'calcite', 'isotherms', 'ion_exchange', 'farea']
                         },
                        'tracer' : 
                            {'from_dir' : 'testing/benchmarking/chemistry/tracer_1d',
                             'dest_dir' : 'doc/user_guide/benchmarking/chemistry/tracer_1d',
                             'index_entry' : 'tracer_1d/amanzi_u-1d-tracer.rst'
                             },
                        'calcite' : 
                            {'from_dir' : 'testing/benchmarking/chemistry/calcite_1d',
                             'dest_dir' : 'doc/user_guide/benchmarking/chemistry/calcite_1d',
                             'index_entry' : 'calcite_1d/amanzi_u-1d-calcite.rst'
                             },
                        'isotherms' : 
                            {'from_dir' : 'testing/benchmarking/chemistry/isotherms_1d',
                             'dest_dir' : 'doc/user_guide/benchmarking/chemistry/isotherms_1d',
                             'index_entry' : 'isotherms_1d/amanzi_u-1d-isotherms.rst'
                             },
                        'ion_exchange' : 
                             {'from_dir' : 'testing/benchmarking/chemistry/ion_exchange_1d',
                              'dest_dir' : 'doc/user_guide/benchmarking/chemistry/ion_exchange_1d',
                              'index_entry' : 'ion_exchange_1d/amanzi_u-1d-ion-exchange.rst'
                              },
                        'tritium':
                            {'from_dir' : 'testing/benchmarking/chemistry/tritium_1d',
                             'dest_dir' : 'doc/user_guide/benchmarking/chemistry/tritium_1d',
                             'index_entry' : 'tritium_1d/amanzi_u-1d-tritium.rst'
                             },
                        'farea':
                            {'from_dir' : 'testing/benchmarking/chemistry/farea_1d',
                             'dest_dir' : 'doc/user_guide/benchmarking/chemistry/farea_1d',
                             'index_entry' : 'farea_1d/amanzi_u-1d-farea.rst'
                             },
                        }


#
#  Benchmarks
#
mycase={}
mycase['index']={'index_title' : 'My Prototype Documentation',
                 'index_file' : 'doc/user_guide/mycase/index.rst',
                 'index_list' : ['newcase'],
             }

mycase['newcase']={'from_dir' : 'testing/verification/flow/saturated/transient/hantush_anisotropic_2d',
                   'dest_dir' : 'doc/user_guide/mycase/hantush_anisotropic_2d',
                   'index_entry' : 'hantush_anisotropic_2d/amanzi_hantush_anisotropic_2d.rst',
               }


# =========================================================================================================================
#
#  Create parser and options
#
p = optparse.OptionParser()
p.add_option('--full-guide', help='Build the full User Guide', default=False, dest='full_guide', action='store_true')
p.add_option('--mycase', help='Build the "mycase" test', default=False, dest='mycase', action='store_true')
p.add_option('--install', default=False, dest='install', action='store_true')
p.add_option('--tutorial', default=False, dest='tutorial', action='store_true')
p.add_option('--verification', default=False, dest='verification', action='store_true')
p.add_option('--benchmarking', default=False, dest='benchmarking', action='store_true')

(opts,args) = p.parse_args()

#
#  Create dictionary for sections
#
sections={}

#
#  Table of Contents (Top Level)
#
toc_user_guide = {'index_list' : [ 'background', 'quickstart', 'capabilities', 'input' ],
                  'background'   : { 'index_entry' : 'background/index.rst' },
                  'quickstart'   : { 'index_entry' : 'install/getting_started.rst'},
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
    
# =========================================================================================================================

# Set the logfile
logfile=sys.stdout

# Set Amanzi source directory
amanzi_home=utils.AmanziHome(logfile)

# Set level counter
level=1

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

