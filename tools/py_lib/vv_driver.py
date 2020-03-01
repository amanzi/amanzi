#!/usr/bin/env python

# ############################################################################ #
#
# Amanzi
#   Verification and Validation Test Suite 
#
# ############################################################################ #

import os, sys
import glob 
from optparse import OptionParser
from xml.etree import ElementTree as ETree


try:
  import amanzi as Amanzi
except ImportError:
  amanzi_python_install_prefix='/usr/local/lib/python2.6/site-packages'
  sys.path.append(amanzi_python_install_prefix)
  import amanzi as Amanzi

try:
  import subprocess as process
except ImportError:
  raise ImportError('Script requires the subprocess module found in Python 2.4 and higher')


# --- Parse command line
parser = OptionParser()

# Binary file
amanzi_dflt_binary='/usr/local/bin/amanzi'
parser.add_option("-b", "--binary", dest="binary", default=amanzi_dflt_binary, 
                  help="Amanzi binary file",metavar="FILE")

# Input file
parser.add_option("-i", "--input", dest="input", 
                  help="Amanzi input file",metavar="FILE")

# Unscramble Viz binary file
unscram_viz_dflt_binary='/usr/local/bin/unscramble_viz'
parser.add_option("--unscramble-viz", dest="unscramble_viz", default=unscram_viz_dflt_binary, 
                  help="Amanzi unscramble viz files binary ",metavar="FILE")

# Unscramble Restart binary file
unscram_restart_dflt_binary='/usr/local/bin/unscramble_restart'
parser.add_option("--unscramble-restart", dest="unscramble_restart", default=unscram_restart_dflt_binary, 
                  help="Amanzi unscramble restart files binary ",metavar="FILE")

# STDOUT file
parser.add_option("-o", "--output", dest="output",
                  help="Redirect STDOUT to file", metavar="FILE")

# STDERR file
parser.add_option("-e", "--error", dest="error",
                  help="Redirect STDERR to file", metavar="FILE")

# Clobber old output files
parser.add_option("--clobber", action="store_true", dest="clobber",
                  default=False, help="Remove old output files")

# Number of processors
parser.add_option("-n", "--nprocs", dest="nprocs", 
                  help="Number of processors",metavar="NUM")

# HDF5 Difference Binary (h5diff)
h5diff_dflt_binary='/local/lpritch/ASCEM/rc2-2012/tpl-new/bin/h5diff'
parser.add_option("--h5diff", dest="h5diff", default=h5diff_dflt_binary,
                  help="HDF5 difference tool", metavar="FILE")

# HDF5 List Binary (h5ls)
h5ls_dflt_binary='/local/lpritch/ASCEM/rc2-2012/tpl-new/bin/h5ls'
parser.add_option("--h5ls", dest="h5ls", default=h5ls_dflt_binary,
                  help="HDF5 ls tool", metavar="FILE")

# HDF5 Dump Binary (h5dump)
h5dump_dflt_binary='/local/lpritch/ASCEM/rc2-2012/tpl-new/bin/h5dump'
parser.add_option("--h5dump", dest="h5dump", default=h5dump_dflt_binary,
                  help="HDF5 data dump tool", metavar="FILE")

# HDF5 Copy Binary (h5copy)
h5copy_dflt_binary='/local/lpritch/ASCEM/rc2-2012/tpl-new/bin/h5copy'
parser.add_option("--h5copy", dest="h5copy", default=h5copy_dflt_binary,
                  help="HDF5 data copy tool", metavar="FILE")

# Extract Dataset
parser.add_option("--extract-data", dest="extract_data",
                  help="Dataset name to extract", metavar="NAME")

# Extract Time
parser.add_option("--extract-time", dest="extract_time", default='last',
                  help="Extract data at time T", metavar="TIME")

# V&V output results file
parser.add_option("--vv-results", dest="vv_results", default='vv_results.h5',
                  help="Output V&V Results to FILE", metavar="FILE")

# Compare file
parser.add_option("--compare-file", dest="compare_file",
                  help="Compare output results to data found in FILE", metavar="FILE")

# Compare object
parser.add_option("--compare-dataset", dest="compare_dataset",
                 help="Compare V&V results against this dataset STRING", metavar="STRING")

# Compare tolerance (relative)
parser.add_option("--compare-tol", dest="compare_tol", default=0.05,
                  help="Run h5diff with relative tolerance set to NUM", metavar="NUM")

# --- Process the arguments
(options,args) = parser.parse_args()

# Check the number of processors
if options.nprocs is None:
  print('Setting number of processors to default value of 1')
  options.nprocs=1

# Check the input file
if options.input is None:
  raise NameError('Must define an input file')

if not os.path.exists(options.input):
  raise ValueError(options.input + ' does not exist')

input_tree=Amanzi.trilinos.InputList(options.input)

# Check h5copy
if options.h5copy is None:
  print('Searching for h5copy')
  search_cmd = Amanzi.command.Command('which','h5copy')
  if search_cmd.exit_code != 0:
    raise RuntimeError('Failed to find h5copy')
  options.h5copy=search_cmd.output
else:
  if not os.path.exists(options.h5copy):
    raise ValueError(options.h5copy + ' does not exist')

# --- Pre-run activities

# - Define the output file patterns

# Grab the sublists in the input file
output_ctrl = input_tree.find_sublist('Output')
viz_ctrl    = output_ctrl.find_sublist('Visualization Data')

# Output file base name 
output_basename=viz_ctrl.find_parameter('File Name Base').get_value()

#  - Remove old output files
if options.clobber:
  print('Remove old checkpoint, plot and V&V result files')
  for old_checkpoint in glob.glob('checkpoint*.h5'):
    try:
      os.remove(old_checkpoint)
    except:
      print('Failed to delete ' + old_checkpoint)
      raise
  plot_regex=output_basename+'_*'
  for old_plot_file in glob.glob(plot_regex):
    try:
      os.remove(old_plot_file)
    except:
      print('Failed to delete ' + old_plot_file)
      raise
  if os.path.exists(options.vv_results):
    try:
      os.remove(options.vv_results)
    except:
      print('Failed to delete ' + options.vv_results)
      raise

# --- Run amanzi
print('>>>>>>>>> LAUNCHING AMANZI <<<<<<<<')
amanzi=Amanzi.interface.AmanziInterface(options.binary)
amanzi.input=options.input
amanzi.output=options.output
amanzi.error=options.error
amanzi.nprocs=options.nprocs
amanzi.run()
print('Amanzi exit code ' + str(amanzi.exit_code))
print('>>>>>>>>> Amanzi RUN COMPLETE <<<<<<<<')

# --- Process output


# - Locate output files

# Mesh file
mesh_file=amanzi.find_mesh_file()

# Viz Files
amanzi.find_data_files(output_basename)
if len(amanzi.data_files) == 0:
  print('No XDMF XML files found')
else:
  print('Found ' + str(len(amanzi.data_files)) + ' XMF files.')

# - Find the correct file to data mine
if options.extract_data != None:
  print('>>>>>>>> Processing Amanzi Output <<<<<<<<')
  print('Searching for dataset "' + options.extract_data + '"')
  search_files = []
  for output in amanzi.data_files:
    if options.extract_data in output.datasets:
      search_files.append(output)
  print('Found ' + str(len(search_files)) + ' files containing dataset "' + options.extract_data + '"')    
  
  source_file=None
  if options.extract_time == 'last':
    sorted_list=amanzi.sort_data_files(files=search_files)

    try:
      source_file=sorted_list[-1]
    except:
      print('Search file list is empty')
  else:
    idx=0
    while source_file == None and idx < len(search_files):
      if float(search_files[idx].time) == float(options.extract_time):
        source_file=search_files[idx]
      idx=idx+1

  if source_file != None:
    if options.nprocs > 1:
      extract_source_file=output_basename+'_data_unscrambled.h5'
      serialize_args=[options.unscramble_viz]
      serialize_args.append(mesh_file)
      serialize_args.append(output_basename+'_data.h5')
      serialize_args.append(extract_source_file)
      try:
        serialize=process.Popen(serialize_args)
      except:
        print('Failed to serialize the output data')
        raise
      else:
        serialize.wait()
        print('Serialize command returned ' + str(serialize.returncode))
    else:
      extract_source_file=output_basename+'_data.h5'
    print('Will extract "' + options.extract_data + '" from ' + source_file.filename)
    group_name=options.extract_data+'/'+str(source_file.cycle)
    h5copy_args = [options.h5copy]
    h5copy_args.append('--input='+extract_source_file)
    h5copy_args.append('--source='+group_name)
    h5copy_args.append('--output='+options.vv_results)
    h5copy_args.append('--destination='+options.extract_data)
    try:
      h5copy=process.Popen(h5copy_args)
    except:
      print('Failed to extract "' +options.extract_data + '"')
    else:
      h5copy.wait()
      print('h5copy returned ' + str(h5copy.returncode))
  
    # Compare data if data file is set
    if options.compare_file != None:
      try:
        print('Will compare dataset ' + options.compare_dataset + ' found in ' + options.compare_file)
      except:
        print('Please define a dataset name to compare')
        raise

      h5diff_args=[options.h5diff]
      h5diff_args.append('--relative='+str(options.compare_tol))
      h5diff_args.append(options.compare_file)
      h5diff_args.append(options.vv_results)
      h5diff_args.append(options.compare_dataset)
      h5diff_args.append(options.extract_data)
      try:
        h5diff=process.Popen(h5diff_args)
      except:
        print('Failed to open h5diff process pipe')
      else:
        h5diff.wait()
        print('h5diff returned ' + str(h5diff.returncode))

      if h5diff.returncode == 0:
        print('Results passed h5diff test with relative tolerance set to ' + str(options.compare_tol))
      elif h5diff.returncode == 1:
        print('Results FAILED h5diff test with relative tolerance set to ' + str(options.compare_tol))
      else:
        print(h5diff_args)
        print('h5diff failed')

  else:
    print('Failed to locate dataset ' + options.extract_data + ' at time ' + options.extract_time)

  print('>>>>>>>> Processing Amanzi Output COMPLETE <<<<<<<<')

sys.exit(amanzi.exit_code)
