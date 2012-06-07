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
  amanzi_python_install_prefix='@AmanziPython_INSTALL_PREFIX@'
  sys.path.append(amanzi_python_install_prefix)
  import amanzi as Amanzi

# --- Parse command line
parser = OptionParser()

# Binary file
amanzi_dflt_binary='@Amanzi_EXECUTABLE@'
parser.add_option("-b", "--binary", dest="binary", default=amanzi_dflt_binary, 
                  help="Amanzi binary file",metavar="FILE")

# Input file
parser.add_option("-i", "--input", dest="input", 
                  help="Amanzi input file",metavar="FILE")

# STDOUT file
parser.add_option("-o", "--output", dest="output",
                  help="Redirect STDOUT to file", metavar="FILE")

# STDERR file
parser.add_option("-e", "--error", dest="error",
                  help="Redirect STDERR to file", metavar="FILE")

# Number of processors
parser.add_option("-n", "--nprocs", dest="nprocs", 
                  help="Number of processors",metavar="NUM")

# HDF5 Difference Binary (h5diff)
h5diff_dflt_binary='@HDF5_H5DIFF_BINARY@'
parser.add_option("--h5diff", dest="h5diff", default=h5diff_dflt_binary,
                  help="HDF5 difference tool", metavar="FILE")

# HDF5 List Binary (h5ls)
h5ls_dflt_binary='@HDF5_H5LS_BINARY@'
parser.add_option("--h5ls", dest="h5ls", default=h5ls_dflt_binary,
                  help="HDF5 ls tool", metavar="FILE")

# HDF5 Dump Binary (h5dump)
h5dump_dflt_binary='@HDF5_H5DUMP_BINARY@'
parser.add_option("--h5dump", dest="h5dump", default=h5dump_dflt_binary,
                  help="HDF5 data dump tool", metavar="FILE")

# HDF5 Copy Binary (h5copy)
h5copy_dflt_binary='@HDF5_H5COPY_BINARY@'
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

# --- Process the arguments
(options,args) = parser.parse_args()

# Check the number of processors
if options.nprocs is None:
  print 'Setting number of processors to default value of 1'
  options.nprocs=1

# Check the input file
if options.input is None:
  raise NameError, 'Must define an input file'

if not os.path.exists(options.input):
  raise ValueError, options.input + ' does not exist'

input_tree=Amanzi.trilinos.InputList(options.input)

# Check h5copy
if options.h5copy is None:
  print 'Searching for h5copy'
  search_cmd = Amanzi.command.Command('which','h5copy')
  if search_cmd.exit_code != 0:
    raise RuntimeError, 'Failed to find h5copy'
  options.h5copy=search_cmd.output
else:
  if not os.path.exists(options.h5copy):
    raise ValueError, options.h5copy + ' does not exist'


# --- Run amanzi
print '>>>>>>>>> LAUNCHING AMANZI <<<<<<<<'
amanzi=Amanzi.interface.AmanziInterface(options.binary)
amanzi.input=options.input
amanzi.output=options.output
amanzi.error=options.error
amanzi.nprocs=options.nprocs
amanzi.run()
print 'Amanzi exit code ' + str(amanzi.exit_code)
print '>>>>>>>>> Amanzi RUN COMPLETE <<<<<<<<'

# --- Process output

# - Define the output file patterns

# Grab the sublists in the input file
output_ctrl = input_tree.find_sublist('Output')
viz_ctrl    = output_ctrl.find_sublist('Visualization Data')

# Output file base name 
output_basename=viz_ctrl.find_parameter('File Name Base').get_value()

# - Locate output files
amanzi.find_data_files(output_basename)
if len(amanzi.data_files) == 0:
  print 'No XDMF XML files found'
else:
  print 'Found ' + str(len(amanzi.data_files)) + ' XMF files.'
  #for data_file in amanzi.data_files:
  #  print 'File:'+data_file.filename
  #  print '\tcycle='+str(data_file.cycle)
  #  print '\ttime='+str(data_file.time)
  #  print '\tdatasets='+str(data_file.datasets)

# - Find the correct file to data mine
if options.extract_data != None:
  print '>>>>>>>> Processing Amanzi Output <<<<<<<<'
  print 'Searching for dataset ' + options.extract_data
  search_files = []
  for output in amanzi.data_files:
    if options.extract_data in output.datasets:
      search_files.append(output)
  
  source_file=None
  if options.extract_time == 'last':
    source_file=search_files[-1]
  else:
    idx=0
    while source_file == None and idx < len(search_files):
      if float(search_files[idx].time) == float(options.extract_time):
        source_file=search_files[idx]
      idx=idx+1

  if source_file != None:
    print 'Will extract ' + options.extract_data + ' from ' + source_file.filename
    root_data_file=output_basename+'_data.h5'
    group_name=options.extract_data+'/'+str(source_file.cycle)
    h5copy_args = ['--parents']
    h5copy_args.append('--input='+root_data_file)
    h5copy_args.append('--source='+group_name)
    h5copy_args.append('--output='+options.vv_results)
    h5copy_args.append('--destination='+options.extract_data)
    h5copy_cmd=Amanzi.command.Command(options.h5copy,h5copy_args)
    if h5copy_cmd.exit_code != 0:
      h5copy_cmd._dump_state()
      raise RuntimeError, 'Failed to create V&V results file:'+options.vv_results
  else:
    print 'Failed to locate dataset ' + options.extract_data + ' at time ' + options.extract_time

  print '>>>>>>>> Processing Amanzi Output COMPLETE <<<<<<<<'

sys.exit(amanzi.exit_code)
