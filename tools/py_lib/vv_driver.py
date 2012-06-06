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
  import amanzi
except ImportError:
  amanzi_python_install_prefix='@AmanziPython_INSTALL_PREFIX@'
  sys.path.append(amanzi_python_install_prefix)
  import amanzi

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
                  help="HDF5 difference tool", metavar="FILE")


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

input_tree=amanzi.trilinos.InputList(options.input)

# Check h5diff
if options.h5diff is None:
  print 'Searching for h5diff'
  search_cmd = amanzi.command.Command('which','h5diff')
  if search_cmd.exit_code != 0:
    raise RuntimeError, 'Failed to find h5diff'
  options.h5diff=search_cmd.output
else:
  if not os.path.exists(options.h5diff):
    raise ValueError, options.h5diff + ' does not exist'


# --- Run amanzi
print '>>>>>>>>> LAUNCHING AMANZI <<<<<<<<'
amanzi=amanzi.interface.AmanziInterface(options.binary)
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

# File base name 
output_basename=viz_ctrl.find_parameter('File Name Base')

# Build the regular expression for globbing
output_regex=output_basename.get_value()
output_regex=output_regex+'_data.h5'
output_xmf_regex=output_basename.get_value()+'_data.h5'+'.'+'[0-9]*'+'.xmf'

# - Create the h5diff command - will add options and files later
#h5diff_cmd=amanzi.command.CommandInterface(options.h5diff)

# - Locate output files
print 'Search for output files matching pattern \'' + output_regex + '\''
output_files=glob.glob(output_regex)
print 'Search for XDMF XML files matching pattern \'' + output_xmf_regex + '\''
xmf_files=glob.glob(output_xmf_regex)
if len(output_files) == 0:
  print 'No output files found'
else:
  for file in output_files:
    print 'Will process output file:' + file

if len(xmf_files) == 0:
  print 'No XDMF XML files found'
else:
  for file in xmf_files:
    #print 'Will process XDMF XML output file:' + file
    xmf_tree=ETree.parse(file)
    time_element=xmf_tree.getiterator(tag='Time')[0]
    time=time_element.get('Value')
    print 'File ' + file + ' dumped out at time=' + time + ' (years)'


# - Call h5diff to compare files    
print '>>>>>>>> Processing Amanzi Output COMPLETE <<<<<<<<'

sys.exit(amanzi.exit_code)
