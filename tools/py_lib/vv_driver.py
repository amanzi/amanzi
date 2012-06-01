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

import amanzi

# --- Parse command line
parser = OptionParser()

# Input file
parser.add_option("-i", "--input-file", dest="input", 
                  help="Amanzi input file",metavar="FILE")

# Number of processors
parser.add_option("-n", "--nprocs", dest="nprocs", 
                  help="Number of processors",metavar="NUM")

# HDF5 Difference Binary (h5diff)
parser.add_option("--h5diff", dest="h5diff",
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

# --- Define the output file patterns

# Grab the sublists in the input file
exec_ctrl = input_tree.find_sublist('Execution Control')
output_ctrl = input_tree.find_sublist('Output')
viz_ctrl = output_ctrl.find_sublist('Visualization Data')

# File base name and number of digits
output_basename=viz_ctrl.find_parameter('File Name Base')
output_ndigits=viz_ctrl.find_parameter('File Name Digits')

# Build the regular expression for globbing
output_regex=output_basename.get('value')
i=1
if i > int(output_ndigits.get('value')):
  print 'yes'
while i <= int(output_ndigits.get('value')):
  output_regex = output_regex + '[0-9]'
  i = i + 1
print 'Search for output files matching \'' + output_regex + '\''

# Search for files
output_files=glob.glob(output_regex)
if len(output_files) == 0:
  print 'No output files found'


sys.exit(0)
