#!/usr/bin/env python
'''
Simple2D XML Generator
'''
import os, sys
from optparse import OptionParser

# Amanzi Python Modules
# set environment variable PYTHONPATH <amnzi source tree location>/tools/py_lib
from amanzi.trilinos import *
from amanzi.mpc import *


# Define and parse the command line options
parser = OptionParser()
parser.add_option("-o", "--output", dest="outfile",
	          help="write XML to FILE", metavar="FILE")

(options,args) = parser.parse_args()

# ################################### #
# Mesh                                #
# ################################### #
# Setup a simple mesh
mesh = ParameterList("Mesh")
mesh.add_parameter("Framework","Simple");
mesh_generator = mesh.add_sublist("Generate")

# Add the number of cells and domain definitions
labels = ['X', 'Y', 'Z']
ncells = [10, 10, 3]
domain_min = [-1.0, -1.0, -1.0]
domain_max = [1.0, -0.075, -0.075]

# Most this will be replaced with a 
# call mesh_generator = SimpleMesh(ncells,domain_min,domain_max)
idx = 0
while idx < 3:
    keyword = 'Number of Cells in ' + labels[idx]
    mesh_generator.add_parameter(keyword,ncells[idx])
    idx = idx + 1

idx = 0
while idx < 3:
    min_keyword = labels[idx] + '_Min'
    max_keyword = labels[idx] + '_Max'
    mesh_generator.add_parameter(min_keyword,domain_min[idx])
    mesh_generator.add_parameter(max_keyword,domain_max[idx])
    idx = idx + 1

# ################################### #
# MPC                                 #
# ################################### #
#  Create the MPC
mpc = MPC()

mpc.set_end_time(86400)
mpc.enable_flow = bool(False)
mpc.enable_transport = bool(True)
mpc.enable_chemistry = bool(True)


# Setup the viz input
mpc.viz.set_file('simple2d.cgns')
mpc.viz.set_dt(86400)

# ################################### #
# State                               #
# ################################### #
# Define the State
state = ParameterList("State")
state.add_parameter("Number of mesh blocks",1)
state.add_parameter("Number of component concentrations",1)

# Water parameters should be under a simple ParameterList
state.add_parameter("Constant water saturation",1.0)
state.add_parameter("Constant water density",1000.0)
state.add_parameter("Constant vicosity",0.001308)

# Need to convert this to a vector not three paramters
state.add_parameter("Gravity x", 0.0)
state.add_parameter("Gravity y", 0.0)
state.add_parameter("Gravity z", -9.8)

# Mesh blocks
mesh_block1 = state.add_sublist("Mesh Block 1")
mesh_block1.set_parameter("Mesh Block ID",0)
mesh_block1.set_parameter("Constant porosity",0.4)
mesh_block1.set_parameter("Constant permeability",1.0e-3)
mesh_block1.set_parameter("Constant component concentration 0",0.0)


# ################################### #
# Flow                                #
# ################################### #
flow = ParameterList("Flow")
flow.add_parameter("Max Iterations", 25)
flow.add_parameter("Error Tolerance", 1.0e-20)

# Define a Darcy Flow Problem
darcy = flow.add_sublist("Darcy Problem")
diff_precon = darcy.add_sublist("Diffusion Preconditioner")
ml = diff_precon.add_sublist("ML Parameters")

# Define ML Parameters
ml.add_parameter("max levels",40)
ml.add_parameter("prec type","MGV")
ml.add_parameter("increasing or decreasing","increasing")
ml.add_parameter("aggregation: type","Uncoupled-MIS")
ml.add_parameter("aggregation: damping factor",1.333)
ml.add_parameter("aggregation: threshold",0.03)
ml.add_parameter("eigen-analysis: type","cg")
ml.add_parameter("eigen-analysis: iterations",20)
ml.add_parameter("smoother: sweeps",5)
ml.add_parameter("smoother: damping factor",1.0)
ml.add_parameter("smoother: pre or post","both")
ml.add_parameter("smoother: type","symmetric Gauss-Seidel")
ml.add_parameter("coarse: type","Amesos-KLU")
ml.add_parameter("coarse: max size", 128)

# Flow boundary conditions
flow_bc = flow.add_sublist("Flow BC")
flow_bc.add_parameter("number of BC's",6)

bc_types  = ['Static Head', 'Static Head', 'Static Head', 'Static Head', 'No Flow', 'No Flow']
bc_values = [0.0, 0.0, 1.0, 0.0, 'skip', 'skip']
bc_ids    = [3, 1, 0, 2, 4, 5]

idx = 0
while idx < 6:
    bc_label = "BC%02d" % (idx)
    bc = flow_bc.add_sublist(bc_label)
    bc.add_parameter('type', bc_types[idx])
    if bc_types[idx] == 'Static Head':
	bc.add_parameter('BC value', bc_values[idx])

    bc.add_parameter('Side set ID', bc_ids[idx])	
    idx = idx + 1

	   



# ################################### #
# Generate Input                      #
# ################################### #

# Finally we assemble the input list
main = ParameterList("Main")
main.add_sublist(mesh)
main.add_sublist(mpc)
main.add_sublist(state)
main.add_sublist(flow)


# And print out
if options.outfile != None:
    ostream = file(options.outfile,mode='w+')
else:
    ostream = sys.stdout

main.dumpXML(file=ostream)
