#!/usr/bin/env python
'''
BC Cribs XML Generator
'''
import os, sys
from optparse import OptionParser

# Amanzi Python Modules
# set environment variable PYTHONPATH <amnzi source tree location>/tools/py_lib
from amanzi.trilinos import *
from amanzi.mpc      import *


class WaterModel(ParameterList):

   def __init__(self,num=None):
        name = "Water retention model %d" % (num)
        ParameterList.__init__(self,name)
        self.add_parameter("Water retention model", "van Genuchten")
        
   def set_id(self,id=None):
       return self.add_parameter('Region ID',id)
   
   def van_genuchten_m(self,m=None):
       return self.add_parameter('van Genuchten m',m)

   def van_genuchten_alpha(self,alpha=None):
       return self.add_parameter('van Genuchten alpha',alpha)

   def van_genuchten_res_saturation(self,sat=None):
       return self.add_parameter('van Genuchten residual saturation',sat)


# Define and parse the command line options
parser = OptionParser()
parser.add_option("-o", "--output", dest="outfile",
                  help="write XML to FILE", metavar="FILE")
parser.add_option("-m", "--exo-mesh-file", dest="meshfile",
                  help="Exodus Mesh File", metavar="FILE")
parser.add_option("--viz-output", dest="vizfile",
                  help="Viz output file (CGNS)", metavar="FILE")
parser.add_option("--stop-time", dest="end_time",
                  help="Final simulation time", metavar="TIME")
(options,args) = parser.parse_args()

if options.meshfile == None:
    raise ValueError("Please define a mesh file")

# ################################### #
# Mesh                                #
# ################################### #
# Setup a STK mesh
# replace with mesh = Mesh(type='stk',file=options.meshfile)
mesh = ParameterList("Mesh")
mesh.add_parameter("Framework","stk::mesh");
mesh.add_parameter("Read",options.meshfile);

# ################################### #
# Observations                        #
# ################################### #
observ = ParameterList("observation")
mass_of_water = observ.add_sublist("mass or water")
mass_of_water.add_parameter("state id", "water")
mass_of_water.add_parameter("region", "all")
mass_of_water.add_parameter("functional", "integral")
mass_of_water.add_parameter("times", [0.0])

# ################################### #
# MPC                                 #
# ################################### #
#  Create the MPC
mpc = MPC()

if options.end_time != None:
    mpc.set_end_time(float(options.end_time)) # End to convert to float to get the type right
else:
    mpc.set_end_time(3.1e+6)

mpc.enable_flow      = bool(True)
mpc.enable_transport = bool(False)
mpc.enable_chemistry = bool(False)

# Setup the viz input
if options.vizfile != None:
    mpc.viz.set_file(options.vizfile)
    mpc.viz.set_dnc(100)

mpc.add_parameter("Flow model", "Richard")
mpc.add_parameter("Gnuplot output", bool(True))
mpc.add_parameter("Restart", bool(False))
mpc.add_parameter("Restart file", "restart-bc-cribs-1D-4PE.h5")


# ################################### #
# State                               #
# ################################### #
# Define the State
state = ParameterList("State")
state.add_parameter("Number of mesh blocks",8)
state.add_parameter("Number of component concentrations",1)

# Water parameters should be under a simple ParameterList
state.add_parameter("Constant water saturation",1.0)
state.add_parameter("Constant water density",998.32)
state.add_parameter("Constant vicosity", 8.9e-4)

# Need to convert this to a vector not three paramters
state.add_parameter("Gravity x", 0.0)
state.add_parameter("Gravity y", 0.0)
state.add_parameter("Gravity z", -9.81)

# Mesh blocks
# Suggested change:
# backfill = state.add_meshblock(id=,permeability_type='constant|functional', porosity_type=...)
# backfill.porosity = 0.158
# backfill.permeability = 5.4344e-13
block_porosity = [0.158, 0.364, 0.388, 0.237, 0.360, 0.237, 0.360, 0.267]
block_permeability = [ 5.4344e-13, 4.83646e-12, 2.0447e-12, 2.9989e-13, 5.0618e-14, 2.9989e-13, 5.0618e-14, 3.7532e-13]
idx = 0
while idx < 8:
    label = "Mesh block ID %d" % (idx+1)
    block = state.add_sublist(label)
    block.add_parameter("Mesh Block ID", idx+1)
    block.add_parameter("Constant porosity", block_porosity[idx])
    block.add_parameter("Constant permeability", block_permeability[idx])
    block.add_parameter("Constant component concentration 0", 0.0)
    idx = idx + 1


# ################################### #
# Flow                                #
# ################################### #
flow = ParameterList("Flow")
richards = flow.add_sublist("Richards Problem")

# Steady State Parameters --- Should be a parameter list
richards.add_parameter("Steady state calculation initial time", 0.0)
richards.add_parameter("Steady state calculation final time", 3e10)
richards.add_parameter("Steady state calculation initial time step", 1.0e-7)
richards.add_parameter("Steady state calculation initial hydrostatic pressure height", 103.2 )
richards.add_parameter("Atmospheric pressure", 101325.0)

# Richards Model Evaluator
richards_eval = richards.add_sublist("Richards model evaluator")
richards_eval.add_parameter("Absolute error tolerance", 1.0)
richards_eval.add_parameter("Relative error tolerance", 1.0e-5)
richards_eval.add_verbose()

# Time Integration
time_integrator = richards.add_sublist("Time integrator")
time_integrator.add_parameter("Nonlinear solver max iterations", 6)
time_integrator.add_parameter("Nonlinear solver tolerance", 0.01)
time_integrator.add_parameter("NKA max vectors", 5)
time_integrator.add_parameter("NKA drop tolerance", 5.0e-2)
time_integrator.add_parameter("Maximum number of BDF tries", 10)
time_integrator.add_verbose('medium')

# Water retention model
H2O_models = richards.add_sublist("Water rentention models")

# Water Model: Backfill
backfill_h2o = WaterModel(0)
backfill_h2o.set_id(1)
backfill_h2o.van_genuchten_m(0.2857)
backfill_h2o.van_genuchten_alpha(1.9401e-4)
backfill_h2o.van_genuchten_res_saturation(0.103)
H2O_models.add_sublist(backfill_h2o)

# Water Model: Hanford Coarse Sand
coarse_sand_h2o = WaterModel(1)
coarse_sand_h2o.set_id(2)
coarse_sand_h2o.van_genuchten_m(0.5115)
coarse_sand_h2o.van_genuchten_alpha(7.3518e-4)
coarse_sand_h2o.van_genuchten_res_saturation(0.074)
H2O_models.add_sublist(coarse_sand_h2o)

# Water Model: Hanford Fine Sand
fine_sand_h2o = WaterModel(2)
fine_sand_h2o.set_id(3)
fine_sand_h2o.van_genuchten_m(0.6011)
fine_sand_h2o.van_genuchten_alpha(2.1443e-4)
fine_sand_h2o.van_genuchten_res_saturation(0.089)
H2O_models.add_sublist(fine_sand_h2o)

# Water Model: Hanford Gravel
hanford_gravel_h2o = WaterModel(3)
hanford_gravel_h2o.set_id(4)
hanford_gravel_h2o.van_genuchten_m(0.6011)
hanford_gravel_h2o.van_genuchten_alpha(2.1443e-4)
hanford_gravel_h2o.van_genuchten_res_saturation(0.089)
H2O_models.add_sublist(hanford_gravel_h2o)

# Water Model: Cold Creek caliche
creek_caliche_h2o = WaterModel(4)
creek_caliche_h2o.set_id(5)
creek_caliche_h2o.van_genuchten_m(0.5554)
creek_caliche_h2o.van_genuchten_alpha(5.1054e-5)
creek_caliche_h2o.van_genuchten_res_saturation(0.097)
H2O_models.add_sublist(creek_caliche_h2o)

# Water Model: Cold Creek gravel
creek_gravel_h2o = WaterModel(5)
creek_gravel_h2o.set_id(6)
creek_gravel_h2o.van_genuchten_m(0.4203)
creek_gravel_h2o.van_genuchten_alpha(1.7358e-4)
creek_gravel_h2o.van_genuchten_res_saturation(0.134)
H2O_models.add_sublist(creek_gravel_h2o)

# Water Model: Ringold Lower Mud 
mud_h2o = WaterModel(6)
mud_h2o.set_id(7)
mud_h2o.van_genuchten_m(0.5554)
mud_h2o.van_genuchten_alpha(5.1054e-5)
mud_h2o.van_genuchten_res_saturation(0.097)
H2O_models.add_sublist(mud_h2o)

# Water Model: Wooded Island
island_h2o = WaterModel(7)
island_h2o.set_id(8)
island_h2o.van_genuchten_m(0.3976)
island_h2o.van_genuchten_alpha(8.1687e-5)
island_h2o.van_genuchten_res_saturation(0.135)
H2O_models.add_sublist(island_h2o)




# Preconditioner Parameters
diff_precon = richards.add_sublist("Diffusion Preconditioner")
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
bc0 = richards.add_sublist("BC00")
bc0.add_parameter("Type", "Time Dependent Pressure Constant")
bc0.add_parameter("Initial BC Value", 1112016.18144)
bc0.add_parameter("BC Value", 101325.0)
bc0.add_parameter("Final Time", 1e10)
bc0.add_parameter("Side set ID", 2)

bc1 = richards.add_sublist("BC01")
bc1.add_parameter("Type", "Darcy Constant")
bc1.add_parameter("BC Value", -1.10984271943176e-10)
bc1.add_parameter("Side set ID", 1)

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
