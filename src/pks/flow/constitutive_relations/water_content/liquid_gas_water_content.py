"""Water content for a two-phase, liquid+water vapor evaluator."""

import sys, os
sys.path.append(os.path.join(os.environ['ATS_SRC_DIR'], "tools", "evaluator_generator"))
from evaluator_generator import generate_evaluator

deps = [("porosity", "phi"),
        ("saturation_liquid", "sl"),
        ("molar_density_liquid", "nl"),
        ("saturation_gas", "sg"),
        ("molar_density_gas", "ng"),
        ("mol_frac_gas", "omega"),
        ("cell_volume", "cv")
        ]
params = []

import sympy
phi, sl, nl, sg, ng, omega, cv = sympy.var("phi,sl,nl,sg,ng,omega,cv")
expression = phi*(sl*nl + sg*ng*omega)*cv;

generate_evaluator("liquid_gas_water_content", "Flow",
                   "liquid+gas water content", "water_content",
                   deps, params, expression=expression, doc=__doc__)
