"""Water content for a three-phase, gas+liquid+ice evaluator."""

import sys, os
sys.path.append(os.path.join(os.environ['ATS_SRC_DIR'], "tools", "evaluator_generator"))
from evaluator_generator import generate_evaluator

deps = [("porosity", "phi"),
        ("saturation_liquid", "sl"),
        ("molar_density_liquid", "nl"),
        ("saturation_ice", "si"),
        ("molar_density_ice", "ni"),
        ("saturation_gas", "sg"),
        ("molar_density_gas", "ng"),
        ("mol_frac_gas", "omega"),
        ("cell_volume", "cv")
        ]
params = []

import sympy
phi, sl, nl, si, ni, sg, ng, omega, cv = sympy.var("phi,sl,nl,si,ni,sg,ng,omega,cv")
expression = phi*(sl*nl + si*ni + sg*ng*omega)*cv;

generate_evaluator("three_phase_water_content", "Flow",
                   "three phase water content", "water_content",
                   deps, params, expression=expression, doc=__doc__)
