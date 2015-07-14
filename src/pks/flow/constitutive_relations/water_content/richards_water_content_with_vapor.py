"""Richards water content with vapor."""

import sys, os
sys.path.append(os.path.join(os.environ['ATS_SRC_DIR'], "tools", "evaluator_generator"))
from evaluator_generator import generate_evaluator

deps = [("porosity", "phi"),
        ("saturation_liquid", "sl"),
        ("saturation_gas", "sg"),
        ("molar_density_liquid", "nl"),
        ("molar_density_gas", "ng"),
        ("mol_fraction_vapor_in_gas", "omega"),
        ("cell_volume", "cv")
        ]
params = []

import sympy
phi, sl, nl, cv = sympy.var("phi,sl,nl,cv")
sg, ng, omega = sympy.var("sg,ng,omega")
expression = cv*phi*(sl*nl + sg*ng*omega)

generate_evaluator("richards_water_content_with_vapor", "Flow",
                   "richards water content with vapor", "water_content",
                   deps, params, expression=expression, doc=__doc__)
