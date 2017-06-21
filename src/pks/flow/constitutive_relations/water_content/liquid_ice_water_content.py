"""Water content for a two-phase, liquid+ice evaluator."""

import sys, os
sys.path.append(os.path.join(os.environ['ATS_SRC_DIR'], "tools", "evaluator_generator"))
from evaluator_generator import generate_evaluator

deps = [("porosity", "phi"),
        ("saturation_liquid", "sl"),
        ("molar_density_liquid", "nl"),
        ("saturation_ice", "si"),
        ("molar_density_ice", "ni"),
        ("cell_volume", "cv")
        ]
params = []

import sympy
phi, sl, nl, si, ni, cv = sympy.var("phi,sl,nl,si,ni,cv")
expression = phi*(sl*nl + si*ni)*cv;

generate_evaluator("liquid_ice_water_content", "Flow",
                   "liquid+ice water content", "water_content",
                   deps, params, expression=expression, doc=__doc__)
