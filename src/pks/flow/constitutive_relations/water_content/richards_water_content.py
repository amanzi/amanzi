"""Richards water content with no vapor."""

import sys, os
sys.path.append(os.path.join(os.environ['ATS_SRC_DIR'], "tools", "evaluator_generator"))
from evaluator_generator import generate_evaluator

deps = [("porosity", "phi"),
        ("saturation_liquid", "sl"),
        ("molar_density_liquid", "nl"),
        ("cell_volume", "cv")
        ]
params = []

import sympy
phi, sl, nl, cv = sympy.var("phi,sl,nl,cv")
expression = phi*sl*nl*cv;

generate_evaluator("richards_water_content", "Flow",
                   "richards water content", "water_content",
                   deps, params, expression=expression, doc=__doc__)
