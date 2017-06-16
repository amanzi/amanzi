"""Richards energy: the standard form as a function of liquid saturation and specific internal energy."""

import sys, os
sys.path.append(os.path.join(os.environ['ATS_SRC_DIR'], "tools", "evaluator_generator"))
from evaluator_generator import generate_evaluator

deps = [("porosity", "phi"),
        ("saturation_liquid", "sl"),
        ("molar_density_liquid", "nl"),
        ("internal_energy_liquid", "ul"),
        ("density_rock", "rho_r"),
        ("internal_energy_rock", "ur"),
        ("cell_volume", "cv")
        ]
params = []

import sympy
phi, sl, nl, ul, rho_r, ur, cv = sympy.var("phi, sl, nl, ul, rho_r, ur, cv")
expression = (phi*(sl*nl*ul) + (1-phi)*rho_r*ur) * cv;

generate_evaluator("richards_energy", "Energy",
                   "richards energy", "energy",
                   deps, params, expression=expression, doc=__doc__)
