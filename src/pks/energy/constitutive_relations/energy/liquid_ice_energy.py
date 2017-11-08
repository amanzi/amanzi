"""Energy for a two-phase, liquid+ice evaluator."""

import sys, os
sys.path.append(os.path.join(os.environ['ATS_SRC_DIR'], "tools", "evaluator_generator"))
from evaluator_generator import generate_evaluator

deps = [("porosity", "phi"),
        ("base_porosity", "phi0"),
        ("saturation_liquid", "sl"),
        ("molar_density_liquid", "nl"),
        ("internal_energy_liquid", "ul"),
        ("saturation_ice", "si"),
        ("molar_density_ice", "ni"),
        ("internal_energy_ice", "ui"),
        ("density_rock", "rho_r"),
        ("internal_energy_rock", "ur"),
        ("cell_volume", "cv")
        ]
params = []

import sympy
phi, phi0, sl, nl, ul, si, ni, ui, rho_r, ur, cv = sympy.var("phi, phi0, sl, nl, ul, si, ni, ui, rho_r, ur, cv")
expression = (phi*(sl*nl*ul + si*ni*ui) + (1-phi0)*rho_r*ur) * cv;

generate_evaluator("liquid_ice_energy", "Energy",
                   "liquid+ice energy", "energy",
                   deps, params, expression=expression, doc=__doc__)
