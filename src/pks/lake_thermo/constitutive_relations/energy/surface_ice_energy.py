"""Energy evaulator for ice+liquid surface water."""

import sys, os
sys.path.append(os.path.join(os.environ['ATS_SRC_DIR'], "tools", "evaluator_generator"))
from evaluator_generator import generate_evaluator

deps = [("ponded_depth", "h"),
        ("unfrozen_fraction", "eta"),
        ("molar_density_liquid", "nl"),
        ("internal_energy_liquid", "ul"),
        ("molar_density_ice", "ni"),
        ("internal_energy_ice", "ui"),
        ("cell_volume", "cv")
        ]
params = []

import sympy
h, eta, nl, ul, ni, ui, cv = sympy.var("h, eta, nl, ul, ni, ui, cv")
expression = h * ( eta*nl*ul + (1-eta)*ni*ui ) * cv;

generate_evaluator("surface_ice_energy", "Energy",
                   "surface ice energy", "energy",
                   deps, params, expression=expression, doc=__doc__)
