"""Interfrost water content portion sl."""

import sys, os
sys.path.append(os.path.join(os.environ['ATS_SRC_DIR'], "tools", "evaluator_generator"))
from evaluator_generator import generate_evaluator
import sympy

# deps = [("porosity", "phi"),
#         ("saturation_liquid", "sl"),
#         ("molar_density_liquid", "nl"),
#         ("molar_density_ice", "ni"),
#         ("cell_volume", "cv")
#         ]
# params = []


# phi, sl, nl, ni, cv = sympy.var("phi,sl,nl,ni,cv")
# expression = phi*(nl-ni)*sl*cv;

# generate_evaluator("interfrost_sl_wc", "Flow",
#                    "interfrost sl water content", "water_content",
#                    deps, params, expression=expression, doc=__doc__)


# deps2 = [("molar_density_liquid", "nl"),
#          ("saturation_liquid", "sl"),
#          ("porosity", "phi")
#          ]
# params = [("beta", "double", "compressibility [1/Pa]"),]

# phi, sl, nl, beta_ = sympy.var("phi,sl,nl,beta_")
# expression2 = phi*sl*nl*beta_

# generate_evaluator("interfrost_dtheta_dpressure", "Flow",
#                    "interfrost dtheta_dpressure", "DThetaDp_coef",
#                    deps2, params, expression=expression2, doc=__doc__)

deps3 = [("porosity", "phi"),
         ("saturation_liquid", "sl"),
         ("molar_density_liquid", "nl"),
         ("saturation_ice", "si"),
         ("molar_density_ice", "ni"),
         ("density_rock", "rhos"),
         ("temperature", "T"),
         ]
params3 = [("W", "double", "W [K]")]

phi,sl,nl,si,ni,rhos,T,W_ = sympy.var("phi,sl,nl,si,ni,rhos,T,W_")

from sympy import Piecewise, exp

expression_part = (1.-.05)*exp(-((T-273.15)/W_)**2) + .05

expression3 = 1.e-6 * (phi * (sl*nl*75.3399846 + si*ni*37.111518) + (1-phi)*rhos*835 + phi * ni * 6017.1102 * Piecewise( (0., T >= 273.15), (expression_part.diff(T), True)))

generate_evaluator("interfrost_denergy_dtemperature", "Flow",
                   "interfrost denergy_dtemperature", "DEnergyDT_coef",
                   deps3, params3, expression=expression3, doc=__doc__)
          
