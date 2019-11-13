"""Interception factor: fraction of incoming precip that is intercepted."""

import sys, os
sys.path.append(os.path.join(os.environ['ATS_SRC_DIR'], "tools", "evaluator_generator"))
import evaluator_generator

# latent heat
deps = [("surface-area_index", "ai"),]
params = [("alpha", "double", "scaling factor [-]", 0.25),] # Lawrence 2007

import sympy
ai, alpha_ = sympy.var("ai, alpha")
expression = alpha_ * (1 - sympy.exp(-0.5 * ai))


evaluator_generator.generate_evaluator("interception_fraction", "SurfaceBalance", "interception fraction",
                                       "interception_fraction", deps, params, expression=expression)



