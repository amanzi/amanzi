import sys, os
sys.path.append(os.path.join(os.environ['ATS_SRC_DIR'], "tools", "evaluator_generator"))
import evaluator_generator

deps = [("temperature", "temp"), ("pressure", "pres")]
params = [("cv", "double", "heat capacity"),
          ("T0", "double", "reference temperature [K]")]

import sympy
cv, T0 = sympy.var("cv_,T0_")
p, T = sympy.var("pres,temp")
expression = cv * (temp-T0)


evaluator_generator.generate_evaluator("eos_ideal_gas", "General", "ideal gas equation of state",
                                       "density", deps, params, expression=expression)
