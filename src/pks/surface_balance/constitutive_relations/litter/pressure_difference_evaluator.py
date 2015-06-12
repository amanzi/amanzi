import sys, os
sys.path.append(os.path.join(os.environ['ATS_SRC_DIR'], "tools", "evaluator_generator"))
import evaluator_generator

# 
deps = [("micropore_pressure", "pm"),
        ("pressure", "pM"),
        ("relative_permeability", "krM"),
        ("micropore_relative_permeability", "krm"),
        ("micropore_absolute_permeability", "K")]
params = [("gamma", "double", "gamma [-]"),
          ("delta", "double", "delta [m]")]

evaluator_generator.generate_evaluator("micropore_macropore_flux", "SurfaceBalance", "micropore-macropore flux",
                                       "micropore_macropore_flux", deps, params)



