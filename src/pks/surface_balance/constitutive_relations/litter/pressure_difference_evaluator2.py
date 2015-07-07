import sys, os
sys.path.append(os.path.join(os.environ['ATS_SRC_DIR'], "tools", "evaluator_generator"))
import evaluator_generator

# 
deps = [("macropore_pressure", "pM"),
        ("surface_pressure", "ps"),
        ("surface_relative_permeability", "krs"),
        ("macropore_relative_permeability", "krM"),
        ("macropore_absolute_permeability", "K")]
params = [("gamma", "double", "gamma [-]"),
          ("delta", "double", "delta [m]")]

evaluator_generator.generate_evaluator("macropore_surface_flux", "SurfaceBalance", "macropore-surface flux",
                                       "macropore_surface_flux", deps, params)



