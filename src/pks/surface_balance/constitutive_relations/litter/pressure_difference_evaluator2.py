import sys, os
sys.path.append(os.path.join(os.environ['ATS_SRC_DIR'], "tools", "evaluator_generator"))
import evaluator_generator

# 
deps = [("macropore_pressure", "pM"),
        ("surface_pressure", "ps"),
        ("surface_relative_permeability", "krs"),
        ("macropore_relative_permeability", "krM"),
        ("macropore_absolute_permeability", "K"),
        ("surface_molar_density_liquid", "nl"),
        ("surface_viscosity", "mu")]
params = [("gamma", "double", "gamma [-]"),
          ("delta", "double", "delta [m]"),
          ("patm", "double", "atmospheric pressure [Pa]")]

evaluator_generator.generate_evaluator("macropore_surface_flux", "SurfaceBalance", "macropore-surface flux",
                                       "macropore_surface_flux", deps, params)



