import sys, os
sys.path.append(os.path.join(os.environ['ATS_SRC_DIR'], "tools", "evaluator_generator"))
import evaluator_generator

# latent heat
deps = [("litter_water_content", "wc"),
        ("surface_molar_density_liquid", "rho"),
        ("litter_thickness", "L")]
params = [("wc_sat", "double", "saturated litter water content [-]"),
          ("tau", "double", "drying time [s]")]

evaluator_generator.generate_evaluator("evaporative_flux_relaxation", "SurfaceBalance", "evaporative flux relaxation",
                                       "evaporative_flux", deps, params)


# latent heat
deps = [("evaporative_flux", "qe")]
params = [("Le", "double", "latent heat of vaporization [MJ/mol]", .0449994810744)]

evaluator_generator.generate_evaluator("latent_heat", "SurfaceBalance", "latent heat from evaporative flux",
                                       "latent_heat", deps, params)


