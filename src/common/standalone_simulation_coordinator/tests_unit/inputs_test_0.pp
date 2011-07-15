#*******************************************************************************
# INPUTS_TEST_0: Single Phase Flow
#*******************************************************************************

# Model name
# ----------
prob.model_name = single-phase


# Geometry
# --------
geometry.coord_sys   = 0    
geometry.prob_lo     = 0. 0. 
geometry.prob_hi     = 8. 8.
geometry.is_periodic = 0  0

# Define regions
geometry.region = left_region right_region obs_region

geometry.left_region
{
  purpose = 0
  type = 1
  param = 0. 0. 4. 8.
}

geometry.right_region
{
  purpose = 0
  type = 1
  param = 4. 0. 8. 8.
}

geometry.obs_region
{
  purpose = 0
  type = 1
  param = 0. 0. 2. 8.
}

# Rock 
# ----
rock.rock = left right
rock.left
{
  density      = 2.8e3
  porosity     = 0.258
  permeability = 336 336
  kr_type      = 3
  kr_param     = 0.4203 0.081 0
  cpl_type     = 3
  cpl_param    = 0.4203 17.329 0.081 0
}
rock.right
{
  density      = 2.8e3
  porosity     = 0.422
  permeability = 2294 2294
  kr_type      = 3
  kr_param     = 0.6011 0.081 0
  cpl_type     = 3
  cpl_param    = 0.6011 21.41 0.081 0
}

# assign rock properties to region
rock.assign = left_region right_region
rock.left_region
{
  type  = left
  phi   = 1
  kappa = 1
}
rock.right_region
{
  type  = right
  phi   = 1
  kappa = 1
}

# save the generated map of the permeability and porosity
rock.permeability_file = test_0/kp
rock.porosity_file     = test_0/pp

# phase and components
# --------------------
phase.phase  = Liquid 
comp.comp = Water Air
comp.Water
{
  phase       = Liquid
  density     = 1.0e3
  viscosity   = 1.0
  diffusivity = 0.0
}
comp.Air
{
  phase       = Liquid
  density     = 1.0e3
  viscosity   = 1.0
  diffusivity = 0.0
}
comp.dominant = Water
comp.lo_bc = 1 4 
comp.hi_bc = 2 4

# initialization of component
comp.init = ALL
comp.ALL
{
   type  = scalar
   Water = 1000
   Air   = 0
}

# boundary condition
comp.inflow = XLOBC 
comp.XLOBC
{
   type = scalar
   Water = 0
   Air   = 1000
}

# tracer
# ------
tracer.ntracer = 0

# chemistry
# ---------
prob.do_chem = -1

# pressure (atm)
# -------------
# note the 1 again to indicate inflow
press.lo_bc = 2 4 
press.hi_bc = 2 4 
press.press_lo = 1. 0.
press.press_hi = 0. 0.
press.inflow_bc_lo  = 0 0 
press.inflow_bc_hi  = 0 0 
press.inflow_vel_lo = 0 0
press.inflow_vel_hi = 0 0
press.water_table_lo = 0
press.water_table_hi = 0

# source term
# -----------
source.do_source = 0

# observation
# -----------
observation.nobs = 1
observation.observation = output1
observation.output_file = test_0/obs.out
observation.output1
{
  var_type  = comp
  var_id    = Air
  region    = obs_region
  obs_type  = integral
  times     = 1000 
}

# flow related
prob.have_capillary = 0
prob.no_corrector = 0
prob.full_cycle = 1
prob.do_reflux  = 1
prob.do_kappa_refine = 1
amr.check_file = test_0/chk
amr.plot_file  = test_0/plt
amr.max_grid_size = 32

amr.n_error_buf = 2 2 2 
amr.grid_eff = 0.75
#prob.fixed_dt  = 20000
#max_step = 1
stop_time = 1000
#*******************************************************************************
# Number of cells in each coordinate direction at the coarsest level
amr.n_cell = 16 16

#*******************************************************************************
# Maximum level (defaults to 0 for single level calculation)
amr.max_level = 0

#*******************************************************************************

#Verbosity controls
amr.v = 1
mac.v = 0
mg.v  = 2
prob.v  = 2
cg.v  = 0
diffuse.v = 0

#Solver controls
#mac.use_cg_solver = 0
prob.visc_tol      = 1.e-12
prob.visc_abs_tol  = 1.e-14
mac.mac_tol      = 1.e-14
mac.mac_sync_tol = 1.e-16
mac.mac_abs_tol  = 1.e-14
#mac.use_fboxlib_mg = true
#diffuse.use_fboxlib_mg = true
mg.cg_solver = 1
mg.maxiter = 100
cg.unstable_criterion = 100
cg.use_jacobi_precond = 1
#mg.usecg = 2
#mg.use_jbb_precond = 1
#mg.maxiter = 100
##mg.nu_b = 4
#mg.nu_0 = 2
#mg.nu_1 = 4
#mg.nu_2 = 4
#mg.nu_f = 10
#mg.rtol_b = 1.e-4
mg.smooth_on_cg_unstable = 1	
#*******************************************************************************
# Interval (in number of level l timesteps) between regridding
amr.regrid_int		= 2 

#*******************************************************************************
# Refinement ratio as a function of level
amr.ref_ratio		= 2 2

#*******************************************************************************
# Interval (in number of coarse timesteps) between plotfiles / checkpoint files
amr.plot_int		= 1000
amr.check_int		= 1000

#*******************************************************************************

# CFL number to be used in calculating the time step : dt = dx / max(velocity)
prob.cfl                  = .8 # CFL number used to set dt

#*******************************************************************************

# Factor by which the first time is shrunk relative to CFL constraint
prob.init_shrink          = .1  # factor which multiplies the very first time step

#*******************************************************************************

# Name of the file which specifies problem-specific parameters (defaults to "probin")
amr.probin_file 	= probin


#*******************************************************************************

# Factor by which grids must be coarsenable.
amr.blocking_factor     = 4

#*******************************************************************************

# Add vorticity to the variables in the plot files.
amr.derive_plot_vars    = gradpx gradpy gradn

#*******************************************************************************
