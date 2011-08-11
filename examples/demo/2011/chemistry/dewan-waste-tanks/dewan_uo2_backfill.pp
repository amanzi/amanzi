#*******************************************************************************
# INPUTS_ASCEM: Two-phase model with tracers
#*******************************************************************************
# Geometry
# --------
geometry.coord_sys   = 0    
geometry.prob_lo     = 0. 0. 
geometry.prob_hi     = 16 16
geometry.is_periodic = 0 0

geometry.region = Backfill_Region Waste_Region 

geometry.Backfill_Region.purpose = 0
geometry.Backfill_Region.type = 1
geometry.Backfill_Region.param = 0 0 16 16

geometry.Waste_Region.purpose = 0
geometry.Waste_Region.type = 1
geometry.Waste_Region.param = 6 10 10 14


# Rock
# ----

rock.rock = Backfill Waste

rock.Backfill.density = 2.63e3
rock.Backfill.porosity = 0.35
rock.Backfill.permeability = 78.58 78.58
rock.Backfill.kr_type = 3
rock.Backfill.kr_param = 0.11238 0.1893 0
rock.Backfill.cpl_type = 3
rock.Backfill.cpl_param = 0.11238 78 0.1893 0

rock.Waste.density = 2.63e3
rock.Waste.porosity = 0.35
rock.Waste.permeability = 100.0 100.0
rock.Waste.kr_type = 3
rock.Waste.kr_param = 0.11238 0.1893 0
rock.Waste.cpl_type = 3
rock.Waste.cpl_param = 0.11238 78 0.1893 0

#rock.Waste
#{
#  density = 2.51e3
#  porosity = 0.266
#  permeability = 0.0372 0.0372
#  kr_type = 3
#  kr_param = 0.48541 0 0
#  cpl_type = 3
#  cpl_param = 0.48541 0.002 0 0
#}


rock.assign = Backfill_Region Waste_Region

rock.Backfill_Region.type = Backfill
rock.Backfill_Region.phi = 1
rock.Backfill_Region.kappa = 1

rock.Waste_Region.type = Waste
rock.Waste_Region.phi = 1
rock.Waste_Region.kappa = 1

rock.permeability_file = test/kp
rock.porosity_file = test/pp

# phase and components
# --------------------
phase.phase  = Liquid Gas
# water must come before air for now.

comp.comp = Water Air

comp.Water.phase = Liquid
comp.Water.density = 1.0e3
comp.Water.viscosity = 1.0
comp.Water.diffusivity = 0.0

comp.Air.phase = Gas
comp.Air.density = 1.2
comp.Air.viscosity = 0.018
comp.Air.diffusivity = 0.0

comp.dominant = Air
comp.lo_bc = 4 1 
comp.hi_bc = 4 1

comp.init = ALL

comp.ALL.type = zero_total_velocity
comp.ALL.inflow = -2.53e-8

# boundary conditions
comp.inflow = YHIBC YLOBC

comp.YHIBC.type = zero_total_velocity
comp.YHIBC.rock = Backfill
comp.YHIBC.inflow = -2.535e-8

comp.YLOBC.type = scalar
comp.YLOBC.Water = 1000
comp.YLOBC.Air = 0.



# tracer
# ------
tracer.ntracer = 13
tracer.group   = Total Solid Absorbed

# these are the correct species names
tracer.tracer  = Al+++ H+ HPO4-- SiO2 UO2++ Koalinite Quartz Uranium SAl+++ SH+ SHPO4-- SSiO2 SUO2++

tracer.Al+++.group = Total
tracer.H+.group = Total
tracer.HPO4--.group = Total
tracer.SiO2.group = Total
tracer.UO2++.group = Total
tracer.Koalinite.group = Solid
tracer.Quartz.group = Solid
tracer.Uranium.group = Solid
tracer.SAl+++.group = Absorbed
tracer.SH+.group = Absorbed
tracer.SHPO4--.group = Absorbed
tracer.SSiO2.group = Absorbed
tracer.SUO2++.group = Absorbed


tracer.init = Backfill_Region Waste_Region

tracer.Backfill_Region.type      = scalar
tracer.Backfill_Region.Al+++     = 6.5874e-9
tracer.Backfill_Region.H+        = -3.1408e-7
tracer.Backfill_Region.HPO4--    = 1.0000e-6
tracer.Backfill_Region.SiO2      = 1.8703e-4
tracer.Backfill_Region.UO2++     = 1.000e-15
tracer.Backfill_Region.Koalinite = 0.15
tracer.Backfill_Region.Quartz    = 0.21
tracer.Backfill_Region.Uranium   = 0.00
tracer.Backfill_Region.SAl+++    = 0.00
tracer.Backfill_Region.SH+       = 0.00
tracer.Backfill_Region.SHPO4--   = 0.00
tracer.Backfill_Region.SSiO2     = 0.00
tracer.Backfill_Region.SUO2++    = 0.00


tracer.Waste_Region.type      = scalar
tracer.Waste_Region.Al+++     = 2.8909e-5
tracer.Waste_Region.H+        = 1.2786e-3
tracer.Waste_Region.HPO4--    = 7.1028e-5
tracer.Waste_Region.SiO2      = 2.5280e-4
tracer.Waste_Region.UO2++     = 3.5414e-5
#tracer.Waste_Region.UO2++     = 3.5414e-2
tracer.Waste_Region.Koalinite = 0.0
tracer.Waste_Region.Quartz    = 0.0
tracer.Waste_Region.Uranium   = 0.0
tracer.Waste_Region.SAl+++    = 0.0
tracer.Waste_Region.SH+       = 0.00      
tracer.Waste_Region.SHPO4--   = 0.00      
tracer.Waste_Region.SSiO2     = 0.00
tracer.Waste_Region.SUO2++    = 0.00

tracer.inflow = YHIBC YLOBC

tracer.YHIBC.type      = scalar
tracer.YHIBC.Al+++     = 6.5874e-9
tracer.YHIBC.H+        = -3.1408e-7
tracer.YHIBC.HPO4--    = 1.0000e-6
tracer.YHIBC.SiO2      = 1.8703e-4
tracer.YHIBC.UO2++     = 1.000e-15
tracer.YHIBC.Koalinite = 0.15
tracer.YHIBC.Quartz    = 0.21
tracer.YHIBC.Uranium   = 0.00
tracer.YHIBC.SAl+++    = 0.00
tracer.YHIBC.SH+       = 0.00
tracer.YHIBC.SHPO4--   = 0.00
tracer.YHIBC.SSiO2     = 0.00
tracer.YHIBC.SUO2++    = 0.00

tracer.YLOBC.type      = scalar
tracer.YLOBC.Al+++     = 6.5874e-9
tracer.YLOBC.H+        = -3.1408e-7
tracer.YLOBC.HPO4--    = 1.0000e-6
tracer.YLOBC.SiO2      = 1.8703e-4
tracer.YLOBC.UO2++     = 1.000e-15
tracer.YLOBC.Koalinite = 0.15
tracer.YLOBC.Quartz    = 0.21
tracer.YLOBC.Uranium   = 0.00
tracer.YLOBC.SAl+++    = 0.00
tracer.YLOBC.SH+       = 0.00
tracer.YLOBC.SHPO4--   = 0.00
tracer.YLOBC.SSiO2     = 0.00
tracer.YLOBC.SUO2++    = 0.00


# chemistry
# ---------
prob.amanzi.file = uo2-5-component.bgd
prob.n_chem_interval = 0

prob.do_chem = 2

# pressure (atm)
# -------------

press.press_lo = 0 0
press.press_hi = 0 0
press.inflow_bc_lo  = 0 0 
press.inflow_bc_hi  = 0 1 
press.inflow_vel_lo = 0 0
press.inflow_vel_hi = 0 0
press.lo_bc = 4 1
press.hi_bc = 4 1
press.water_table_lo = 0
press.water_table_hi = 0

# source term
# -----------
source.do_source = 0


source.source = Infiltration 
#Discharge Tracer_Discharge

source.Infiltration.var_type  = comp
source.Infiltration.var_id    = Water
source.Infiltration.region    = Top_Region
source.Infiltration.dist_type = constant
source.Infiltration.val       = 9.524e-6


# capillary pressure
prob.have_capillary = 1

# flow related
prob.model  = 2
prob.no_corrector = 0
prob.full_cycle = 1
prob.do_reflux  = 1
#amr.restart = temp/chk10000
amr.check_file = temp/chk
amr.plot_file = temp/plt
amr.max_grid_size = 64

amr.n_error_buf = 2 2 2 
amr.grid_eff = 0.75

#prob.fixed_dt  = 8000
max_step = 20000
#stop_time = 500000

#*******************************************************************************
# Number of cells in each coordinate direction at the coarsest level
amr.n_cell = 32 32

#*******************************************************************************
# Maximum level (defaults to 0 for single level calculation)
amr.max_level = 0

#*******************************************************************************

#Verbosity controls
amr.v = 2
mac.v = 3
mg.v  = 1
prob.v  = 3
cg.v  = 0
diffuse.v = 2

#Solver controls
#mac.use_cg_solver = 0
prob.visc_tol      = 1.e-12
prob.visc_abs_tol  = 1.e-14
mac.mac_tol      = 1.e-12
mac.mac_sync_tol = 1.e-12
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
amr.plot_int		= 100
amr.check_int		= 100

#*******************************************************************************

# CFL number to be used in calculating the time step : dt = dx / max(velocity)
prob.cfl                  = 0.8 # CFL number used to set dt

#*******************************************************************************

# Factor by which the first time is shrunk relative to CFL constraint
prob.init_shrink          = 1  # factor which multiplies the very first time step

#*******************************************************************************

# Name of the file which specifies problem-specific parameters (defaults to "probin")
amr.probin_file 	= probin


#*******************************************************************************

# Factor by which grids must be coarsenable.
amr.blocking_factor     = 16

#*******************************************************************************

# Add vorticity to the variables in the plot files.
amr.derive_plot_vars    = gradpx gradpy gradn

#*******************************************************************************
