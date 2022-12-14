/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include <winstd.H>
#include <ParmParse.H>
#include <Interpolater.H>
#include <MultiGrid.H>
#include <ArrayLim.H>
#include <TagBox.H>
#include <DataServices.H>
#include <AmrData.H>
#include <Utility.H>
#include <time.h>

#include <PorousMedia.H>
#include <PMAMR_Labels.H>
#include <RegType.H>
#include <Prob_PM_F.H>
#include <PMAMR_Labels.H>
#include <PMAmr.H>
#include <ChemConstraintEval.H>
#include <DiffDomRelSrc.H>

#ifdef _OPENMP
#include "omp.h"
#endif

#define SHOWVALARR(val)                        \
{                                              \
    std::cout << #val << " = ";                \
    for (int i=0;i<val.size();++i)             \
    {                                          \
        std::cout << val[i] << " " ;           \
    }                                          \
    std::cout << std::endl;                    \
}
#define SHOWVALARRA(val) { SHOWVALARR(val); BoxLib::Abort();}
#define SHOWVAL(val) { std::cout << #val << " = " << val << std::endl;}
#define SHOWVALA(val) { SHOWVAL(val); BoxLib::Abort();}


#include <AlquimiaHelper_Structured.H>
#include <AmanziChemHelper_Structured.H>

#include "ActivityModelFactory.hh"
#include "SimpleThermoDatabase.hh"

#include <TabularFunction.H>

std::ostream& operator<< (std::ostream& os, const Array<std::string>& rhs)
{
    for (int i=0; i<rhs.size(); ++i) {
        os << rhs[i] << " ";
    }
    return os;
}

std::ostream& operator<< (std::ostream& os, const Array<int>& rhs)
{
    for (int i=0; i<rhs.size(); ++i) {
        os << rhs[i] << " ";
    }
    return os;
}

std::ostream& operator<< (std::ostream& os, const Array<Real>& rhs)
{
    for (int i=0; i<rhs.size(); ++i) {
        os << rhs[i] << " ";
    }
    return os;
}

namespace
{
  const std::string solid("Solid");
  const std::string absorbed("Absorbed");
}

//
//**********************************************************************
//
// Set all default values for static variables in InitializeStaticVariables()!!!
//
//**********************************************************************
//

int PorousMedia::echo_inputs;
//
// The num_state_type actually varies with model.
//
// Add 2 if do_tracer_chemistry>0 later.
//
int PorousMedia::num_state_type;
//
// Region.
//
std::string      PorousMedia::surf_file;
//
// Source.
//
bool          PorousMedia::do_source_term;
Array<Real>   PorousMedia::source_volume;

//
// Phases and components.
//
Array<std::string>  PorousMedia::pNames;
Array<std::string>  PorousMedia::cNames;
Array<int >         PorousMedia::pType;
Array<Real>         PorousMedia::density;
PArray<RegionData>  PorousMedia::ic_array;
PArray<RegionData>  PorousMedia::bc_array;
PArray<RegionData>  PorousMedia::source_array;
Array<Real>         PorousMedia::muval;
int                 PorousMedia::nphases;
int                 PorousMedia::ncomps;
int                 PorousMedia::ndiff;
//
// Tracers.
//
Array<std::string>  PorousMedia::qNames;
Array<std::string>  PorousMedia::tNames;
int                 PorousMedia::ntracers;
Array<int>          PorousMedia::tType;
Array<Real>         PorousMedia::tDen;
Array<PArray<IdxRegionData> > PorousMedia::tic_array;
Array<PArray<IdxRegionData> > PorousMedia::tbc_array;
Array<PArray<RegionData> > PorousMedia::tsource_array;
std::map<std::string,Array<int> > PorousMedia::group_map;
RockManager::ChemICMap       PorousMedia::solute_chem_ics; // sc[icname][solute][property] = val
RockManager::ICLabelParmPair PorousMedia::sorption_chem_ics; // sc[icname][property] = val

//
// Minerals and Sorption sites
//
double             PorousMedia::uninitialized_data;
int                PorousMedia::nminerals;
Array<std::string> PorousMedia::minerals;
int                PorousMedia::nsorption_sites;
Array<std::string> PorousMedia::sorption_sites;
int                PorousMedia::ncation_exchange;
int                PorousMedia::nsorption_isotherms;
bool               PorousMedia::using_sorption;

// Pressure.
//
Real        PorousMedia::wt_lo;
Real        PorousMedia::wt_hi;
Array<Real> PorousMedia::press_lo;
Array<Real> PorousMedia::press_hi;
Array<int>  PorousMedia::inflow_bc_lo;
Array<int>  PorousMedia::inflow_bc_hi;
Array<int>  PorousMedia::rinflow_bc_lo;
Array<int>  PorousMedia::rinflow_bc_hi;
//
// Temperature.
//
Real  PorousMedia::temperature;
//
// Observations
//
int PorousMedia::verbose_observation_processing;

std::map<std::string,Array<std::string> > PorousMedia::write_region_sets;

//
// Flow.
//
int  PorousMedia::verbose;
Real PorousMedia::cfl;
Real PorousMedia::init_shrink;
Real PorousMedia::dt_grow_max;
Real PorousMedia::dt_shrink_max;
Real PorousMedia::fixed_dt;
Real PorousMedia::steady_max_dt;
Real PorousMedia::transient_max_dt;
Real PorousMedia::dt_cutoff;
Real PorousMedia::gravity;
int  PorousMedia::gravity_dir;
Real PorousMedia::z_location;
Real PorousMedia::domain_thickness;
int  PorousMedia::initial_step;
int  PorousMedia::initial_iter;
int  PorousMedia::sum_interval;
int  PorousMedia::NUM_SCALARS;
int  PorousMedia::NUM_STATE;
Real PorousMedia::dt_init;
int  PorousMedia::max_n_subcycle_transport;
int  PorousMedia::max_dt_iters_flow;
bool PorousMedia::abort_on_chem_fail;
bool PorousMedia::show_selected_runtimes;

//
// Viscosity parameters.
//
Real PorousMedia::be_cn_theta;
Real PorousMedia::visc_tol;
Real PorousMedia::visc_abs_tol;
bool PorousMedia::def_harm_avg_cen2edge;
//
// Capillary pressure flag.
//
Real PorousMedia::atmospheric_pressure_atm;
std::map<std::string,bool> PorousMedia::use_gauge_pressure;

//
// Molecular diffusion flag.
//
bool PorousMedia::variable_scal_diff;

Array<int>  PorousMedia::is_diffusive;
Array<Real> PorousMedia::visc_coef;
Array<Real> PorousMedia::molecular_diffusivity;
//
// Transport flags
//
bool PorousMedia::do_tracer_advection;
bool PorousMedia::do_tracer_diffusion;
bool PorousMedia::setup_tracer_transport;
bool PorousMedia::advect_tracers;
bool PorousMedia::diffuse_tracers;
bool PorousMedia::tensor_tracer_diffusion;
bool PorousMedia::solute_transport_limits_dt;

//
// Chemistry flag.
//
bool PorousMedia::do_tracer_chemistry;
bool PorousMedia::react_tracers;
bool PorousMedia::do_full_strang;
int  PorousMedia::n_chem_interval;
int  PorousMedia::it_chem;
Real PorousMedia::dt_chem;
int  PorousMedia::max_grid_size_chem;
bool PorousMedia::no_initial_values;
bool PorousMedia::use_funccount;
Real PorousMedia::max_chemistry_time_step;
Array<Real> PorousMedia::first_order_decay_constant;
//
// Lists.
//
std::map<std::string, int> PorousMedia::phase_list;
std::map<std::string, int> PorousMedia::comp_list;
std::map<std::string, int> PorousMedia::tracer_list;
Array<std::string> PorousMedia::user_derive_list;
PorousMedia::FLOW_MODEL_ID PorousMedia::flow_model;
PorousMedia::FLOW_EVAL_ID PorousMedia::flow_eval;
bool PorousMedia::flow_is_static;
//
// AMANZI flags.
//
#ifdef AMANZI

#ifdef ALQUIMIA_ENABLED
Amanzi::AmanziChemistry::ChemistryEngine* PorousMedia::chemistry_engine;
#endif
std::string PorousMedia::chemistry_model_name;
std::string PorousMedia::chemistry_engine_name;
ChemistryHelper_Structured* PorousMedia::chemistry_helper;

std::string PorousMedia::amanzi_database_file;
std::string PorousMedia::amanzi_activity_model;

#endif
//
// Internal switches.
//
bool PorousMedia::do_reflux;
Real PorousMedia::ic_chem_relax_dt;
int  PorousMedia::nGrowHYP;
int  PorousMedia::nGrowMG;
int  PorousMedia::nGrowEIGEST;
bool PorousMedia::do_constant_vel;
Real PorousMedia::be_cn_theta_trac;
bool PorousMedia::do_output_flow_time_in_years;
bool PorousMedia::do_output_chemistry_time_in_years;
bool PorousMedia::do_output_transport_time_in_years;
//
// Init to steady
//
int  PorousMedia::steady_min_iterations;
int  PorousMedia::steady_min_iterations_2;
int  PorousMedia::steady_max_iterations;
int  PorousMedia::steady_limit_iterations;
Real PorousMedia::steady_time_step_reduction_factor;
Real PorousMedia::steady_time_step_increase_factor;
Real PorousMedia::steady_time_step_increase_factor_2;
Real PorousMedia::steady_time_step_retry_factor_1;
Real PorousMedia::steady_time_step_retry_factor_2;
Real PorousMedia::steady_time_step_retry_factor_f;
int  PorousMedia::steady_max_consecutive_failures_1;
int  PorousMedia::steady_max_consecutive_failures_2;
int  PorousMedia::steady_max_time_steps;
Real PorousMedia::steady_max_time_step_size;
int  PorousMedia::steady_max_num_consecutive_success;
Real PorousMedia::steady_extra_time_step_increase_factor;
int  PorousMedia::steady_max_num_consecutive_increases;
Real PorousMedia::steady_consecutive_increase_reduction_factor;
int  PorousMedia::steady_limit_function_evals;
Real PorousMedia::steady_abs_tolerance;
Real PorousMedia::steady_rel_tolerance;
Real PorousMedia::steady_abs_update_tolerance;
Real PorousMedia::steady_rel_update_tolerance;
bool PorousMedia::steady_do_grid_sequence;
Array<Real> PorousMedia::steady_grid_sequence_new_level_dt_factor;
std::string PorousMedia::steady_record_file;

int  PorousMedia::richard_max_ls_iterations;
Real PorousMedia::richard_min_ls_factor;
Real PorousMedia::richard_ls_acceptance_factor;
Real PorousMedia::richard_ls_reduction_factor;
int  PorousMedia::richard_monitor_linear_solve;
int  PorousMedia::richard_monitor_line_search;
Real PorousMedia::richard_perturbation_scale_for_J;
bool PorousMedia::richard_use_fd_jac;
bool PorousMedia::richard_use_dense_Jacobian;
std::string PorousMedia::richard_rel_perm_method;
int  PorousMedia::richard_pressure_maxorder;
bool PorousMedia::richard_scale_solution_before_solve;
bool PorousMedia::richard_semi_analytic_J;
bool PorousMedia::richard_centered_diff_J;
Real PorousMedia::richard_variable_switch_saturation_threshold;
Real PorousMedia::richard_dt_thresh_pure_steady;

RichardSolver* PorousMedia::richard_solver;
NLScontrol* PorousMedia::richard_solver_control;
RSdata* PorousMedia::richard_solver_data;

PorousMedia::ExecutionMode PorousMedia::execution_mode;

namespace
{
    static void PM_Setup_CleanUpStatics()
    {
      ChemistryHelper_Structured *chemistry_helper = PorousMedia::GetChemistryHelper();
      delete chemistry_helper; chemistry_helper = 0;

#ifdef ALQUIMIA_ENABLED
      Amanzi::AmanziChemistry::ChemistryEngine *chemistry_engine = PorousMedia::GetChemistryEngine();
      delete chemistry_engine; chemistry_engine = 0;
#endif
    }
}


//
// Components are  Interior, Inflow, Outflow, Symmetry, SlipWall, NoSlipWall.
//

static int scalar_bc[] =
  {
    INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_EVEN, FOEXTRAP, SEEPAGE
    //INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_EVEN, REFLECT_ODD, SEEPAGE
  };

static int tracer_bc[] =
  {
    //INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_EVEN, REFLECT_EVEN, SEEPAGE
    INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_EVEN, FOEXTRAP, SEEPAGE
  };

static int press_bc[] =
  {
    // Interior, Inflow,   Outflow, Symmetry,     SlipWall, NoSlipWall.
    INT_DIR,     FOEXTRAP, EXT_DIR, REFLECT_EVEN, FOEXTRAP, FOEXTRAP
  };

static int norm_vel_bc[] =
  {
    INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_ODD, EXT_DIR, EXT_DIR
  };

static int tang_vel_bc[] =
  {
    INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_EVEN, HOEXTRAP, EXT_DIR
  };

static int aux_bc[] =
  {
    INT_DIR, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN
  };

static BCRec trac_bc; // Set in read_trac, used in variableSetUp

static
void
set_scalar_bc (BCRec&       bc,
               const BCRec& phys_bc)
{
  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();
  for (int i = 0; i < BL_SPACEDIM; i++)
    {
      bc.setLo(i,scalar_bc[lo_bc[i]]);
      bc.setHi(i,scalar_bc[hi_bc[i]]);
    }
}

static
void
set_tracer_bc (BCRec&       bc,
               const BCRec& phys_bc)
{
  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();
  for (int i = 0; i < BL_SPACEDIM; i++)
    {
      bc.setLo(i,tracer_bc[lo_bc[i]]);
      bc.setHi(i,tracer_bc[hi_bc[i]]);
    }
}

static
void
set_pressure_bc (BCRec&       bc,
                 const BCRec& phys_bc)
{
  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();
  for (int i = 0; i < BL_SPACEDIM; i++)
    {
      bc.setLo(i,press_bc[lo_bc[i]]);
      bc.setHi(i,press_bc[hi_bc[i]]);
    }
}

static
void
set_aux_bc (BCRec&       bc,
	    const BCRec& phys_bc)
{
  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();
  for (int i = 0; i < BL_SPACEDIM; i++)
    {
      bc.setLo(i,aux_bc[lo_bc[i]]);
      bc.setHi(i,aux_bc[hi_bc[i]]);
    }
}

static
void
set_x_vel_bc (BCRec&       bc,
              const BCRec& phys_bc)
{
  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();
  bc.setLo(0,norm_vel_bc[lo_bc[0]]);
  bc.setHi(0,norm_vel_bc[hi_bc[0]]);
  bc.setLo(1,tang_vel_bc[lo_bc[1]]);
  bc.setHi(1,tang_vel_bc[hi_bc[1]]);
#if (BL_SPACEDIM == 3)
  bc.setLo(2,tang_vel_bc[lo_bc[2]]);
  bc.setHi(2,tang_vel_bc[hi_bc[2]]);
#endif
}

static
void
set_y_vel_bc (BCRec&       bc,
              const BCRec& phys_bc)
{
  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();
  bc.setLo(0,tang_vel_bc[lo_bc[0]]);
  bc.setHi(0,tang_vel_bc[hi_bc[0]]);
  bc.setLo(1,norm_vel_bc[lo_bc[1]]);
  bc.setHi(1,norm_vel_bc[hi_bc[1]]);
#if (BL_SPACEDIM == 3)
  bc.setLo(2,tang_vel_bc[lo_bc[2]]);
  bc.setHi(2,tang_vel_bc[hi_bc[2]]);
#endif
}

#if (BL_SPACEDIM == 3)
static
void
set_z_vel_bc (BCRec&       bc,
              const BCRec& phys_bc)
{
  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();
  bc.setLo(0,tang_vel_bc[lo_bc[0]]);
  bc.setHi(0,tang_vel_bc[hi_bc[0]]);
  bc.setLo(1,tang_vel_bc[lo_bc[1]]);
  bc.setHi(1,tang_vel_bc[hi_bc[1]]);
  bc.setLo(2,norm_vel_bc[lo_bc[2]]);
  bc.setHi(2,norm_vel_bc[hi_bc[2]]);
}
#endif

typedef StateDescriptor::BndryFunc BndryFunc;
typedef ErrorRec::ErrorFunc ErrorFunc;

void
PorousMedia::InitializeStaticVariables ()
{
  //
  // Set all default values for static variables here!!!
  //
  PorousMedia::num_state_type = -1;
  PorousMedia::do_source_term = false;

  PorousMedia::flow_model     = PorousMedia::PM_FLOW_MODEL_INVALID;
  PorousMedia::flow_eval      = PorousMedia::PM_FLOW_EVAL_INVALID;
  PorousMedia::flow_is_static = true;
  PorousMedia::nphases        = 0;
  PorousMedia::ncomps         = 0;
  PorousMedia::ndiff          = 0;

  PorousMedia::ntracers = 0;
  PorousMedia::uninitialized_data = 1.0e30;
  PorousMedia::nminerals = 0;
  PorousMedia::minerals.clear();
  PorousMedia::nsorption_sites = 0;
  PorousMedia::sorption_sites.clear();
  PorousMedia::ncation_exchange = 0;
  PorousMedia::nsorption_isotherms = 0;
  PorousMedia::using_sorption = false;

  PorousMedia::wt_lo = 0;
  PorousMedia::wt_hi = 0;

  PorousMedia::temperature = 0;

  PorousMedia::verbose_observation_processing = 0;

  PorousMedia::verbose      = 0;
  PorousMedia::cfl          = 0.8;
  PorousMedia::init_shrink  = 1.0;
  PorousMedia::dt_init      = -1.0; // Ignore if < 0
  PorousMedia::dt_grow_max  = -1;
  PorousMedia::dt_shrink_max  = 10;
  PorousMedia::fixed_dt     = -1.0;
  PorousMedia::steady_max_dt = -1; // Ignore if < 0
  PorousMedia::transient_max_dt = -1; // Ignore if < 0
  PorousMedia::dt_cutoff    = 0.0;
  PorousMedia::gravity      = 9.807 / BL_ONEATM;
  PorousMedia::gravity_dir  = BL_SPACEDIM-1;
  PorousMedia::z_location   = 0;
  PorousMedia::domain_thickness = 1.0; // Not used in 3D
  PorousMedia::initial_step = false;
  PorousMedia::initial_iter = false;
  PorousMedia::sum_interval = 1;
  PorousMedia::NUM_SCALARS  = 0;
  PorousMedia::NUM_STATE    = 0;

  PorousMedia::be_cn_theta           = 0.5;
  PorousMedia::visc_tol              = 1.0e-10;
  PorousMedia::visc_abs_tol          = 1.0e-10;
  PorousMedia::def_harm_avg_cen2edge = true;

  PorousMedia::atmospheric_pressure_atm = 1;

  PorousMedia::variable_scal_diff = true;

  PorousMedia::do_tracer_chemistry = false;
  PorousMedia::do_tracer_advection = false;
  PorousMedia::do_tracer_diffusion = false;
  PorousMedia::setup_tracer_transport = false;
  PorousMedia::advect_tracers     = false;
  PorousMedia::diffuse_tracers    = false;
  PorousMedia::tensor_tracer_diffusion = false;
  PorousMedia::do_full_strang     = false;
  PorousMedia::n_chem_interval    = 0;
  PorousMedia::it_chem            = 0;
  PorousMedia::dt_chem            = 0;
  PorousMedia::max_grid_size_chem = 16;
  PorousMedia::no_initial_values  = true;
  PorousMedia::use_funccount      = false;
  PorousMedia::max_chemistry_time_step = -1;

  PorousMedia::do_reflux           = true;
  PorousMedia::execution_mode      = PorousMedia::INVALID;
  PorousMedia::ic_chem_relax_dt    = -1; // < 0 implies not done
  PorousMedia::solute_transport_limits_dt = false;
  PorousMedia::do_constant_vel = false;
  PorousMedia::nGrowHYP = 3;
  PorousMedia::nGrowMG = 1;
  PorousMedia::nGrowEIGEST = 1;
  PorousMedia::max_n_subcycle_transport = 10;
  PorousMedia::max_dt_iters_flow = 20;
  PorousMedia::abort_on_chem_fail = true;
  PorousMedia::show_selected_runtimes = false;
  PorousMedia::be_cn_theta_trac = 1.0;
  //PorousMedia::do_output_flow_time_in_years = true;
  PorousMedia::do_output_flow_time_in_years = false;
  PorousMedia::do_output_chemistry_time_in_years = false;
  PorousMedia::do_output_transport_time_in_years = false;

  PorousMedia::steady_min_iterations = 10;
  PorousMedia::steady_min_iterations_2 = 2;
  PorousMedia::steady_max_iterations = 15;
  PorousMedia::steady_limit_iterations = 20;
  PorousMedia::steady_time_step_reduction_factor = 0.8;
  PorousMedia::steady_time_step_increase_factor = 1.6;
  PorousMedia::steady_time_step_increase_factor_2 = 10;
  PorousMedia::steady_time_step_retry_factor_1 = 0.2;
  PorousMedia::steady_time_step_retry_factor_2 = 0.01;
  PorousMedia::steady_time_step_retry_factor_f = 0.001;
  PorousMedia::steady_max_consecutive_failures_1 = 3;
  PorousMedia::steady_max_consecutive_failures_2 = 4;
  PorousMedia::steady_max_time_steps = 8000;
  PorousMedia::steady_max_time_step_size = 1.e20;
  PorousMedia::steady_max_num_consecutive_success = 0;
  PorousMedia::steady_extra_time_step_increase_factor = 10.;
  PorousMedia::steady_max_num_consecutive_increases = 3;
  PorousMedia::steady_consecutive_increase_reduction_factor = 0.4;
  PorousMedia::steady_limit_function_evals = 1e8;
  PorousMedia::steady_abs_tolerance = 1.e-10;
  PorousMedia::steady_rel_tolerance = 1.e-20;
  PorousMedia::steady_abs_update_tolerance = 1.e-12;
  PorousMedia::steady_rel_update_tolerance = -1;
  PorousMedia::steady_do_grid_sequence = true;
  PorousMedia::steady_grid_sequence_new_level_dt_factor.resize(1,1);
  PorousMedia::steady_record_file.clear();

  PorousMedia::richard_max_ls_iterations = 10;
  PorousMedia::richard_min_ls_factor = 1.e-8;
  PorousMedia::richard_ls_acceptance_factor = 1.4;
  PorousMedia::richard_ls_reduction_factor = 0.1;
  PorousMedia::richard_monitor_linear_solve = 0;
  PorousMedia::richard_monitor_line_search = 0;
  PorousMedia::richard_perturbation_scale_for_J = 1.e-8;
  PorousMedia::richard_use_fd_jac = true;
  PorousMedia::richard_use_dense_Jacobian = false;
  PorousMedia::richard_rel_perm_method = "upwind-darcy_velocity";
  PorousMedia::richard_pressure_maxorder = 3;
  PorousMedia::richard_scale_solution_before_solve = true;
  PorousMedia::richard_semi_analytic_J = false;
  PorousMedia::richard_centered_diff_J = true;
  PorousMedia::richard_variable_switch_saturation_threshold = -1;
  PorousMedia::richard_dt_thresh_pure_steady = -1;

  PorousMedia::echo_inputs    = 0;
  PorousMedia::richard_solver = 0;
  PorousMedia::richard_solver_control = 0;
  PorousMedia::richard_solver_data = 0;

#if ALQUIMIA_ENABLED
  PorousMedia::chemistry_engine = 0;
#endif
  PorousMedia::chemistry_engine_name = "";
  PorousMedia::chemistry_helper = 0;
}

std::pair<std::string,std::string>
SplitDirAndName(const std::string& orig)
{
    if (orig[orig.length()-1] == '/') {
        BoxLib::Abort(std::string("Invalid filename:" + orig).c_str());
    }
    vector<std::string> tokens = BoxLib::Tokenize(orig,std::string("/"));
    std::pair<std::string,std::string> result;
    int size = tokens.size();
    BL_ASSERT(tokens.size()>0);
    if (size>1) {
        for (int i=0; i<size-2; ++i) {
            result.first += tokens[i] + "/";
        }
        result.first += tokens[size-2];
    }
    else {
        result.first = ".";
    }
    result.second = tokens[size-1];
    return result;
}

void
PorousMedia::variableSetUp ()
{

  InitializeStaticVariables();
  ParmParse pproot;
  pproot.query("echo_inputs",echo_inputs);

  BL_ASSERT(desc_lst.size() == 0);

  for (int dir = 0; dir < BL_SPACEDIM; dir++)
  {
    phys_bc.setLo(dir,SlipWall);
    phys_bc.setHi(dir,SlipWall);
  }

  std::string pp_dump_file = "";
  if (pproot.countval("dump_parmparse_table")) {
      pproot.get("dump_parmparse_table",pp_dump_file);
      std::ofstream ofs;
      std::pair<std::string,std::string> df = SplitDirAndName(pp_dump_file);
      if (ParallelDescriptor::IOProcessor()) {
          if (!BoxLib::UtilCreateDirectory(df.first, 0755)) {
              BoxLib::CreateDirectoryFailed(df.first);
          }

          ofs.open(pp_dump_file.c_str());
          if (ofs.fail()) {
              BoxLib::Abort(std::string("Cannot open pp dump file: "+pp_dump_file).c_str());
          }
          if (verbose>1)
          {
              std::cout << "\nDumping ParmParse table:\n";
          }

          // NOTE: Formatting useless since all data are strings at this point
          //
          // std::ios::fmtflags oflags = ofs.flags();
          // ofs.setf(std::ios::floatfield, std::ios::scientific);
          // int old_prec = ofs.precision(15);

          bool prettyPrint = false;
          ParmParse::dumpTable(ofs,prettyPrint);

          // ofs.flags(oflags);
          // ofs.precision(old_prec);

          if (!ofs.good())
              BoxLib::Error("Write of pp dump file failed");


          if (verbose>1)
          {
              std::cout << "... done dumping ParmParse table.\n" << '\n';
          }
          ofs.close();
      }
      ParallelDescriptor::Barrier();
  }

  read_params();
  BCRec bc;

  //
  // Set state variables Ids.
  //
  // NUM_SCALARS = ncomps + 2; // Currently unused last 2 components
  NUM_SCALARS = ncomps;

  if (ntracers > 0)
    NUM_SCALARS = NUM_SCALARS + ntracers;

  // No longer supported
  // if (flow_model == PM_FLOW_MODEL_POLYMER)
  // {
  //   NUM_SCALARS = NUM_SCALARS + 2;
  // }

  // add velocity and correction velocity
  NUM_STATE = NUM_SCALARS + BL_SPACEDIM + BL_SPACEDIM ;

  //
  // **************  DEFINE SCALAR VARIABLES  ********************
  //

  Array<BCRec>       bcs(ncomps);
  Array<std::string> names(ncomps);

  desc_lst.addDescriptor(State_Type,IndexType::TheCellType(),
			 StateDescriptor::Point,1,NUM_SCALARS,
			 &cell_cons_interp);

  set_scalar_bc(bc,phys_bc);
  for (int i = 0; i < ncomps; i++)
  {
    bcs[i] = bc;
    names[i] = cNames[i];
  }

  desc_lst.setComponent(State_Type,
			0,
			names,
			bcs,
			BndryFunc(FORT_ONE_N_FILL,FORT_ALL_N_FILL));

  if (ntracers > 0)
  {
    Array<BCRec>       tbcs(ntracers);
    Array<std::string> tnames(ntracers);

    for (int i = 0; i < ntracers; i++)
    {
      tbcs[i]   = trac_bc;
      tnames[i] = tNames[i];
    }

    desc_lst.setComponent(State_Type,
			  ncomps,
			  tnames,
			  tbcs,
			  BndryFunc(FORT_ONE_N_FILL,FORT_ALL_T_FILL));
  }

#if 0
  // Currently unused
  desc_lst.setComponent(State_Type,ncomps+ntracers,"Aux1",
			bc,BndryFunc(FORT_ENTHFILL));
  desc_lst.setComponent(State_Type,ncomps+ntracers+1,"Aux2",
			bc,BndryFunc(FORT_ADVFILL));
#endif

  // FIXME: No longer supported
  // if (flow_model == PM_FLOW_MODEL_POLYMER) {
  //   desc_lst.setComponent(State_Type,ncomps+2,"s",
  // 			  bc,BndryFunc(FORT_ONE_N_FILL));
  //   desc_lst.setComponent(State_Type,ncomps+3,"c",
  // 			  bc,BndryFunc(FORT_ONE_N_FILL));
  // }

  if (chemistry_helper != 0) {
    const std::map<std::string,int>& aux_chem_variables_map = chemistry_helper->AuxChemVariablesMap();
    int num_aux_chem_variables = aux_chem_variables_map.size();
    if (num_aux_chem_variables > 0)
    {
      Array<BCRec> cbcs(num_aux_chem_variables);
      set_aux_bc(bc,phys_bc);
      Array<std::string> tmp_aux(num_aux_chem_variables);
      for (std::map<std::string,int>::const_iterator it=aux_chem_variables_map.begin();
	   it!=aux_chem_variables_map.end(); ++it)
      {
	int i = it->second;
	tmp_aux[i] = it->first;
	cbcs[i] = bc;
      }

      FORT_AUXPARAMS(&num_aux_chem_variables);

      desc_lst.addDescriptor(Aux_Chem_Type,IndexType::TheCellType(),
                             StateDescriptor::Point,0,num_aux_chem_variables,
                             &cell_cons_interp);
      desc_lst.setComponent(Aux_Chem_Type,0,tmp_aux,cbcs,
                            BndryFunc(FORT_ONE_A_FILL,FORT_ALL_A_FILL));
    }
  }

  //
  // **************  DEFINE PRESSURE VARIABLE  ********************
  //

  desc_lst.addDescriptor(Press_Type,IndexType::TheCellType(),
			 StateDescriptor::Point,1,1,
			 &cell_cons_interp);
  set_pressure_bc(bc,pres_bc);
  desc_lst.setComponent(Press_Type,Pressure,"pressure",
			bc,BndryFunc(FORT_PRESFILL));

  //
  // **************  DEFINE VELOCITY VARIABLES  ********************

  desc_lst.addDescriptor(Vcr_Type,IndexType::TheCellType(),
			 StateDescriptor::Point,1,BL_SPACEDIM,
			 &cell_cons_interp);
  set_x_vel_bc(bc,phys_bc);
  desc_lst.setComponent(Vcr_Type,Xvcr,"x_vcorr",
			bc,BndryFunc(FORT_XVELFILL));
  set_y_vel_bc(bc,phys_bc);
  desc_lst.setComponent(Vcr_Type,Yvcr,"y_vcorr",
			bc,BndryFunc(FORT_YVELFILL));
#if (BL_SPACEDIM == 3)
  set_z_vel_bc(bc,phys_bc);
  desc_lst.setComponent(Vcr_Type,Zvcr,"z_vcorr",
			bc,BndryFunc(FORT_ZVELFILL));
#endif

#if defined(AMANZI)
  if (do_tracer_chemistry>0)
    {
      // add function count
      int nfunccountghost = 0;
      if (do_full_strang) nfunccountghost=1;
      desc_lst.addDescriptor(FuncCount_Type, IndexType::TheCellType(),
			     StateDescriptor::Point,nfunccountghost,1, &cell_cons_interp);
      desc_lst.setComponent(FuncCount_Type, 0, "FuncCount",
			    bc, BndryFunc(FORT_ONE_A_FILL));
    }
#endif

  // "User defined" - although these must correspond to those in PorousMedia::derive
  IndexType regionIDtype(IndexType::TheCellType());
  std::string amr_prefix = "amr";
  ParmParse pp(amr_prefix);
  int num_user_derives = pp.countval("user_derive_list");
  if (num_user_derives != 0) {
    pp.getarr("user_derive_list",user_derive_list,0,num_user_derives);
  }

  // Make list of internally required or otherwise desired variables
  Array<std::string> intern_reqd_der;
  for (int i=0; i<pNames.size(); ++i) {
    intern_reqd_der.push_back(pNames[i] + "_Pressure");
  }
  for (int i=0; i<cNames.size(); ++i) {
    intern_reqd_der.push_back(cNames[i] + "_Saturation");
    intern_reqd_der.push_back("Volumetric_" + cNames[i] + "_Content");
  }
  intern_reqd_der.push_back("Capillary_Pressure");
  for (int i=0; i<tNames.size(); ++i) {
    for (int j=0; j<pNames.size(); ++j) {
      for (int k=0; k<cNames.size(); ++k) {
	if (pNames.size() == 1 &&  cNames.size() == 1) {
	  intern_reqd_der.push_back(tNames[i] + "_" + pNames[k] + "_Concentration");
	}
	else {
	  intern_reqd_der.push_back(tNames[i] + "_Concentration_in_" + pNames[k] + "_" + cNames[j]);
	}
      }
    }
    intern_reqd_der.push_back("Volumetric_" + tNames[i] + "_Content");
  }

  intern_reqd_der.push_back("Porosity");
  std::string dirStr[3] = {"X", "Y", "Z"};
  for (int i=0; i<cNames.size(); ++i) {
    for (int d=0; d<BL_SPACEDIM; ++d) {
      intern_reqd_der.push_back("Volumetric_"+cNames[i]+"_Flux_"+dirStr[d]);
    }
  }

  if (do_tracer_diffusion) {
    for (int d=0; d<BL_SPACEDIM; ++d) {
      intern_reqd_der.push_back("Tortuosity_"+dirStr[d]);
    }
  }

  if (flow_model==PM_FLOW_MODEL_SATURATED) {
    intern_reqd_der.push_back("Specific_Storage");
    intern_reqd_der.push_back("Specific_Yield");
    intern_reqd_der.push_back("Particle_Density");
  }

  for (int d=0; d<BL_SPACEDIM; ++d) {
    intern_reqd_der.push_back("Intrinsic_Permeability_"+dirStr[d]);
  }

  if (rock_manager!=0  && do_tracer_chemistry!=0)
  {
    // primary species concentration required
    for (int i=0; i<tNames.size(); ++i) {
      intern_reqd_der.push_back(tNames[i]+"_Free_Ion_Guess");
    }

    // SM: use PM_setup static variables in lieu of RockManager
    //     this fixes incomplete translation from 2.x > parmparse format
    //     and disconnect between RockManager and output -- with the
    //     advantage that chemistry does not need to be in native parmparse

    // sorption
    //  if (rock_manager->UsingSorption())
    if (using_sorption)
    {
      for (int i=0; i<tNames.size(); ++i) {
	intern_reqd_der.push_back(tNames[i]+"_Sorbed_Concentration");
      }

      // cation exchange
      //if (rock_manager->CationExchangeCapacityICs().size()>0) {
      if (ncation_exchange>0) {
	intern_reqd_der.push_back("Cation_Exchange_Capacity");
      }

      // surface complexation
      //const Array<std::string>& sorptionSiteNames = rock_manager->SorptionSiteNames();
      const Array<std::string>& sorptionSiteNames = sorption_sites;
      if (sorptionSiteNames.size()>0) {
	for (int i=0; i<sorptionSiteNames.size(); ++i) {
	  intern_reqd_der.push_back(sorptionSiteNames[i]+"_Surface_Site_Density");
	}
      }
    } //sorption

    // minerals
    //const Array<std::string>& mineralNames = rock_manager->MineralNames();
    const Array<std::string>& mineralNames = minerals;
    for (int i=0; i<mineralNames.size(); ++i) {
      intern_reqd_der.push_back(mineralNames[i]+"_Volume_Fraction");
    }
    for (int i=0; i<mineralNames.size(); ++i) {
      intern_reqd_der.push_back(mineralNames[i]+"_Specific_Surface_Area");
    }
  }

  intern_reqd_der.push_back("Material_ID");
  intern_reqd_der.push_back("Grid_ID");
  intern_reqd_der.push_back("Core_ID");
  intern_reqd_der.push_back("Cell_ID");
  if (gravity != 0) {
    intern_reqd_der.push_back("Hydraulic_Head");
  }

  for (int i=0; i<intern_reqd_der.size(); ++i) {
    int j = -1;
    for (int k=0; k<user_derive_list.size() && j<0; ++k) {
      j = intern_reqd_der[i] == user_derive_list[k] ? k : -1;
    }
    if (j<0) {
      user_derive_list.push_back(intern_reqd_der[i]);
    }
  }

  // Add variables in list above to derivable list if not already added
  for (int i=0; i<user_derive_list.size(); ++i) {
    int nCompThis = (user_derive_list[i] == "Dispersivity" ? 2 : 1);
    derive_lst.add(user_derive_list[i], regionIDtype, nCompThis);
  }

  int nwr = pp.countval("write_regions");
  if (nwr > 0) {
    Array<std::string> wrNames;
    pp.getarr("write_regions",wrNames,0,nwr);
    std::string prefix("amr.write_region");
    ParmParse ppwr(prefix.c_str());

    for (int i = 0; i<nwr; i++) {
      const std::string& wrname = wrNames[i];
      int nwrst = ppwr.countval(wrname.c_str());
      if (nwrst == 0) {
	BoxLib::Abort(std::string(prefix+"."+wrname+" = <region names> required").c_str());
      }
      Array<std::string> wrRegions;
      ppwr.getarr(wrname.c_str(),wrRegions,0,nwrst);
      write_region_sets[wrname] = wrRegions;
    }
  }


  //
  // **************  DEFINE ERROR ESTIMATION QUANTITIES  *************
  //
  if (region_manager == 0) {
    BoxLib::Abort("static Region manager must be set up prior to reading AMR refinement indicators");
  }

  Array<std::string> refinement_indicators;
  pp.queryarr("refinement_indicators",refinement_indicators,0,pp.countval("refinement_indicators"));
  for (int i=0; i<refinement_indicators.size(); ++i)
  {
      std::string ref_prefix = amr_prefix + "." + refinement_indicators[i];
      ParmParse ppr(ref_prefix);
      Real min_time = 0; ppr.query("start_time",min_time);
      Real max_time = -1; ppr.query("end_time",max_time);
      int max_level = -1;  ppr.query("max_level",max_level);
      Array<std::string> region_names(1,"All");
      int nreg = ppr.countval("regions");
      if (nreg) {
          ppr.getarr("regions",region_names,0,nreg);
      }
      Array<const Region*> regions = region_manager->RegionPtrArray(region_names);
      if (ppr.countval("value_greater")) {
          Real value; ppr.get("value_greater",value);
          std::string field; ppr.get("field_name",field);
          err_list.add(field.c_str(),0,ErrorRec::Special,
                       PM_Error_Value(FORT_VALGTERROR,value,min_time,max_time,max_level,regions));
      }
      else if (ppr.countval("value_less")) {
          Real value; ppr.get("value_less",value);
          std::string field; ppr.get("field_name",field);
          err_list.add(field.c_str(),0,ErrorRec::Special,
                       PM_Error_Value(FORT_VALLTERROR,value,min_time,max_time,max_level,regions));
      }
      else if (ppr.countval("adjacent_difference_greater")) {
          Real value; ppr.get("adjacent_difference_greater",value);
          std::string field; ppr.get("field_name",field);
          err_list.add(field.c_str(),1,ErrorRec::Special,
                       PM_Error_Value(FORT_DIFFGTERROR,value,min_time,max_time,max_level,regions));
      }
      else if (ppr.countval("inside_region")) {
	  //Real value; ppr.get("inside_region",value);
          err_list.add("PMAMR_DUMMY",1,ErrorRec::Special,
                       PM_Error_Value(min_time,max_time,max_level,regions));
      }
      else {
          BoxLib::Abort(std::string("Unrecognized refinement indicator for " + refinement_indicators[i]).c_str());
      }
  }

  num_state_type = desc_lst.size();

  BoxLib::ExecOnFinalize(PM_Setup_CleanUpStatics);

}

void PorousMedia::read_prob()
{
  ParmParse pp;

  pp.query("do_output_flow_time_in_years;",do_output_flow_time_in_years);
  pp.query("do_output_transport_time_in_years;",do_output_transport_time_in_years);
  pp.query("do_output_chemistry_time_in_years;",do_output_chemistry_time_in_years);

  ParmParse pb("prob");
  std::string flow_state_str;
  pb.get("flow_state",flow_state_str);

  if (flow_state_str == "on") {
    std::string flow_model_str;
    pb.get("flow_model",flow_model_str);
    if (flow_model_str == "richards") {
      flow_model = PM_FLOW_MODEL_RICHARDS;
      flow_eval  = PM_FLOW_EVAL_EVOLVE;
      flow_is_static = false;
    } else if (flow_model_str == "saturated") {
      flow_model = PM_FLOW_MODEL_SATURATED;
      flow_eval  = PM_FLOW_EVAL_EVOLVE;
      flow_is_static = false;
    } else if (flow_model_str == "constant") {
      flow_model = PM_FLOW_MODEL_SATURATED;

      // FIXME: How to decide here?  overwrite later when we get proper info
      flow_eval  = PM_FLOW_EVAL_SOLVE_GIVEN_PBC;
      flow_is_static = true;
    } else {
      flow_model = PM_FLOW_MODEL_INVALID;
      flow_eval  = PM_FLOW_EVAL_INVALID;
      flow_is_static = true;
    }
  }
  else {
    flow_model = PM_FLOW_MODEL_OFF; // FIXME: Richards or saturated??
    flow_eval  = PM_FLOW_EVAL_CONSTANT;
    flow_is_static = true;
  }

  if (flow_model==PM_FLOW_MODEL_SATURATED || flow_model==PM_FLOW_MODEL_OFF) {
    solute_transport_limits_dt = true;
  }

  pb.query("do_tracer_advection",do_tracer_advection);
  pb.query("do_tracer_diffusion",do_tracer_diffusion);
  if (do_tracer_advection || do_tracer_diffusion) {
      setup_tracer_transport = true; // NOTE: May want these data structures regardless...
  }

  if (setup_tracer_transport &&
      ( flow_model==PM_FLOW_MODEL_SATURATED
	|| flow_model == PM_FLOW_MODEL_RICHARDS) )
  {
      advect_tracers = do_tracer_advection;
      diffuse_tracers = do_tracer_diffusion;
      react_tracers = do_tracer_chemistry;
  }

  // Verbosity
  pb.query("v",verbose);

  // Get timestepping parameters.  Some will be used to default values for int-to-steady solver
  pb.get("cfl",cfl);
  pb.query("init_shrink",init_shrink);
  pb.query("dt_init",dt_init);
  pb.query("dt_cutoff",dt_cutoff);
  pb.query("dt_grow_max",dt_grow_max);
  pb.query("dt_shrink_max",dt_shrink_max);
  pb.query("fixed_dt",fixed_dt);
  pb.query("steady_max_dt",steady_max_dt);
  pb.query("transient_max_dt",transient_max_dt);
  pb.query("sum_interval",sum_interval);
  pb.query("max_n_subcycle_transport",max_n_subcycle_transport);

  pb.query("max_dt_iters_flow",max_dt_iters_flow);
  pb.query("show_selected_runtimes",show_selected_runtimes);
  pb.query("abort_on_chem_fail",abort_on_chem_fail);

  pb.query("steady_record_file",steady_record_file);
  pb.query("steady_min_iterations",steady_min_iterations);
  pb.query("steady_min_iterations_2",steady_min_iterations_2);
  pb.query("steady_max_iterations",steady_max_iterations);
  pb.query("steady_limit_iterations",steady_limit_iterations);
  pb.query("steady_time_step_reduction_factor",steady_time_step_reduction_factor);
  pb.query("steady_time_step_increase_factor_2",steady_time_step_increase_factor_2);
  pb.query("steady_time_step_increase_factor",steady_time_step_increase_factor);
  pb.query("steady_time_step_retry_factor_1",steady_time_step_retry_factor_1);
  pb.query("steady_time_step_retry_factor_2",steady_time_step_retry_factor_2);
  pb.query("steady_time_step_retry_factor_f",steady_time_step_retry_factor_f);
  pb.query("steady_max_consecutive_failures_1",steady_max_consecutive_failures_1);
  pb.query("steady_max_consecutive_failures_2",steady_max_consecutive_failures_2);
  pb.query("steady_max_time_steps",steady_max_time_steps);
  pb.query("steady_max_time_step_size",steady_max_time_step_size);
  pb.query("steady_max_num_consecutive_success",steady_max_num_consecutive_success);
  pb.query("steady_extra_time_step_increase_factor",steady_extra_time_step_increase_factor);
  pb.query("steady_max_num_consecutive_increases",steady_max_num_consecutive_increases);
  pb.query("steady_consecutive_increase_reduction_factor",steady_consecutive_increase_reduction_factor);
  pb.query("steady_limit_function_evals",steady_limit_function_evals);
  pb.query("steady_abs_tolerance",steady_abs_tolerance);
  pb.query("steady_rel_tolerance",steady_rel_tolerance);
  pb.query("steady_abs_update_tolerance",steady_abs_update_tolerance);
  pb.query("steady_rel_update_tolerance",steady_rel_update_tolerance);
  pb.query("steady_do_grid_sequence",steady_do_grid_sequence);
  int ndt = pb.countval("steady_grid_sequence_new_level_dt_factor");
  if (ndt > 0) {
      pb.getarr("steady_grid_sequence_new_level_dt_factor",steady_grid_sequence_new_level_dt_factor,0,ndt);
  }

  pb.query("richard_max_ls_iterations",richard_max_ls_iterations);
  pb.query("richard_min_ls_factor",richard_min_ls_factor);
  pb.query("richard_ls_acceptance_factor",richard_ls_acceptance_factor);
  pb.query("richard_ls_reduction_factor",richard_ls_reduction_factor);
  pb.query("richard_monitor_linear_solve",richard_monitor_linear_solve);
  pb.query("richard_monitor_line_search",richard_monitor_line_search);
  pb.query("richard_perturbation_scale_for_J",richard_perturbation_scale_for_J);
  pb.query("richard_use_fd_jac",richard_use_fd_jac);
  pb.query("richard_use_dense_Jacobian",richard_use_dense_Jacobian);
  if (flow_model == PM_FLOW_MODEL_SATURATED) {
    richard_rel_perm_method = "other-harmonic_average";
  }
  pb.query("richard_rel_perm_method",richard_rel_perm_method);
  pb.query("richard_pressure_maxorder",richard_pressure_maxorder);
  pb.query("richard_scale_solution_before_solve",richard_scale_solution_before_solve);
  pb.query("richard_semi_analytic_J",richard_semi_analytic_J);
  pb.query("richard_centered_diff_J",richard_centered_diff_J);
  pb.query("richard_variable_switch_saturation_threshold",richard_variable_switch_saturation_threshold);
  pb.query("richard_dt_thresh_pure_steady",richard_dt_thresh_pure_steady);

  // Gravity are specified as m/s^2 in the input file
  // This is converted to the unit that is used in the code.
  if (pb.contains("gravity")) {
    pb.get("gravity",gravity);
    gravity /= BL_ONEATM;
  }
  pb.query("gravity_dir",gravity_dir);
  BL_ASSERT(gravity_dir>=0 && gravity_dir<3); // Note: can set this to 2 for a 2D problem
  if (BL_SPACEDIM<3 && gravity_dir>BL_SPACEDIM-1) {
    pb.query("z_location",z_location);
  }
  if (BL_SPACEDIM<3) {
    pb.query("domain_thickness",domain_thickness);
    if (domain_thickness <= 0) {
      BoxLib::Abort("domain_thickness, if specified, must be > 0");
    }
  }
  if (pb.countval("richard_atmospheric_pressure")) {
    pb.get("richard_atmospheric_pressure",atmospheric_pressure_atm);
    atmospheric_pressure_atm *= 1 / BL_ONEATM;
  }

  // Get algorithmic flags and options
  pb.query("do_reflux",  do_reflux );

  // Get solver tolerances
  pb.query("visc_tol",visc_tol);
  pb.query("visc_abs_tol",visc_abs_tol);
  pb.query("be_cn_theta",be_cn_theta);
  if (be_cn_theta > 1.0 || be_cn_theta < .5)
    BoxLib::Abort("PorousMedia::Must have be_cn_theta <= 1.0 && >= .5");
  pb.query("be_cn_theta_trac",be_cn_theta_trac);
  if (be_cn_theta > 1.0 || be_cn_theta < 0)
    BoxLib::Abort("PorousMedia::Must have be_cn_theta_trac <= 1.0 && >= 0");
  pb.query("harm_avg_cen2edge", def_harm_avg_cen2edge);

  pb.query("verbose_observation_processing",verbose_observation_processing);
}

//
// Construct bc functions
//

struct PressToRhoSat
    : public ArrayTransform
{
    PressToRhoSat() {}
    virtual ArrayTransform* clone() const {return new PressToRhoSat(*this);}
    virtual Array<Real> transform(Real inData) const;
protected:
};

Array<Real>
PressToRhoSat::transform(Real aqueous_pressure) const
{
    // FIXME: Requires Water
    const Array<std::string>& cNames = PorousMedia::componentNames();
    const Array<Real>& density = PorousMedia::Density();

    int ncomps = cNames.size();
    int idx = -1;
    for (int j=0; j<ncomps; ++j) {
        if (cNames[j] == "Water" || cNames[j] == "water") {
            idx = j;
        }
    }
    BL_ASSERT(idx>=0);

    Array<double> rhoSat(ncomps,0);
    rhoSat[idx] = density[idx] * 1; // Fully saturated...an assumption
    return rhoSat;
}

void  PorousMedia::read_comp()
{
  //
  // Read in parameters for phases and components
  //
  ParmParse pp("phase");

  // Get number and names of phases
  nphases = pp.countval("phases");
  pp.getarr("phases",pNames,0,nphases);
  for (int i = 0; i<nphases; i++) phase_list[pNames[i]] = i;

  // Build flattened list of components
  ndiff = 0;
  for (int i = 0; i<nphases; i++) {
    const std::string prefix("phase." + pNames[i]);
    ParmParse ppr(prefix.c_str());
    int p_nc = ppr.countval("comps");
    BL_ASSERT(p_nc==1); // An assumption all over the place...
    ncomps += p_nc;
    Array<std::string> p_cNames; ppr.getarr("comps",p_cNames,0,p_nc);
    for (int j=0; j<p_cNames.size(); ++j) {
      cNames.push_back(p_cNames[j]);
    }
    Real p_rho; ppr.get("density",p_rho); density.push_back(p_rho);
    // viscosity in units of kg/(m.s)
    Real p_visc; ppr.get("viscosity",p_visc); muval.push_back(p_visc);
    Real p_diff; ppr.get("diffusivity",p_diff); visc_coef.push_back(p_diff);

    // Tracer diffusion handled during tracer read
    if (visc_coef.back() > 0)
    {
      is_diffusive[visc_coef.size()-1] = 1;
    }
    else {
      variable_scal_diff = false;
    }
    ++ndiff;

    pType.push_back(phase_list[pNames[i]]);
  }

  ParmParse cp("comp");
  for (int i = 0; i<ncomps; i++) comp_list[cNames[i]] = i;

  //
  // Initial condition and boundary condition
  //
  // Component ics, bcs will be set all at once
  int n_ics = cp.countval("ic_labels");
  if (n_ics > 0)
  {
    Array<std::string> ic_names;
    cp.getarr("ic_labels",ic_names,0,n_ics);
    ic_array.resize(n_ics,PArrayManage);
    do_constant_vel = false;
    if (region_manager == 0) {
      BoxLib::Abort("static Region manager must be set up prior to reading tracer ICs");
    }
    for (int i = 0; i<n_ics; i++)
    {
      const std::string& icname = ic_names[i];
      const std::string prefix("comp.ics." + icname);
      ParmParse ppr(prefix.c_str());

      int n_ic_regions = ppr.countval("regions");
      Array<std::string> region_names;
      ppr.getarr("regions",region_names,0,n_ic_regions);
      Array<const Region*> ic_regions = region_manager->RegionPtrArray(region_names);

      std::string ic_type; ppr.get("type",ic_type);
      BL_ASSERT(!do_constant_vel); // If this is ever set, it must be the only IC so we should never see this true here
      if (ic_type == "uniform_pressure")
      {
	int nPhase = pNames.size();
	Array<Real> vals(nPhase);

	int num_phases_reqd = nPhase;
	std::map<std::string,bool> phases_set;
	for (int j = 0; j<pNames.size(); j++)
	{
	  std::string val_name = "val";
	  ppr.get(val_name.c_str(),vals[0]);
	  phases_set[pNames[j]] = true;
	}

	// convert to atm
	for (int j=0; j<vals.size(); ++j) {
	  vals[j] = vals[j] / BL_ONEATM;
	}

	int num_phases = phases_set.size();
	if (num_phases != num_phases_reqd) {
	  std::cerr << icname << ": Insufficient number of phases specified" << std::endl;
	  std::cerr << " ngiven, nreqd: " << num_phases << ", " << num_phases_reqd << std::endl;
	  std::cerr << " current flow_model: " << flow_model << std::endl;
	  BoxLib::Abort();
	}

	Array<Real> times(1,0);
	Array<std::string> forms(0);
	ic_array.set(i, new ArrayRegionData(icname,times,vals,forms,ic_regions,ic_type,1));
	if (flow_is_static) flow_eval  = PM_FLOW_EVAL_COMPUTE_GIVEN_P;
      }
      else if (ic_type == "linear_pressure")
      {
	int nPhase = pNames.size();
	if (nPhase!=1) {
	  std::cerr << "Multiphase not currently surrported" << std::endl;
	  BoxLib::Abort();
	}

	Real press_val;
	std::string val_name = "val";
	ppr.get(val_name.c_str(),press_val);
	press_val = press_val / BL_ONEATM;

	int ngrad = ppr.countval("grad");
	if (ngrad<BL_SPACEDIM) {
	  std::cerr << "Insufficient number of components given for pressure gradient" << std::endl;
	  BoxLib::Abort();
	}
	Array<Real> pgrad(BL_SPACEDIM);
	ppr.getarr("grad",pgrad,0,ngrad);
	for (int j=0; j<pgrad.size(); ++j) {
	  pgrad[j] = pgrad[j] / BL_ONEATM;
	}

	int nref = ppr.countval("loc");
	if (nref<BL_SPACEDIM) {
	  if (ParallelDescriptor::IOProcessor()) {
	    std::cerr << "Insufficient number of components given for pressure reference location" << std::endl;
	  }
	  BoxLib::Abort();
	}
	Array<Real> pref(BL_SPACEDIM);
	ppr.getarr("loc",pref,0,nref);

	int ntmp = 2*BL_SPACEDIM+1;
	Array<Real> tmp(ntmp);
	tmp[0] = press_val;
	for (int j=0; j<BL_SPACEDIM; ++j) {
	  tmp[1+j] = pgrad[j];
	  tmp[1+j+BL_SPACEDIM] = pref[j];
	}
	ic_array.set(i, new RegionData(icname,ic_regions,ic_type,tmp));
	if (flow_is_static) flow_eval  = PM_FLOW_EVAL_COMPUTE_GIVEN_P;
      }
      else if (ic_type == "uniform_saturation")
      {
	Array<Real> vals(ncomps);
	for (int j = 0; j<cNames.size(); j++) {
	  ppr.get(cNames[j].c_str(),vals[j]);
	  vals[j] *= density[j];
	}
	std::string generic_type = "scalar";
	ic_array.set(i, new RegionData(icname,ic_regions,generic_type,vals));
      }
      else if (ic_type == "velocity")
      {
	if ((flow_eval != PM_FLOW_EVAL_CONSTANT)
	    || (flow_model != PM_FLOW_MODEL_SATURATED)
	    || (flow_is_static != true)) {
	  if (ParallelDescriptor::IOProcessor()) {
	    std::cerr << "constant-velocity settings may only be used with steady-saturated flow" << std::endl;
	    BoxLib::Abort();
	  }
	}
	Array<Real> vals(BL_SPACEDIM);
	ppr.getarr("vel",vals,0,BL_SPACEDIM);
	std::string generic_type = "velocity";
	do_constant_vel = true;

	flow_eval  = PM_FLOW_EVAL_CONSTANT;

	ic_array.set(i, new RegionData(icname,ic_regions,generic_type,vals));
      }
      else if (ic_type == "hydrostatic")
      {
	Array<Real> water_table_height(1); ppr.get("water_table_height",water_table_height[0]);
	Array<Real> times(1,0);
	Array<std::string> forms;
	ic_array.set(i, new ArrayRegionData(icname,times,water_table_height,
					    forms,ic_regions,ic_type,1));
	if (flow_is_static) flow_eval  = PM_FLOW_EVAL_COMPUTE_GIVEN_P;
      }
      else if (ic_type == "zero_total_velocity")
      {
	Array<Real> vals(4);
	Array<Real> times(1,0);
	Array<std::string> forms;
	ppr.get("aqueous_vol_flux",vals[0]);
	ppr.get("water_table_height",vals[1]);
	ppr.query("val",vals[2]);
	ppr.query("grad",vals[3]);
	ic_array.set(i,new RegionData(icname,ic_regions,ic_type,vals));
	if (flow_is_static) flow_eval  = PM_FLOW_EVAL_COMPUTE_GIVEN_P;
      }
      else {
	BoxLib::Abort(std::string("Unsupported comp ic: \""+ic_type+"\"").c_str());
      }
    }
  }

  // default to no flow first.
  for (int j=0;j<BL_SPACEDIM;j++) {
    phys_bc.setLo(j,Symmetry);
    pres_bc.setLo(j,Symmetry);
    phys_bc.setHi(j,Symmetry);
    pres_bc.setHi(j,Symmetry);
  }
  rinflow_bc_lo.resize(BL_SPACEDIM,0);
  rinflow_bc_hi.resize(BL_SPACEDIM,0);
  inflow_bc_lo.resize(BL_SPACEDIM,0);
  inflow_bc_hi.resize(BL_SPACEDIM,0);

  int n_bcs = cp.countval("bc_labels");
  if (n_bcs > 0)
  {
    bc_array.resize(n_bcs,PArrayManage);
    Array<std::string> bc_names;
    cp.getarr("bc_labels",bc_names,0,n_bcs);

    if (region_manager == 0) {
      BoxLib::Abort("static Region manager must be set up prior to reading tracer BCs");
    }
    for (int i = 0; i<n_bcs; i++)
    {
      int ibc = i;
      const std::string& bcname = bc_names[i];
      const std::string prefix("comp.bcs." + bcname);
      ParmParse ppr(prefix.c_str());

      int n_bc_regions = ppr.countval("regions");
      Array<std::string> region_names;
      ppr.getarr("regions",region_names,0,n_bc_regions);
      Array<const Region*> bc_regions = region_manager->RegionPtrArray(region_names);
      std::string bc_type; ppr.get("type",bc_type);

      bool is_inflow = false;
      int component_bc = 1;
      int pressure_bc  = 1;

      use_gauge_pressure[bcname] = false; // Default value

      /*
	Supported types
	string bc_type_labels[10] = {"inward_mass_flux",
	"outward_mass_flux",
	"inward_volumetric_flux",
	"outward_volumetric_flux",
	"uniform_pressure",
	"linear_pressure",
	"seepage_face",
	"hydrostatic",
	"linear_hydrostatic",
	"no_flow"};
      */
      if (bc_type == "uniform_pressure")
      {
	int nPhase = pNames.size();
	BL_ASSERT(nPhase==1); // FIXME
	Array<Real> vals, times;
	Array<std::string> forms;

	std::string val_name = "vals";
	int nv = ppr.countval(val_name.c_str());
	if (nv) {
	  ppr.getarr(val_name.c_str(),vals,0,nv);
	  times.resize(nv,0);
	  if (nv>1) {
	    ppr.getarr("times",times,0,nv);
	    ppr.getarr("forms",forms,0,nv-1);
	  }
	}

	// convert to atm
	for (int j=0; j<vals.size(); ++j) {
	  vals[j] = vals[j] / BL_ONEATM;
	}

	is_inflow = false;
	if (flow_model == PM_FLOW_MODEL_SATURATED) {
	  component_bc = Outflow;
	} else {
	  component_bc = Inflow;
	}
	pressure_bc = Outflow;
	bc_array.set(ibc, new RegionData(bcname,bc_regions,bc_type,vals));
	// Override flow_eval mode from IC if pbc given
	if (flow_is_static && (flow_eval!=PM_FLOW_EVAL_CONSTANT))
	  flow_eval  = PM_FLOW_EVAL_SOLVE_GIVEN_PBC;
      }
      else if (bc_type == "hydraulic_head"
	       || bc_type == "hydrostatic")
      {
	Array<Real> vals, times;
	Array<std::string> forms;
	std::string val_name = "vals";
	int nv = ppr.countval(val_name.c_str());
	BL_ASSERT(nv>0);
	ppr.getarr(val_name.c_str(),vals,0,nv);

	times.resize(nv,0);
	if (nv>1) {
	  ppr.getarr("times",times,0,nv);
	  ppr.getarr("forms",forms,0,nv-1);
	}

	if (pp.countval("normalization")>0) {
	  std::string norm_str; pp.get("normalization",norm_str);
	  if (norm_str == "Absolute") {
	    use_gauge_pressure[bcname] = false;
	  } else if (norm_str == "Relative") {
	    use_gauge_pressure[bcname] = true;
	  } else {
	    BoxLib::Abort("hydraulic_head BC normalization must be \"Absolute\" or \"Relative\"");
	  }
	}

	is_inflow = false;
	if (flow_model == PM_FLOW_MODEL_SATURATED)
	{
	  component_bc = Outflow;
	} else {
	  component_bc = Inflow;
	}
	pressure_bc = Outflow;

	bc_array.set(ibc,new ArrayRegionData(bcname,times,vals,forms,bc_regions,bc_type,vals.size()));
	// Override flow_eval mode from IC if pbc given
	if (flow_is_static && (flow_eval!=PM_FLOW_EVAL_CONSTANT))
	  flow_eval  = PM_FLOW_EVAL_SOLVE_GIVEN_PBC;
      }
      else if (bc_type == "linear_pressure")
      {
	Real val; ppr.get("val",val);
	int ng = ppr.countval("grad");
	BL_ASSERT(ng>=BL_SPACEDIM);
	Array<Real> grad(BL_SPACEDIM); ppr.getarr("grad",grad,0,BL_SPACEDIM);

	int nl = ppr.countval("loc");
	BL_ASSERT(nl>=BL_SPACEDIM);
	Array<Real> loc(BL_SPACEDIM); ppr.getarr("loc",loc,0,BL_SPACEDIM);

	Array<Real> vals(2*BL_SPACEDIM+1);
	vals[0] = val / BL_ONEATM;
	for (int d=0; d<BL_SPACEDIM; ++d) {
	  vals[1+d] = grad[d] / BL_ONEATM;
	  vals[1+d+BL_SPACEDIM] = loc[d];
	}

	is_inflow = false;
	if (flow_model == PM_FLOW_MODEL_SATURATED) {
	  component_bc = Outflow;
	} else {
	  component_bc = Inflow;
	}
	pressure_bc = Outflow;

	Array<Array<Real> > values(vals.size(),Array<Real>(1,0));
	for (int j=0; j<vals.size(); ++j) {
	  values[j][0] = vals[j];
	}
	Array<Array<Real> > times(vals.size(),Array<Real>(1,0));
	Array<Array<std::string> > forms(vals.size(),Array<std::string>(0));
	bc_array.set(ibc,new ArrayRegionData(bcname,times,values,forms,bc_regions,bc_type));
	// Override flow_eval mode from IC if pbc given
	if (flow_is_static && (flow_eval!=PM_FLOW_EVAL_CONSTANT))
	  flow_eval  = PM_FLOW_EVAL_SOLVE_GIVEN_PBC;
      }
      else if (bc_type == "inward_volumetric_flux"
	       || bc_type == "outward_volumetric_flux"
	       || bc_type == "inward_mass_flux"
	       || bc_type == "outward_mass_flux" )
      {
	Array<Real> vals, times;
	Array<std::string> forms;

	int nv = ppr.countval("vals");
	if (nv) {
	  ppr.getarr("vals",vals,0,nv);
	  times.resize(nv,0);
	  if (nv>1) {
	    ppr.getarr("times",times,0,nv);
	    ppr.getarr("forms",forms,0,nv-1);
	  }
	}
	else {
	  vals.resize(1,0);
	  times.resize(1,0);
	  forms.resize(0);
	}

	// If mass flux, convert to volumetric flux.  Assume that we have already decided what the density is
	// Note: we should delay this conversion until bc applied if density is not constant
	if (bc_type == "inward_mass_flux"
	    || bc_type == "outward_mass_flux" ) {
	  for (int k=0; k<vals.size(); ++k) {
	    vals[k] *= 1/density[0]; // Note: this ASSUMES fluxes are of comp[0] mass
	  }
	}

	// Work out sign of flux for this boundary
	int is_hi = -1;
	for (int j=0; j<bc_regions.size(); ++j)
	{
	  const std::string purpose = bc_regions[j]->purpose;
	  for (int k=0; k<7; ++k) {
	    if (purpose == PMAMR::RpurposeDEF[k]) {
	      if (k == 6) {
		BoxLib::Abort(std::string("BC \""+bcname+"\" must be applied on a face region").c_str());
	      }
	      bool this_is_hi = (k>3);
	      if (is_hi < 0) {
		is_hi = this_is_hi;
	      }
	      else {
		if (this_is_hi != is_hi) {
		  BoxLib::Abort("BC must apply to a single face only");
		}
	      }
	    }
	  }
	}
	if (is_hi) {
	  for (int k=0; k<vals.size(); ++k) {
	    vals[k] = -vals[k];
	  }
	}

	// Now that we have converted them all, set to generic flux name
	bc_type = "volumetric_flux";

	is_inflow = true;
	component_bc = Inflow;
	pressure_bc = Inflow;
	bc_array.set(ibc,new ArrayRegionData(bcname,times,vals,forms,bc_regions,bc_type,1));
	// Override flow_eval mode from IC if pbc given
	if (flow_is_static && (flow_eval!=PM_FLOW_EVAL_CONSTANT))
	  flow_eval  = PM_FLOW_EVAL_SOLVE_GIVEN_PBC;
      }
      else if (bc_type == "no_flow")
      {
	Array<Real> vals(1,0), times(1,0);
	Array<std::string> forms(0);
	is_inflow = true;
	component_bc = Inflow;
	pressure_bc = Inflow;
	bc_array.set(ibc,new ArrayRegionData(bcname,times,vals,forms,bc_regions,bc_type,1));
	// Override flow_eval mode from IC if pbc given
	if (flow_is_static && (flow_eval!=PM_FLOW_EVAL_CONSTANT))
	  flow_eval  = PM_FLOW_EVAL_SOLVE_GIVEN_PBC;
      }
      else
      {
	if (bc_type == "seepage_face"
	    || bc_type == "linear_hydrostatic") {
	  BoxLib::Abort(std::string(bc_type+" not yet implemented").c_str());
	}
	BoxLib::Abort(std::string(bc_type+" not a valid bc type").c_str());
      }

      // Some clean up
      std::set<std::string> o_set;

      for (int j=0; j<bc_regions.size(); ++j)
      {
	const std::string purpose = bc_regions[j]->purpose;
	int dir = -1, is_hi;
	for (int k=0; k<7; ++k) {
	  if (purpose == PMAMR::RpurposeDEF[k]) {
	    BL_ASSERT(k != 6);
	    dir = k%3;
	    is_hi = k>=3;
	  }
	}
	if (dir<0 || dir > BL_SPACEDIM) {
	  std::cout << "Bad region for boundary: \n" << bc_regions[j] << std::endl;
	  BoxLib::Abort();
	}

	if (o_set.find(purpose) == o_set.end())
	{
	  o_set.insert(purpose);

	  if (is_hi) {
	    rinflow_bc_hi[dir] = (is_inflow ? 1 : 0);
	    phys_bc.setHi(dir,component_bc);
	    pres_bc.setHi(dir,pressure_bc);
	  }
	  else {
	    rinflow_bc_lo[dir] = (is_inflow ? 1 : 0);
	    phys_bc.setLo(dir,component_bc);
	    pres_bc.setLo(dir,pressure_bc);
	  }
	}
	else {

	  bool is_consistent = true;
	  if (is_hi) {

	    is_consistent = ( (rinflow_bc_hi[dir] == is_inflow)
			      && (phys_bc.hi()[dir] == component_bc)
			      && (pres_bc.hi()[dir] == pressure_bc) );
	  }
	  else {

	    is_consistent = ( (rinflow_bc_lo[dir] == is_inflow)
			      && (phys_bc.lo()[dir] == component_bc)
			      && (pres_bc.lo()[dir] == pressure_bc) );
	  }

	  if (!is_consistent) {
	    BoxLib::Abort("Inconconsistent type for boundary ");
	  }
	}
      }
    }
  }
}

using PMAMR::RlabelDEF;
void  PorousMedia::read_tracer()
{
  //
  // Read in parameters for tracers
  //
  ParmParse pp("tracer");

  // Get number of tracers
  ntracers = pp.countval("tracers");
  if (ntracers > 0)
  {
    tic_array.resize(ntracers);
    tbc_array.resize(ntracers);
    pp.getarr("tracers",tNames,0,ntracers);
  }

  // Now that we have solute names, set up materials
  if (region_manager == 0) {
    BoxLib::Abort("static Region manager must be set up prior to building Rock Manager");
  }

  BL_ASSERT(density.size()>0);
  BL_ASSERT(muval.size()>0);
  Real gravity_value_mks = gravity * BL_ONEATM;
  Real liquid_density_mks = density[0];
  Real liquid_viscosity_mks = muval[0];
  rock_manager = new RockManager(region_manager,&tNames,liquid_density_mks,liquid_viscosity_mks,gravity_value_mks);

  if (do_tracer_diffusion) {
    tensor_tracer_diffusion = rock_manager->DoTensorDiffusion();
  }

  ParmParse ppp("prob");

  // get Chemistry stuff
  const std::string Chemistry_Model_stru = "chemistry_model";
  ppp.query(Chemistry_Model_stru.c_str(),chemistry_model_name);
  do_tracer_chemistry = chemistry_model_name != "Off";
  ppp.query("do_full_strang",do_full_strang);
  ppp.query("n_chem_interval",n_chem_interval);
  ppp.query("ic_chem_relax_dt",ic_chem_relax_dt);
  if (n_chem_interval > 0) {
    do_full_strang = false;
  }

#if ALQUIMIA_ENABLED
  chemistry_engine = 0;
#endif
  chemistry_helper = 0;

  first_order_decay_constant.resize(ntracers,0);
  bool have_tracer_chemistry = chemistry_model_name == "Amanzi" || chemistry_model_name == "Alquimia";
  if (have_tracer_chemistry) {
    const std::string chemistry_str = "Chemistry";

    ParmParse ppc(chemistry_str.c_str());

    if (chemistry_model_name != "Off") {

      const std::string Chemistry_Max_Time_Step_str = "Max_Time_Step";
      max_chemistry_time_step = -1;
      if (int nmts = ppc.countval(Chemistry_Max_Time_Step_str.c_str())) {
        ppc.get(Chemistry_Max_Time_Step_str.c_str(),max_chemistry_time_step);
      }

      for (int i = 0; i<ntracers; i++) {
        const std::string prefix("tracer." + tNames[i]);
        ParmParse ppr(prefix.c_str());
        if (ppr.countval("firstOrderDecayConstant") > 0) {
          BoxLib::Abort("Radioactive decay constants cannot yet be specified in Amanzi input");
        }
        ppr.query("firstOrderDecayConstant",first_order_decay_constant[i]);
      }

      if (chemistry_model_name == "Amanzi") {

        Teuchos::ParameterList plist;
        ParmParse pb("prob.amanzi");
        std::string verbose_chemistry_init = "silent"; ppc.query("verbose_chemistry_init",verbose_chemistry_init);

        const std::string thermo_str = "Thermodynamic_Database";
        const std::string thermo_fmt_str = thermo_str + "_Format";
        const std::string thermo_file_str = thermo_str + "_File";

        std::string amanzi_thermo_fmt, amanzi_thermo_file;
        ppc.get(thermo_fmt_str.c_str(),amanzi_thermo_fmt);
        ppc.get(thermo_file_str.c_str(),amanzi_thermo_file);

        const std::string& activity_model_u  = Amanzi::AmanziChemistry::ActivityModelFactory::unit;
        std::string activity_model = activity_model_u;
        ppc.query("Activity_Model",activity_model);

        Real tolerance=1.5e-12; ppc.query("Tolerance",tolerance);
        int max_num_Newton_iters = 150; ppc.query("Maximum_Newton_Iterations",max_num_Newton_iters);
        std::string outfile=""; ppc.query("Output_File_Name",outfile);
        bool use_stdout = true; ppc.query("Use_Standard_Out",use_stdout);
        //int num_aux = ppc.countval("Auxiliary_Data");
        // if (num_aux>0) {
        //   Array<std::string> tmpaux(num_aux);
        //   aux_chem_variables.clear();
        //   ppc.getarr("Auxiliary_Data",tmpaux,0,num_aux);
        //   for (int i=0;i<num_aux;i++)
        //     aux_chem_variables[tmpaux[i]] = i;
        // }

        nminerals = rock_manager->NumMinerals();
        minerals.resize(nminerals);
        for (int i=0; i<nminerals; ++i) {
          minerals[i] = rock_manager->MineralNames()[i];
        }

        nsorption_sites = rock_manager->NumSorptionSites();
        sorption_sites.resize(nsorption_sites);
        for (int i=0; i<nsorption_sites; ++i) {
          sorption_sites[i] = rock_manager->SorptionSiteNames()[i];
        }
        using_sorption = rock_manager->UsingSorption();
        ncation_exchange = rock_manager->NumCationExchange();
        bool hasCationExchangeCapacity = ncation_exchange > 0;

        Array<std::string> sorbedPrimarySpecies;
        if (using_sorption) {
          sorbedPrimarySpecies.resize(ntracers);
          for (int i=0; i<ntracers; ++i) {
            sorbedPrimarySpecies[i] = tNames[i];
          }
        }

        int nisotherms = rock_manager->NumSorptionIsotherms();
        Array<std::string> isothermNames;
        if (nisotherms > 0) {
          if (nisotherms != ntracers) {
            BoxLib::Abort("Disallowed number of isotherms");
          }
          isothermNames.resize(ntracers);
          for (int i=0; i<ntracers; ++i) {
            isothermNames[i] = tNames[i];
          }
        }

        chemistry_helper = new AmanziChemHelper_Structured(tNames,sorbedPrimarySpecies,minerals,sorption_sites,hasCationExchangeCapacity,
                                                           isothermNames,tNames,amanzi_thermo_file,amanzi_thermo_fmt,activity_model);
#if ALQUIMIA_ENABLED
      } else {
        BL_ASSERT(chemistry_model_name == "Alquimia");
        const std::string Chemistry_Engine_stru = "Engine";
        const std::string Chemistry_Engine_Input_stru = "Engine_Input_File";
        ppc.query(Chemistry_Engine_stru.c_str(),chemistry_engine_name);
        std::string chem_engine_input_filename; ppc.get(Chemistry_Engine_Input_stru.c_str(),chem_engine_input_filename);
        chemistry_engine = new Amanzi::AmanziChemistry::ChemistryEngine(chemistry_engine_name,chem_engine_input_filename);

	//SM: alquimia helper chemistry_engine gets Alquimia engine
	//    but overwriting amanzi's stuff with alquimia's (below)
	//    does have any effect on output unless in PorousMedia::variableSetUp
	//    PM_setup static variables set here are used there
        chemistry_helper = new AlquimiaHelper_Structured(chemistry_engine);

        //
        // FIXME: THIS WILL OVERWRITE THE LIST OF AMANZI TRACERS
	// SM: FIX not needed for now: overwriting is good - see above and PorousMedia::variableSetUp
        //
        // convert arrays to those of PM internals
        std::vector<std::string> primarySpeciesNames;
        chemistry_engine->GetPrimarySpeciesNames(primarySpeciesNames);
        ntracers = primarySpeciesNames.size();
        tNames.resize(ntracers);
        for (int i=0; i<ntracers; ++i) {
          tNames[i] = primarySpeciesNames[i];
        }

        std::vector<std::string> mineralNames;
        chemistry_engine->GetMineralNames(mineralNames);
        nminerals = mineralNames.size();
        minerals.resize(nminerals);
        for (int i=0; i<nminerals; ++i) {
          minerals[i] = mineralNames[i];
        }

	//SM: add isotherms, surface complexation, ion exchange below

	//ion exchange
        ncation_exchange = chemistry_engine->NumIonExchangeSites();

	//surface complexation
	std::vector<std::string> surfaceNames;
        chemistry_engine->GetSurfaceSiteNames(surfaceNames);
        nsorption_sites = surfaceNames.size();
        sorption_sites.resize(nsorption_sites);
        for (int i=0; i<nsorption_sites; ++i) {
          sorption_sites[i] = surfaceNames[i];
        }

	//isotherms
	nsorption_isotherms = chemistry_engine->NumIsothermSpecies();

	//using_sorption
	if (ncation_exchange ||
	    nsorption_sites ||
	    nsorption_isotherms) {
	  using_sorption=true;
	}
#endif
      }
    }
  }

  ppp.query("use_funccount",use_funccount);
  ppp.query("max_grid_size_chem",max_grid_size_chem);
  BL_ASSERT(max_grid_size_chem > 0);

  molecular_diffusivity.resize(ntracers,0);
  if (ntracers > 0)
  {
    int Nmobile = ntracers;
#if ALQUIMIA_ENABLED
    if (chemistry_engine != 0) {
      // FIXME: Amanzi input MUST be consistent with chemistry class
      std::vector<std::string> primarySpeciesNames; chemistry_engine->GetPrimarySpeciesNames(primarySpeciesNames);
      BL_ASSERT(ntracers==chemistry_engine->NumPrimarySpecies());
      Nmobile = chemistry_engine->NumPrimarySpecies();
      for (int i = 0; i<Nmobile; i++) {
        BL_ASSERT(primarySpeciesNames[i] == tNames[i]);
      }
    }
#endif
    if (do_tracer_chemistry>0  ||  do_tracer_advection  ||  do_tracer_diffusion) {
      setup_tracer_transport = true;
      for (int i = 0; i<Nmobile; i++) {
        group_map["Total"].push_back(i+ncomps);
      }
    }
    else {
      setup_tracer_transport = false;
    }

    for (int i = 0; i<ntracers; i++)
    {
      const std::string prefix("tracer." + tNames[i]);
      ParmParse ppr(prefix.c_str());

      // Initial condition and boundary condition
      Array<std::string> tic_names;
      int n_ic = ppr.countval("tinits");
      if (n_ic <= 0)
      {
        BoxLib::Abort("each tracer must be initialized");
      }
      ppr.getarr("tinits",tic_names,0,n_ic);
      tic_array[i].resize(n_ic,PArrayManage);

      for (int n = 0; n<n_ic; n++)
      {
        const std::string prefixIC(prefix + ".ic." + tic_names[n]);
        ParmParse ppri(prefixIC.c_str());
        int n_ic_region = ppri.countval("regions");
        Array<std::string> region_names;
        ppri.getarr("regions",region_names,0,n_ic_region);
        Array<const Region*> tic_regions = region_manager->RegionPtrArray(region_names);
        std::string tic_type; ppri.get("type",tic_type);

        if (tic_type == "concentration")
        {
          if (ppri.countval("geochemical_condition")) {

            if ( !(chemistry_model_name == "Alquimia") ) {
              BoxLib::Abort("Cannot use geochemical conditions if chemistry model not Alquimia");
            }
            std::string geocond; ppri.get("geochemical_condition",geocond);
            tic_array[i].set(n, new ChemConstraint(tNames[i],tic_regions,tic_type,
                                                   ChemConstraintEval(geocond,i,rock_manager,chemistry_helper)));
          }
          else {
            int nv = ppri.countval("val");
            if (nv != 1) {
              std::string m = "Solute IC: \"" + tic_names[n]
                + "\" for \"" + tNames[i] + "\" requires a single value";
              BoxLib::Abort(m.c_str());
            }
            Real val = 0; ppri.query("val",val);
            tic_array[i].set(n, new IdxRegionData(tNames[i],tic_regions,tic_type,val));
          }

          // Check for "Free_Ion_Guess", load structure used to set aux_chem components
          const std::string FIG_str = "Free_Ion_Guess";
          int nfig = ppri.countval(FIG_str.c_str());
          if (nfig > 0) {
            Real val(0.0);
            if (nfig == 1) {
              ppri.get(FIG_str.c_str(), val);
            }
            else if (nfig > 1) {
              std::string m = "Solute IC: \"" + tic_names[n]
                + "\": Free Ion Guess parameter for \"" + tNames[i] + "\" requires a single value";
              BoxLib::Abort(m.c_str());
            }
            solute_chem_ics[tic_names[n]][tNames[i]][FIG_str] = val; // sc[rockname][solute][property] = val
          }

          // Check for "Activity_Coefficient", load structure used to set aux_chem components
          const std::string AC_str = "Activity_Coefficient";
          int nac = ppri.countval(AC_str.c_str());
          if (nac > 0) {
            Real valac(0.0);
            if (nac == 1) {
              ppri.get(AC_str.c_str(), valac);
            }
            else if (nac > 1) {
              std::string m = "Solute IC: \"" + tic_names[n]
                + "\": Activity Coefficient parameter for \"" + tNames[i] + "\" requires a single value";
              BoxLib::Abort(m.c_str());
            }
            solute_chem_ics[tic_names[n]][tNames[i]][AC_str] = valac; // sc[rockname][solute][property] = val
          }
        }
        else {
          std::string m = "Solute IC: \"" + tic_names[n]
            + "\": Unsupported Solute IC type: \"" + tic_type + "\"";
          BoxLib::Abort(m.c_str());
        }
      }

      if (setup_tracer_transport)
      {
        Array<std::string> tbc_names;
        int n_tbc = ppr.countval("tbcs");
        ppr.getarr("tbcs",tbc_names,0,n_tbc);
        //tbc_array[i].resize(n_tbc+2*BL_SPACEDIM,PArrayManage);
        tbc_array[i].resize(n_tbc,PArrayManage);

        // Explicitly build default BCs
        int tbc_cnt = 0;

        // FIXME:
        // When these are used, we pick up a cross derivative term that can be seen when a front that is
        // perpendicular to the boundary moves tangential to that boundary, even when the normal velocity
        // across that boundary is identically zero.  This is an error and should be fixed since information
        // should not propagate through a zero velocity wall.  For the time being, we have set it up so that
        // the default BC is instead FOEXTRAP, minimizing this effect.  However we should go more carefully
        // through the advection code to find why the cross terms are not correctly dealt with.
        //
        // for (int n=0; n<BL_SPACEDIM; ++n) {
        //   tbc_array[i].set(tbc_cnt++,
        //                    new RegionData(RlabelDEF[n] + "_DEFAULT",
        //                                   region_manager->RegionPtrArray(Array<std::string>(1,RlabelDEF[n])),
        //                                   std::string("noflow"),0));
        //   tbc_array[i].set(tbc_cnt++,
        //                    new RegionData(RlabelDEF[n+3] + "_DEFAULT",
        //                                   region_manager->RegionPtrArray(Array<std::string>(1,RlabelDEF[n+3])),
        //                                   std::string("noflow"),0));
        // }

        Array<int> orient_types(6,-1);
        for (int n = 0; n<n_tbc; n++)
        {
          const std::string prefixTBC(prefix + ".bc." + tbc_names[n]);
          ParmParse ppri(prefixTBC.c_str());

          int n_tbc_region = ppri.countval("regions");
          Array<std::string> tbc_region_names;
          ppri.getarr("regions",tbc_region_names,0,n_tbc_region);

          Array<const Region*> tbc_regions = region_manager->RegionPtrArray(tbc_region_names);
          std::string tbc_type; ppri.get("type",tbc_type);

          // When we get the BCs, we need to translate to AMR-standardized type id.  By
          // convention, components are  Interior, Inflow, Outflow, Symmetry, SlipWall, NoSlipWall.
          int AMR_BC_tID = -1;

          if (tbc_type == "concentration")
          {
            Array<Real> times, vals;
            Array<std::string> forms;

            if (ppri.countval("geochemical_conditions")) {
              if ( !(chemistry_model_name == "Alquimia") ) {
                BoxLib::Abort("Cannot use geochemical conditions if chemistry model not Alquimia");
              }
              int nv = ppri.countval("geochemical_conditions");
              Array<std::string> geoconds; ppri.getarr("geochemical_conditions",geoconds,0,nv);
              if (nv > 1) {
                ppri.getarr("times",times,0,nv);
              }
              else {
                times.resize(1,0);
              }
              tbc_array[i].set(tbc_cnt++, new ChemConstraint(tbc_names[n],tbc_regions,tbc_type,
							     ChemConstraintEval(geoconds,times,i,rock_manager,chemistry_helper)));
            }
            else {
              int nv = ppri.countval("vals");
              if (nv) {
                ppri.getarr("vals",vals,0,nv);
                if (nv>1) {
                  ppri.getarr("times",times,0,nv);
                  ppri.getarr("forms",forms,0,nv-1);
                }
                else {
                  times.resize(1,0);
                }
              }
              else {
                vals.resize(1,0); // Default tracers to zero for all time
                times.resize(1,0);
                forms.resize(0);
              }
	      tbc_array[i].set(tbc_cnt++, new IdxRegionData(tbc_names[n],tbc_regions,tbc_type,vals,times,forms));
            }
            AMR_BC_tID = 1; // Inflow
          }
          else if (tbc_type == "noflow")
          {
            Real val = 0;
            tbc_array[i].set(tbc_cnt++, new IdxRegionData(tbc_names[n],tbc_regions,tbc_type,val));
            AMR_BC_tID = 2;
          }
          else if (tbc_type == "outflow")
          {
            Real val=0;
            tbc_array[i].set(tbc_cnt++, new IdxRegionData(tbc_names[n],tbc_regions,tbc_type,val));
            AMR_BC_tID = 3; // Outflow
          }
          else {
            std::string m = "Tracer BC: \"" + tbc_names[n]
              + "\": Unsupported tracer BC type: \"" + tbc_type + "\"";
            BoxLib::Abort(m.c_str());
          }


          for (int j=0; j<tbc_regions.size(); ++j)
          {
            const std::string purpose = tbc_regions[j]->purpose;
            int dir = -1, is_hi, k;
            for (int kt=0; kt<7 && dir<0; ++kt) {
              if (purpose == PMAMR::RpurposeDEF[kt]) {
                BL_ASSERT(kt != 6);
                dir = kt%3;
                is_hi = kt>=3;
                k = kt;
              }
            }
            if (dir<0 || dir > BL_SPACEDIM) {
              std::cout << "Bad region for boundary: \n" << tbc_regions[j] << std::endl;
              BoxLib::Abort();
            }

            if (orient_types[k] < 0) {
              orient_types[k] = AMR_BC_tID;
            } else {
              if (orient_types[k] != AMR_BC_tID) {
                BoxLib::Abort("BC for tracers must all be of same type on each side");
              }
            }
          }
        }
        // Set the default BC type
        for (int k=0; k<orient_types.size(); ++k) {
          if (orient_types[k] < 0) orient_types[k] = 2;
        }

        BCRec phys_bc_trac;
        for (int i = 0; i < BL_SPACEDIM; i++) {
          phys_bc_trac.setLo(i,orient_types[i]);
          phys_bc_trac.setHi(i,orient_types[i+3]);
        }
        set_tracer_bc(trac_bc,phys_bc_trac);
      }

      ppr.query("molecularDiffusivity",molecular_diffusivity[i]);
    }
    ndiff += ntracers;
  }
}

static
int loc_in_array(const std::string& val,const Array<std::string>& arr)
{
  int location = -1;
  for (int i=0; i<arr.size() && location<0; ++i) {
    if (val == arr[i]) location = i;
  }
  return location;
}


void  PorousMedia::read_source()
{
  //
  // Read in parameters for sources
  //
  ParmParse pp("source");
  ParmParse ppb("prob");
  ppb.query("do_source_term",do_source_term);

  int nsources = pp.countval("sources");
  if (nsources>0) {
    if (region_manager == 0) {
      BoxLib::Abort("static Region manager must be set up prior to reading sources");
    }
    source_array.resize(nsources,PArrayManage);
    tsource_array.resize(nsources);
    Array<std::string> source_names(nsources);
    pp.getarr("sources",source_names,0,nsources);
    for (int i=0; i<nsources; ++i) {
      const std::string& source_name = source_names[i];
      const std::string prefix("source." + source_name);
      ParmParse pps(prefix.c_str());

      int n_src_regions = pps.countval("regions");
      Array<std::string> src_region_names;
      pps.getarr("regions",src_region_names,0,n_src_regions);
      const Array<const Region*> source_regions = region_manager->RegionPtrArray(src_region_names);

      if (pps.countval("type")) {
	std::string source_type; pps.get("type",source_type);
	if (source_type == "uniform"
	    || source_type == "volume_weighted"
	    || source_type == "permeability_weighted"
	    || source_type == "point")
	  {
	    int nvars = pps.countval("vals");
	    BL_ASSERT(nvars>0);
	    Array<Real> vals; pps.getarr("vals",vals,0,nvars);

            if (source_type == "point") {
              BL_ASSERT(source_regions.size() == 1);
              BL_ASSERT(source_regions[0]->type=="point");
            }

	    source_array.set(i, new RegionData(source_name,source_regions,source_type,vals));
	  }
	else {
	  std::string m = "Source: \"" + source_names[i]
	    + "\": Unsupported source type: \"" + source_type + "\"";
	  BoxLib::Abort(m.c_str());
	}
      }
      else {
	std::string m = "Source: \"" + source_names[i]
	  + "\": Requires \"type\" specifier";
	BoxLib::Abort(m.c_str());
      }

      for (int ip=0; ip<pNames.size(); ++ip) {
	const std::string& pName = pNames[ip];
	const std::string p_prefix(prefix+"."+pName);
	ParmParse pps_p(p_prefix.c_str());

	for (int ic=0; ic<cNames.size(); ++ic) {
	  const std::string& cName = cNames[ic];
	  const std::string c_prefix(p_prefix+"."+cName);
	  ParmParse pps_c(c_prefix.c_str());

	  int ntracers_with_sources = pps_c.countval("tracers_with_sources");
	  if (ntracers_with_sources>0) {
	    Array<std::string> tracers_with_sources;
	    pps_c.getarr("tracers_with_sources",tracers_with_sources,0,ntracers_with_sources);
	    tsource_array[i].resize(ntracers, PArrayManage);

	    for (int it=0; it<tracers_with_sources.size(); ++it) {
	      const std::string& tName = tracers_with_sources[it];
	      int t_pos = loc_in_array(tName,tNames);
	      if (t_pos>=0) {
		const std::string c_t_prefix(c_prefix+"."+tName);
		ParmParse pps_c_t(c_t_prefix.c_str());

		if (pps_c_t.countval("type")) {
		  std::string tsource_type; pps_c_t.get("type",tsource_type);
		  if (tsource_type == "uniform"
		      || tsource_type == "flow_weighted"
		      || tsource_type == "point")
		    {
		      int ntvars = pps_c_t.countval("vals");
		      BL_ASSERT(ntvars>0);
		      Array<Real> tvals; pps_c_t.getarr("vals",tvals,0,ntvars);
		      tsource_array[i].set(t_pos, new RegionData(source_name,source_regions,tsource_type,tvals));
		    }
                  else if (tsource_type == "diffusion_dominated_release_model") {
                    Real total_inventory; pps_c_t.get("total_inventory",total_inventory);
                    Real mixing_length; pps_c_t.get("mixing_length",mixing_length);
                    Real D_eff; pps_c_t.get("effective_diffusion_coef",D_eff);
                    Real start_time; pps_c_t.get("start_time",start_time);
                    Real end_time; pps_c_t.get("end_time",end_time);
                    Real time_scale; pps_c_t.get("time_scale",time_scale);
                    tsource_array[i].set(t_pos, new DiffDomRelSrc(source_name,source_regions,tsource_type,mixing_length,D_eff,total_inventory,start_time,end_time,time_scale));
                  }
		  else {
		    BoxLib::Abort(std::string("Source: \"" + source_names[i] +
					      "\", Comp: \"" + cName + "\", Solute SOURCE: \"" + tName
					      + "\": Unsupported source type: \"" + tsource_type + "\"").c_str());
		  }
		} else {
		  BoxLib::Abort(std::string("Source: \"" + source_names[i]
					    + "\": Requires \"type\" specifier for solute \""+tName+"\"").c_str());
		}
		if (pps_c_t.countval("Concentration_Units")) {
		  // FIXME: We do not currently do anything with this parameter
		}
	      }
	      else {
		BoxLib::Abort(std::string("Source: \"" + source_names[i]
					  + "\" contains unknown tracer: \""+tName+"\"").c_str());
	      }
	    }

	    // Set default source (uniform=0) for all tracers not set explicitly
	    const std::string default_tsource_type = "uniform";
	    const Array<Real> default_tsource_tvals(1,0);
	    for (int it=0; it<ntracers; ++it) {
	      if ( !(tsource_array[i].defined(it)) ) {
		tsource_array[i].set(it, new RegionData(source_name,source_regions,default_tsource_type,default_tsource_tvals));
	      }
	    }
	  }
	}
      }
    }
  }
}

void PorousMedia::read_params()
{
  // problem-specific
  read_prob();

  // Require regions prior to setting up phases/comps
  if (verbose > 1 && ParallelDescriptor::IOProcessor())
    std::cout << "Reading geometry" << std::endl;
  region_manager = new RegionManager();
  if (echo_inputs && ParallelDescriptor::IOProcessor()) {
      std::cout << "The Regions: " << std::endl;
      const Array<const Region*> regions = region_manager->RegionPtrArray();
      for (int i=0; i<regions.size(); ++i) {
	std::cout << *(regions[i]) << std::endl;
      }
  }

  // components and phases
  if (verbose > 1 && ParallelDescriptor::IOProcessor())
    std::cout << "Reading phases/components."<< std::endl;
  read_comp();

  // tracers and chemistry
  if (verbose > 1 && ParallelDescriptor::IOProcessor())
    std::cout << "Read tracers/chemistry"<< std::endl;
  read_tracer();

  // source
  if (verbose > 1 && ParallelDescriptor::IOProcessor())
    std::cout << "Reading sources."<< std::endl;
  read_source();

  int model_int = FlowModel();
  FORT_INITPARAMS(&ncomps,&nphases,&model_int,density.dataPtr(),
		  muval.dataPtr(),pType.dataPtr(),
		  &gravity,&gravity_dir);

  if (ntracers > 0)
    FORT_TCRPARAMS(&ntracers);
}
