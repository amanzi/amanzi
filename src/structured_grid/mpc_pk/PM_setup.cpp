#include <winstd.H>
#include <ParmParse.H>
#include <Interpolater.H>
#include <MultiGrid.H>
#include <ArrayLim.H>
#include <Profiler.H>
#include <TagBox.H>
#include <DataServices.H>
#include <AmrData.H>
#include <time.h> 

#include <PorousMedia.H>
#include <PMAMR_Labels.H>
#include <RegType.H> 
#include <PROB_PM_F.H>
#include <DERIVE_F.H>
#include <PMAMR_Labels.H>
#include <PMAmr.H> 

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


#ifdef AMANZI
#include "simple_thermo_database.hh"
#include "activity_model_factory.hh"
#endif

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

//
// The num_state_type actually varies with model.
//
// Add 2 if do_chem>0 later.
//
int PorousMedia::num_state_type;
//
// Region.
//
std::string    PorousMedia::surf_file;
PArray<Region> PorousMedia::regions;
//
// Rock
//
std::string PorousMedia::gsfile;
MultiFab*   PorousMedia::kappadata;
MultiFab*   PorousMedia::phidata;
bool        PorousMedia::porosity_from_fine;
bool        PorousMedia::permeability_from_fine;
PArray<Rock> PorousMedia::rocks;
bool        PorousMedia::material_is_layered;
Real        PorousMedia::saturation_threshold_for_vg_Kr;
int         PorousMedia::use_shifted_Kr_eval;
//
// Source.
//
bool          PorousMedia::do_source_term;
Array<Source> PorousMedia::source_array;
//
// Observation.
//
std::string        PorousMedia::obs_outputfile;
PArray<Observation> PorousMedia::observations;
Array<std::string>  PorousMedia::vis_cycle_macros;
Array<std::string>  PorousMedia::vis_time_macros;
Array<std::string>  PorousMedia::chk_cycle_macros;
Array<std::string>  PorousMedia::chk_time_macros;
int                 PorousMedia::echo_inputs;
//
// Phases and components.
//
Array<std::string>  PorousMedia::pNames;
Array<std::string>  PorousMedia::cNames;
Array<int >         PorousMedia::pType;
Array<Real>         PorousMedia::density;
PArray<RegionData>  PorousMedia::ic_array;
PArray<RegionData>  PorousMedia::bc_array;
Array<Real>         PorousMedia::muval;
std::string         PorousMedia::model_name;
int                 PorousMedia::model;
int                 PorousMedia::nphases;
int                 PorousMedia::ncomps;
int                 PorousMedia::ndiff;
int                 PorousMedia::idx_dominant;
//
// Tracers.
//
Array<std::string>  PorousMedia::qNames;
Array<std::string>  PorousMedia::tNames;
int                 PorousMedia::ntracers;
Array<int>          PorousMedia::tType; 
Array<Real>         PorousMedia::tDen;
Array<PArray<RegionData> > PorousMedia::tic_array;
Array<PArray<RegionData> > PorousMedia::tbc_array;
std::map<std::string,Array<int> > PorousMedia::group_map;

//
// Minerals and Sorption sites
//
int                PorousMedia::nminerals;
Array<std::string> PorousMedia::minerals;
int                PorousMedia::nsorption_sites;
Array<std::string> PorousMedia::sorption_sites;
int                PorousMedia::ncation_exchange;
Array<std::string> PorousMedia::aux_chem_variables;
PorousMedia::ChemICMap PorousMedia::sorption_isotherm_ics;
PorousMedia::ChemICMap PorousMedia::mineralogy_ics;
PorousMedia::ChemICMap PorousMedia::surface_complexation_ics;
PorousMedia::ICLabelParmPair PorousMedia::cation_exchange_ics;
PorousMedia::LabelIdx PorousMedia::mineralogy_label_map;
PorousMedia::LabelIdx PorousMedia::sorption_isotherm_label_map;
PorousMedia::LabelIdx PorousMedia::surface_complexation_label_map;
std::map<std::string,int> PorousMedia::cation_exchange_label_map;

// Pressure.
//
#ifdef MG_USE_FBOXLIB
int         PorousMedia::richard_iter;
#endif
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
// Flow.
//
int  PorousMedia::verbose;
Real PorousMedia::cfl;
Real PorousMedia::init_shrink;
Real PorousMedia::change_max;
Real PorousMedia::fixed_dt;
Real PorousMedia::richard_max_dt;
Real PorousMedia::dt_cutoff;
Real PorousMedia::gravity;
int  PorousMedia::initial_step;
int  PorousMedia::initial_iter;
int  PorousMedia::sum_interval;
int  PorousMedia::NUM_SCALARS;
int  PorousMedia::NUM_STATE;
int  PorousMedia::full_cycle;
int  PorousMedia::max_step;
Real PorousMedia::stop_time;
Real PorousMedia::dt_init;

Array<AdvectionForm> PorousMedia::advectionType;
Array<DiffusionForm> PorousMedia::diffusionType;
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
int  PorousMedia::have_capillary;
//
// Molecular diffusion flag.
//
int  PorousMedia::variable_scal_diff;

Array<int>  PorousMedia::is_diffusive;
Array<Real> PorousMedia::visc_coef;
//
// Transport flame
//
int  PorousMedia::do_tracer_transport;
//
// Chemistry flag.
//
int  PorousMedia::do_chem;
int  PorousMedia::do_full_strang;
int  PorousMedia::n_chem_interval;
int  PorousMedia::it_chem;
Real PorousMedia::dt_chem;
int  PorousMedia::max_grid_size_chem;
bool PorousMedia::no_initial_values;
bool PorousMedia::use_funccount;
bool PorousMedia::do_richard_sat_solve;
//
// Lists.
//
std::map<std::string, int> PorousMedia::model_list;
std::map<std::string, int> PorousMedia::phase_list;
std::map<std::string, int> PorousMedia::comp_list;
std::map<std::string, int> PorousMedia::tracer_list;
Array<std::string> PorousMedia::user_derive_list;
//
// AMANZI flags.
//
#ifdef AMANZI
int         PorousMedia::n_total;
int         PorousMedia::n_minerals;
int         PorousMedia::n_sorbed;
std::string PorousMedia::amanzi_input_file;
std::string PorousMedia::amanzi_activity_model;

PArray<amanzi::chemistry::SimpleThermoDatabase>    PorousMedia::chemSolve(PArrayManage);
Array<amanzi::chemistry::Beaker::BeakerComponents> PorousMedia::components;
Array<amanzi::chemistry::Beaker::BeakerParameters> PorousMedia::parameters;
#endif
//
// Internal switches.
//
int  PorousMedia::do_simple;
int  PorousMedia::do_multilevel_full;
int  PorousMedia::do_reflux;
int  PorousMedia::do_correct;
int  PorousMedia::no_corrector;
int  PorousMedia::do_kappa_refine;
int  PorousMedia::n_pressure_interval;
int  PorousMedia::it_pressure;
bool PorousMedia::do_any_diffuse;
int  PorousMedia::do_cpl_advect;

int  PorousMedia::richard_solver_verbose;
//
// Init to steady
//
bool PorousMedia::do_richard_init_to_steady;
int  PorousMedia::richard_init_to_steady_verbose;
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
Real PorousMedia::steady_tolerance;
Real PorousMedia::steady_init_time_step;
Real PorousMedia::steady_max_time_steps;
Real PorousMedia::steady_max_psuedo_time;
int  PorousMedia::steady_max_num_consecutive_success;
Real PorousMedia::steady_extra_time_step_increase_factor;
int  PorousMedia::steady_max_num_consecutive_increases;
Real PorousMedia::steady_consecutive_increase_reduction_factor;
int  PorousMedia::richard_max_ls_iterations;
Real PorousMedia::richard_min_ls_factor;
Real PorousMedia::richard_ls_acceptance_factor;
Real PorousMedia::richard_ls_reduction_factor;
int PorousMedia::richard_monitor_linear_solve;
int PorousMedia::richard_monitor_line_search;

static std::map<std::string,EventCoord::Event*> defined_events; // accumulate all defined, register as needed
namespace
{
    static void PM_Setup_CleanUpStatics() 
    {
        for (std::map<std::string,EventCoord::Event*>::iterator it=defined_events.begin(); it!=defined_events.end(); ++it) {
            delete it->second;
        }
    }
}

static Box grow_box_by_one (const Box& b) { return BoxLib::grow(b,1); }

//
// Components are  Interior, Inflow, Outflow, Symmetry, SlipWall, NoSlipWall.
//

static int scalar_bc[] =
  {
    //    INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_EVEN, FOEXTRAP, SEEPAGE
    INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_EVEN, REFLECT_ODD, SEEPAGE
  };

static int press_bc[] =
  {
    INT_DIR, FOEXTRAP, EXT_DIR, REFLECT_EVEN, FOEXTRAP, FOEXTRAP
  };

static int norm_vel_bc[] =
  {
    INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_ODD, EXT_DIR, EXT_DIR
  };

static int tang_vel_bc[] =
  {
    INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_EVEN, HOEXTRAP, EXT_DIR
  };

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
PorousMedia::setup_list()
{
  // model list
  model_list["single-phase"] = 0;
  model_list["single-phase-solid"] = 1;
  model_list["two-phase"] = 2;
  model_list["polymer"] = 3;
  model_list["richard"] = 4;

#if 0
  // bc_list
  bc_list["file"] = 0;
  bc_list["scalar"] = 1;
  bc_list["hydrostatic"] = 2;
  bc_list["rockhold"] = 3;
  bc_list["zero_total_velocity"] = 4;
#endif
}

void
PorousMedia::InitializeStaticVariables ()
{
  //
  // Set all default values for static variables here!!!
  //
  PorousMedia::num_state_type = 4;

  PorousMedia::kappadata = 0;
  PorousMedia::phidata   = 0;

  PorousMedia::porosity_from_fine     = false;
  PorousMedia::permeability_from_fine = false;

  PorousMedia::do_source_term = false;

  PorousMedia::model        = 0;
  PorousMedia::nphases      = 0;
  PorousMedia::ncomps       = 0; 
  PorousMedia::ndiff        = 0;
  PorousMedia::idx_dominant = -1;

  PorousMedia::ntracers = 0; 
  PorousMedia::nminerals = 0; 
  PorousMedia::minerals.clear();
  PorousMedia::nsorption_sites = 0; 
  PorousMedia::sorption_sites.clear();
  PorousMedia::aux_chem_variables.clear();
  PorousMedia::sorption_isotherm_ics.clear();
  PorousMedia::mineralogy_ics.clear();
  PorousMedia::surface_complexation_ics.clear();
  PorousMedia::cation_exchange_ics.clear();
  PorousMedia::mineralogy_label_map.clear();
  PorousMedia::sorption_isotherm_label_map.clear();
  PorousMedia::surface_complexation_label_map.clear();
  PorousMedia::cation_exchange_label_map.clear();
  
#ifdef MG_USE_FBOXLIB
  PorousMedia::richard_iter = 100;
#endif
  PorousMedia::wt_lo = 0;
  PorousMedia::wt_hi = 0;

  PorousMedia::temperature = 300;

  PorousMedia::verbose      = 0;
  PorousMedia::cfl          = 0.8;
  PorousMedia::init_shrink  = 1.0;
  PorousMedia::dt_init      = -1.0; // Ignore if < 0
  PorousMedia::change_max   = 1.1;
  PorousMedia::fixed_dt     = -1.0;
  PorousMedia::richard_max_dt = 5.e5;
  PorousMedia::dt_cutoff    = 0.0;
  PorousMedia::gravity      = 9.70297e-5; // 9.81/1.01e5
  PorousMedia::initial_step = false;
  PorousMedia::initial_iter = false;
  PorousMedia::sum_interval = -1;
  PorousMedia::NUM_SCALARS  = 0;
  PorousMedia::NUM_STATE    = 0;
  PorousMedia::full_cycle   = 0;
  PorousMedia::max_step     = 0;
  PorousMedia::stop_time    = 0;

  PorousMedia::be_cn_theta           = 0.5;
  PorousMedia::visc_tol              = 1.0e-10;  
  PorousMedia::visc_abs_tol          = 1.0e-10;  
  PorousMedia::def_harm_avg_cen2edge = false;

  PorousMedia::have_capillary = 0;

  PorousMedia::saturation_threshold_for_vg_Kr = -1; // <0 bypasses smoothing
  PorousMedia::use_shifted_Kr_eval = 0; //

  PorousMedia::variable_scal_diff = 1; 

  PorousMedia::do_chem            = 0;
  PorousMedia::do_tracer_transport          = -1;
  PorousMedia::do_full_strang     = 1;
  PorousMedia::n_chem_interval    = 0;
  PorousMedia::it_chem            = 0;
  PorousMedia::dt_chem            = 0;
  PorousMedia::max_grid_size_chem = 16;
  PorousMedia::no_initial_values  = true;
  PorousMedia::use_funccount      = false;

#ifdef AMANZI
  PorousMedia::n_total    = 0;
  PorousMedia::n_minerals = 0;
  PorousMedia::n_sorbed   = 0;
#endif

  PorousMedia::do_simple           = 0;
  PorousMedia::do_multilevel_full  = 0;
  PorousMedia::do_reflux           = 1;
  PorousMedia::do_correct          = 0;
  PorousMedia::no_corrector        = 0;
  PorousMedia::do_kappa_refine     = 0;
  PorousMedia::n_pressure_interval = 0;
  PorousMedia::it_pressure         = 0;  
  PorousMedia::do_any_diffuse      = false;
  PorousMedia::do_cpl_advect       = 0;
  PorousMedia::do_richard_sat_solve = false;

  PorousMedia::richard_solver_verbose = 1;

  PorousMedia::do_richard_init_to_steady = false;
  PorousMedia::richard_init_to_steady_verbose = 1;
  PorousMedia::steady_min_iterations = 10;
  PorousMedia::steady_max_iterations = 15;
  PorousMedia::steady_limit_iterations = 20;
  PorousMedia::steady_time_step_reduction_factor = 0.8;
  PorousMedia::steady_time_step_increase_factor = 1.25;
  PorousMedia::steady_time_step_retry_factor_1 = 0.5;
  PorousMedia::steady_time_step_retry_factor_2 = 0.1;
  PorousMedia::steady_time_step_retry_factor_f = 0.01;
  PorousMedia::steady_max_consecutive_failures_1 = 3;
  PorousMedia::steady_max_consecutive_failures_2 = 4;
  PorousMedia::steady_tolerance = 1.e-8;
  PorousMedia::steady_init_time_step = 1.e2;
  PorousMedia::steady_max_time_steps = 8000;
  PorousMedia::steady_max_psuedo_time = 1.e10;
  PorousMedia::steady_max_num_consecutive_success = 3;
  PorousMedia::steady_extra_time_step_increase_factor = 10.;
  PorousMedia::steady_max_num_consecutive_increases = 3;
  PorousMedia::steady_consecutive_increase_reduction_factor = 0.15;

  PorousMedia::richard_max_ls_iterations = 10;
  PorousMedia::richard_min_ls_factor = 1.e-8;
  PorousMedia::richard_ls_acceptance_factor = 1.4;
  PorousMedia::richard_ls_reduction_factor = 0.1;
  PorousMedia::richard_monitor_linear_solve = 0;
  PorousMedia::richard_monitor_line_search = 0;

  PorousMedia::echo_inputs         = 0;
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

  setup_list();
  std::string pp_dump_file = ""; 
  if (pproot.countval("dump_parmparse_table")) {
      pproot.get("dump_parmparse_table",pp_dump_file);
      std::ofstream ofs;
      std::pair<std::string,std::string> df = SplitDirAndName(pp_dump_file);
      if (ParallelDescriptor::IOProcessor()) {
          if (!BoxLib::UtilCreateDirectory(df.first, 0755)) {
              BoxLib::CreateDirectoryFailed(df.first);
          }
      }
      ParallelDescriptor::Barrier();

      ofs.open(pp_dump_file.c_str());
      if (ofs.fail()) {
          BoxLib::Abort(std::string("Cannot open pp dump file: "+pp_dump_file).c_str());
      }
      if (verbose>1 && ParallelDescriptor::IOProcessor())
      {
          std::cout << "\nDumping ParmParse table:\n";
      }
      ParmParse::dumpTable(ofs);
      if (verbose>1 && ParallelDescriptor::IOProcessor())
      {
          std::cout << "... done dumping ParmParse table.\n" << '\n';
      }
      ofs.close();
  }

  read_params(); 
  BCRec bc;

  //
  // Set state variables Ids.
  //
  int num_gradn = ncomps;
  NUM_SCALARS   = ncomps + 2;

  if (ntracers > 0)
    NUM_SCALARS = NUM_SCALARS + ntracers;

  if (model == model_list["polymer"])
  {
    NUM_SCALARS = NUM_SCALARS + 2;
  }

  // add velocity and correction velocity
  NUM_STATE = NUM_SCALARS + BL_SPACEDIM + BL_SPACEDIM ;

  // add COREREACT stuff
#if defined(COREREACT)
    NUM_STATE = NUM_STATE + sz_corereact;
#endif

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
      tbcs[i]   = bc;
      tnames[i] = tNames[i];
    }

    desc_lst.setComponent(State_Type,
			  ncomps,
			  tnames,
			  tbcs,
			  BndryFunc(FORT_ONE_N_FILL,FORT_ALL_T_FILL));
  }

  desc_lst.setComponent(State_Type,ncomps+ntracers,"Aux1",
			bc,BndryFunc(FORT_ENTHFILL));
  desc_lst.setComponent(State_Type,ncomps+ntracers+1,"Aux2",
			bc,BndryFunc(FORT_ADVFILL));

  if (model == model_list["polymer"]) {
    desc_lst.setComponent(State_Type,ncomps+2,"s",
			  bc,BndryFunc(FORT_ONE_N_FILL));
    desc_lst.setComponent(State_Type,ncomps+3,"c",
			  bc,BndryFunc(FORT_ONE_N_FILL));
  }

  is_diffusive.resize(NUM_SCALARS);
  advectionType.resize(NUM_SCALARS);
  diffusionType.resize(NUM_SCALARS);

  // For components.
  for (int i=0; i<ncomps; i++) 
    {
      advectionType[i] = Conservative;
      diffusionType[i] = Laplacian_S;
      is_diffusive[i] = false;
      if (visc_coef[i] > 0.0 && solid.compare(pNames[pType[i]])!=0)
	is_diffusive[i] = true;
    }

  for (int i = ncomps; i < NUM_SCALARS; i++)
    {
      advectionType[i] = NonConservative;
      diffusionType[i] = Laplacian_S;
      is_diffusive[i] = false;
    }

  bool needs_aux_chem_data = aux_chem_variables.size() || 
      ( (nminerals>0)  ||  (nsorption_sites>0)  ||   (ncation_exchange>0) );

  if (needs_aux_chem_data)
  {
      if (ntracers>0)
      {
          Array<std::string> solute_data_names;
          solute_data_names.push_back("Total_Sorbed");
          solute_data_names.push_back("Kd");
          solute_data_names.push_back("Freundlich_n");
          solute_data_names.push_back("Langmuir_b");
          if (do_chem>0) {
              solute_data_names.push_back("Free_Ion_Guess");
              solute_data_names.push_back("Activity_Coefficient");
          }
          for (int i=0; i<tNames.size(); ++i) {
              const std::string& name = tNames[i];
              for (int j=0; j<solute_data_names.size(); ++j) {
                  const std::string label = solute_data_names[j]+"_"+name;
                  sorption_isotherm_label_map[name][solute_data_names[j]] = aux_chem_variables.size();
                  aux_chem_variables.push_back(label);
              }
          }
      }      

      if (nminerals>0)
      {
          Array<std::string> mineral_data_names;
          mineral_data_names.push_back("Volume_Fraction");
          mineral_data_names.push_back("Specific_Surface_Area");

          for (int i=0; i<minerals.size(); ++i) {
              const std::string& name = minerals[i];
              for (int j=0; j<mineral_data_names.size(); ++j) {
                  const std::string label = mineral_data_names[j]+"_"+name;
                  mineralogy_label_map[name][mineral_data_names[j]] = aux_chem_variables.size();
                  aux_chem_variables.push_back(label);
              }
          }
      }      

      if (nsorption_sites>0)
      {
          Array<std::string> sorption_site_data_names;
          sorption_site_data_names.push_back("Site_Density");

          for (int i=0; i<sorption_sites.size(); ++i) {
              const std::string& name = sorption_sites[i];
              for (int j=0; j<sorption_site_data_names.size(); ++j) {
                  const std::string label = sorption_site_data_names[j]+"_"+name;
                  surface_complexation_label_map[name][sorption_site_data_names[j]] = aux_chem_variables.size();
                  aux_chem_variables.push_back(label);
              }
          }
      }

      if (ncation_exchange>0)
      {
          Array<std::string> cation_exchange_data_names;
          cation_exchange_data_names.push_back("Cation_Exchange_Capacity");

          for (int j=0; j<cation_exchange_data_names.size(); ++j) {
              cation_exchange_label_map[cation_exchange_data_names[j]] = aux_chem_variables.size();
              aux_chem_variables.push_back(cation_exchange_data_names[j]);
          }
      }

      int num_aux_chem_variables = aux_chem_variables.size();
      Array<BCRec> cbcs(num_aux_chem_variables);
      for (int i = 0; i < num_aux_chem_variables; i++) 
      {
          cbcs[i] = bc;
      }

      desc_lst.addDescriptor(Aux_Chem_Type,IndexType::TheCellType(),
                             StateDescriptor::Point,0,num_aux_chem_variables,
                             &cell_cons_interp);
      desc_lst.setComponent(Aux_Chem_Type,0,aux_chem_variables,cbcs,
                            BndryFunc(FORT_ONE_N_FILL,FORT_ALL_T_FILL));

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

  desc_lst.addDescriptor(Vel_Type,IndexType::TheCellType(),
			 StateDescriptor::Point,1,BL_SPACEDIM,
			 &cell_cons_interp);
  set_x_vel_bc(bc,phys_bc);
  desc_lst.setComponent(Vel_Type,Xvel,"x_velocity",
			bc,BndryFunc(FORT_XVELFILL));
  set_y_vel_bc(bc,phys_bc);
  desc_lst.setComponent(Vel_Type,Yvel,"y_velocity",
			bc,BndryFunc(FORT_YVELFILL));
#if (BL_SPACEDIM == 3)
  set_z_vel_bc(bc,phys_bc);
  desc_lst.setComponent(Vel_Type,Zvel,"z_velocity",
			bc,BndryFunc(FORT_ZVELFILL));
#endif

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

  //
  // **************  DEFINE DERIVED QUANTITIES ********************
  //
  derive_lst.add("gradn",IndexType::TheCellType(),1,
		 FORT_DERGRDN,grow_box_by_one);
  derive_lst.addComponent("gradn",desc_lst,State_Type,0,num_gradn);
  //if (ntracers > 0)
  {
    //derive_lst.add("gradt",IndexType::TheCellType(),1,
	//	   FORT_DERGRDN,grow_box_by_one);
    //derive_lst.addComponent("gradn",desc_lst,State_Type,ncomps,1);
  }
  derive_lst.add("gradpx",IndexType::TheCellType(),1,
		 FORT_DERGRDPX,grow_box_by_one);
  derive_lst.addComponent("gradpx",desc_lst,Press_Type,Pressure,1);
  derive_lst.add("gradpy",IndexType::TheCellType(),1,
		 FORT_DERGRDPY,grow_box_by_one);
  derive_lst.addComponent("gradpy",desc_lst,Press_Type,Pressure,1);
#if (BL_SPACEDIM == 3)
  derive_lst.add("gradpz",IndexType::TheCellType(),1,
		 FORT_DERGRDPZ,grow_box_by_one);
  derive_lst.addComponent("gradpz",desc_lst,Press_Type,Pressure,1);
#endif

#if defined(AMANZI)
  if (do_chem>0)
    {
      // add function count
      desc_lst.addDescriptor(FuncCount_Type, IndexType::TheCellType(),
			     StateDescriptor::Point,1,1, &cell_cons_interp);
      desc_lst.setComponent(FuncCount_Type, 0, "FuncCount", 
			    bc, BndryFunc(FORT_ONE_N_FILL));
    }
#endif

  // "User defined" - atthough these must correspond to those in PorousMedia::derive
  IndexType regionIDtype(IndexType::TheCellType());
  int nCompRegion = 1;
  std::string amr_prefix = "amr";
  ParmParse pp(amr_prefix);
  int num_user_derives = pp.countval("user_derive_list");
  Array<std::string> user_derive_list(num_user_derives);
  pp.getarr("user_derive_list",user_derive_list,0,num_user_derives);
  for (int i=0; i<num_user_derives; ++i) {
      derive_lst.add(user_derive_list[i], regionIDtype, nCompRegion);
  }

  //
  // **************  DEFINE ERROR ESTIMATION QUANTITIES  *************
  //
  //err_list.add("gradn",1,ErrorRec::Special,ErrorFunc(FORT_ADVERROR));

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
      PArray<Region> regions = build_region_PArray(region_names);
      if (ppr.countval("val_greater_than")) {
          Real value; ppr.get("val_greater_than",value);
          std::string field; ppr.get("field",field);
          err_list.add(field.c_str(),0,ErrorRec::Special,
                       PM_Error_Value(FORT_VALGTERROR,value,min_time,max_time,max_level,regions));
      }
      else if (ppr.countval("val_less_than")) {
          Real value; ppr.get("val_less_than",value);
          std::string field; ppr.get("field",field);
          err_list.add(field.c_str(),0,ErrorRec::Special,
                       PM_Error_Value(FORT_VALLTERROR,value,min_time,max_time,max_level,regions));
      }
      else if (ppr.countval("diff_greater_than")) {
          BoxLib::Abort("Difference refinement not yet supported");
          Real value; ppr.get("diff_greater_than",value);
          std::string field; ppr.get("field",field);
          err_list.add(field.c_str(),1,ErrorRec::Special,
                       PM_Error_Value(FORT_DIFFGTERROR,value,min_time,max_time,max_level,regions));
      }
      else if (ppr.countval("in_region")) {
          Real value; ppr.get("in_region",value);
          err_list.add("PMAMR_DUMMY",1,ErrorRec::Special,
                       PM_Error_Value(min_time,max_time,max_level,regions));
      }
      else {
          BoxLib::Abort(std::string("Unrecognized refinement indicator for " + refinement_indicators[i]).c_str());
      }
  }

  BoxLib::ExecOnFinalize(PM_Setup_CleanUpStatics);

}

//
//  Read input file
//
void PorousMedia::read_geometry()
{
  //
  // Get geometry-related parameters.  
  // Note: 1. The domain size and periodity information are read in 
  //          automatically.  This function deals primarily with region
  //          definition.
  //       2. regions defined in PorousMedia.H as PArray<Region>
  //
  ParmParse pp("geometry");

  Array<Real> problo, probhi;
  pp.getarr("prob_lo",problo,0,BL_SPACEDIM);
  pp.getarr("prob_hi",probhi,0,BL_SPACEDIM);
  Region::domlo = problo;
  Region::domhi = probhi;
  
  Real geometry_eps = -1; pp.get("geometry_eps",geometry_eps);
  Region::geometry_eps = geometry_eps;

  // set up  1+2*BL_SPACEDIM default regions
  bool generate_default_regions = true; pp.query("generate_default_regions",generate_default_regions);
  int nregion_DEF = 0;
  regions.clear();
  if (generate_default_regions) {
      nregion_DEF = 1 + 2*BL_SPACEDIM;
      regions.resize(nregion_DEF,PArrayManage);
      regions.set(0, new   allRegion(problo,probhi));
      regions.set(1, new allBCRegion(0,0,problo,probhi));
      regions.set(2, new allBCRegion(0,1,problo,probhi));
      regions.set(3, new allBCRegion(1,0,problo,probhi));
      regions.set(4, new allBCRegion(1,1,problo,probhi));
#if BL_SPACEDIM == 3
      regions.set(5, new allBCRegion(2,0,problo,probhi));
      regions.set(6, new allBCRegion(2,1,problo,probhi));
#endif
  }

  // Get parameters for each user defined region 
  int nregion = nregion_DEF;

  int nregion_user = pp.countval("regions");

  if (!generate_default_regions  && nregion_user==0) {
      BoxLib::Abort("Default regions not generated and none provided.  Perhaps omitted regions list?");
  }
  if (nregion_user)
    {
      std::string r_purpose, r_type;
      Array<std::string> r_name;
      pp.getarr("regions",r_name,0,nregion_user);
      nregion += nregion_user;
      regions.resize(nregion,PArrayManage);

      for (int j=0; j<nregion_user; ++j)
      {
          const std::string prefix("geometry." + r_name[j]);
          ParmParse ppr(prefix.c_str());
	  ppr.get("purpose",r_purpose);
	  ppr.get("type",r_type);      

	  if (r_type == "point")
          {
	      Array<Real> coor;
	      ppr.getarr("coordinate",coor,0,BL_SPACEDIM);
              regions.set(nregion_DEF+j, new pointRegion(r_name[j],r_purpose,r_type,coor));
	    }
	  else if (r_type == "box" || r_type == "surface")
	    {
	      Array<Real> lo_coor,hi_coor;
	      ppr.getarr("lo_coordinate",lo_coor,0,BL_SPACEDIM);
	      ppr.getarr("hi_coordinate",hi_coor,0,BL_SPACEDIM);
              regions.set(nregion_DEF+j, new boxRegion(r_name[j],r_purpose,r_type,lo_coor,hi_coor));
	    }
	  else if (r_type == "color_function")
          {
              int color_value; ppr.get("color_value",color_value);
              std::string color_file; ppr.get("color_file",color_file);
              colorFunctionRegion* cfr = new colorFunctionRegion(r_name[j],r_purpose,r_type,color_file,color_value);
	      regions.set(nregion_DEF+j, cfr);
          }
          else {
              std::string m = "region type not supported \"" + r_type + "\"";
              BoxLib::Abort(m.c_str());
          }
	}
      pp.query("surf_file",surf_file);
    }
}

bool check_if_layered(const Array<std::string>& material_region_names,
                      const PArray<Region>&     regions,
                      const Array<Real>&        plo,
                      const Array<Real>&        phi)
{
    for (int i=0; i<material_region_names.size(); ++i)
    {
        const std::string& name = material_region_names[i];
        for (int j=0; j<regions.size(); ++j)
        {
            const Region* r = &(regions[j]);
            if (r->name == name) {
                const boxRegion* testp = dynamic_cast<const boxRegion*>(r);
                if (testp) {
                    for (int d=0; d<BL_SPACEDIM-1; ++d) {
                        if (testp->lo[d] > plo[d]  ||  testp->hi[d] < phi[d]) {
                            return false;
                        }
                    }
                }
                else {
                    std::cout << *r << std::endl;
                    return false; // cant be layered if region is not a box
                }
            }
        }
    }
    return true;
}

void
PorousMedia::read_rock()
{
    //
    // Get parameters related to rock
    //
    ParmParse pp("rock");
    int nrock = pp.countval("rock");
    if (nrock <= 0) {
        BoxLib::Abort("At least one rock type must be defined.");
    }
    Array<std::string> r_names;  pp.getarr("rock",r_names,0,nrock);
    rocks.clear();
    rocks.resize(nrock,PArrayManage);

    Array<std::string> material_regions;
    for (int i = 0; i<nrock; i++)
    {
        const std::string& rname = r_names[i];
        const std::string prefix("rock." + rname);
        ParmParse ppr(prefix.c_str());
        
        Real rdensity = -1; // ppr.get("density",rdensity); // not actually used anywhere
        Real rporosity; ppr.get("porosity",rporosity);    
        Array<Real> rpermeability; ppr.getarr("permeability",rpermeability,0,ppr.countval("permeability"));
        BL_ASSERT(rpermeability.size() == 2); // Horizontal, Vertical

        // The permeability is specified in mDa.  
        // This needs to be multiplied with 1e-7 to be consistent 
        // with the other units in the code
        for (int j=0; j<rpermeability.size(); ++j) {
            rpermeability[j] *= 1.e-7;
        }

        // relative permeability: include kr_coef, sat_residual
        int rkrType = 0;  ppr.query("kr_type",rkrType);
        Array<Real> rkrParam;
        if (rkrType > 0) {
            ppr.getarr("kr_param",rkrParam,0,ppr.countval("kr_param"));
        }

        // capillary pressure: include cpl_coef, sat_residual, sigma
        int rcplType = 0;  ppr.query("cpl_type", rcplType);
        Array<Real> rcplParam;
        if (rcplType > 0) {
            ppr.getarr("cpl_param",rcplParam,0,ppr.countval("cpl_param"));
        }

        Array<std::string> region_names;
        ppr.getarr("regions",region_names,0,ppr.countval("regions"));
        PArray<Region> rregions = build_region_PArray(region_names);
        for (int j=0; j<region_names.size(); ++j) {
            material_regions.push_back(region_names[j]);
        }

        std::string porosity_dist; ppr.get("porosity_dist",porosity_dist);
        int rporosity_dist_type = Rock::rock_dist_map[porosity_dist];
        Array<Real> rporosity_dist_param;
        if (rporosity_dist_type > 1) {
            ppr.getarr("porosity_dist_param",rporosity_dist_param,
                       0,ppr.countval("porosity_dist_param"));
        }
        
        std::string permeability_dist; ppr.get("permeability_dist",permeability_dist);
        int rpermeability_dist_type = Rock::rock_dist_map[permeability_dist];
        Array<Real> rpermeability_dist_param;
        if (rpermeability_dist_type > 1)
        {
            ppr.getarr("permeability_dist_param",rpermeability_dist_param,
                       0,ppr.countval("permeability_dist_param"));
        }

        rocks.set(i, new Rock(rname,rdensity,rporosity,rporosity_dist_type,rporosity_dist_param,
                              rpermeability,rpermeability_dist_type,rpermeability_dist_param,
                              rkrType,rkrParam,rcplType,rcplParam,rregions));
                  
        if (ntracers>0 && do_chem>0) {
            Array<std::pair<std::string,Real> > parameters;
            parameters.push_back(std::make_pair<std::string,Real>(        "Total_Sorbed", 0));
            parameters.push_back(std::make_pair<std::string,Real>(                  "Kd", 0));
            parameters.push_back(std::make_pair<std::string,Real>(          "Langmuir_b", 0));
            parameters.push_back(std::make_pair<std::string,Real>(        "Freundlich_n", 1));
            parameters.push_back(std::make_pair<std::string,Real>(      "Free_Ion_Guess", 0));
            parameters.push_back(std::make_pair<std::string,Real>("Activity_Coefficient", 0));

            for (int k=0; k<tNames.size(); ++k) {
                for (int j=0; j<parameters.size(); ++j) {
                    const std::string& str = parameters[j].first;
                    const std::string prefixM(prefix+".Sorption_Isotherms."+tNames[k]);
                    ParmParse pprm(prefixM.c_str());
                    sorption_isotherm_ics[rname][tNames[k]][str] = parameters[j].second;
                    pprm.query(str.c_str(),sorption_isotherm_ics[rname][tNames[k]][str]);
                }
            }
        }

        if (nminerals>0 && do_chem>0) {
            Array<std::pair<std::string,Real> > parameters;
            parameters.push_back(std::make_pair<std::string,Real>(       "Volume_Fraction", 0));
            parameters.push_back(std::make_pair<std::string,Real>("Specific_Surface_Area", 0));

            for (int k=0; k<minerals.size(); ++k) {
                for (int j=0; j<parameters.size(); ++j) {
                    const std::string& str = parameters[j].first;
                    const std::string prefixM(prefix+".Mineralogy."+minerals[k]);
                    ParmParse pprm(prefixM.c_str());
                    mineralogy_ics[rname][minerals[k]][str] = parameters[j].second;
                    pprm.query(str.c_str(),mineralogy_ics[rname][minerals[k]][str]);
                }
            }
        }
        
        if (nsorption_sites>0 && do_chem>0) {
            Array<std::pair<std::string,Real> > parameters;
            parameters.push_back(std::make_pair<std::string,Real>("Site_Density", 0));

            for (int k=0; k<sorption_sites.size(); ++k) {
                for (int j=0; j<parameters.size(); ++j) {
                    const std::string& str = parameters[j].first;
                    const std::string prefixM(prefix+".Surface_Complexation_Sites."+sorption_sites[k]);
                    ParmParse pprm(prefixM.c_str());
                    surface_complexation_ics[rname][sorption_sites[k]][str] = parameters[j].second;
                    pprm.query(str.c_str(),surface_complexation_ics[rname][sorption_sites[k]][str]);
                }
            }
        }

        if (ncation_exchange>0 && do_chem>0) {
            Array<std::pair<std::string,Real> > parameters;
            parameters.push_back(std::make_pair<std::string,Real>("Cation_Exchange_Capacity", 0));
            for (int j=0; j<parameters.size(); ++j) {
                cation_exchange_ics[rname][parameters[j].first] = parameters[j].second;
            }
        }
        
    }
    
    bool read_full_pmap  = false;
    bool read_full_kmap  = false;
    bool build_full_pmap = true;
    bool build_full_kmap = true;
    permeability_from_fine = true;
    porosity_from_fine = true;
    
    std::string kfile, pfile,gsfile;
    pp.query("permeability_file", kfile);
    pp.query("porosity_file", pfile);
    pp.query("gslib_file",gsfile);

    pp.query("Use_Shifted_Kr_Eval",use_shifted_Kr_eval);
    pp.query("Saturation_Threshold_For_vg_Kr",saturation_threshold_for_vg_Kr);
    FORT_KR_INIT(&saturation_threshold_for_vg_Kr,
		 &use_shifted_Kr_eval);
    
    // The I/O processor makes the directory if it doesn't already exist.
    if (ParallelDescriptor::IOProcessor())
        if (!BoxLib::UtilCreateDirectory(pfile, 0755))
            BoxLib::CreateDirectoryFailed(pfile);
    ParallelDescriptor::Barrier();
    if (ParallelDescriptor::IOProcessor())
        if (!BoxLib::UtilCreateDirectory(kfile, 0755))
            BoxLib::CreateDirectoryFailed(kfile);
    ParallelDescriptor::Barrier();
    
    pfile += "/pp";
    kfile += "/kp";
    
    if (read_full_kmap)
    {
        std::cout << "Current code only allows reading in "
                  << "the distribution in full." << std::endl;
        
        if (kappadata == 0)
            kappadata = new MultiFab;
        
        VisMF::Read(*kappadata,kfile);
    }
    
    if (read_full_pmap)
    {
        std::cout << "Current code only allows reading in "
                  << "the distribution in full.\n";
        
        if (phidata == 0)
            phidata = new MultiFab;
        
        VisMF::Read(*phidata,pfile);
    }
    
    // determine parameters needed to build kappadata and phidata
    int max_level;
    Array<int> n_cell, fratio;
    Array<Real> problo, probhi;
    if (build_full_kmap || build_full_pmap)
    { 
        ParmParse am("amr");
        am.query("max_level",max_level);
        am.getarr("n_cell",n_cell,0,BL_SPACEDIM);
        if (max_level>0) {
            am.getarr("ref_ratio",fratio,0,max_level);
        }
        
        ParmParse gm("geometry");
        gm.getarr("prob_lo",problo,0,BL_SPACEDIM);
        gm.getarr("prob_hi",probhi,0,BL_SPACEDIM);
    }
    
    // construct permeability field based on the specified parameters
    if (build_full_kmap)
    {
        
        if (verbose > 1 && ParallelDescriptor::IOProcessor()) 
            std::cout << "Building kmap on finest level..." << std::endl;

        BoxArray ba = Rock::build_finest_data(max_level, n_cell, fratio);
        
        if (kappadata == 0) {

            kappadata = new MultiFab(ba,BL_SPACEDIM,0,Fab_allocate);
        
            for (int i=0; i<rocks.size(); ++i) 
            {
                // these are temporary work around.   
                // Should utilizes region to determine size.
                Rock& r = rocks[i];
                r.max_level = max_level;
                r.n_cell = n_cell;
                r.fratio = fratio;
                r.problo = problo;
                r.probhi = probhi;
                r.build_kmap(*kappadata, gsfile);
            }
            
            if (verbose > 1 && ParallelDescriptor::IOProcessor()) 
                std::cout << "   Finished building kmap on finest level.  Writing..." << std::endl;
            
            VisMF::SetNOutFiles(10); // FIXME: Should not be hardwired here
            VisMF::Write(*kappadata,kfile);
            
            if (verbose > 1 && ParallelDescriptor::IOProcessor()) 
                std::cout << "   Finished writing kmap..." << std::endl;
        }
    }
    
    if (build_full_pmap)
    {
        if (verbose > 1 && ParallelDescriptor::IOProcessor()) 
            std::cout << "Building pmap on finest level..." << std::endl;

        BoxArray ba = Rock::build_finest_data(max_level, n_cell, fratio);
        
        if (phidata == 0) {

            phidata = new MultiFab(ba,1,0,Fab_allocate);

            for (int i=0; i<rocks.size(); ++i)
            {
                // these are temporary work around.   
                // Should utilizes region to determine size.
                Rock& r = rocks[i];
                r.max_level = max_level;
                r.n_cell = n_cell;
                r.fratio = fratio;
                r.problo = problo;
                r.probhi = probhi;
                r.build_pmap(*phidata, gsfile);
            }
            
            if (verbose > 1 && ParallelDescriptor::IOProcessor()) 
                std::cout << "   Finished building pmap on finest level.  Writing..." << std::endl;
            
            VisMF::SetNOutFiles(10); // FIXME: Should not be hardwired here
            VisMF::Write(*phidata,pfile);
            
            if (verbose > 1 && ParallelDescriptor::IOProcessor()) 
                std::cout << "   Finished writing pmap..." << std::endl;
        }
    }

    if (ParallelDescriptor::IOProcessor()) {
        std::cout << "Rock name mapping in output: " << std::endl;
    }
    for (int i=0; i<rocks.size(); ++i) 
    {
        if (ParallelDescriptor::IOProcessor()) {
            std::cout << "Rock: " << rocks[i].name << " -> " << i << std::endl;
        }
    }

    // Check if material is actually layered in the vertical coordinate
    material_is_layered = check_if_layered(material_regions,regions,problo,probhi);

    // Check that at the coarse level material properties have been set over the entire domain
    BoxArray ba(Box(IntVect(D_DECL(0,0,0)),
                    IntVect(D_DECL(n_cell[0]-1,n_cell[1]-1,n_cell[2]-1))));
    ba.maxSize(32);
    MultiFab rockTest(ba,1,0);
    rockTest.setVal(-1);
    Array<Real> dx0(BL_SPACEDIM);
    for (int i=0; i<BL_SPACEDIM; ++i) {
        dx0[i] = (probhi[i] - problo[i])/n_cell[i];
    }
    for (MFIter mfi(rockTest); mfi.isValid(); ++mfi)
    {
        FArrayBox& fab = rockTest[mfi];
        for (int i=0; i<rocks.size(); ++i)
        {
            Real val = (Real)i;
            const PArray<Region>& rock_regions = rocks[i].regions;
            for (int j=0; j<rock_regions.size(); ++j) {
                rock_regions[j].setVal(fab,val,0,dx0.dataPtr(),0);
            }
        }
    }
    if (rockTest.min(0) < 0) {
        BoxLib::Abort("Material has not been defined over the entire domain");
    }
}

void PorousMedia::read_prob()
{
  ParmParse pp;

  max_step  = -1;    
  stop_time = -1.0;  

  pp.query("max_step",max_step);
  pp.query("stop_time",stop_time);

  std::string event_name = "Stop_Time";
  defined_events[event_name] = new EventCoord::TimeEvent(Array<Real>(1,stop_time));
  PMAmr::eventCoord().Register(event_name,defined_events[event_name]);

  //
  // Get run options.
  //  
  ParmParse pb("prob");

  // determine the model based on model_name
  pb.query("model_name",model_name);
  model = model_list[model_name];

  // if model is specified, this supersedes model_name
  pb.query("model",model);

  pb.query("do_tracer_transport",do_tracer_transport);

  // Verbosity
  pb.query("v",verbose);
  pb.query("richard_solver_verbose",richard_solver_verbose);
  pb.query("do_richard_sat_solve",do_richard_sat_solve);

  pb.query("richard_init_to_steady_verbose",richard_init_to_steady_verbose);
  pb.query("do_richard_init_to_steady",do_richard_init_to_steady);
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
  pb.query("steady_tolerance",steady_tolerance);
  pb.query("steady_init_time_step",steady_init_time_step);
  pb.query("steady_max_time_steps",steady_max_time_steps);
  pb.query("steady_max_psuedo_time",steady_max_psuedo_time);
  pb.query("steady_max_num_consecutive_success",steady_max_num_consecutive_success);
  pb.query("steady_extra_time_step_increase_factor",steady_extra_time_step_increase_factor);
  pb.query("steady_max_num_consecutive_increases",steady_max_num_consecutive_increases);
  pb.query("consecutive_increase_reduction_factor",steady_consecutive_increase_reduction_factor);
  pb.query("richard_monitor_linear_solve",richard_monitor_linear_solve);
  pb.query("richard_monitor_linear_solve",richard_monitor_line_search);
  
  // Get timestepping parameters.
  pb.get("cfl",cfl);
  pb.query("init_shrink",init_shrink);
  pb.query("dt_init",dt_init);
  pb.query("dt_cutoff",dt_cutoff);
  pb.query("change_max",change_max);
  pb.query("fixed_dt",fixed_dt);
  pb.query("max_dt",richard_max_dt);
  pb.query("sum_interval",sum_interval);

  // Gravity are specified as m/s^2 in the input file
  // This is converted to the unit that is used in the code.
  if (pb.contains("gravity"))
    {
      pb.get("gravity",gravity);
      gravity /= 1.01e5;
    }

  // Get algorithmic flags and options
  pb.query("full_cycle", full_cycle);
  //pb.query("algorithm", algorithm);
  pb.query("do_multilevel_full",  do_multilevel_full );
  pb.query("do_simple",  do_simple );
  pb.query("do_reflux",  do_reflux );
  pb.query("do_correct", do_correct);
  pb.query("do_cpl_advect", do_cpl_advect);
  pb.query("no_corrector",no_corrector);
  pb.query("do_kappa_refine",do_kappa_refine);
  pb.query("n_pressure_interval",n_pressure_interval);

  // Get solver tolerances
  pb.query("visc_tol",visc_tol);
  pb.query("visc_abs_tol",visc_abs_tol);
  pb.query("be_cn_theta",be_cn_theta);
  if (be_cn_theta > 1.0 || be_cn_theta < .5)
    BoxLib::Abort("PorousMedia::read_prob():Must have be_cn_theta <= 1.0 && >= .5");   
  pb.query("harm_avg_cen2edge", def_harm_avg_cen2edge);

  // if capillary pressure flag is true, then we make sure 
  // the problem can handle capillary pressure.
  pb.query("have_capillary",have_capillary);
  if (have_capillary == 1) 
    {
      if (nphases != 2 && ncomps !=nphases) 
	{
	  if (ParallelDescriptor::IOProcessor())
	    {
	      std::cerr << "PorousMedia::read_prob: nphases != 2 && ncomps !=nphases "
			<< "although have_capillary == 1.\n ";
	      BoxLib::Abort("PorousMedia::read_prob()");
	    }
	}
    }
}

//
// Construct bc functions
//

PArray<Region>
PorousMedia::build_region_PArray(const Array<std::string>& region_names)
{
    PArray<Region> ret(region_names.size(), PArrayNoManage);
    for (int i=0; i<region_names.size(); ++i)
    {
        const std::string& name = region_names[i];
        bool found = false;
        for (int j=0; j<regions.size() && !found; ++j)
        {
            Region& r = regions[j];
            if (regions[j].name == name) {
                found = true;
                ret.set(i,&r);
            }
        }
        if (!found) {
            std::string m = "Named region not found " + name;
            BoxLib::Error(m.c_str());
        }
    }
    return ret;
}

const Rock&
PorousMedia::find_rock(const std::string& name)
{
    bool found=false;
    int iRock = -1;
    for (int i=0; i<rocks.size() && !found; ++i)
    {
        const Rock& rock = rocks[i];
        if (name == rock.name) {
            found = true;
            iRock = i;
        }
    } 
    if (iRock < 0) {
        std::string m = "Named rock not found " + name;
        BoxLib::Abort(m.c_str());
    }
    return rocks[iRock];
}


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
        if (cNames[j] == "Water") {
            idx = j;
        }
    }
    BL_ASSERT(idx>=0);

    Array<double> rhoSat(ncomps,0);
    rhoSat[idx] = density[idx] * 1; // Fully saturated...an assumption
    return rhoSat;
}

struct FluxToRhoSat
    : public ArrayTransform
{
    FluxToRhoSat(const Rock& rock)
        : rock(rock) {}
    virtual ArrayTransform* clone() const {return new FluxToRhoSat(*this);}
    virtual Array<Real> transform(Real inData) const;
protected:
    const Rock& rock;
};

Array<Real>
FluxToRhoSat::transform(Real aqueous_Darcy_flux) const
{
    Real gravity = PorousMedia::getGravity();
    const Array<Real>& density = PorousMedia::Density(); // Assumes 1 component per phase
    int ncomps = density.size();
    BL_ASSERT(ncomps>0 && ncomps<=2);
        
    Real lkappa = rock.permeability[0];
    Real gstar;
    if (density.size() > 1)
        gstar = -lkappa*(density[0]-density[1])*gravity;
    else
        gstar = -lkappa*(density[0])*gravity;
    
    // Compute saturation given Aqueous flow rate
    int lkrtype = rock.krType;
    Real lkrcoef = rock.krParam[0];
    Real lsatres = rock.krParam[1];            
    int nc = 1;
    Real vtot = 0.; // Zero total velocity
    Real sol;
    const Array<Real>& visc = PorousMedia::Viscosity();
    
    Array<Real> rhoSat(ncomps);

    //std::cout << "aqueous_Darcy_flux: " << aqueous_Darcy_flux << std::endl;
    //std::cout << "gstar: " << gstar << std::endl;
    //std::cout << "visc: " << visc[0] << std::endl;

    if (ncomps > 1) 
      FORT_FIND_INV_FLUX(&sol, &aqueous_Darcy_flux, &nc, &vtot, &gstar,
			 visc.dataPtr(),&ncomps,&lkrtype,&lkrcoef);
    else
      FORT_FIND_INV_RFLUX(&sol,&aqueous_Darcy_flux, &gstar,visc.dataPtr(),
			  &ncomps,&lkrtype,&lkrcoef);
    
    rhoSat[0] = density[0]*(sol*(1.0-lsatres)+lsatres);
    if (ncomps > 1) {
        rhoSat[1] = density[1]*(1.0-rhoSat[0]/density[0]);
    }
    //std::cout << "rhoSat: " << rhoSat[0] << std::endl;
    //BoxLib::Abort();
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
      // FIXME: Assume visc given in mass units, George hardwaired rho top 
      Real p_visc; ppr.get("viscosity",p_visc); p_visc *= 1.e3; muval.push_back(p_visc);
      Real p_diff; ppr.get("diffusivity",p_diff); visc_coef.push_back(p_diff);

      // Only components have diffusion at the moment.  
      if (visc_coef.back() > 0)
      {
	  do_any_diffuse = true;
	  is_diffusive[visc_coef.size()-1] = 1;
      }
      else {
          variable_scal_diff = 0;
      }
      ++ndiff;

      pType.push_back(phase_list[pNames[i]]);
  }

  ParmParse cp("comp");
  for (int i = 0; i<ncomps; i++) comp_list[cNames[i]] = i;
#if 0

  // Get the dominant component
  std::string domName;
  cp.query("dominant",domName);
  if (!domName.empty())
    idx_dominant = comp_list[domName];

  // Get the boundary conditions for the components
  Array<int> lo_bc(BL_SPACEDIM), hi_bc(BL_SPACEDIM);
  cp.getarr("lo_bc",lo_bc,0,BL_SPACEDIM);
  cp.getarr("hi_bc",hi_bc,0,BL_SPACEDIM);
  for (int i = 0; i < BL_SPACEDIM; i++)
    {
      phys_bc.setLo(i,lo_bc[i]);
      phys_bc.setHi(i,hi_bc[i]);
    }
    
  // Check phys_bc against possible periodic geometry: 
  //  if periodic, that boundary must be internal BC.
  if (Geometry::isAnyPeriodic())
    {      
      // Do idiot check.  Periodic means interior in those directions.
      for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
	  if (Geometry::isPeriodic(dir))
	    {
	      if (lo_bc[dir] != Interior)
		{
		  std::cerr << "PorousMedia::variableSetUp:periodic in direction "
			    << dir
			    << " but low BC is not Interior\n";
		  BoxLib::Abort("PorousMedia::read_comp()");
		}
	      if (hi_bc[dir] != Interior)
		{
		  std::cerr << "PorousMedia::variableSetUp:periodic in direction "
			    << dir
			    << " but high BC is not Interior\n";
		  BoxLib::Abort("PorousMedia::read_comp()");
		}
	    } 
        }
    }
  else
    {
      
      // Do idiot check.  If not periodic, should be no interior.
      for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
	  if (!Geometry::isPeriodic(dir))
	    {
	      if (lo_bc[dir] == Interior)
		{
		  std::cerr << "PorousMedia::variableSetUp:Interior bc in direction "
			    << dir
			    << " but not defined as periodic\n";
		  BoxLib::Abort("PorousMedia::read_comp()");
		}
	      if (hi_bc[dir] == Interior)
		{
		  std::cerr << "PorousMedia::variableSetUp:Interior bc in direction "
			    << dir
			    << " but not defined as periodic\n";
		  BoxLib::Abort("PorousMedia::read_comp()");
		}
	    }
        }
    }
#endif

  // Initial condition and boundary condition
  //
  // Component ics, bcs will be set all at once
  int n_ics = cp.countval("ic_labels");
  if (n_ics > 0)
  {
      Array<std::string> ic_names;
      cp.getarr("ic_labels",ic_names,0,n_ics);
      ic_array.resize(n_ics,PArrayManage);
      for (int i = 0; i<n_ics; i++)
      {
          const std::string& icname = ic_names[i];
	  const std::string prefix("comp.ics." + icname);
	  ParmParse ppr(prefix.c_str());
          
	  int n_ic_regions = ppr.countval("regions");
          Array<std::string> region_names;
	  ppr.getarr("regions",region_names,0,n_ic_regions);
          PArray<Region> ic_regions = build_region_PArray(region_names);

          std::string ic_type; ppr.get("type",ic_type);          
          if (ic_type == "pressure")
          {
              int nPhase = pNames.size();
              Array<Real> vals(nPhase);
              
              int num_phases_reqd = nPhase;
              std::map<std::string,bool> phases_set;
              for (int j = 0; j<pNames.size(); j++)
              {
                  std::string val_name = "val";                  
#if 1
                  ParmParse& pps = ppr;
#else
                  ParmParse pps(prefix + "." + pNames[j]);
#endif
                  ppr.get(val_name.c_str(),vals[0]);          
                  phases_set[pNames[j]] = true;
              }

              int num_phases = phases_set.size();
              if (num_phases != num_phases_reqd) {
                  std::cerr << icname << ": Insufficient number of phases specified" << std::endl;
                  std::cerr << " ngiven, nreqd: " << num_phases << ", " << num_phases_reqd << std::endl;
                  std::cerr << " current model: " << model << std::endl;
                  BoxLib::Abort();
              }
              
              Array<Real> times(1,0);
              Array<std::string> forms(0);
              ic_array.set(i, new ArrayRegionData(icname,times,vals,forms,ic_regions,ic_type,1));
          }
          else if (ic_type == "saturation")
          {
              Array<Real> vals(ncomps);
              for (int j = 0; j<cNames.size(); j++) {
                  ppr.get(cNames[j].c_str(),vals[j]);
                  vals[j] *= density[j];
              }
              std::string generic_type = "scalar";
              ic_array.set(i, new RegionData(icname,ic_regions,generic_type,vals));
          }
          else if (ic_type == "hydrostatic")
          {
              Array<Real> water_table_height(1); ppr.get("water_table_height",water_table_height[0]);
              Array<Real> times(1,0);
              Array<std::string> forms;
              ic_array.set(i, new ArrayRegionData(icname,times,water_table_height,
                                                  forms,ic_regions,ic_type,1));
          }
          else if (ic_type == "zero_total_velocity")
          {
	      Array<Real> vals(4);
              Array<Real> times(1,0);
              Array<std::string> forms;
	      ppr.get("aqueous_vol_flux",vals[0]);
              ppr.get("water_table_height",vals[1]);
              Real aqueous_ref_pres = 0; ppr.query("val",vals[2]);
              Real aqueous_pres_grad = 0; ppr.query("grad",vals[3]);
	      ic_array.set(i,new RegionData(icname,ic_regions,ic_type,vals));
          }
          else {
              BoxLib::Abort("Unsupported comp ic");
          }
      }
  }

  int n_bcs = cp.countval("bc_labels");
  if (n_bcs > 0)
  {
      rinflow_bc_lo.resize(BL_SPACEDIM,0); 
      rinflow_bc_hi.resize(BL_SPACEDIM,0); 
      inflow_bc_lo.resize(BL_SPACEDIM,0); 
      inflow_bc_hi.resize(BL_SPACEDIM,0); 

      bc_array.resize(n_bcs,PArrayManage);
      Array<std::string> bc_names;
      cp.getarr("bc_labels",bc_names,0,n_bcs);

      // default to no flow first.
      for (int j=0;j<BL_SPACEDIM;j++) {
	phys_bc.setLo(j,4);
	pres_bc.setLo(j,4);
	phys_bc.setHi(j,4);
	pres_bc.setHi(j,4);
      }	  

      EventCoord& event_coord = PMAmr::eventCoord();
      for (int i = 0; i<n_bcs; i++)
      {
          const std::string& bcname = bc_names[i];
	  const std::string prefix("comp.bcs." + bcname);
	  ParmParse ppr(prefix.c_str());
          
	  int n_bc_regions = ppr.countval("regions");
          Array<std::string> region_names;
	  ppr.getarr("regions",region_names,0,n_bc_regions);
          const PArray<Region> bc_regions = build_region_PArray(region_names);
          std::string bc_type; ppr.get("type",bc_type);

          bool is_inflow = false;
          int component_bc = 4;
	  int pressure_bc  = 4;

          if (bc_type == "pressure")
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
              
              // convert to atm with datum at atmospheric pressure
              for (int j=0; j<vals.size(); ++j) {
                  vals[j] = vals[j] / 1.01325e5 - 1.e0;
              }
              
              is_inflow = false;
              component_bc = 1;
              pressure_bc = 2;

              PressToRhoSat p_to_sat;
              bc_array.set(i, new Transform_S_AR_For_BC(bcname,times,vals,forms,bc_regions,
                                                        bc_type,ncomps,p_to_sat));
              defined_events[bcname] = new EventCoord::TimeEvent(times);
              event_coord.Register(bcname,defined_events[bcname]);
          }
          else if (bc_type == "zero_total_velocity")
          {
              Array<Real> vals, times;
              Array<std::string> forms;

              int nv = ppr.countval("aqueous_vol_flux");
              if (nv) {
                  ppr.getarr("aqueous_vol_flux",vals,0,nv); // "inward" flux
                  times.resize(nv,0);
                  if (nv>1) {
                      ppr.getarr("inflowtimes",times,0,nv);
                      ppr.getarr("inflowfncs",forms,0,nv-1);
                  }
              }
              else {
                  vals.resize(1,0);
                  times.resize(1,0);
                  forms.resize(0);
              }        

              // Work out sign of flux for this boundary
              int is_hi = -1;
              for (int j=0; j<bc_regions.size(); ++j)
              {
                  const std::string purpose = bc_regions[j].purpose;
                  for (int k=0; k<7; ++k) {
                      if (purpose == PMAMR::RpurposeDEF[k]) {
                          BL_ASSERT(k != 6);
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

              is_inflow = true;
              component_bc = 1;
              pressure_bc = 1;
	      bc_array.set(i,new ArrayRegionData(bcname,times,vals,forms,bc_regions,bc_type,1));
              defined_events[bcname] = new EventCoord::TimeEvent(times);
              event_coord.Register(bcname,defined_events[bcname]);
          }
          else if (bc_type == "noflow")
          {
              is_inflow = false;
              component_bc = 4;
              pressure_bc = 4;

              Array<Real> val(1,0);
              bc_array.set(i, new RegionData(bcname,bc_regions,bc_type,val));
          }
          else
          {
	    std::cout << bc_type << " not a valid bc_type " << std::endl;
	    BoxLib::Abort();
          }

          // Some clean up 
          std::set<std::string> o_set;

          for (int j=0; j<bc_regions.size(); ++j)
          {
              const std::string purpose = bc_regions[j].purpose;
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

                  if (is_consistent) {
                      BoxLib::Abort("Inconconsistent type for boundary ");
                  }
              }
          }
      }
  }
}


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

      for (int i = 0; i<ntracers; i++)
      {
          const std::string prefix("tracer." + tNames[i]);
	  ParmParse ppr(prefix.c_str());
          if (do_chem>0  ||  do_tracer_transport == 1) {
              std::string g="Total"; ppr.query("group",g); // FIXME: is this relevant anymore?
              group_map[g].push_back(i+ncomps);
          }

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
              const std::string prefixIC(prefix + "." + tic_names[n]);
              ParmParse ppri(prefixIC.c_str());
              int n_ic_region = ppri.countval("regions");
              Array<std::string> region_names;
              ppri.getarr("regions",region_names,0,n_ic_region);
              const PArray<Region> tic_regions = build_region_PArray(region_names);
              std::string tic_type; ppri.get("type",tic_type);
              
              if (tic_type == "concentration")
              {
                  Real val = 0; ppri.query("val",val);
                  tic_array[i].set(n, new RegionData(tNames[i],tic_regions,tic_type,val));
              }
              else {
                  std::string m = "Tracer IC: \"" + tic_names[n] 
                      + "\": Unsupported tracer IC type: \"" + tic_type + "\"";
                  BoxLib::Abort(m.c_str());
              }
          }
              
          if (do_tracer_transport)
          {
              Array<std::string> tbc_names;
              int n_tbc = ppr.countval("tbcs");
              if (n_tbc <= 0)
              {
                  BoxLib::Abort("each tracer requires boundary conditions");
              }
              ppr.getarr("tbcs",tbc_names,0,n_tbc);
              tbc_array[i].resize(n_tbc,PArrayManage);
              
              for (int n = 0; n<n_tbc; n++)
              {
                  const std::string prefixTBC(prefix + "." + tbc_names[n]);
                  ParmParse ppri(prefixTBC.c_str());
                  
                  int n_tbc_region = ppri.countval("regions");
                  Array<std::string> tbc_region_names;
                  ppri.getarr("regions",tbc_region_names,0,n_tbc_region);
                  const PArray<Region> tbc_regions = build_region_PArray(tbc_region_names);
                  std::string tbc_type; ppri.get("type",tbc_type);
                  
                  if (tbc_type == "concentration")
                  {
                      Array<Real> times, vals;
                      Array<std::string> forms;
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
                      int nComp = 1;
                      tbc_array[i].set(n, new ArrayRegionData(tbc_names[n],times,vals,forms,tbc_regions,tbc_type,nComp));
                  }
                  else if (tbc_type == "noflow"  ||  tbc_type == "outflow")
                  {
                      Array<Real> val(1,0);
                      tbc_array[i].set(n, new RegionData(tbc_names[n],tbc_regions,tbc_type,val));
                  }
                  else {
                      std::string m = "Tracer BC: \"" + tbc_names[n] 
                          + "\": Unsupported tracer BC type: \"" + tbc_type + "\"";
                      BoxLib::Abort(m.c_str());
                  }
              }
          }
      }
  }
}

void  PorousMedia::read_source()
{
  //
  // Read in parameters for sources
  //
  ParmParse pp("source");

  // determine number of sources 
  pp.query("do_source",do_source_term);

  if (do_source_term) {
      BoxLib::Abort("Sources no longer supported");
  }

#if 0
  int nsource = pp.countval("source");
  pp.query("nsource",nsource);
  if (pp.countval("source") != nsource) 
    {
      std::cerr << "Number of sources specified and listed "
		<< "do not match.\n";
      BoxLib::Abort("read_source()");
    }
  if (do_source_term > 0 && nsource > 0)
    {
      source_array.resize(nsource);
      Array<std::string> sname(nsource);
      pp.getarr("source",sname,0,nsource);

      // Get parameters for each source
      // influence function:0=constant,1=linear,2=spherical,3=exponential
      std::string buffer;
      for (int i=0; i<nsource; i++)
	{
          const std::string prefix("source." + sname[i]);
	  ParmParse ppr(prefix.c_str());
	  source_array[i].name = sname[i];
	  ppr.get("var_type",source_array[i].var_type);
	  ppr.get("var_id",source_array[i].var_id);
	  if (source_array[i].var_type == "comp")
	    {
	      source_array[i].id.resize(1);
	      for (int j=0; j<cNames.size();j++)
		{
                    if (source_array[i].var_id == cNames[j]) {
                        source_array[i].id[0] = j;
                    }
		}
	    }
	  else if (source_array[i].var_type == "tracer")
	    {
	      if (source_array[i].var_id == "ALL")
		{
		  source_array[i].id.resize(ntracers);
		  for (int j=0;j<ntracers;j++)
		    source_array[i].id[j] = j ;
		}
	      else
		{
		  source_array[i].id.resize(1);
		  for (int j=0; j<ntracers;j++)
		    {
                      if (source_array[i].var_id == tNames[j]) {
                          source_array[i].id[0] = j;
                      }
		    }
		}
	    }

	  ppr.get("regions",buffer);
          bool region_set=false;
	  for (int j=0; j<region_array.size();j++)
          {
	      if (buffer==region_array[j]->name)
              {
                  source_array[i].region = j;
                  region_set = true;
              }
          }
          BL_ASSERT(region_set);
	  ppr.get("dist_type",buffer);
	  if (!buffer.compare("constant"))
	      source_array[i].dist_type = 0;
	  else if (!buffer.compare("linear"))
	      source_array[i].dist_type = 1;
	  else if (!buffer.compare("quadratic"))
	      source_array[i].dist_type = 2;
	  else if (!buffer.compare("exponential"))
	      source_array[i].dist_type = 3;

	  ppr.getarr("val",source_array[i].val_param,
		      0,ppr.countval("val"));
	  if (ppr.countval("val")< source_array[i].id.size())
	    std::cout << "Number of values does not match the number of var_id.\n" ;
	  if (source_array[i].dist_type != 0)
	    ppr.getarr("dist_param",source_array[i].dist_param,
			0,ppr.countval("dist_param"));
	 
	}
    }
#endif
}

void PorousMedia::read_observation()
{
  Real time_eps = 1.e-6; // FIXME: needs to be computed

  // Build time macros
  ParmParse ppa("amr");

  EventCoord& event_coord = PMAmr::eventCoord();

  int n_cmac = ppa.countval("cycle_macros");
  Array<std::string> cmacroNames;
  ppa.getarr("cycle_macros",cmacroNames,0,n_cmac);
  std::map<std::string,int> cmacro_map;
  for (int i=0; i<n_cmac; ++i) {
      std::string prefix = "amr.cycle_macro." + cmacroNames[i];
      ParmParse ppc(prefix);
      std::string type; ppc.get("type",type);
      if (type == "period") {
          int start, period, stop;
          ppc.get("start",start);
          ppc.get("period",period);
          ppc.get("stop",stop);
          defined_events[cmacroNames[i]] = new EventCoord::CycleEvent(start,period,stop);
      }
      else if (type == "cycles" ){
          Array<int> cycles; ppc.getarr("cycles",cycles,0,ppc.countval("cycles"));
          defined_events[cmacroNames[i]] = new EventCoord::CycleEvent(cycles);
      }
      else {
          BoxLib::Abort("Unrecognized cycle macros type");
      }
      cmacro_map[cmacroNames[i]] = i;
  }

  int n_tmac = ppa.countval("time_macros");
  Array<std::string> tmacroNames;
  ppa.getarr("time_macros",tmacroNames,0,n_tmac);
  std::map<std::string,int> tmacro_map;
  for (int i=0; i<n_tmac; ++i) {
      std::string prefix = "amr.time_macro." + tmacroNames[i];
      ParmParse ppt(prefix);
      std::string type; ppt.get("type",type);
      if (type == "period") {
          Real start, period, stop;
          ppt.get("start",start);
          ppt.get("period",period);
          ppt.get("stop",stop);
          defined_events[cmacroNames[i]] = new EventCoord::TimeEvent(start,period,stop);
      }
      else if (type == "times" ){
          Array<Real> times; ppt.getarr("times",times,0,ppt.countval("times"));
          defined_events[cmacroNames[i]] = new EventCoord::TimeEvent(times);
      }
      else {
          BoxLib::Abort("Unrecognized time macros type");
      }
      tmacro_map[tmacroNames[i]] = i;
  }

  ParmParse pp("observation");
  
  // determine number of observation
  int n_obs = pp.countval("observation");
  if (n_obs > 0)
  {
      observations.resize(n_obs,PArrayManage);
      Array<std::string> obs_names;
      pp.getarr("observation",obs_names,0,n_obs);

      // Get time and cycle macros

      // Get parameters for each observation
      // observation type:0=production,1=mass_fraction,2=mole_fraction,3=saturation
      for (int i=0; i<n_obs; i++)
      {
          const std::string prefix("observation." + obs_names[i]);
	  ParmParse ppr(prefix.c_str());

          std::string obs_type; ppr.get("obs_type",obs_type);
          std::string obs_field; ppr.get("field",obs_field);
          Array<std::string> region_names(1); ppr.get("region",region_names[0]);
          const PArray<Region> obs_regions = build_region_PArray(region_names);
          
          std::string obs_time_macro; ppr.get("time_macro",obs_time_macro);

          observations.set(i, new Observation(obs_names[i],obs_field,obs_regions[0],obs_type,obs_time_macro));
	}
      
      // filename for output
      pp.query("output_file",obs_outputfile);
      
    }

  ppa.queryarr("vis_cycle_macros",vis_cycle_macros,0,ppa.countval("vis_cycle_macros"));
  ppa.queryarr("vis_time_macros",vis_time_macros,0,ppa.countval("vis_time_macros"));
  ppa.queryarr("chk_cycle_macros",chk_cycle_macros,0,ppa.countval("chk_cycle_macros"));
  ppa.queryarr("chk_time_macros",chk_time_macros,0,ppa.countval("chk_time_macros"));

  std::map<std::string,EventCoord::Event*>::const_iterator eit;

  for (int i=0; i<vis_cycle_macros.size(); ++i)
  {
      eit = defined_events.find(vis_cycle_macros[i]);
      if (eit != defined_events.end()  && eit->second->IsCycle() ) {
          event_coord.Register(eit->first,eit->second);
      }
      else {
          std::string m = "vis_cycle_macros contains unrecognized macro name \"" + vis_cycle_macros[i] + "\"";
          BoxLib::Abort(m.c_str());
      }      
  }

  for (int i=0; i<vis_time_macros.size(); ++i)
  {
      eit = defined_events.find(vis_time_macros[i]);
      if (eit != defined_events.end()  && eit->second->IsTime() ) {
          event_coord.Register(eit->first,eit->second);
      }
      else {
          std::string m = "vis_time_macros contains unrecognized macro name \"" + vis_cycle_macros[i] + "\"";
          BoxLib::Abort(m.c_str());
      }      
  }
  for (int i=0; i<chk_cycle_macros.size(); ++i)
  {
      eit = defined_events.find(chk_cycle_macros[i]);
      if (eit != defined_events.end()  && eit->second->IsCycle() ) {
          event_coord.Register(eit->first,eit->second);
      }
      else {
          std::string m = "chk_cycle_macros contains unrecognized macro name \"" + chk_cycle_macros[i] + "\"";
          BoxLib::Abort(m.c_str());
      }      
  }
  for (int i=0; i<chk_time_macros.size(); ++i)
  {
      eit = defined_events.find(chk_time_macros[i]);
      if (eit != defined_events.end()  && eit->second->IsTime() ) {
          event_coord.Register(eit->first,eit->second);
      }
      else {
          std::string m = "chk_time_macros contains unrecognized macro name \"" + chk_cycle_macros[i] + "\"";
          BoxLib::Abort(m.c_str());
      }      
  }
}

void  PorousMedia::read_chem()
{

  ParmParse pp("prob");

  // get Chemistry stuff
  pp.query("do_chem",do_chem);
  pp.query("do_full_strang",do_full_strang);
  pp.query("n_chem_interval",n_chem_interval);
  if (n_chem_interval > 0) 
    {
      do_full_strang = 0;
    }
      
#ifdef AMANZI


#ifdef AMANZI
  amanzi::chemistry::SetupDefaultChemistryOutput();
  amanzi::chemistry::chem_out->AddLevel("silent");
#endif

  // get input file name, create SimpleThermoDatabase, process
  if (do_chem>0)
    {

        ParmParse ppm("mineral");
        nminerals = ppm.countval("minerals");
        minerals.resize(nminerals);
        if (nminerals>0) {
            ppm.getarr("minerals",minerals,0,nminerals);
        }

        ParmParse pps("sorption_site");
        nsorption_sites = pps.countval("sorption_sites");
        sorption_sites.resize(nsorption_sites);
        if (nsorption_sites>0) {
            pps.getarr("sorption_sites",sorption_sites,0,nsorption_sites);
        }

        ncation_exchange = 1;

        int tnum = 1;
#ifdef _OPENMP
	tnum = omp_get_max_threads();
#endif
        ParmParse pb("prob.amanzi");

        std::string verbose_chemistry_init = "silent"; pb.query("verbose_chemistry_init",verbose_chemistry_init);
        std::string fmt = "simple"; pb.query("Thermodynamic_Database_Format",fmt);
	pb.query("file", amanzi_input_file);

        const std::string& activity_model_dh = amanzi::chemistry::ActivityModelFactory::debye_huckel;
        const std::string& activity_model_ph = amanzi::chemistry::ActivityModelFactory::pitzer_hwm;
        const std::string& activity_model_u  = amanzi::chemistry::ActivityModelFactory::unit;
        std::string activity_model = activity_model_dh; pp.query("Activity_Model",activity_model);

        Real tolerance=1.5e-12; pp.query("Tolerance",tolerance);
        int max_num_Newton_iters = 150; pp.query("Maximum_Newton_Iterations",max_num_Newton_iters);
        std::string outfile=""; pp.query("Output_File_Name",outfile);
        bool use_stdout = true; pp.query("user_stdout",use_stdout);
        int num_aux = pp.countval("Auxiliary_Data");
        if (num_aux>0) {
            aux_chem_variables.resize(num_aux);
            pp.getarr("Auxiliary_Data",aux_chem_variables,0,num_aux);
        }
        //
	// In order to thread the AMANZI chemistry, we had to give each thread 
	// its own chemSolve and components object.
        //
	chemSolve.resize(tnum);
        components.resize(tnum);
	parameters.resize(tnum);


	for (int ithread = 0; ithread < tnum; ithread++)
        {
            chemSolve.set(ithread, new amanzi::chemistry::SimpleThermoDatabase());
	  
            parameters[ithread] = chemSolve[ithread].GetDefaultParameters();
            parameters[ithread].thermo_database_file = amanzi_input_file;
            parameters[ithread].activity_model_name = activity_model;
            if (ParallelDescriptor::IOProcessor() && ithread == 0) {
                BoxLib::Warning("PM_setup::read_chem: Translate ntracers, nminerals, nsorption_sites, etc into amanzi chem parlance");
            }
#if 0
            n_minerals = 0;
            n_sorbed   = 0;
            n_total    = 0;
            for (int icmp = 0; icmp < ntracers; icmp++)
	    {
                if (solid.compare(qNames[tType[icmp]]) == 0)
                    n_minerals += 1;
                else if (absorbed.compare(qNames[tType[icmp]]) == 0)
                    n_sorbed += 1;
                else
                    n_total += 1;
	    }
            BL_ASSERT(n_total == n_sorbed);
	  
            components[ithread].minerals.resize(n_minerals,0);
            components[ithread].total.resize(n_total,0);
            components[ithread].free_ion.resize(n_total,1.0e-9);
            components[ithread].total_sorbed.resize(n_sorbed,0);
            components[ithread].ion_exchange_sites.resize(0);
	  
            parameters[ithread].porosity   = 1.0; 
            parameters[ithread].saturation = 1.0;
            parameters[ithread].volume     = 1.0;
	  
            chemSolve[ithread].verbosity(amanzi::chemistry::kTerse);
	  
            chemSolve[ithread].Setup(components[ithread],parameters[ithread]);

            if (verbose_chemistry_init!="silent" && ParallelDescriptor::IOProcessor() && ithread == 0) {
                chemSolve[ithread].Display();
                chemSolve[ithread].DisplayComponents(components[ithread]);
            }
#endif
	}
    }
#endif

  pp.query("use_funccount",use_funccount);
  pp.query("max_grid_size_chem",max_grid_size_chem);
  BL_ASSERT(max_grid_size_chem > 0);
}


void PorousMedia::read_params()
{
  // problem-specific
  read_prob();

  // geometry
  read_geometry();
  if (verbose > 1 && ParallelDescriptor::IOProcessor()) 
    std::cout << "Read geometry." << std::endl;

  if (echo_inputs && ParallelDescriptor::IOProcessor()) {
      std::cout << "The Regions: " << std::endl;
      for (int i=0; i<regions.size(); ++i) {
          std::cout << regions[i] << std::endl;
      }
  }

  // chemistry
  if (verbose > 1 && ParallelDescriptor::IOProcessor()) 
    std::cout << "Read chemistry."<< std::endl;
  read_chem();

  // components and phases
  read_comp();
  if (verbose > 1 && ParallelDescriptor::IOProcessor()) 
    std::cout << "Read components."<< std::endl;
  
  // tracers
  read_tracer();
  if (verbose > 1 && ParallelDescriptor::IOProcessor()) 
    std::cout << "Read tracers."<< std::endl;

  // rock
  read_rock();
  if (verbose > 1 && ParallelDescriptor::IOProcessor()) 
    std::cout << "Read rock."<< std::endl;

  if (echo_inputs && ParallelDescriptor::IOProcessor()) {
      std::cout << "The Materials: " << std::endl;
      for (int i=0; i<rocks.size(); ++i) {
          std::cout << rocks[i] << std::endl;
      }
  }

  // source
  //if (verbose > 1 && ParallelDescriptor::IOProcessor()) 
  //  std::cout << "Read sources."<< std::endl;
  //read_source();

  // observation
  if (verbose > 1 && ParallelDescriptor::IOProcessor()) 
    std::cout << "Read observation."<< std::endl;
  read_observation();

  FORT_INITPARAMS(&ncomps,&nphases,&model,density.dataPtr(),
		  muval.dataPtr(),pType.dataPtr(),
		  &gravity);
    
  if (ntracers > 0)
    FORT_TCRPARAMS(&ntracers);

}

