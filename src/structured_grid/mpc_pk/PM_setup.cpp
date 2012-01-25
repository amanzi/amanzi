#include <winstd.H>
#include <ParmParse.H>
#include <ErrorList.H>
#include <Interpolater.H>
#include <MultiGrid.H>
#include <ArrayLim.H>
#include <Profiler.H>
#include <TagBox.H>
#include <DataServices.H>
#include <AmrData.H>
#include <time.h> 

#include <PorousMedia.H>
#include <RegType.H> 
#include <PROB_PM_F.H>
#include <DERIVE_F.H>

#ifdef BL_USE_OMP
#include "omp.h"
#endif

#ifdef AMANZI
#include "simple_thermo_database.hh"
#include "activity_model_factory.hh"
#endif

#include <TabularFunction.H>

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
// Add 2 if do_chem == 1 later.
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
//
// Phases and components.
//
Array<std::string>  PorousMedia::pNames;
Array<std::string>  PorousMedia::cNames;
Array<int >         PorousMedia::pType;
Array<Real>         PorousMedia::density;
PArray<RegionData>  PorousMedia::ics;
PArray<RegionData>  PorousMedia::bcs;
Array<Real>         PorousMedia::muval;
Array<std::string>  PorousMedia::PM_specific_derives;
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
Array<PArray<RegionData> > PorousMedia::tics;
Array<PArray<RegionData> > PorousMedia::tbcs;
//
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
Array<Real> PorousMedia::inflow_vel_lo;
Array<Real> PorousMedia::inflow_vel_hi;
Array<int>  PorousMedia::rinflow_bc_lo;
Array<int>  PorousMedia::rinflow_bc_hi;
Array<Real> PorousMedia::rinflow_vel_lo;
Array<Real> PorousMedia::rinflow_vel_hi;
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
//
// Lists.
//
std::map<std::string, int> PorousMedia::model_list;
std::map<std::string, int> PorousMedia::phase_list;
std::map<std::string, int> PorousMedia::comp_list;
std::map<std::string, int> PorousMedia::tracer_list;

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
  bc_list["richard"] = 5;
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

#ifdef MG_USE_FBOXLIB
  PorousMedia::richard_iter = 100;
#endif
  PorousMedia::wt_lo = 0;
  PorousMedia::wt_hi = 0;

  PorousMedia::temperature = 300;

  PorousMedia::verbose      = 0;
  PorousMedia::cfl          = 0.8;
  PorousMedia::init_shrink  = 1.0;
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

  PorousMedia::variable_scal_diff = 1; 

  PorousMedia::do_chem            = -1;
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
}

void
PorousMedia::variableSetUp ()
{

  InitializeStaticVariables();


  BL_ASSERT(desc_lst.size() == 0);

  for (int dir = 0; dir < BL_SPACEDIM; dir++)
  {
    phys_bc.setLo(dir,SlipWall);
    phys_bc.setHi(dir,SlipWall);
  }

  setup_list();
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
  if (do_chem>-1)
    {
      // add function count
      desc_lst.addDescriptor(FuncCount_Type, IndexType::TheCellType(),
			     StateDescriptor::Point,0,1, &cell_cons_interp);
      desc_lst.setComponent(FuncCount_Type, 0, "FuncCount", 
			    bc, BndryFunc(FORT_ONE_N_FILL));
    }
#endif
  IndexType regionIDtype(IndexType::TheCellType());
  int nCompRegion = 1;
  int cnt = 0;
  PM_specific_derives.resize(0);
  PM_specific_derives.push_back("MaterialID");
  PM_specific_derives.push_back("Capillary_Pressure");
  PM_specific_derives.push_back("Volumetric_Water_Content");
  PM_specific_derives.push_back("Porosity");
  PM_specific_derives.push_back("Aqueous_Saturation");
  PM_specific_derives.push_back("Aqueous_Pressure");
  for (int i=0; i<PM_specific_derives.size(); ++i) {
      derive_lst.add(PM_specific_derives[i], regionIDtype, nCompRegion);
  }

  //
  // **************  DEFINE ERROR ESTIMATION QUANTITIES  *************
  //
  err_list.add("gradn",1,ErrorRec::Special,FORT_ADVERROR);

  if (verbose && ParallelDescriptor::IOProcessor())
  {
      std::cout << "\nDumping ParmParse table:\n \n";
      ParmParse::dumpTable(std::cout);
      std::cout << "\n... done dumping ParmParse table.\n" << '\n';
  }
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

  // Get number of regions
  int nregion_user = pp.countval("regions");
  
  // set up  1+2*BL_SPACEDIM default regions
  int nregion_DEF = 1 + 2*BL_SPACEDIM;
  int nregion = nregion_user + nregion_DEF;
  regions.resize(nregion);
  regions.set(0, new   allRegion(problo,probhi));
  regions.set(1, new allBCRegion(0,0,problo,probhi));
  regions.set(2, new allBCRegion(0,1,problo,probhi));
  regions.set(3, new allBCRegion(1,0,problo,probhi));
  regions.set(4, new allBCRegion(1,1,problo,probhi));
#if BL_SPACEDIM == 3
  regions.set(5, new allBCRegion(2,0,problo,probhi));
  regions.set(6, new allBCRegion(2,1,problo,probhi));
#endif

  // Get parameters for each user defined region 
  if (nregion_user)
    {
      std::string r_purpose, r_type;
      Array<std::string> r_name;
      pp.getarr("regions",r_name,0,nregion_user);

      Real geometry_eps = 1.e-3;
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
	  else if (r_type == "box")
	    {
	      Array<Real> lo_coor,hi_coor;
	      ppr.getarr("lo_coordinate",lo_coor,0,BL_SPACEDIM);
	      ppr.getarr("hi_coordinate",hi_coor,0,BL_SPACEDIM);
	      
	      // check if it is at the boundary.  If yes, then include boundary.
	      for (int dir=0;dir<BL_SPACEDIM; dir++)
		{
		  if (lo_coor[dir] <= problo[dir]+geometry_eps) 
		    lo_coor[dir] = -2e20;
		  if (hi_coor[dir] >= probhi[dir]-geometry_eps)
		    hi_coor[dir] = 2e20;
		}
              regions.set(nregion_DEF+j, new boxRegion(r_name[j],r_purpose,r_type,lo_coor,hi_coor));
	    }
	  else if (r_type == "color_function")
          {
              int color_value; ppr.get("color_value",color_value);
              std::string color_file; ppr.get("color_file",color_file);
              colorFunctionRegion* cfr = new colorFunctionRegion(r_name[j],r_purpose,r_type,color_file,color_value);

	      // check if it is at the boundary.  If yes, then include boundary.
              Array<Real>& lo = cfr->lo;
              Array<Real>& hi = cfr->hi;
              
	      for (int dir=0;dir<BL_SPACEDIM; dir++)
              {
		  if (lo[dir] <= problo[dir]+geometry_eps) 
                      lo[dir] = -2e20;
		  if (hi[dir] >= probhi[dir]-geometry_eps)
                      hi[dir] = 2e20;
              }
	      regions.set(nregion_DEF+j, cfr);
          }
          else BoxLib::Abort("region type not supported");
	}
      pp.query("surf_file",surf_file);
    }
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
    rocks.resize(nrock,PArrayManage);

    for (int i = 0; i<nrock; i++)
    {
        const std::string& rname = r_names[i];
        const std::string prefix("rock." + rname);
        ParmParse ppr(prefix.c_str());
        
        Real rdensity; ppr.get("density",rdensity);
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
        am.getarr("ref_ratio",fratio,0,max_level);
        
        ParmParse gm("geometry");
        gm.getarr("prob_lo",problo,0,BL_SPACEDIM);
        gm.getarr("prob_hi",probhi,0,BL_SPACEDIM);
    }
    
    // construct permeability field based on the specified parameters
    if (build_full_kmap)
    {
        
        BoxArray ba = Rock::build_finest_data(max_level, n_cell, fratio);
        
        if (kappadata == 0)
            kappadata = new MultiFab;
        
        kappadata->define(ba,BL_SPACEDIM,0,Fab_allocate);
        
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
        
      VisMF::SetNOutFiles(10); // FIXME: Should not be hardwired here
      VisMF::Write(*kappadata,kfile);
    }
    
    if (build_full_pmap)
    {
        BoxArray ba = Rock::build_finest_data(max_level, n_cell, fratio);
        
        if (phidata == 0)
            phidata = new MultiFab;
        
        phidata->define(ba,1,0,Fab_allocate);
        
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
        
        VisMF::SetNOutFiles(10); // FIXME: Should not be hardwired here
        VisMF::Write(*phidata,pfile);
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
}

void PorousMedia::read_prob()
{
  ParmParse pp;

  max_step  = -1;    
  stop_time = -1.0;  

  pp.query("max_step",max_step);
  pp.query("stop_time",stop_time);

  //
  // Get run options.
  //  
  ParmParse pb("prob");

  // determine the model based on model_name
  pb.query("model_name",model_name);
  model = model_list[model_name];

  // if model is specified, this supersedes model_name
  pb.query("model",model);

  // Verbosity
  pb.query("v",verbose);
  
  // Get timestepping parameters.
  pb.get("cfl",cfl);
  pb.query("init_shrink",init_shrink);
  pb.query("dt_cutoff",dt_cutoff);
  pb.query("change_max",change_max);
  pb.query("fixed_dt",fixed_dt);
  pb.query("richard_max_dt",richard_max_dt);
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
    BoxLib::Abort("PorousMedia::read_params():Must have be_cn_theta <= 1.0 && >= .5");   
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
	      std::cerr << "PorousMedia::read_params: nphases != 2 && ncomps !=nphases "
			<< "although have_capillary == 1.\n ";
	      BoxLib::Abort("PorousMedia::read_params()");
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
            BoxLib::Error("Named region not found");
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
        BoxLib::Abort("Named rock not found");
    }
    return rocks[iRock];
}

void  PorousMedia::read_comp()
{
  //
  // Read in parameters for phases and components
  //
  ParmParse pp("phase");

  // Get number and names of phases
  nphases = pp.countval("phase");
  pp.getarr("phase",pNames,0,nphases);
  for (int i = 0; i<nphases; i++) phase_list[pNames[i]] = i;

  // Build flattened list of components
  int ndiff = 0;
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
      Real p_visc; ppr.get("viscosity",p_visc); muval.push_back(p_visc);
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
		  BoxLib::Abort("PorousMedia::read_params()");
		}
	      if (hi_bc[dir] != Interior)
		{
		  std::cerr << "PorousMedia::variableSetUp:periodic in direction "
			    << dir
			    << " but high BC is not Interior\n";
		  BoxLib::Abort("PorousMedia::read_params()");
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
		  BoxLib::Abort("PorousMedia::read_params()");
		}
	      if (hi_bc[dir] == Interior)
		{
		  std::cerr << "PorousMedia::variableSetUp:Interior bc in direction "
			    << dir
			    << " but not defined as periodic\n";
		  BoxLib::Abort("PorousMedia::read_params()");
		}
	    }
        }
    }

  // Initial condition and boundary condition
  //
  // Component ics, bcs will be set all at once
  int n_ics = cp.countval("inits");
  if (n_ics > 0)
  {
      Array<std::string> ic_names;
      cp.getarr("inits",ic_names,0,n_ics);
      ics.resize(n_ics);
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
          if (ic_type == "scalar")
          {
              std::string func; ppr.get("function",func);
              if (func=="constant") {
                  const std::string prefixC(prefix + ".constant");
                  Array<Real> vals(ncomps);
                  ParmParse pprc(prefixC.c_str());
                  for (int j = 0; j<cNames.size(); j++) {
                      pprc.get(cNames[j].c_str(),vals[j]);
                  }
                  ics.set(i, new RegionData(ic_regions,ic_type,vals));
              }
              else {
                  BoxLib::Abort("Unsupported comp ic");
              }
          }
          else if (ic_type == "hydrostatic")
          {
              Array<Real> water_table(1); ppr.get("water_table",water_table[0]);
              ics.set(i, new RegionData(ic_regions,ic_type,water_table));
          }
          else if (ic_type == "zero_total_velocity")
          {
              BoxLib::Abort("comp ic zero_total_velocity not implemented yet");
          }
          else if (ic_type == "richard")
          {
              BoxLib::Abort("comp ic richard not implemented yet");
#if 0
              ppr.getarr("inflow_bc_lo",rinflow_bc_lo,0,BL_SPACEDIM);
              ppr.getarr("inflow_bc_hi",rinflow_bc_hi,0,BL_SPACEDIM);
              ppr.getarr("inflow_vel_lo",rinflow_vel_lo,0,BL_SPACEDIM);
              ppr.getarr("inflow_vel_hi",rinflow_vel_hi,0,BL_SPACEDIM);
#endif
          }          
      }
  }

  int n_bcs = cp.countval("bcs");
  if (n_bcs > 0)
  {
      Array<std::string> bc_names;
      cp.getarr("bcs",bc_names,0,n_bcs);
      for (int i = 0; i<n_bcs; i++)
      {
          const std::string& bcname = bc_names[i];
	  const std::string prefix("comp.bcs" + bcname);
	  ParmParse ppr(prefix.c_str());
          
	  int n_bc_regions = ppr.countval("regions");
          Array<std::string> region_names;
	  ppr.getarr("regions",region_names,0,n_bc_regions);
          const PArray<Region> bc_regions = build_region_PArray(region_names);
          std::string bc_type; ppr.get("type",bc_type);
          
          if (bc_type == "scalar")
          {
              int nComp = cNames.size();
              Array<Array<Real> > vals(nComp), times(nComp);
              Array<Array<std::string> > forms(nComp);
              for (int j = 0; j<cNames.size(); j++)
              {
                  ParmParse pps(prefix + "." + cNames[j]);
                  int nv = pps.countval("vals");
                  if (nv) {
                      pps.getarr("vals",vals[j],0,nv);
                      pps.getarr("times",times[j],0,nv);
                      if (nv>1) {
                          pps.getarr("forms",forms[j],0,nv-1);
                      }
                  }
                  else {
                      vals[j].resize(1,0);
                      times[j].resize(1,0);
                      forms[j].resize(0);
                  }            
              }
              
              bcs.set(i, new ArrayRegionData(times,vals,forms,bc_regions,bc_type));
              
          }
          else if (bc_type == "hydrostatic")
          {
              Array<Real> water_table(1); ppr.get("water_table",water_table[0]);
              bcs.set(i, new RegionData(bc_regions,bc_type,water_table));
          }
          else if (bc_type == "zero_total_velocity")
          {
              std::string rocklabel=""; ppr.get("rock",rocklabel);
              const Rock& rock = find_rock(rocklabel);

              Array<Real> vals, times;
              Array<std::string> forms;

              int nv = ppr.countval("aqueous_vol_flux");
              if (nv) {
                  ppr.getarr("aqueous_vol_flux",vals,0,nv);
                  ppr.getarr("inflowtimes",times,0,nv);
                  if (nv>1) {
                      ppr.getarr("inflowfncs",forms,0,nv-1);
                  }
              }
              else {
                  vals.resize(1,0);
                  times.resize(1,0);
                  forms.resize(0);
              }        
              bcs.set(i, new FluxToArrayBC(times,vals,forms,bc_regions,bc_type,ncomps,rock));
          }
          else if (bc_type == "richard")
          {
              BoxLib::Abort("comp bc richard not implemented yet");
#if 0
              ppr.getarr("inflow_bc_lo",rinflow_bc_lo,0,BL_SPACEDIM);
              ppr.getarr("inflow_bc_hi",rinflow_bc_hi,0,BL_SPACEDIM);
              ppr.getarr("inflow_vel_lo",rinflow_vel_lo,0,BL_SPACEDIM);
              ppr.getarr("inflow_vel_hi",rinflow_vel_hi,0,BL_SPACEDIM);
#endif
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
  ntracers = pp.countval("tracer");
  std::map<std::string,Array<std::string> > group_map;
  if (ntracers > 0)
  {
      tics.resize(ntracers);
      pp.getarr("tracer",tNames,0,ntracers);

      for (int i = 0; i<ntracers; i++)
      {
          const std::string prefix("tracer." + tNames[i]);
	  ParmParse ppr(prefix.c_str());
          std::string g; ppr.get("group",g);
          group_map[g].push_back(tNames[i]);
      
          // Initial condition and boundary condition  
          Array<std::string> ic_names;
          int n_ic = ppr.countval("inits");
          if (n_ic <= 0)
          {
              BoxLib::Abort("each tracer must be initialized");
          }
          ppr.getarr("inits",ic_names,0,n_ic);
          
          Array<std::string> region_names;
          for (int n = 0; n<n_ic; n++)
          {
              const std::string prefixIC(prefix + "." + ic_names[n]);
              ParmParse ppri(prefix.c_str());
              
              int n_ic_region = ppri.countval("regions");
              ppri.getarr("regions",region_names,0,n_ic_region);
              const PArray<Region>& tic_regions = build_region_PArray(region_names);
              std::string tic_type; ppri.get("type",tic_type);
              
              if (tic_type == "scalar")
              {
                  Real val = 0; ppri.query("val",val);
                  tics[i].set(i, new RegionData(tic_regions,tic_type,val));
              }
              else {
                  BoxLib::Abort("Unsupported tracer IC type");
              }
          }


          Array<std::string> tbc_names;
          int n_tbc = ppr.countval("tbcs");
          if (n_tbc <= 0)
          {
              BoxLib::Abort("each tracer requires boundary conditions");
          }
          ppr.getarr("bcs",tbc_names,0,n_tbc);
          
          Array<std::string> tbc_region_names;
          for (int n = 0; n<n_tbc; n++)
          {
              const std::string prefixTBC(prefix + "." + tbc_names[n]);
              ParmParse ppri(prefix.c_str());
              
              int n_tbc_region = ppri.countval("regions");
              ppri.getarr("regions",tbc_region_names,0,n_tbc_region);
              const PArray<Region>& tbc_regions = build_region_PArray(tbc_region_names);
              std::string tbc_type; ppri.get("type",tbc_type);
              
              if (tbc_type == "scalar")
              {
                  Array<Real> times, vals;
                  Array<std::string> forms;
                  int nv = ppri.countval("vals");
                  if (nv) {
                      ppri.getarr("vals",vals,0,nv);
                      ppri.getarr("times",times,0,nv);
                      ppri.getarr("forms",forms,0,nv-1);
                  }
                  else {
                      vals.resize(1,0); // Default tracers to zero for all time
                      times.resize(1,0);
                      forms.resize(0);
                  }
                  int nComp = 1;
                  tbcs[i].set(i, new ArrayRegionData(times,vals,forms,tbc_regions,tbc_type,nComp));
              
              }
              else {
                  BoxLib::Abort("Unsupported tracer BC type");
              }
          }


      }
  }
}
  
void  PorousMedia::read_pressure()
{
  //
  // Read in parameters for pressure
  //
  ParmParse pp("press");
  press_lo.resize(BL_SPACEDIM);
  press_hi.resize(BL_SPACEDIM);
  inflow_bc_lo.resize(BL_SPACEDIM);
  inflow_bc_hi.resize(BL_SPACEDIM);
  inflow_vel_lo.resize(BL_SPACEDIM);
  inflow_vel_hi.resize(BL_SPACEDIM);
  for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
      press_lo[dir] = 0;
      press_hi[dir] = 0;
      inflow_bc_lo[dir] = 0;
      inflow_bc_hi[dir] = 0;
      inflow_vel_lo[dir] = 0;
      inflow_vel_hi[dir] = 0;
    }
  pp.query("water_table_lo",wt_lo);
  pp.query("water_table_hi",wt_hi);
  if (pp.countval("press_lo") == BL_SPACEDIM)
      pp.getarr("press_lo",press_lo,0,BL_SPACEDIM);
  if (pp.countval("press_hi") == BL_SPACEDIM)
    pp.getarr("press_hi",press_hi,0,BL_SPACEDIM);
  if (pp.countval("inflow_bc_lo") == BL_SPACEDIM)
    pp.getarr("inflow_bc_lo",inflow_bc_lo,0,BL_SPACEDIM);
  if (pp.countval("inflow_bc_hi") == BL_SPACEDIM)
    pp.getarr("inflow_bc_hi",inflow_bc_hi,0,BL_SPACEDIM);
  if (pp.countval("inflow_vel_lo") == BL_SPACEDIM) 
    pp.getarr("inflow_vel_lo",inflow_vel_lo,0,BL_SPACEDIM);
  if (pp.countval("inflow_vel_hi") == BL_SPACEDIM) 
    pp.getarr("inflow_vel_hi",inflow_vel_hi,0,BL_SPACEDIM);

#ifdef MG_USE_FBOXLIB
  if (model == model_list["richard"])
    {
      rinflow_bc_lo = inflow_bc_lo;
      rinflow_bc_hi = inflow_bc_hi;
      rinflow_vel_lo = inflow_vel_lo;
      rinflow_vel_hi = inflow_vel_hi;
    }
#endif

  Array<int> plo_bc(BL_SPACEDIM), phi_bc(BL_SPACEDIM);
  pp.getarr("lo_bc",plo_bc,0,BL_SPACEDIM);
  pp.getarr("hi_bc",phi_bc,0,BL_SPACEDIM);
  for (int i = 0; i < BL_SPACEDIM; i++)
    {
      pres_bc.setLo(i,plo_bc[i]);
      pres_bc.setHi(i,phi_bc[i]);
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
  //
  // Read in parameters for sources
  //
  ParmParse pp("observation");

  // determine number of observation
  int n_obs = pp.countval("observation");
  if (n_obs > 0)
  {
      observations.resize(n_obs);
      Array<std::string> obs_names;
      pp.getarr("observation",obs_names,0,n_obs);

      // Get parameters for each observation
      // observation type:0=production,1=mass_fraction,2=mole_fraction,3=saturation
      for (int i=0; i<n_obs; i++)
      {
          const std::string prefix("observation." + obs_names[i]);
	  ParmParse ppr(prefix.c_str());

          std::string obs_type; ppr.get("obs_type",obs_type);
          std::string obs_op_type; ppr.get("obs_op_type",obs_op_type);
          std::string obs_field; ppr.get("field",obs_field);
          Array<std::string> region_names(1); ppr.get("region",region_names[0]);
          const PArray<Region>& obs_regions = build_region_PArray(region_names);
          
	  Array<Real> obs_times; ppr.getarr("times",obs_times,0,ppr.countval("times"));

          observations.set(i, new Observation(obs_names[i],obs_field,obs_regions[0],obs_type,obs_op_type,obs_times));
	}
      
      // filename for output
      pp.query("output_file",obs_outputfile);
      
    }
}

void  PorousMedia::read_chem()
{

  ParmParse pp("chem");

  // get Chemistry stuff
  pp.query("do_chem",do_chem);
  pp.query("do_full_strang",do_full_strang);
  pp.query("n_chem_interval",n_chem_interval);
  if (n_chem_interval > 0) 
    {
      do_full_strang = 0;
    }
      
#ifdef AMANZI
  // get input file name, create SimpleThermoDatabase, process
  if (do_chem > -1)
    {
int tnum = 1;

#ifdef BL_USE_OMP
	tnum = omp_get_max_threads();
#endif
        ParmParse pb("prob");
	pb.query("amanzi.file", amanzi_input_file);
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
            parameters[ithread].activity_model_name = amanzi::chemistry::ActivityModelFactory::debye_huckel;    
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

  if (ParallelDescriptor::IOProcessor()) {
      std::cout << "The Regions: " << std::endl;
      for (int i=0; i<regions.size(); ++i) {
          std::cout << regions[i] << std::endl;
      }
  }

  // rock
  read_rock();
  if (verbose > 1 && ParallelDescriptor::IOProcessor()) 
    std::cout << "Read rock."<< std::endl;

  // components and phases
  read_comp();
  if (verbose > 1 && ParallelDescriptor::IOProcessor()) 
    std::cout << "Read components."<< std::endl;
  
  if (ParallelDescriptor::IOProcessor()) {
      std::cout << "The Components: " << std::endl;
      for (int i=0; i<regions.size(); ++i) {
          std::cout << regions[i] << std::endl;
      }
  }

  // tracers
  read_tracer();
  if (verbose > 1 && ParallelDescriptor::IOProcessor()) 
    std::cout << "Read tracers."<< std::endl;

  std::cout << "read_tracer DONE." << std::endl;  
  BoxLib::Abort("MARC");

  // pressure
  read_pressure();
  if (verbose > 1 && ParallelDescriptor::IOProcessor()) 
    std::cout << "Read sources."<< std::endl;

  // source
  read_source();

  // chemistry
  read_chem();

  // read amr
  //read_amr();

  // source
  read_observation();

  FORT_INITPARAMS(&ncomps,&nphases,&model,density.dataPtr(),
		  muval.dataPtr(),pType.dataPtr(),
		  &gravity);
    
  if (ntracers > 0)
    FORT_TCRPARAMS(&ntracers);

}
