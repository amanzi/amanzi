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
Array<Region*> PorousMedia::region_array;
std::string    PorousMedia::surf_file;
//
// Rock
//
Array<Rock> PorousMedia::rock_array;
std::string PorousMedia::gsfile;
MultiFab*   PorousMedia::kappadata;
MultiFab*   PorousMedia::phidata;
bool        PorousMedia::porosity_from_fine;
bool        PorousMedia::permeability_from_fine;
//
// Source.
//
bool          PorousMedia::do_source_term;
Array<Source> PorousMedia::source_array;
//
// Observation.
//
std::string        PorousMedia::obs_outputfile;
Array<Observation> PorousMedia::observation_array;
//
// Phases and components.
//
Array<std::string>  PorousMedia::pNames;
Array<std::string>  PorousMedia::cNames;
Array<int >         PorousMedia::pType;
Array<Real>         PorousMedia::density;
Array<BCData>       PorousMedia::ic_array;
Array<BCData>       PorousMedia::bc_array;
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
Array<BCData>       PorousMedia::tic_array;
Array<BCData>       PorousMedia::tbc_array;
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
std::map<std::string, int> PorousMedia::bc_list;
std::map<std::string, int> PorousMedia::obs_list;
std::map<std::string, int> PorousMedia::phase_list;
std::map<std::string, int> PorousMedia::comp_list;
std::map<std::string, int> PorousMedia::tracer_list;
std::map<std::string, int> PorousMedia::region_list;
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

  // bc_list
  bc_list["file"] = 0;
  bc_list["scalar"] = 1;
  bc_list["hydrostatic"] = 2;
  bc_list["rockhold"] = 3;
  bc_list["zero_total_velocity"] = 4;
  bc_list["richard"] = 5;

  // obs_list

  // default region_list
  region_list["ALL"] = 0;
  region_list["XLOBC"] = 1;
  region_list["XHIBC"] = 2;
  region_list["YLOBC"] = 3;
  region_list["YHIBC"] = 4;
#if BL_SPACEDIM == 3
  region_list["ZLOBC"] = 5;
  region_list["ZHIBC"] = 6;
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
  //       2. region_array defined in PorousMedia.H as Array<*Region>
  //
  ParmParse pp("geometry");

  // Get number of regions
  int nregion = pp.countval("region");
  pp.query("nregions",nregion);
  if (pp.countval("region") != nregion) 
    {
      std::cerr << "Number of regions specified and listed "
		<< "do not match.\n";
      BoxLib::Abort("PorousMedia::read_geometry()");
    }
  
  // set up  1+2*BL_SPACEDIM default regions
  int nRegion_DEF = 1 + 2*BL_SPACEDIM;
  nregion = nregion + nRegion_DEF;

  region_array.resize(nregion);
  int cnt=0;
  region_array[cnt++] = new allRegion();
  for (int dir = 0; dir < BL_SPACEDIM; dir++)
  {
      region_array[cnt++] = new allBCRegion(dir,0);
      region_array[cnt++] = new allBCRegion(dir,1);
  }
  Array<Real> problo, probhi;
  pp.getarr("prob_lo",problo,0,BL_SPACEDIM);
  pp.getarr("prob_hi",probhi,0,BL_SPACEDIM);
  Array<Real> temp(2*BL_SPACEDIM);
  for (int i=0;i<BL_SPACEDIM;i++)
    {
      temp[i] = problo[i];
      temp[i+BL_SPACEDIM] = probhi[i];
    }
  for (int i = 0; i < nRegion_DEF; i++)
  {
      BL_ASSERT(region_array[i]!=0);
      region_array[i]->set(temp);
  }

  // Get parameters for each user defined region 
  if (nregion > nRegion_DEF)
    {
      std::string r_purpose, r_type;
      int nRegion_Read = nregion-nRegion_DEF;
      Array<std::string> r_name(nRegion_Read);

      pp.getarr("region",r_name,0,nRegion_Read);
      
      for (int i = nRegion_DEF; i<nregion; i++)
	{
	  int idx = i-nRegion_DEF;      
	  region_list[r_name[idx]] = i;
          const std::string prefix("geometry." + r_name[idx]);
          ParmParse ppr(prefix.c_str());
	  ppr.get("purpose",r_purpose);
	  ppr.get("type",r_type);      
	  if (Region::region_map[r_type] == Region::region_map["point"])
	    {
	      Array<Real> coor;
	      ppr.getarr("coordinate",coor,0,BL_SPACEDIM);
	      region_array[i] = new pointRegion(r_name[idx],r_purpose,r_type);
	      region_array[i]->set(coor);
	    }
	  else if (Region::region_map[r_type] == Region::region_map["box"])
	    {
	      Array<Real> lo_coor,hi_coor;
	      ppr.getarr("lo_coordinate",lo_coor,0,BL_SPACEDIM);
	      ppr.getarr("hi_coordinate",hi_coor,0,BL_SPACEDIM);
	      region_array[i] = new boxRegion(r_name[idx],r_purpose,r_type);
	      
	      // check if it is at the boundary.  If yes, then include boundary.
	      Array<Real>r_param(2*BL_SPACEDIM);
	      for (int dir=0;dir<BL_SPACEDIM; dir++)
		{
		  if (lo_coor[dir] == problo[dir]) 
		    lo_coor[dir] = -2e20;
		  if (hi_coor[dir] == probhi[dir])
		    hi_coor[dir] = 2e20;
		  r_param[dir]=lo_coor[dir];
		  r_param[dir+BL_SPACEDIM] = hi_coor[dir];
		}
	      region_array[i]->set(r_param);
	    }
          else BoxLib::Abort("type not supported");
	}
      pp.query("surf_file",surf_file);
    }

}

void PorousMedia::read_rock()
{
  //
  // Get parameters related to rock
  //
  ParmParse pp("rock");
  int nrock = pp.countval("rock");
  if (nrock <= 0)
    {
      std::cerr << "At least one rock type must be defined.\n";
      BoxLib::Abort("read_rock()");
    }
  rock_array.resize(nrock);
  Array<std::string> r_names;
  pp.getarr("rock",r_names,0,nrock);
  for (int i = 0; i<nrock; i++)
    {
      const std::string prefix("rock." + r_names[i]);
      ParmParse ppr(prefix.c_str());
      rock_array[i].name = r_names[i];
      ppr.get("density",rock_array[i].density);
      ppr.get("porosity",rock_array[i].porosity);    
      ppr.getarr("permeability",rock_array[i].permeability,
		  0,ppr.countval("permeability"));
      BL_ASSERT(rock_array[i].permeability.size() == BL_SPACEDIM);
      // The permeability is specified in mDa.  
      // This needs to be multiplied with 1e-7 to be consistent 
      // with the other units in the code
      for (int j=0; j<rock_array[i].permeability.size();j++)
	rock_array[i].permeability[j] *= 1.e-7;
      // relative permeability: include kr_coef, sat_residual
      rock_array[i].krType = 0;
      ppr.query("kr_type",rock_array[i].krType);
      if (rock_array[i].krType > 0)
	ppr.getarr("kr_param",rock_array[i].krParam,
		    0,ppr.countval("kr_param"));

      // capillary pressure: include cpl_coef, sat_residual, sigma
      rock_array[i].cplType = 0;
      ppr.query("cpl_type", rock_array[i].cplType);
      if (rock_array[i].cplType > 0)
	ppr.getarr("cpl_param",rock_array[i].cplParam,
		    0,ppr.countval("cpl_param"));
      
      Array<std::string> assign_region_name;
      std::string permeability_dist, porosity_dist;
      // assigned rock to regions      
      ppr.getarr("region",assign_region_name,0,ppr.countval("region"));
      ppr.get("permeability_dist",permeability_dist);
      ppr.get("porosity_dist",porosity_dist);

      for (Array<std::string>::iterator it=assign_region_name.begin(); 
	   it!=assign_region_name.end(); it++)
	rock_array[i].region.push_back(region_list[*it]);	

      rock_array[i].porosity_dist_type = Rock::rock_dist_map[porosity_dist];
      if (rock_array[i].porosity_dist_type > 1)
	{
	  Array<Real> param;   
	  ppr.getarr("porosity_dist_param",param,0,ppr.countval("porosity_dist_param"));
	  rock_array[i].porosity_dist_param = param;
	}

      rock_array[i].permeability_dist_type = Rock::rock_dist_map[permeability_dist];
      if (rock_array[i].permeability_dist_type > 1)
	{
	  Array<Real> param;
	  ppr.getarr("permeability_dist_param",param,0,ppr.countval("permeability_dist_param"));
	  rock_array[i].permeability_dist_param = param;
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

      for (int i=0; i<nrock; i++)
	{
	  // these are temporary work around.   
	  // Should utilizes region to determine size.
	  rock_array[i].max_level = max_level;
	  rock_array[i].n_cell = n_cell;
	  rock_array[i].fratio = fratio;
	  rock_array[i].problo = problo;
	  rock_array[i].probhi = probhi;
	  rock_array[i].build_kmap(*kappadata, region_array, gsfile);
	}

      VisMF::SetNOutFiles(10);
      VisMF::Write(*kappadata,kfile);
    }

  if (build_full_pmap)
    {
      BoxArray ba = Rock::build_finest_data(max_level, n_cell, fratio);

      if (phidata == 0)
          phidata = new MultiFab;

      phidata->define(ba,1,0,Fab_allocate);
      
      for (int i=0; i<nrock; i++)
	{
	  // these are temporary work around.   
	  // Should utilizes region to determine size.
	  rock_array[i].max_level = max_level;
	  rock_array[i].n_cell = n_cell;
	  rock_array[i].fratio = fratio;
	  rock_array[i].problo = problo;
	  rock_array[i].probhi = probhi;

	  rock_array[i].build_pmap(*phidata, region_array, gsfile);
	}

      VisMF::SetNOutFiles(10);
      VisMF::Write(*phidata,pfile);
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
// assign_bc_coef assigns the type and associated coefficients of 
// initial and boundary condtiions.
//

void PorousMedia::assign_bc_coef(int type_id, 
				 ParmParse& ppr,
				 Array<std::string>& names, 
				 Array<Real>& coef)
{
  int ncoef = names.size();
  std::string rtype;
  Real vel;
  if (type_id == bc_list["scalar"])
    {
      coef.resize(ncoef);
      for (int j = 0; j<ncoef; j++)
	{
	  coef[j] = 0;
	  ppr.get(names[j].c_str(),coef[j]);
	}
    }
  else if (type_id == bc_list["hydrostatic"])
    {
      coef.resize(1);
      ppr.get("water_table",coef[0]);
    }
  else if (type_id == bc_list["zero_total_velocity"])
    {
      coef.resize(ncomps+1);
      ppr.query("rock",rtype);
      ppr.get("inflow",vel);
      for (int i=0;i<ncomps;i++)
	coef[i] = density[i];
      coef[ncomps] = vel;
      if (!rtype.empty())
	{
	  int lkrtype;
	  Real lkappa, lkrcoef, lsatres;
	  for (int i = 0; i<rock_array.size(); i++)
	    {
	      if (!rtype.compare(rock_array[i].name))
		{
		  lkappa = rock_array[i].permeability[0];
		  lkrtype = rock_array[i].krType;
		  lkrcoef = rock_array[i].krParam[0];
		  lsatres = rock_array[i].krParam[1];
		}
	    }
	  int nc = 1;
	  Real vtot = 0.;
	  Real gstar;
	  if (ncomps > 1)
	    gstar = -lkappa*(density[0]-density[1])*gravity;
	  else
	    gstar = -lkappa*(density[0])*gravity;
	      
	  Real sol;
	  FORT_FIND_INV_FLUX(&sol, &vel, &nc, &vtot,&gstar,muval.dataPtr(),&ncomps,&lkrtype,&lkrcoef);
	  coef[0] = density[0]*(sol*(1.0-lsatres)+lsatres);
	  if (ncomps > 1)
	    coef[1] = density[1]*(1.0-coef[0]/density[0]);
	}
      else if (type_id == bc_list["richard"])
	{
	  coef.resize(1);
	  ppr.getarr("inflow_bc_lo",rinflow_bc_lo,0,BL_SPACEDIM);
	  ppr.getarr("inflow_bc_hi",rinflow_bc_hi,0,BL_SPACEDIM);
	  ppr.getarr("inflow_vel_lo",rinflow_vel_lo,0,BL_SPACEDIM);
	  ppr.getarr("inflow_vel_hi",rinflow_vel_hi,0,BL_SPACEDIM);
	}
    }
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

  ParmParse cp("comp");

  // Get number of components
  ncomps = cp.countval("comp");
  cp.query("ncomp",ncomps);
  if (cp.countval("comp") != ncomps) 
    {
      std::cerr << "Number of components specified and listed "
		<< "do not match.\n";
      BoxLib::Abort("read_comp()");
    }

  // Get parameters for each component
  cp.getarr("comp",cNames,0,ncomps);
  for (int i = 0; i<ncomps; i++) comp_list[cNames[i]] = i;
  pType.resize(ncomps);
  density.resize(ncomps);
  muval.resize(ncomps);
  visc_coef.resize(ncomps);
  is_diffusive.resize(ncomps);
  std::string buffer;
  for (int i = 0; i<ncomps; i++)
    {
      const std::string prefix("comp." + cNames[i]);
      ParmParse ppr(prefix.c_str());
      ppr.get("phase",buffer);
      ppr.get("density",density[i]);
      ppr.get("viscosity",muval[i]);
      ppr.get("diffusivity",visc_coef[i]);
      pType[i] = phase_list[buffer];

      // Only components have diffusion at the moment.  
      if (visc_coef[i] > 0)
	{
	  do_any_diffuse = true;
	  is_diffusive[i]   = 1;
	}
      else
	variable_scal_diff = 0;
      ++ndiff;
    }

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
  Array<std::string> ic_array_name;
  int n_ic_array = cp.countval("init");
  if (n_ic_array <= 0)
    {
      std::cerr << "the domain must be initialized.\n";
      BoxLib::Abort("read_comp()");
    }
  cp.getarr("init",ic_array_name,0,n_ic_array);
  ic_array.resize(n_ic_array);
  std::string type_name;
  Array<std::string> region_name;
  for (int i = 0; i<n_ic_array; i++)
    {
      const std::string prefix("comp." + ic_array_name[i]);
      ParmParse ppr(prefix.c_str());
      ppr.get("type",type_name);
      ic_array[i].type = bc_list[type_name];
      int n_ic_region = ppr.countval("region");
      ppr.getarr("region",region_name,0,n_ic_region);
      ic_array[i].region.resize(n_ic_region);
      for (int j=0;j<region_name.size();j++)
	ic_array[i].region[j] = region_list[region_name[j]];
      assign_bc_coef(ic_array[i].type,ppr,cNames,ic_array[i].param);
    }

  Array<std::string> bc_array_name;
  int n_bc_array = cp.countval("inflow");
  if (n_bc_array > 0)
    {
      cp.getarr("inflow",bc_array_name,0,n_bc_array);
      bc_array.resize(n_bc_array);
      for (int i = 0; i<n_bc_array; i++)
	{
	  const std::string prefix("comp." + bc_array_name[i]);
	  ParmParse ppr(prefix.c_str());
	  ppr.get("type",type_name);
	  bc_array[i].type = bc_list[type_name];
	  int n_bc_region = ppr.countval("region");
	  ppr.getarr("region",region_name,0,n_bc_region);
	  bc_array[i].region.resize(n_bc_region);
	  for (int j=0;j<region_name.size();j++)
	    bc_array[i].region[j] = region_list[region_name[j]];
	  assign_bc_coef(bc_array[i].type,ppr,cNames,bc_array[i].param);
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
  pp.query("ntracer",ntracers);
  if (pp.countval("tracer") != ntracers) 
    {
      std::cerr << "Number of tracers specified and listed "
		<< "do not match.\n";
      BoxLib::Abort("read_tracer()");
    }

  // Get parameters for each component
  if (ntracers > 0)
    {
      // Get names of tracer groups
      pp.getarr("group",qNames,0,pp.countval("group"));
      pp.getarr("tracer",tNames,0,ntracers);
      for (int i = 0; i< ntracers; i++) tracer_list[tNames[i]] = i+ncomps;
      tType.resize(ntracers);
      std::string buffer;
      for (int i = 0; i<ntracers; i++)
	{
          const std::string prefix("tracer." + tNames[i]);
	  ParmParse ppr(prefix.c_str());
	  ppr.get("group",buffer);
	  for (int j=0;j<qNames.size(); j++) 
	    if (buffer.compare(qNames[j])==0) tType[i] = j;	   
	}
      
      // Initial condition and boundary condition  
      Array<std::string> ic_array_name;
      int n_ic_array = pp.countval("init");
      if (n_ic_array <= 0)
	{
	  std::cerr << "the domain must be initialized.\n";
	  BoxLib::Abort("read_tracer()");
	}
      pp.getarr("init",ic_array_name,0,n_ic_array);
      tic_array.resize(n_ic_array);
      std::string type_name;
      Array<std::string> region_name;
      for (int i = 0; i<n_ic_array; i++)
	{
	  const std::string prefix("tracer." + ic_array_name[i]);
	  ParmParse ppr(prefix.c_str());
	  ppr.get("type",type_name);
	  tic_array[i].type = bc_list[type_name];
	  int n_ic_region = ppr.countval("region");
	  ppr.getarr("region",region_name,0,n_ic_region);
	  tic_array[i].region.resize(n_ic_region);
	  for (int j=0;j<region_name.size();j++)
	    tic_array[i].region[j] = region_list[region_name[j]];
	  assign_bc_coef(tic_array[i].type,ppr,tNames,tic_array[i].param);
	}

      Array<std::string> bc_array_name;
      int n_bc_array = pp.countval("inflow");
      if (n_bc_array > 0)
	{
	  pp.getarr("inflow",bc_array_name,0,n_bc_array);
	  tbc_array.resize(n_bc_array);
	  for (int i = 0; i<n_bc_array; i++)
	    {
	      const std::string prefix("tracer." + bc_array_name[i]);
	      ParmParse ppr(prefix.c_str());
	      ppr.get("type",type_name);
	      tbc_array[i].type = bc_list[type_name];
	      int n_bc_region = ppr.countval("region");
	      ppr.getarr("region",region_name,0,n_bc_region);
	      tbc_array[i].region.resize(n_bc_region);
	      for (int j=0;j<region_name.size();j++)
		tbc_array[i].region[j] = region_list[region_name[j]];
	      assign_bc_coef(tbc_array[i].type,ppr,tNames,tbc_array[i].param);
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

	  ppr.get("region",buffer);
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
}

void PorousMedia::read_observation()
{
  //
  // Read in parameters for sources
  //
  ParmParse pp("observation");

  // determine number of observation
  int nobs = pp.countval("observation");
  pp.query("nobs",nobs);
  if (pp.countval("observation") != nobs) 
    {
      std::cerr << "Number of observations specified and listed "
		<< "do not match.\n";
      BoxLib::Abort("read_observation()");
    }
  if (nobs > 0)
    {
      observation_array.resize(nobs);
      Array<std::string> oname;
      pp.getarr("observation",oname,0,nobs);

      // Get parameters for each observation
      // observation type:0=production,1=mass_fraction,2=mole_fraction,3=saturation
      std::string buffer;
      for (int i=0; i<nobs; i++)
	{
          const std::string prefix("observation." + oname[i]);
	  ParmParse ppr(prefix.c_str());
	  observation_array[i].name = oname[i];
	  obs_list[oname[i]] = i;
	  ppr.get("obs_type",observation_array[i].obs_type);
	  ppr.get("var_type",observation_array[i].var_type);
	  ppr.get("var_id",observation_array[i].var_id);

	  if (!observation_array[i].var_type.compare("comp"))
	    observation_array[i].id = comp_list[observation_array[i].var_id];
	  else if (!observation_array[i].var_type.compare("tracer"))
	    observation_array[i].id = tracer_list[observation_array[i].var_id];

	  ppr.get("region",buffer);
          observation_array[i].region = region_list[buffer];
          
	  ppr.getarr("times",observation_array[i].times,
		      0,ppr.countval("times"));

	  //observation_array[i].vals.resize(observation_array[i].times.size());	 
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

  // rock
  read_rock();
  if (verbose > 1 && ParallelDescriptor::IOProcessor()) 
    std::cout << "Read rock."<< std::endl;

  // components and phases
  read_comp();
  if (verbose > 1 && ParallelDescriptor::IOProcessor()) 
    std::cout << "Read components."<< std::endl;
  
  // tracers
  read_tracer();
  if (verbose > 1 && ParallelDescriptor::IOProcessor()) 
    std::cout << "Read tracers."<< std::endl;

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
