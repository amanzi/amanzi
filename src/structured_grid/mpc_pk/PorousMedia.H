#ifndef _PorousMedia_H_
#define _PorousMedia_H_
//
// $Id: PorousMedia.H,v 1.97 2011-07-08 08:09:47 marc Exp $
//
#include <Geometry.H>
#include <FArrayBox.H>
#include <BoxDomain.H>
#include <LevelBld.H>
#include <BC_TYPES.H>
#include <BC_TYPES_PM.H>
#include <Utility.H>
#include <AmrLevel.H>
#include <ErrorList.H>
#include <FluxRegister.H>
#include <ViscBndry.H>
#include <Godunov.H>
#include <Diffusion.H>
#include <MacProj.H> 
#include <AuxBoundaryData.H>
#ifdef COREREACT
  #include <PChemDriver.H>
#endif
#include <Region.H>
#include <Rock.H>
#include <Source.H>
#include <Observation.H>

#ifdef AMANZI
#include "simple_thermo_database.hh"
#include "beaker.hh"
#endif

enum StateType {State_Type=0, Press_Type, Vel_Type, Vcr_Type, FuncCount_Type};

enum PressureNames { Pressure=0 };
enum VelNames { Xvel=0, Yvel, Zvel};
enum VcrNames { Xvcr=0, Yvcr, Zvcr};

class PMBld
  :
  public LevelBld
{
  virtual void variableSetUp ();
  virtual void variableCleanUp ();
  virtual AmrLevel *operator() ();
  virtual AmrLevel *operator() (Amr&            papa,
				int             lev,
				const Geometry& level_geom,
				const BoxArray& ba,
				Real            time);
};

class PorousMedia
  :
  public AmrLevel 
{
public:
  friend class MacProj;
  friend class Diffusion;
  friend class Observation;

  PorousMedia ();

  PorousMedia (Amr&            papa,
	       int             lev,
	       const Geometry& level_geom,
	       const BoxArray& bl,
	       Real            time);

  virtual ~PorousMedia ();

  virtual void restart (Amr&          papa,
			std::istream& is,
			bool          bReadSpecial = false);
  //
  // Read problem input file. GP 08/06/2007.
  //
  void probDefine(char * filename);
  //
  // Define data descriptors.
  //
  static void variableSetUp ();
  //
  // Cleanup data descriptors at end of run.
  //
  static void variableCleanUp ();
  //
  // Init grid data at problem start-up.
  //
  virtual void initData ();
  //
  // Reset time levels for the initial iterations.
  //
  virtual void resetState (Real time,
			   Real dt_old,
			   Real dt_new);
  //
  // Set time levels of state data.
  //
  virtual void setTimeLevel (Real time,
			     Real dt_old,
			     Real dt_new);
  //
  // This function fills a new level n with the best
  // level n and coarser data available (usually during regrid).
  //
  virtual void init (AmrLevel& old);
  //
  // Init data on this level after regridding if old level
  // did not exist previously.
  //
  virtual void init ();

  virtual void allocOldData ();

  virtual void removeOldData ();
  //
  // Error estimation for regridding.
  //
  virtual void errorEst (TagBoxArray& tb,
			 int          clearval,
			 int          tagval,
			 Real         time,
			 int          n_error_buf = 0,
			 int          ngrow = 0);
  //
  // Integrate derived quantities over domain.
  //
  Real sumDerive (const std::string& name,
		  Real           time);

  Real volWgtSum (const std::string& name,
		  Real           time);
  //
  // A string written as the first item in writePlotFile() at level zero.
  // It is so we can distinguish between different types of plot files.
  //
  // For PorousMedia it has the form: PorousMedia-Vnnn
  //
  virtual std::string thePlotFileType () const;

  //
  // Write plot file stuff to specified Edirectory.
  //
  virtual void writePlotFile (const std::string& dir,
			      std::ostream&  os,
			      VisMF::How     how);
  //
  // Modify list of variables to be plotted
  //
  virtual void setPlotVariables();

  //
  // Timestep estimation functions follow ...
  //
  virtual Real estTimeStep (MultiFab* u_mac);

  void predictDT(MultiFab* u_mac);

  Real initialTimeStep (MultiFab* u_mac);

  virtual void computeInitialDt (int                   finest_level,
				 int                   sub_cycle,
				 Array<int>&           n_cycle,
				 const Array<IntVect>& ref_ratio,
				 Array<Real>&          dt_level,
				 Real                  stop_time);

  virtual void computeNewDt (int                   finest_level,
			     int                   sub_cycle,
			     Array<int>&           n_cycle,
			     const Array<IntVect>& ref_ratio,
			     Array<Real>&          dt_min,
			     Array<Real>&          dt_level,
			     Real                  stop_time,
			     int                   post_regrid_flag);
  //
  // This function estimates the initial timesteping used by the model.
  //
  void post_init_estDT (Real&        dt_init,
			Array<int>&  nc_save,
			Array<Real>& dt_save,
			Real         stop_time);
  //
  // Estimate the end of the simulation for amrLevel.
  //
  virtual int okToContinue ();
  //
  // Sync state at the end of a composite timestep.
  //
  virtual void post_timestep (int iteration);
  //
  // Build any additional data structures after regrid.
  //
  virtual void post_regrid (int lbase, int new_finest);
  //
  // Build any additional data structures after restart.
  //
  virtual void post_restart ();
  //
  // Insure state, and pressure are consistent.
  //
  virtual void post_init (Real stop_time);
  //
  // Advance grids at this level in time.
  //
  virtual Real advance (Real time,
			Real dt,
			int  iteration,
			int  ncycle);

  void advance_incompressible  (Real time,
				Real dt);

  void advance_simple  (Real time,
			Real dt);

#ifdef MG_USE_FBOXLIB
  void advance_richard  (Real time,
			Real dt);
#endif

  void advance_tracer  (Real time,
			Real dt);
  //
  // Public so MacProj can access it.
  //
  virtual void SetGodunov ();
  //
  // Returns the value of "gravity" for use in the projection outflow bcs.
  //
  Real getGravity () 
  {
    return gravity;
  }
  //
  // Fill the 'dcomp' component of MultiFab 'mf' from the 'varname'
  // component of plotfile 'pltfile'.
  //
  // This code assumes that the physical domain size of the current
  // problem is the same as that of the one that generated the pltfile.
  // It also assumes that the pltfile has at least as many levels as does
  // the current problem, and that the reference ratios between levels are
  // the same for the current problem.
  //
  void fill_from_plotfile (MultiFab&          mf,
			   int                dcomp,
			   const std::string& pltfile,
			   const std::string& varname);

  static const Array<Observation>& TheObservationArray() {return observation_array;}

protected:

  DistributionMapping getFuncCountDM (const BoxArray& bxba, int ngrow = 0);
    
  void init_rock_properties (); // Initialize the rock properties
  void buildMetrics ();         // 1-D metric arrays for RZ
  static void setup_list ();
  static void read_params ();   // Read input file
  static void read_geometry ();
  static void read_rock ();
  static void read_prob ();
  static void read_comp (); 
  static void read_tracer (); 
  static void read_pressure (); 
  static void read_source (); 
  static void read_observation (); 
  static void read_chem (); 
  static void read_amr ();
  static void assign_bc_coef(int type_id, 
			     ParmParse::Record& ppr,
			     Array<std::string>& names, 
			     Array<Real>& coef);
  virtual void sum_integrated_quantities ();

  //
  // Setup, clean up for a level timestep.
  //
  virtual void advance_setup (Real time,
			      Real dt,
			      int  iteration,
			      int  ncycle);

  void advance_cleanup (Real dt,
			int  iteration,
			int  ncycle);
  //
  // Compute divergent mac velocities, estimate best timestep
  // from t^n data, and compile cfl number.
  //
  Real predict_velocity (Real  dt,
			 Real& comp_cfl);

  //
  // Create the coefficient for the MAC projection.
  //
  void create_lambda (Real time);

  //
  // Create the coefficient for the MAC projection.
  //
  void initialize_umac (MultiFab* u_mac, 
			MultiFab& RhoG, 
			MultiFab* RhoD, 
			Real      time);

  void compute_vel_phase(MultiFab* u_phase,
			 MultiFab* u_mac,
			 Real      time) ;

  void compute_vel_phase(MultiFab* u_phase,
                         int       nc,
                         Real      time) ;  
  //
  // Impose divergence constraint on MAC velocities.
  //
  void mac_project (MultiFab* u_mac, MultiFab* RhoD, Real time);


  //
  // Do an initial mac projection to get a good first dt.
  //
  void initial_mac_project (MultiFab* u_mac, MultiFab* RhoD, Real time);
  //
  // Grow by 1 and fillpatch the MAC-projected velocities: AT LEVEL 0.
  //
  void create_umac_grown (MultiFab* u_mac, MultiFab* u_macG);
  //
  // Grow by 1 and fillpatch the MAC-projected velocities: AT LEVEL > 0.
  //
  void create_umac_grown (MultiFab* u_mac, PArray<MultiFab>& umac_crse, MultiFab* u_macG) ;
  //
  // Get the crse u_mac to be used in create_umac_grown
  //
  void GetCrseUmac (PArray<MultiFab>& u_mac_crse, Real time) const;
  //
  // Get the crse pressure to define the boundary conditions for the fine solve
  //
  void GetCrsePressure(MultiFab& phi_crse,
		       Real      time      ) const;

  //
  // Advect scalars.
  //
  void scalar_advection (MultiFab* u_mac,
			 Real dt,
			 int  first_scalar,
			 int  last_scalar,
			 bool reflux_on_this_call);


  void tracer_advection (MultiFab* u_mac,
			 Real dt,
			 int  first_scalar,
			 int  last_scalar,
			 bool reflux_on_this_call);

  // Update scalars, (viscous solve in scalar_update).
  //
  void scalar_update (Real dt,
		      int  first_scalar,
		      int  last_scalar,
		      int corrector,
		      MultiFab* u_mac = 0);

  void scalar_advection_update (Real dt,
				int  first_scalar,
				int  last_scalar,
				int  corrector);
  void tracer_advection_update (Real dt,
				int  first_scalar,
				int  last_scalar,
				int  corrector);


  void scalar_diffusion_update (Real dt,
				int  first_scalar,
				int  last_scalar,
				int  corrector);

  void scalar_capillary_update (Real dt,
				int  corrector,
				MultiFab* u_mac = 0);

  void diff_capillary_update   (Real dt,
				int  corrector,
				MultiFab* u_mac = 0);

#ifdef MG_USE_FBOXLIB
  void richard_eqb_update (MultiFab* u_mac = 0);

  void richard_scalar_update (Real dt, 
			      int& total_nwt_iter,
			      MultiFab* u_mac = 0);
#endif

  void scalar_adjust_constraint (int  first_scalar,
				 int  last_scalar);

  void diffuse_adjust_dominant(MultiFab&              Phi_new,
			       int                    sigma,
			       Real                   dt,	
			       MultiFab**             fluxn,
			       MultiFab**             fluxnp1=0);


  // Reaction-related functions
  void set_preferred_boundary_values (MultiFab& S,
				      int       state_index,
				      int       src_comp,
				      int       dst_comp,
				      int       num_comp,
				      Real      time) const;
  void set_overdetermined_boundary_cells (Real time);

  //
  // ngrow == 0 -> we're working on the valid region of state.
  //
  // ngrow > 0  -> we're working on aux_boundary_data_old with that many grow cells.
  //
  void strang_chem (MultiFab&  state,
		    Real       dt,
                    int        ngrow = 0);
  //
  // Ensure state is consistent, i.e. velocity field is nondivergent,
  // coarse level data are averages of fine level data, pressure is zero.
  //
  virtual void post_init_state ();
  //
  // Interpolate cell-centered cync correction from coarse to fine.
  //
  enum SyncInterpType
    {
      PC_T,
      CellCons_T,
      CellConsLin_T,
      CellConsProt_T
    };

  void SyncInterp (MultiFab&      CrseSync,
		   int            c_lev,
		   MultiFab&      FineSync,
		   int            f_lev,
		   IntVect&       ratio,
		   int            src_comp,
		   int            dest_comp,
		   int            num_comp,
		   int            increment,
		   Real           dt_clev,
		   int**          bc_orig_qty,
		   SyncInterpType which_interp = CellCons_T,
		   int            state_comp   = -1);
  //
  // Bilinear interpolate nodal pressures from coarse to fine.
  //
  void SyncProjInterp (MultiFab& phi,
		       int       c_lev,
		       MultiFab& P_new,
		       MultiFab& P_old,
		       int       f_lev,
		       IntVect&  ratio,
		       bool      first_crse_step_after_initial_iters,
		       Real      cur_crse_pres_time,
		       Real      prev_crse_pres_time);

  //
  // Average edge from fine to coarse
  //
  void SyncEAvgDown (PArray<MultiFab> u_mac_crse,
		     PArray<MultiFab> u_mac_fine, 
		     int c_lev); 

  void SyncEAvgDown (MultiFab* u_mac,
		     int f_lev,
		     MultiFab* u_mac_crse, 
		     int c_lev); 

  void SyncEAvgDown (MultiFab* u_mac[],
		     int f_lev,
		     MultiFab* u_mac_crse[], 
		     int c_lev); 
  //
  // Average a fine multifab down onto a coarse one.
  //
  void avgDown (const BoxArray& grids,
		const BoxArray& fgrids,
		MultiFab&       S_crse,
		MultiFab&       S_fine,
		MultiFab&       volume,
		MultiFab&       fvolume,
		int             c_level,
		int             f_level,
		int             strt_comp,
		int             num_comp,
		const IntVect&  fratio);
  //
  // Average fine down to coarse in the ovlp intersection.
  //
  void avgDown (const FArrayBox& fine_fab,
		const FArrayBox& crse_fab, 
		const FArrayBox& fine_vol,
		const FArrayBox& crse_vol,
		int              f_level,
		int              c_level,
		const Box&       ovlp,
		int              start_comp,
		int              num_comp,
		const IntVect&   fratio);
  //
  // Average fine down to coarse in the ovlp intersection.
  //
  static void avgDown_doit (const FArrayBox& fine_fab,
			    const FArrayBox& crse_fab, 
			    const FArrayBox& fine_vol,
			    const FArrayBox& crse_vol,
			    int              f_level,
			    int              c_level,
			    const Box&       ovlp,
			    int              start_comp,
			    int              num_comp,
			    const IntVect&   fratio);
  //
  // Compute the mac sync correction.
  //
  virtual void mac_sync ();
#ifdef MG_USE_FBOXLIB
  virtual void richard_sync ();
#endif

  //
  // Reflux function.
  //
  virtual void reflux ();

  virtual void avgDown (); // Average down for all the state types.
  void avgDown (int comp); // Average down for a single StateType scalar

  //
  // Virtual function get fluxes out of the advection routines.
  //
  virtual void pullFluxes (int        gridno,
			   int        start_ind,
			   int        ncomp,
			   FArrayBox& xflux,
			   FArrayBox& yflux,
#if (BL_SPACEDIM == 3)
			   FArrayBox& zflux,
#endif
			   Real       dt);
  //
  // Compute viscous terms.
  //
  virtual void getViscTerms (MultiFab& visc_terms,
			     int       src_comp, 
			     int       num_comp,
			     Real      time);
  //
  // Get the forcing term.
  //
  void getForce (FArrayBox&       force,
		 int              gridno,
		 int              ngrow,
		 int              strt_comp,
		 int              num_comp,
		 const Real       time,
		 int              do_rho_scale = 0);

  void getForce_Tracer (FArrayBox&       force,
			int              gridno,
			int              ngrow,
			int              strt_comp,
			int              num_comp,
			const Real       time);
  //
  // Fill ghost cells of state.
  //
  void FillStateBndry (Real time,
		       int  state_indx,
		       int  src_comp, 
		       int  num_comp); 
  //
  // Calculate nonuniform diffusivity
  //
  void calcDiffusivity (const Real time,
			const int  src_comp = 0,
			const int  num_comp = 1);  
  //
  // Calculate diffusive term due to capillary pressure
  //
  void calcDiffusivity_CPL (MultiFab*  diffusivity[BL_SPACEDIM],
			    const Real time);
  void calcDiffusivity_CPL (MultiFab*  diffusivity[BL_SPACEDIM],
			    const MultiFab* lbd);
  //
  // Calculate diffusive term due to capillary pressure in jacobian
  //
  void calcDiffusivity_CPL_dp (MultiFab*  diffusivity[BL_SPACEDIM],
			       const MultiFab* lbd,
			       const Real time,
			       const int  num_comp = 1);
#ifdef MG_USE_FBOXLIB
  //
  // Caclulate diffusive term in the richard's equation
  //
  void calc_richard_coef (MultiFab*  diffusivity[BL_SPACEDIM],
			  const MultiFab* lbd,
			  const int  nc = 0);
  //
  // Calculate the diffusive term in the jacobian of richard's equation
  //
  void calc_richard_jac (MultiFab*  diffusivity[BL_SPACEDIM],
			 const MultiFab* lbd,
			 const Real time,
			 const int  nc = 1);
  //
  // Calculate divergence of velocity boundary condition 
  //
  void calc_richard_velbc(MultiFab& res,
			  const Real dt = 1.0);
#endif

  //
  // Calculate capillary pressure
  //
  void calcCapillary (const Real time);
  void calcCapillary (MultiFab* pc, const MultiFab& S);
  //
  // Calculate inverse capillary pressure
  //
  void calcInvCapillary (const Real time);
  //
  // Calculate lambda = relative permeability/mu
  //
  void calcLambda (const Real time, MultiFab* lbd_cc = 0);
  void calcLambda (MultiFab* lbd, const MultiFab &S);
  //
  // Calculate first derivative of lambda
  //
  void calcDLambda (const Real time,MultiFab* dlbd_cc = 0);
  
  void getDiffusivity (MultiFab*  diffusivity[BL_SPACEDIM],
			       const Real time,
			       const int  state_comp,
			       const int  dst_comp,
			       const int  num_comp);

  void center_to_edge_plain (const FArrayBox& ccfab,
			     FArrayBox&       ecfab,
			     int              sComp,
			     int              dComp,
			     int              nComp);

  void umac_edge_to_cen(MultiFab* u_mac, 
			int       idx_type);

  void umac_cpy_edge_to_cen(MultiFab* u_mac, 
			    int       idx_type,
			    int       ishift = 0);

  //
  // Boundary conditions
  //  
  void dirichletStateBC (const Real time);
  void dirichletStateBC (FArrayBox& fab, const int ngrow, const Real time);
  void dirichletTracerBC (FArrayBox& fab, const int ngrow, const Real time);
  //
  // Called in grid_places after other tagging routines to modify
  //   the list of tagged points
  //
  virtual void manual_tags_placement (TagBoxArray&    tags, 
				      Array<IntVect>& bf_lev);
  PorousMedia& getLevel (int lev)
  {
    return *(PorousMedia*) &parent->getLevel(lev);
  }

  Diffusion& getDiffusion ()
  {
    BL_ASSERT(diffusion);
    return *diffusion;
  }

  FluxRegister& getAdvFluxReg ()
  {
    BL_ASSERT(advflux_reg);
    return *advflux_reg;
  }

  FluxRegister& getAdvFluxReg (int lev)
  {
    return getLevel(lev).getAdvFluxReg();
  }

  FluxRegister& getViscFluxReg ()
  {
    BL_ASSERT(viscflux_reg);
    return *viscflux_reg;
  }

  FluxRegister& getViscFluxReg (int lev)
  {
    return getLevel(lev).getViscFluxReg();
  }

  //
  // diagnostics
  //
  void check_sum();
  void check_minmax();
  void check_minmax(int fscalar, int lscalar);
  void check_minmax(MultiFab& mf);
  void check_minmax(MultiFab* umac);
  void compute_divu (MultiFab& soln, MultiFab* umac);

  /*
    Flow model.  Options are
    0: single-phase flow with no precipitate.
    1: single-phase flow with precipitate
    2: 2-phase 2-component flow
    3: polymer flow (2-phase 3-component)
  */
  static std::map<std::string, int> model_list;
  static std::map<std::string, int> bc_list;
  static std::string model_name;
  static int model;

protected:

  // Deals with the ghost cells 
  AuxBoundaryData aux_boundary_data_old;
  bool            FillPatchedOldState_ok;


  // Region: use of pointer to allow derived class in array.
  static Array<Region*> region_array;
  static std::string surf_file;

  // Rock
  static Array<Rock> rock_array;
  static std::string gsfile;
  static MultiFab    kappadata;
  static MultiFab    phidata;
  static bool        porosity_from_fine;
  static bool        permeability_from_fine;

  // Phases and components
  static int nphases;
  static int ncomps; 
  static int ndiff;
  static int idx_dominant;
  static Array<int> pType;
  static Array<std::string> pNames;
  static Array<std::string> cNames;
  static Array<Real> density;
  static Array< Array<Real> > rhoinit;
  static Array< Array<int> > rhoinit_param;
  static Array< Array<Real> > rhoinflow;
  static Array< Array<int> > rhoinflow_param;
  static Array<Real> muval;
  
  // Tracers
  static Array<std::string> qNames;
  static Array<std::string> tNames;
  static int  ntracers; 
  static Array<int> tType; 
  static Array<Real> tDen;
  static Array<Array<Real> > tinit;
  static Array< Array<int> > tinit_param;  
  static Array< Array<Real> > tinflow;
  static Array< Array<int> > tinflow_param;
  static Array<Real> tInflow;

  // Pressure
  static Real wt_lo;
  static Real wt_hi;
  static Array<Real> press_lo;
  static Array<Real> press_hi;
  static Array<int>  inflow_bc_lo;
  static Array<int>  inflow_bc_hi;
  static Array<Real> inflow_vel_lo;
  static Array<Real> inflow_vel_hi;

  // Richard boundary conditions		
#ifdef MG_USE_FBOXLIB
  static int richard_iter;
#endif
  static Array<int>  rinflow_bc_lo;
  static Array<int>  rinflow_bc_hi;
  static Array<Real> rinflow_vel_lo;
  static Array<Real> rinflow_vel_hi;

  // temperature
  static Real temperature;

  // source term
  static bool do_source_term;
  static Array<Source> source_array;

  // observations
  static Array<Observation> observation_array;

#ifdef AMANZI
  static PArray<amanzi::chemistry::SimpleThermoDatabase> chemSolve;
  static std::string amanzi_input_file;
  static std::string amanzi_activity_model; 
  static Array<amanzi::chemistry::Beaker::BeakerComponents> components;
  static Array<amanzi::chemistry::Beaker::BeakerParameters> parameters;
  static int n_total;
  static int n_minerals;
  static int n_sorbed;
#endif
  static int  max_grid_size_chem;
  static bool no_initial_values;
  static bool use_funccount;

  // 
  // MAC edge velocities.
  //
  MultiFab* u_mac_prev;
  MultiFab* u_mac_curr;
  MultiFab* u_macG_prev;
  MultiFab* u_macG_curr;
  MultiFab* u_macG_trac;
  MultiFab* u_corr;
  //
  // Advective update terms.
  //
  MultiFab* aofs;
  //
  // Rock descriptors
  //
  MultiFab* kappa;
  MultiFab* kpedge;
  MultiFab* rock_phi;
  MultiFab* kr_coef;
  MultiFab* cpl_coef;
  //
  // lambda values and its derivative
  //
  MultiFab* lambda;       
  MultiFab* lambda_cc;    
  MultiFab* lambdap1_cc;
  MultiFab* dlambda_cc;   
  //
  // RHS for the MAC projection.
  //
  MultiFab* rhs_RhoD;
  //
  //Diffusion and capillary pressure
  //
  Diffusion* diffusion;
  MultiFab *diffn_cc, *diffnp1_cc;
  MultiFab *pcn_cc,   *pcnp1_cc;
  //
  // Volume and area fractions.
  //
  MultiFab volume;
  MultiFab area[BL_SPACEDIM];
  //
  // Scalar sync update storage
  //
  MultiFab* Ssync;    
  //
  // Data structures to store advective and viscous refluxing quantities 
  // on coarse-fine grid interfaces.
  //
  FluxRegister* advflux_reg;
  FluxRegister* viscflux_reg;
  //
  // Flag for doing first step after regrid
  //
  bool is_first_step_after_regrid;
  bool is_grid_changed_after_regrid;
  //
  // Intersection of current BoxArray with that from before regrid.
  //
  BoxArray old_intersect_new;
  //
  // Static objects.
  //
  static Godunov*    godunov;
  static ErrorList   err_list;
  static BCRec       phys_bc;
  static BCRec       pres_bc;
  static MacProj*    mac_projector;
  //
  // Internal parameters for timestepping.
  //
  static Real init_shrink;   // reduction factor of first esimated timestep
  static Real cfl;           // desired maximum cfl
  static Real change_max;    // maximum change in dt over a timestep
  static Real fixed_dt;      // set > 0 to specify dt
  static int  initial_iter;  // flag for initial pressure iterations
  static int  initial_step;  // flag for initial iterations
  static Real dt_cutoff;     // minimum dt allowed
  static int  sum_interval;  // number of timesteps for conservation stats
  static int  full_cycle;
  //
  // Internal parameters for options.
  //
  static int        verbose;
  static Real       gravity;
  static int        NUM_SCALARS;      // number of non-velocity components
  static int        NUM_STATE;        // total number of state variables
  //
  // Controls over how the convective and diffusive terms are calculated
  // for each state variable.
  //
  static Array<AdvectionForm> advectionType;
  static Array<DiffusionForm> diffusionType;
  //
  // Viscosity parameters.
  //
  static Array<int>  is_diffusive;    // does variable diffuse?
  static Array<Real> visc_coef;       // const coef viscosity terms
  static Real        visc_tol;
  static Real        visc_abs_tol;
  static Real        be_cn_theta;
  static int         variable_scal_diff;
  static int         have_capillary;
  //  
  // Chemistry parameters 
  //
  static int  do_chem;
  static int  do_full_strang;
  static int  n_chem_interval;
  static int  it_chem;
  static Real dt_chem;
  //
  // Internal switches.
  //
  static int  do_simple;
  static int  do_sync_proj;
  static int  do_reflux;
  static int  do_correct;
  static int  no_corrector;
  static int  do_kappa_refine;
  static int  n_pressure_interval;
  static int  it_pressure;
  static bool do_any_diffuse;
  static int  do_cpl_advect;

  //
  // Number of states
  //
  static int  num_state_type;
  //
  // Control velocity vs momentum update
  //
  static bool def_harm_avg_cen2edge;

  // Maximum eigenvalue at each level
  Real dt_eig;

  //
  //@ManDoc: Write current state to checkpoint file.
  //
  virtual void checkPoint (const std::string& dir,
			   std::ostream&  os,
			   VisMF::How     how,
			   bool           dump_old);
};

#endif /*_PorousMedia_H_*/
