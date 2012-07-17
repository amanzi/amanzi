#ifndef _PorousMedia_H_
#define _PorousMedia_H_

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
#include <ParmParse.H>
#include <Region.H>
#include <Rock.H>
#include <Source.H>
#include <Observation.H>
#include <RegionData.H>

#include <PMAmr.H>
class RichardSolver;
#include <RichardSolver.H>

#ifdef AMANZI
#include "simple_thermo_database.hh"
#include "beaker.hh"
#endif

enum StateType {State_Type=0, Press_Type, Vel_Type, Vcr_Type, FuncCount_Type, Aux_Chem_Type};

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

extern "C"
{
//
// Function called by BCRec for user-supplied boundary data.
//
    typedef void (*PMEF)(int* tag, ARLIM_P(tlo), ARLIM_P(thi),
                         const int* tagval, const int* clearval,
                         const Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                         const Real* mask, ARLIM_P(mask_lo), ARLIM_P(mask_hi),
                         const int* lo, const int* hi, const int* nvar,
                         const int* domain_lo, const int* domain_hi,
                         const Real* dx, const Real* xlo,
                         const Real* prob_lo, const Real* time,
                         const int* level, const Real* value);
}

class PM_Error_Value
    :
    public ErrorRec::ErrorFunc
{
public:
    PM_Error_Value()
	:
        pmef(0), value(), min_time(), max_time(), max_level() {}
    
    PM_Error_Value (Real min_time, Real max_time, int max_level, 
                    const PArray<Region>& regions);

    PM_Error_Value (PMEF pmef, Real value, Real min_time,
                    Real max_time, int max_level, const PArray<Region>& regions);

    virtual ~PM_Error_Value () {}
    
    virtual ErrorRec::ErrorFunc* clone () const
    {
        return new PM_Error_Value(pmef,value,min_time,max_time,max_level,regions);
    }

    bool regionOnly() const {return pmef==0;}
    void tagCells(int* tag, ARLIM_P(tlo), ARLIM_P(thi),
                  const int* tagval, const int* clearval,
                  const Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                  const Real* mask, ARLIM_P(mask_lo), ARLIM_P(mask_hi),
                  const int* lo, const int* hi, const int* nvar,
                  const int* domain_lo, const int* domain_hi,
                  const Real* dx, const Real* xlo,
                  const Real* prob_lo, const Real* time,
                  const int* level) const;
    int MaxLevel() const {return max_level;}
    Real MinTime() const {return min_time;}
    Real MaxTime() const {return max_time;}
    Real Value() const {return value;}
    const PArray<Region>& Regions () const {return regions;}

protected:
    void set_regions(const PArray<Region>& regions);
    PMEF pmef;
    Real value, min_time, max_time;
    int max_level;
    PArray<Region> regions;
};

struct RichardNLSdata {
  enum Reason
  {
      RICHARD_SUCCESS, 
      RICHARD_NONLINEAR_FAIL,
      RICHARD_LINEAR_FAIL
  };

  RichardNLSdata(int slev, int nlevs, PMAmr* amrp);
  void SetMaxNewtonIterations(int max_iter);
  void SetMaxJacobianReuse(int max_num_reuse);
  void SetMaxConsecutiveFails(int max_num);
  void SetDtRetryFactor(Real factor);
  void SetMaxConsecutiveFails2(int max_num);
  void SetDtRetryFactor2(Real factor);
  void SetDtRetryFactorF(Real factor);
  void SetMaxConsecutiveErrIncrease(int max_num);
  void SetConsecutiveErrIncreaseDtReduction(Real factor);

  void SetMaxNewtonIterationsForDt(int max_iter);
  void SetMaxConsecutiveSuccess(int max_num);
  void SetDtIncreaseFactor(Real factor);
  void SetMinNewtonIterationsForDt(int max_iter);
  void SetDtIncreaseFactor2(Real factor);
  void SetMinNewtonIterationsForDt2(int max_iter);
  void SetDtReductionFactor(Real factor);

  bool UpdateJacobian(int lev);
  void ResetJacobianCounter(int lev);
  void ResetJacobianCounter();
  int NLIterationsTaken() {return nl_iterations_taken;}
  int MaxNLIterations() {return max_nl_iterations;}
  void operator++() {nl_iterations_taken++;}
  void operator++(int) {nl_iterations_taken++;}
  void SetNLIterationsTaken(int iters);
  void ResetCounters();

  bool AdjustDt(Real                dt, 
		Reason              nl_solver_status, 
		Real                abs_err, 
		Real                rel_err, 
		Real&               new_dt); // Note: return bool for whether run should stop
      
  bool AdjustDt(Real                dt, 
		Reason              nl_solver_status, 
		Real&               new_dt); // Note: return bool for whether run should stop

  Array < PArray <MultiFab> > Hcoeffs;
  Array < PArray <MultiFab> > Jacobian;
  PArray <MultiFab>  DAlpha;
  Array <MultiFab*> velPhase;
  PArray <MultiFab>  initialState;

private:
  void Build();
  int start_level, end_level;
  int max_num_Jacobian_reuse;
  Array<int> num_Jacobian_reuses_remaining;
  int max_nl_iterations, nl_iterations_taken;
  Real max_nl_residual_norm, nl_residual_norm; 

  int max_num_consecutive_success, num_consecutive_success;
  int max_num_consecutive_failures_1, num_consecutive_failures_1;
  int max_num_consecutive_failures_2, num_consecutive_failures_2;
  int max_num_consecutive_increases, num_consecutive_increases;
  Real consecutive_increase_reduction_factor;
  int min_nl_iterations_for_dt, min_nl_iterations_for_dt_2, max_nl_iterations_for_dt;
  Real time_step_increase_factor, time_step_increase_factor_2, time_step_reduction_factor;
  Real time_step_retry_factor, time_step_retry_factor_2, time_step_retry_factor_f;
  bool last_chance;
  Real prev_abs_err;
  bool first;

  Array<PArray<MultiFab> > ktmp; // Currently maintains its own copy of this, should be fixed
  MFTower *Rhs, *RhoSatOld, *RhoSatNew, *Lambda, *Porosity, *Pnew;
  PArray<MFTower> Kappa, DarcyVelocity;
  Vec RhsV, SolnV;
  PMAmr* pm_amr;

  SNES           snes;
  ISColoring     iscoloring;
  MatFDColoring  matfdcoloring;

  Array<Real> gravity;
  Array<Real> density;
  MFTFillPatch* mftfp;

};

class PorousMedia
  :
  public AmrLevel 
{
public:
  friend class MacProj;
  friend class Diffusion;
  friend class Observation;

    RichardNLSdata BuildInitNLS();

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

  static void CleanupStatics ();
  //
  // Read problem input file. 
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

  void multilevel_advance (Real time,
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
 
  void advance_multilevel_richard  (Real time,
				    Real dt);
#endif
  void advance_multilevel_saturated (Real time,
				     Real dt);

  void advance_tracer  (Real time,
			Real dt);
  //
  // Public so MacProj can access it.
  //
  virtual void SetGodunov ();
  //
  // Returns the value of "gravity" for use in the projection outflow bcs.
  //
  static Real getGravity () 
  {
    return gravity;
  }

  static void coarseTimeStep (Real stop_time);


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

  static PArray<Observation>& TheObservationArray() {return observations;}

  static void getDirichletFaces (Array<Orientation>& Faces,
				 const int           comp_Type,
				 const BCRec&        _bc);

  static bool grids_on_side_of_domain (const BoxArray&    _grids,
				       const Box&         _domain,
				       const Orientation& _Face); 

  static const Array<std::string>& UserDerives() {return user_derive_list;}
  static PArray<Region>& Regions() {return regions;}
  static PArray<Rock>& Rocks() {return rocks;}
  static const Array<Real>& Density() {return density;}
  static const Array<Real>& Viscosity() {return muval;}
  static int numComponents() {return ncomps;}
  static const Array<std::string>& componentNames() {return cNames;}
  static const Array<std::string>& phaseNames() {return pNames;}
  static const Array<std::string>& soluteNames() {return tNames;}

  static PArray<RegionData>& BCs() {return bc_array;}
  static PArray<RegionData>& TBCs(int iTrac) {return tbc_array[iTrac];}

  static Array<std::string> vis_cycle_macros, vis_time_macros;
  static Array<std::string> chk_cycle_macros, chk_time_macros;
  static int echo_inputs;

  void BuildRichardNLSdata(RichardNLSdata& nld, int slev, int nlevs);

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
  static void read_source (); 
  static void read_observation (); 
  static void read_chem (); 
  static void assign_bc_coef(int type_id, 
			     ParmParse& ppr,
			     Array<std::string>& names, 
			     Array<Real>& coef);

  static PArray<Region> build_region_PArray(const Array<std::string>& region_names);
  static const Rock& find_rock(const std::string& name);

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
#ifdef BL_USE_PETSC
  RichardNLSdata::Reason richard_PETSc_pressure_solve (Real dt, RichardSolver& rs, RichardNLSdata& nl_data);
#endif

#ifdef MG_USE_FBOXLIB
  void richard_eqb_update (MultiFab* u_mac = 0);

  RichardNLSdata::Reason richard_scalar_update (Real dt, 
						int& total_nwt_iter,
						MultiFab* u_mac = 0);

  RichardNLSdata::Reason richard_composite_update (Real dt, 
						   RichardNLSdata& nl_data);

  void richard_init_to_steady();

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
  void avgDown (MultiFab* s_crse,
		int c_lev,
		MultiFab* s_fine, 
		int f_lev);

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
			  const MultiFab* u_mac,
			  const int  nc = 0,
			  const int do_upwind = 0,
			  double time = 0.);
  //
  // Calculate the diffusive term in the jacobian of richard's equation
  //
  void calc_richard_jac (MultiFab*  diffusivity[BL_SPACEDIM],
                         MultiFab*  dalpha,
			 const MultiFab* lbd,
			 const MultiFab* umac,
			 Real time,
                         Real dt,
			 int  nc,
			 int do_upwind,
			 bool do_n);

  //
  // Calculate the dt term in the jacobian for p-formulation
  //
  void calc_richard_alpha (MultiFab*     alpha,
			   const Real    time);

  void calc_richard_alpha (MultiFab*       alpha,
			   const MultiFab& S) const;
  //
  // Calculate divergence of velocity boundary condition 
  //
  void calc_richard_velbc(MultiFab& res,
			  MultiFab* uphase,
			  const Real dt = 1.0);
#endif

  //
  // Calculate capillary pressure
  //
  void calcCapillary (const Real time);
  void calcCapillary (MultiFab* pc, const MultiFab& S);
  void calcCapillary (FArrayBox&       pc,
		      const FArrayBox& S,
		      const FArrayBox& rockphi,
		      const FArrayBox& rockkappa,
		      const FArrayBox& cplcoef,
		      const Box&       grids,
		      const Array<int>& s_bc);
  void smooth_pc (MultiFab* pc);
  //
  // Calculate inverse capillary pressure
  //
  void calcInvCapillary (MultiFab& S,const MultiFab& pc);
  //
  // Calculate inverse pressure
  //
  void calcInvPressure (MultiFab& S,const MultiFab& P);
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
  void setPhysBoundaryValues (FArrayBox& dest,
                              int        state_indx,
                              Real       time,
                              int        dest_comp,
                              int        src_comp,
                              int        num_comp);
  void dirichletStateBC (const Real time);
  void dirichletStateBC (FArrayBox& fab, Real time);
  void dirichletTracerBC(FArrayBox& fab, Real time);
  void dirichletPressBC (FArrayBox& fab, Real time);
  bool get_inflow_velocity(const Orientation& face,
                           FArrayBox&         ccBndFab,
                           Real               time);
  void get_inflow_density(const Orientation& face,
			  const RegionData&  face_bc,
			  FArrayBox&         ccBndFab,
			  Real               time);
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

  void derive (const std::string& name,
               Real           time,
               MultiFab&      mf,
               int            dcomp);

  MultiFab* derive (const std::string& name,
                    Real           time,
                    int            ngrow);

  //
  // diagnostics
  //
  static void check_minmax(MultiFab& mf);
  static void check_minmax(MultiFab* umac);

  void check_sum();
  void check_minmax();
  void check_minmax(int fscalar, int lscalar);
  void compute_divu (MultiFab& soln, MultiFab* umac);

  void setup_bound_desc();

  // list
  static std::map<std::string, int> model_list;
  static std::map<std::string, int> phase_list;
  static std::map<std::string, int> comp_list;
  static std::map<std::string, int> tracer_list;

  /*
    Flow model.  Options are
    0: single-phase flow with no precipitate.
    1: single-phase flow with precipitate
    2: 2-phase 2-component flow
    3: polymer flow (2-phase 3-component)
  */
  static std::string model_name;
  static int model;

  PMAmr* PMParent();

  MultiFab* UMac_Prev() {return u_mac_prev;}
  MultiFab* UMac_Curr() {return u_mac_curr;}
  MultiFab* LambdaCC_Curr() {return lambdap1_cc;}
  MultiFab* Porosity() {return rock_phi;}
  MultiFab* KappaEC() {return kpedge;}

  void set_vel_from_bcs(Real      time,
			MultiFab* vel);

protected:

  static void InitializeStaticVariables();

  // Deals with the ghost cells 
  AuxBoundaryData aux_boundary_data_old;
  bool            FillPatchedOldState_ok;

  static Real geometry_eps;
  static PArray<Region> regions;
  static std::string surf_file;
  static bool material_is_layered;

  // Rock
  static PArray<Rock> rocks;
  static std::map<std::string, int> rock_to_int; // hack to support user output
  static std::string gsfile;
  static MultiFab*   kappadata;
  static MultiFab*   phidata;
  static bool        porosity_from_fine;
  static bool        permeability_from_fine;
  static Real        saturation_threshold_for_vg_Kr;
  static int         use_shifted_Kr_eval;


  // Phases and components
  static int nphases;
  static int ncomps; 
  static int ndiff;
  static int idx_dominant;
  static Array<int> pType;
  static Array<std::string> pNames;
  static Array<std::string> cNames;
  static Array<Real> density;
  static PArray<RegionData> ic_array;
  static PArray<RegionData> bc_array;
  static Array<Real> muval;
  static Array<std::string> user_derive_list;;
  
  // Tracers
  static Array<std::string> qNames;
  static Array<std::string> tNames;
  static int  ntracers; 
  static Array<int> tType; 
  static Array<Real> tDen;
  static Array<PArray<RegionData> > tic_array;
  static Array<PArray<RegionData> > tbc_array;
  static std::map<std::string,Array<int> > group_map;

  // Minerals, Sorption sites
  static int nminerals; 
  static Array<std::string> minerals;
  static int nsorption_sites; 
  static Array<std::string> sorption_sites;
  static int ncation_exchange; 
  static int nsorption_isotherms;
  static Array<std::string> aux_chem_variables;
  static bool using_sorption;

  typedef std::map<std::string,Real> ICParmPair;
  typedef std::map<std::string, ICParmPair > ICLabelParmPair;
  typedef std::map<std::string, ICLabelParmPair> ChemICMap;
  static ChemICMap sorption_isotherm_ics;
  static ChemICMap mineralogy_ics;
  static ChemICMap surface_complexation_ics;
  static ICLabelParmPair cation_exchange_ics;
  static Array<Real> tInflow;

  typedef std::map<std::string, std::map<std::string,int> > LabelIdx;
  static LabelIdx mineralogy_label_map;
  static LabelIdx sorption_isotherm_label_map;
  static LabelIdx surface_complexation_label_map;
  static std::map<std::string,int> cation_exchange_label_map;

  // Pressure
  static Real wt_lo;
  static Real wt_hi;
  static Array<Real> press_lo;
  static Array<Real> press_hi;
  static Array<int>  inflow_bc_lo;
  static Array<int>  inflow_bc_hi;

  // Richard boundary conditions		
#ifdef MG_USE_FBOXLIB
  static int richard_iter;
#endif
  static Array<int>  rinflow_bc_lo;
  static Array<int>  rinflow_bc_hi;

  typedef std::pair<Box,Array<int> > BCDesc;
  std::map<Orientation,BCDesc> bc_descriptor_map;
  std::map<Orientation,BCDesc> pbc_descriptor_map;
  Array<std::map<Orientation,BCDesc> > tbc_descriptor_map;

  // temperature
  static Real temperature;

  // source term
  static bool do_source_term;
  static Array<Source> source_array;

  // observations
  static std::string obs_outputfile;
  static PArray<Observation> observations;

#ifdef AMANZI
  static PArray<amanzi::chemistry::SimpleThermoDatabase> chemSolve;
  static std::string amanzi_database_file;
  static std::string amanzi_activity_model; 
  static Array<amanzi::chemistry::Beaker::BeakerComponents> components;
  static Array<amanzi::chemistry::Beaker::BeakerParameters> parameters;
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
  static Real richard_max_dt; // define maximum time step for the richard solver
  static int  initial_iter;  // flag for initial pressure iterations
  static int  initial_step;  // flag for initial iterations
  static Real dt_cutoff;     // minimum dt allowed
  static int  sum_interval;  // number of timesteps for conservation stats
  static int  full_cycle;
  static Real dt_init; // if >0, initial time step size

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
  static int  do_tracer_transport;
  static int  do_full_strang;
  static int  n_chem_interval;
  static int  it_chem;
  static Real dt_chem;
  //
  // Internal switches.
  //
  static int  do_simple;
  static int  do_sync_proj;
  static int  do_multilevel_full;
  static int  do_reflux;
  static int  do_correct;
  static int  no_corrector;
  static int  do_kappa_refine;
  static int  n_pressure_interval;
  static int  it_pressure;
  static bool do_any_diffuse;
  static int  do_cpl_advect;
  static int  max_step;
  static Real stop_time;
  static bool do_richard_sat_solve;

  static int  richard_solver_verbose;

  //
  // Init-to-steady parameters
  //
  static bool do_richard_init_to_steady;
  static int  richard_init_to_steady_verbose;
  static int  steady_min_iterations;
  static int  steady_min_iterations_2;
  static int  steady_max_iterations;
  static int  steady_limit_iterations;
  static Real steady_time_step_reduction_factor;
  static Real steady_time_step_increase_factor;
  static Real steady_time_step_increase_factor_2;
  static Real steady_time_step_retry_factor_1;
  static Real steady_time_step_retry_factor_2;
  static Real steady_time_step_retry_factor_f;
  static int  steady_max_consecutive_failures_1;
  static int  steady_max_consecutive_failures_2;
  static Real steady_tolerance;
  static Real steady_init_time_step;
  static Real steady_max_time_steps;
  static int  steady_max_num_consecutive_success;
  static Real steady_extra_time_step_increase_factor;
  static int  steady_max_num_consecutive_increases;
  static Real steady_consecutive_increase_reduction_factor;
  static Real steady_max_psuedo_time;
  //
  // Newton step defaults
  //
  static int  richard_max_ls_iterations;
  static Real richard_min_ls_factor;
  static Real richard_ls_acceptance_factor;
  static Real richard_ls_reduction_factor;
  static int richard_monitor_linear_solve;
  static int richard_monitor_line_search;

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
