/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

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
#include <TensorDiffusionBndry.H>
#include <AuxBoundaryData.H>
#include <ParmParse.H>
#include <Region.H>
#include <Material.H>
#include <Source.H>
#include <Observation.H>
#include <RegionData.H>
#include <DataServices.H>

#include <PMAmr.H>
#include <RichardSolver.H>
#include <RockManager.H>
#include <IdxRegionData.H>

#ifdef ALQUIMIA_ENABLED
#include <ChemistryEngine.hh>
#endif
#include <ChemistryHelper_Structured.H>

#include "Teuchos_ParameterList.hpp"

static Real BL_ONEATM = 101325.0;

enum StateType {State_Type=0, Press_Type, Vcr_Type, FuncCount_Type, Aux_Chem_Type, Num_State_Type};

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
                    const Array<const Region*>& regions);

    PM_Error_Value (PMEF pmef, Real value, Real min_time,
                    Real max_time, int max_level, const Array<const Region*>& regions);

    virtual ~PM_Error_Value () {}

    virtual ErrorRec::ErrorFunc* clone () const {
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
    const Array<const Region*>& Regions () const {return regions;}

protected:
    void set_regions(const Array<const Region*>& regions);
    PMEF pmef;
    Real value, min_time, max_time;
    int max_level;
    Array<const Region*> regions;
};

class PorousMedia
  :
  public AmrLevel
{
public:
  friend class Observation;
  friend class PMFillPatchIterator;

  enum ExecutionMode {INIT_TO_STEADY, STEADY, TRANSIENT, INVALID};

  PorousMedia ();

  PorousMedia (Amr&            papa,
	       int             lev,
	       const Geometry& level_geom,
	       const BoxArray& bl,
	       Real            time);

  virtual ~PorousMedia ();

  static const Teuchos::ParameterList&InputParameterList() {return input_parameter_list;}
  static void SetInputParameterList(const Teuchos::ParameterList& pl) {input_parameter_list=pl;}

  virtual void restart (Amr&          papa,
			std::istream& is,
			bool          bReadSpecial = false);

  static void CleanupStatics ();
  void RegisterPhysicsBasedEvents();
  void FinalizeRockManager(bool restart);
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

  virtual void removeOldData () {}
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

  void predictDT(MultiFab* u_mac, Real t_eval);

  Real predictDT_diffusion_explicit (Real t_eval,
                                     MultiFab* saturation = 0);

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
			int  ncycle) { return 0.0; }

  bool multilevel_advance (Real  time,
			   Real  dt,
			   int   iteration,
			   int   ncycle,
                           Real& dt_new);

  bool advance_multilevel_richards_flow (Real  time,
                                         Real  dt,
                                         Real& dt_new);

  bool ml_step_driver(Real  time,
		      int   amr_iteration,
		      int   amr_ncycle,
                      Real  dt_try,
                      Real& dt_taken,
                      Real& dt_suggest,
		      bool  attempt_to_recover_failed_step);

  void advance_richards_transport_dt(Real      t,
                                     MultiFab* saturation = 0);

  void advance_flow_nochange(Real time, Real dt);

  void advance_saturated_transport_dt(Real time);

  bool advance_richards_transport_chemistry(Real  t,
					    Real  dt,
					    int   iteration,
					    Real& dt_new,
					    bool  do_subcycle,
					    bool  do_recursive,
                                            bool  use_cached_sat);

  void set_saturated_velocity();

  MultiFab* AllocateUMacG() const;

  //
  // Returns the value of "gravity" for use in the projection outflow bcs.
  //
  static Real getGravity ()
  {
    return gravity;
  }
  static int getGravityDir ()
  {
    return gravity_dir;
  }

  //
  // Fill the 'dcomp' component of MultiFab 'mf' from the 'varname'
  // component of plotfile 'pltfile'.
  //
  // This code assumes that the physical domain size of the current
  // problem is the same as that of the one that generated the pltfile.
  // It also assumes that the pltfile has at least as many levels as does
  // the current problem, and that the refinement ratios between levels are
  // the same for the current problem.
  //
  void fill_from_plotfile (MultiFab&          mf,
			   int                dcomp,
			   const std::string& pltfile,
			   const std::string& varname);

  static void getDirichletFaces (Array<Orientation>& Faces,
				 const int           comp_Type,
				 const BCRec&        _bc);

  static bool grids_on_side_of_domain (const BoxArray&    _grids,
				       const Box&         _domain,
				       const Orientation& _Face);

  static const Array<std::string>& UserDerives() {return user_derive_list;}

  static Real& Temperature() {return temperature;}
  static Array<Real>& Density() {return density;}
  static Array<Real>& Viscosity() {return muval;}
  static int numComponents() {return ncomps;}
  static const Array<std::string>& componentNames() {return cNames;}
  static const Array<std::string>& phaseNames() {return pNames;}
  static const Array<std::string>& soluteNames() {return tNames;}
  static bool DoSoluteChemistry() {return do_tracer_chemistry;}

  static PArray<RegionData>& BCs() {return bc_array;}
  static PArray<IdxRegionData>& TBCs(int iTrac) {return tbc_array[iTrac];}

  static int echo_inputs;

  void BuildNLScontrolData(NLScontrol&        nlsc,
                           RSdata&            rs_data,
                           const std::string& IDstring);

  DistributionMapping getFuncCountDM (const BoxArray& bxba, int ngrow = 0);

  void init_rock_properties (bool restart); // Initialize the rock properties
  void buildMetrics ();         // 1-D metric arrays for RZ
  static void read_params ();   // Read input file
  static void read_geometry ();
  static void read_prob ();
  static void read_comp ();
  static void read_tracer ();
  static void read_source ();
  static void read_observation ();

  virtual void sum_integrated_quantities ();

  //
  // Setup, clean up for a level timestep.
  //
  virtual void advance_setup (Real time,
			      Real dt,
			      int  iteration,
			      int  ncycle) {}
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

  void tracer_advection (MultiFab* u_mac,
			 bool reflux_on_this_call,
                         bool use_cached_sat,
                         MultiFab* F = 0);

  void tracer_diffusion (bool reflux_on_this_call,
                         bool use_cached_sat,
                         const MultiFab& F);

  void cache_component_saturations(int nGrow);
  void reinstate_component_saturations();
  void get_fillpatched_rhosat(Real t_eval, MultiFab& RhoSat, int nGrow);

  Real solve_for_initial_flow_field();

  bool advance_chemistry (Real time,
			  Real dt,
			  int  ngrow=0);
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

  //
  // Reflux function.
  //
  virtual void reflux ();
  void reflux(FluxRegister& fr, int sComp, int nComp);

  virtual void avgDown (); // Average down for all the state types.
  void avgDown (int comp); // Average down for a single StateType scalar

  //
  // Virtual function get fluxes out of the advection routines.
  //
  /*
  virtual void pullFluxes (int        gridno,
			   int        start_ind,
			   int        ncomp,
			   FArrayBox& xflux,
			   FArrayBox& yflux,
#if (BL_SPACEDIM == 3)
			   FArrayBox& zflux,
#endif
			   Real       dt);
  */
  //
  // Compute viscous terms.
  //
  virtual void getViscTerms (MultiFab& visc_terms,
			     int       src_comp,
			     int       num_comp,
			     Real      time) {}
  //
  // Get the forcing term.
  //
  void getForce (MultiFab& force,
		 int       nGrow,
		 int       strt_comp,
		 int       num_comp,
		 Real      time,
		 Real      dt,
		 bool      do_rho_scale = false);
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
			const int  num_comp = 1,
                        MultiFab*  saturation = 0);
  //
  // Calculate the dt term in the jacobian for p-formulation
  //  (alpha = dN /dPc)
  //
  void calc_richard_alpha (MultiFab&       alpha,
                           const MultiFab& N,
                           Real            time,
                           int             sComp,
                           int             dComp,
                           int             nGrow) const;
  //
  // Calculate divergence of velocity boundary condition
  //
  void calc_richard_velbc(MultiFab& res,
			  MultiFab* uphase,
			  const Real dt = 1.0);
  //
  // Calculate capillary pressure
  //
  void calcCapillary (MultiFab&       pc,
                      const MultiFab& N,
                      Real            time,
                      int             sComp,
                      int             dComp,
                      int             nGrow) const;
  //
  // Calculate inverse capillary pressure
  //
  void calcInvCapillary (MultiFab&       N,
                         const MultiFab& pc,
                         Real            time,
                         int             sComp,
                         int             dComp,
                         int             nGrow) const;
  //
  // Calculate inverse pressure
  //
  void calcInvPressure (MultiFab&       N,
                        const MultiFab& P,
                        Real            time,
                        int             sComp,
                        int             dComp,
                        int             nGrow) const;
  //
  // Calculate lambda = relative permeability/mu
  //
  void calcLambda (MultiFab&       lambda,
                   const MultiFab& N,
                   Real            time,
                   int             sComp,
                   int             dComp,
                   int             nGrow) const;
  void calcLambda (const Real time);

  //
  // Calculate first derivative of lambda
  //
  void calcDLambda (const Real time,MultiFab* dlbd_cc = 0);

  void getDiffusivity (MultiFab*  diffusivity[BL_SPACEDIM],
                       const Real time,
                       const int  state_comp,
                       const int  dst_comp,
                       const int  num_comp);

  void getTensorDiffusivity (MultiFab*  diagonal_diffusivity[BL_SPACEDIM],
                             MultiFab*  off_diagonal_diffusivity[BL_SPACEDIM],
                             const Real time);

  void center_to_edge_plain (const FArrayBox& ccfab,
			     FArrayBox&       ecfab,
			     int              sComp,
			     int              dComp,
			     int              nComp);

  void umac_edge_to_cen(MultiFab* u_mac,
			MultiFab& U_cc,
                        bool      do_upwind = false);

  void umac_cpy_edge_to_cen(MultiFab* u_mac,
			    int       idx_type,
			    int       ishift);
  //
  // Boundary conditions
  //
  void setPhysBoundaryValues (FArrayBox& dest,
                              int        state_indx,
                              Real       time,
                              int        dest_comp,
                              int        src_comp,
                              int        num_comp);
  void PMsetPhysBoundaryValues (FArrayBox&       dest,
				const IArrayBox& matID,
				const FArrayBox& aux,
				int              state_indx,
				Real             time,
				int              dest_comp,
				int              src_comp,
				int              num_comp);
  void dirichletStateBC (FArrayBox& fab, const IArrayBox& matID, const FArrayBox& aux, Real time, int sComp, int dComp, int nComp);
  void dirichletTracerBC(FArrayBox& fab, const IArrayBox& matID, const FArrayBox& aux, Real time, int sComp, int dComp, int nComp);
  void dirichletPressBC (FArrayBox& fab, const IArrayBox& matID, const FArrayBox& aux, Real time);
  void dirichletDefaultBC (FArrayBox& fab, const IArrayBox& matID, Real time);
  bool get_inflow_velocity(const Orientation& face,
                           FArrayBox&         ccBndFab,
                           FArrayBox&         mask,
                           Real               t);

  Real AdjustBCevalTime(int  state_idx,
                        Real time,
                        bool tadj_verbose);

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
  static int SumInterval() {return sum_interval;}

  void check_sum();
  void check_minmax();
  void check_minmax(int fscalar, int lscalar);

  // list
  static std::map<std::string, int> phase_list;
  static std::map<std::string, int> comp_list;
  static std::map<std::string, int> tracer_list;

  enum FLOW_MODEL_ID
  {
    PM_FLOW_MODEL_INVALID,
    PM_FLOW_MODEL_SATURATED,
    PM_FLOW_MODEL_RICHARDS,
    PM_FLOW_MODEL_OFF
  };
  static FLOW_MODEL_ID flow_model;
  static int FlowModel() {return (int)flow_model;}

  enum FLOW_EVAL_ID
  {
    PM_FLOW_EVAL_INVALID,
    PM_FLOW_EVAL_CONSTANT,
    PM_FLOW_EVAL_COMPUTE_GIVEN_P,
    PM_FLOW_EVAL_SOLVE_GIVEN_PBC,
    PM_FLOW_EVAL_EVOLVE,
    PM_FLOW_EVAL_OFF
  };
  static FLOW_EVAL_ID flow_eval;
  static bool flow_is_static;
  static int FlowEval() {return (int)flow_eval;}

  PMAmr* PMParent();
  const PMAmr* PMParent() const;

#ifdef ALQUIMIA_ENABLED
  static Amanzi::AmanziChemistry::ChemistryEngine* GetChemistryEngine() {return chemistry_engine;}
#endif
  static ChemistryHelper_Structured* GetChemistryHelper() {return chemistry_helper;}
  static const std::string& GetChemistryModel() {return chemistry_model_name;}
  static const std::string& GetChemistryEngineName() {return chemistry_engine_name;}
  static RockManager* GetRockManager() {return rock_manager;}
  static RegionManager* GetRegionManager() {return region_manager;}


  MultiFab* UMac_Prev() {return u_mac_prev;}
  MultiFab* UMac_Curr() {return u_mac_curr;}
  MultiFab* LambdaCC_Curr() {return lambdap1_cc;}
  MultiFab* KappaCCavg() {return kappa;};
  MultiFab* Porosity() {return rock_phi;}
  MultiFab* SpecificStorage() {return specific_storage;}
  MultiFab* SpecificYield() {return specific_yield;}
  MultiFab* ParticleDensity() {return particle_density;}
  MultiFab* KappaEC() {return kpedge;}
  MultiFab* PCapParams() {return cpl_coef;}
  MultiFab* KrParams() {return kr_coef;}
  MultiFab* Source() {return source;}
  Array<int>& rinflowBCLo() {return rinflow_bc_lo;}
  Array<int>& rinflowBCHi() {return rinflow_bc_hi;}

  void set_vel_from_bcs(Real      time,
			MultiFab* vel);

  static int NGrowHYP() {return nGrowHYP;}

  void derive_Material_ID(Real      time,
                          MultiFab& mf,
                          int       dcomp);
  void derive_Grid_ID(Real      time,
                      MultiFab& mf,
                      int       dcomp);

  void derive_Core_ID(Real      time,
                      MultiFab& mf,
                      int       dcomp);

  void derive_Cell_ID(Real      time,
                      MultiFab& mf,
                      int       dcomp);

  void derive_Volumetric_Water_Content(Real      time,
                                       MultiFab& mf,
                                       int       dcomp,
                                       int       ntrac = -1);

  void derive_Aqueous_Saturation(Real      time,
                                 MultiFab& mf,
                                 int       dcomp);

  void derive_Aqueous_Pressure(Real      time,
                               MultiFab& mf,
                               int       dcomp);

  void derive_Capillary_Pressure(Real      time,
				 MultiFab& mf,
				 int       dcomp);

  void derive_Hydraulic_Head(Real      time,
                             MultiFab& mf,
                             int       dcomp);

  void derive_Aqueous_Volumetric_Flux(Real      time,
                                      MultiFab& mf,
                                      int       dcomp,
                                      int       dir);

  void derive_Porosity(Real      time,
                       MultiFab& mf,
                       int       dcomp);

  void derive_Intrinsic_Permeability(Real      time,
                                     MultiFab& mf,
                                     int       dcomp,
                                     int       dir);

  void derive_Tortuosity(Real      time,
                         MultiFab& mf,
                         int       dcomp,
                         int       dir);

  void derive_SpecificStorage(Real      time,
                              MultiFab& mf,
                              int       dcomp);

  void derive_SpecificYield(Real      time,
                            MultiFab& mf,
                            int       dcomp);

  void derive_ParticleDensity(Real      time,
                              MultiFab& mf,
                              int       dcomp);

  void derive_CationExchangeCapacity(Real      time,
                                     MultiFab& mf,
                                     int       dcomp);

  void derive_Dispersivity(Real      time,
                           MultiFab& mf,
                           int       dcomp,
                           int       dir);

  void derive_WriteRegion(Real                      time,
			  MultiFab&                 mf,
			  int                       dcomp,
			  const Array<std::string>& region_set);

  void BuildNLScontrol(NLScontrol&        nlsc,
		       const std::string& IDstring);

  const iMultiFab& MaterialID() const {return *materialID;}

protected:

  static void InitializeStaticVariables();

  void setup_bound_desc();
  void ComputeSourceVolume ();
  MultiFab& SourceCellVolume() {return source_cell_volume;}

  static RichardSolver* richard_solver;
  static NLScontrol* richard_solver_control;
  static RSdata* richard_solver_data;

  static Real geometry_eps;
  static std::string surf_file;

  // Rock
  static std::map<std::string, int> rock_to_int; // hack to support user output

  // Phases and components
  static int nphases;
  static int ncomps;
  static int ndiff;
  static Array<int> pType;
  static Array<std::string> pNames;
  static Array<std::string> cNames;
  static Array<Real> density;
  static PArray<RegionData> ic_array;
  static PArray<RegionData> bc_array;
  static PArray<RegionData> source_array;
  static Array<Real> muval;
  static Array<std::string> user_derive_list;
  static Array<Real> source_volume;

  // Tracers
  static Array<std::string> qNames;
  static Array<std::string> tNames;
  static int  ntracers;
  static Array<int> tType;
  static Array<Real> tDen;
  static Array<PArray<IdxRegionData> > tic_array;
  static Array<PArray<IdxRegionData> > tbc_array;
  static Array<PArray<RegionData> > tsource_array;
  static std::map<std::string,Array<int> > group_map;

  // Minerals, Sorption sites
  static double uninitialized_data;
  static int nminerals;
  static Array<std::string> minerals;
  static int nsorption_sites;
  static Array<std::string> sorption_sites;
  static int ncation_exchange;
  static int nsorption_isotherms;
  static bool using_sorption;

  static Array<Real> tInflow;

  static RockManager::ChemICMap solute_chem_ics; // sc[rockname][solute][property] = val
  static RockManager::ICLabelParmPair sorption_chem_ics; // sc[solute][property] = val

  // Pressure
  static Real wt_lo;
  static Real wt_hi;
  static Array<Real> press_lo;
  static Array<Real> press_hi;
  static Array<int>  inflow_bc_lo;
  static Array<int>  inflow_bc_hi;

  // Richard boundary conditions
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

  // observation processing
  static int verbose_observation_processing;

  static std::map<std::string,Array<std::string> > write_region_sets;

#ifdef AMANZI

#ifdef ALQUIMIA_ENABLED
  static Amanzi::AmanziChemistry::ChemistryEngine* chemistry_engine;
#endif
  static std::string chemistry_model_name;
  static std::string chemistry_engine_name;
  static ChemistryHelper_Structured* chemistry_helper;
  static std::string amanzi_database_file;
  static std::string amanzi_activity_model;

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
  MultiFab* specific_storage;
  MultiFab* specific_yield;
  MultiFab* particle_density;
  MultiFab* kr_coef;
  MultiFab* cpl_coef;
  iMultiFab* materialID;
  //
  // lambda values and its derivative
  //
  MultiFab* lambda;
  MultiFab* lambda_cc;
  MultiFab* lambdap1_cc;
  MultiFab* dlambda_cc;
  //
  // component sources
  //
  MultiFab* source;
  //
  //Diffusion and capillary pressure
  //
  MultiFab *diffn_cc, *diffnp1_cc;
  //
  // Volume and area fractions.
  //
  MultiFab volume;
  MultiFab area[BL_SPACEDIM];
  MultiFab source_cell_volume;
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
  bool is_first_flow_step_after_regrid;
  bool is_grid_changed_after_regrid;
  //
  // Intersection of current BoxArray with that from before regrid.
  //
  BoxArray old_intersect_new;
  //
  // Static objects.
  //
  static ErrorList   err_list;
  static BCRec       phys_bc;
  static BCRec       pres_bc;
  static RegionManager* region_manager;
  static RockManager* rock_manager;
  //
  // Internal parameters for timestepping.
  //
  static Real init_shrink;   // reduction factor of first esimated timestep
  static Real cfl;           // desired maximum cfl
  static Real dt_grow_max;    // maximum change in dt over a timestep
  static Real dt_shrink_max;    // maximum change in dt over a timestep
  static Real fixed_dt;      // set > 0 to specify dt
  static Real steady_max_dt; // define maximum time step taken in steady solver mode
  static Real transient_max_dt; // define maximum time step taken in the transient solver mode
  static int  initial_iter;  // flag for initial pressure iterations
  static int  initial_step;  // flag for initial iterations
  static Real dt_cutoff;     // minimum dt allowed
  static int  sum_interval;  // number of timesteps for conservation stats
  static Real dt_init; // if >0, initial time step size
  static int max_n_subcycle_transport; // maximum number of subcycled transport steps
  static int max_dt_iters_flow; // maximum number of tries to reduce dt and retry to solve Richards flow
  static bool abort_on_chem_fail; // do not recover from chemistry failure
  static bool show_selected_runtimes; // dump wall-clock times for selected sections
  static Real max_chemistry_time_step; // Largest dt for chem, triggers subcycling otherwise (-1 means not active)
  //
  // Internal parameters for options.
  //
  static int  verbose;
  static Real gravity;
  static int  gravity_dir;
  static int  NUM_SCALARS;      // number of non-velocity components
  static int  NUM_STATE;        // total number of state variables
  //
  // Viscosity parameters.
  //
  static Array<int>  is_diffusive;    // does variable diffuse?
  static Array<Real> visc_coef;       // const coef viscosity terms
  static Array<Real> molecular_diffusivity; // for tracers
  static Real        visc_tol;
  static Real        visc_abs_tol;
  static Real        be_cn_theta;
  static bool        variable_scal_diff;
  static Real        atmospheric_pressure_atm;
  static std::map<std::string,bool> use_gauge_pressure; // On bc_tag string
  static Real        z_location;
  static Real        domain_thickness;
  //
  // Chemistry/solute/tracer parameters
  //
  static bool do_tracer_chemistry; // Does the user want to react tracers
  static bool do_tracer_advection; // Does the user want to advect tracers
  static bool do_tracer_diffusion; // Does the user want to diffuse tracers
  static bool setup_tracer_transport; // Should the memory for tracers be allocated/initialized
  static bool advect_tracers; // Should the tracers be transported now
  static bool react_tracers; // Should the tracers be reacted now
  static bool diffuse_tracers; // Should the tracers be diffused now
  static bool tensor_tracer_diffusion; // Use tensor operator (includes dispersion) rather than diffusion-only ABec
  static bool do_full_strang;
  static int  n_chem_interval; // Number of transport steps taken per chemistry solve.  Unused if <= 0
  static int  it_chem;
  static Real dt_chem;
  static Real be_cn_theta_trac;
  static Array<Real> first_order_decay_constant;
  //
  // Internal switches.
  //
  static bool do_reflux;
  static ExecutionMode execution_mode;
  static Real ic_chem_relax_dt;
  static bool solute_transport_limits_dt;
  static int nGrowHYP;
  static int nGrowMG;
  static int nGrowEIGEST;
  static bool do_constant_vel;
  static bool do_output_flow_time_in_years;
  static bool do_output_chemistry_time_in_years;
  static bool do_output_transport_time_in_years;

  //
  // Init-to-steady parameters
  //
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
  static int  steady_max_time_steps;
  static Real steady_max_time_step_size;
  static int  steady_max_num_consecutive_success;
  static Real steady_extra_time_step_increase_factor;
  static int  steady_max_num_consecutive_increases;
  static Real steady_consecutive_increase_reduction_factor;
  static int  steady_limit_function_evals;
  static Real steady_abs_tolerance;
  static Real steady_rel_tolerance;
  static Real steady_abs_update_tolerance;
  static Real steady_rel_update_tolerance;
  static bool steady_do_grid_sequence;
  static Array<Real> steady_grid_sequence_new_level_dt_factor;
  static std::string steady_record_file;

  //
  // Newton step defaults
  //
  static int  richard_max_ls_iterations;
  static Real richard_min_ls_factor;
  static Real richard_ls_acceptance_factor;
  static Real richard_ls_reduction_factor;
  static int  richard_monitor_linear_solve;
  static int  richard_monitor_line_search;
  static Real richard_perturbation_scale_for_J;
  static bool richard_use_fd_jac;
  static bool richard_use_dense_Jacobian;
  static string richard_rel_perm_method;
  static int  richard_pressure_maxorder;
  static bool richard_scale_solution_before_solve;
  static bool richard_semi_analytic_J;
  static bool richard_centered_diff_J;
  static Real richard_variable_switch_saturation_threshold;
  static Real richard_dt_thresh_pure_steady;
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

  bool component_saturations_cached;
  MultiFab *sat_old_cached, *sat_new_cached;
  Real t_sat_old_cached, t_sat_new_cached;

  static Teuchos::ParameterList input_parameter_list;
};

/*
class PMFillPatchIterator
  :
  public MFIter
{
public:

  friend class PorousMedia;

  PMFillPatchIterator (AmrLevel& amrlevel,
                       MultiFab& leveldata);

  PMFillPatchIterator (AmrLevel&  amrlevel,
                       MultiFab&  leveldata,
                       int        boxGrow,
                       Real       time,
                       int        state_indx,
                       int        scomp,
                       int        ncomp);

  void Initialize (int             boxGrow,
                   const BoxArray& grids,
                   AmrLevel&       amrlevel);

  ~PMFillPatchIterator ();

  FArrayBox& operator() ();
  const IArrayBox& MatID() const { return m_matID[m_fpi.index()]; }

  const Box& UngrownBox () const { return MFIter::validbox(); }

  // From MFIter
  //
  // Returns the valid Box that current tile resides.
  //
  const Box& validbox () const { return m_fpi.validbox(); }
  //
  // Returns the Box of the FAB at which we currently point.
  //
  const Box fabbox () const { return m_fpi.fabbox(); }
  //
  // Increments iterator to the next tile we own.
  //
  void operator++ () { ++m_fpi;}
  //
  // Is the iterator valid i.e. is it associated with a FAB?
  //
  //bool isValid () { return m_fpi.isValid(); }
  bool isValid ();
  //
  // The index into the underlying BoxArray of the current FAB.
  //
  int index () const { return m_fpi.index(); }

private:
  //
  // Disallowed.
  //
  PMFillPatchIterator ();
  PMFillPatchIterator (const PMFillPatchIterator& rhs);
  PMFillPatchIterator& operator= (const PMFillPatchIterator& rhs);
  //
  // The data.
  //
  FillPatchIterator m_fpi;
  PorousMedia& m_pmlevel;
  iMultiFab m_matID;

  Real m_time;
  int m_stateIndex,m_scomp,m_ncomp;
};

*/

class PMFillPatchIterator
  :
  public FillPatchIterator
{
public:

  friend class PorousMedia;

  PMFillPatchIterator (AmrLevel& amrlevel,
                       MultiFab& leveldata);

  PMFillPatchIterator (AmrLevel&  amrlevel,
                       MultiFab&  leveldata,
                       int        boxGrow,
                       Real       time,
                       int        state_indx,
                       int        scomp,
                       int        ncomp);

  ~PMFillPatchIterator ();

  const IArrayBox& MatID() const { return m_matID[FillPatchIterator::index()]; }
  const FArrayBox& Aux() const { return m_aux[FillPatchIterator::index()]; }

  void Initialize (int             boxGrow,
                   const BoxArray& grids,
                   AmrLevel&       amrlevel);

  FArrayBox& operator() ();

private:
  //
  // Disallowed.
  //
  PMFillPatchIterator ();
  PMFillPatchIterator (const PMFillPatchIterator& rhs);
  PMFillPatchIterator& operator= (const PMFillPatchIterator& rhs);
  //
  // The data.
  //
  //FillPatchIterator m_fpi;
  PorousMedia& m_pmlevel;
  iMultiFab m_matID;
  MultiFab m_aux;

  Real m_time;
  int m_stateIndex,m_scomp,m_ncomp;
};

#endif /*_PorousMedia_H_*/
