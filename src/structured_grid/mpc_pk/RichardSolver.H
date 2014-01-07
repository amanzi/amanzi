#ifndef _RichardSolver_H_
#define _RichardSolver_H_

#include <PorousMedia.H>
#include <PMAmr.H>
#include <MFTower.H>

class NLSstatus
{
public:
  enum NLSstatusVal {
    NLS_SUCCESS           = 0,
    NLS_NONLINEAR_FAIL    = 1,
    NLS_LINEAR_FAIL       = 2,
    NLS_CATASTROPHIC_FAIL = 3,
    NLS_INVALID           = 4};

  NLSstatus(NLSstatusVal _val = NLS_INVALID)
    : val(_val)
  {
    BL_ASSERT(val>=0 && val<=4);
  }
  bool operator==(const NLSstatus& rhs) {return val == rhs.val;}
  bool operator<(const NLSstatus& rhs) {return val < rhs.val;}
  bool operator!=(const NLSstatus& rhs) {return val != rhs.val;}
  NLSstatusVal val;
};


class RSdata; // Forward ref
struct NLScontrol
{
  /*
    Control structure for generic time-stepper algorithm with Newton-based nonlinear solver
   */
  NLScontrol();

  void SetMaxNewtonIterations(int max_iter);
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
  void SetMaxDt(Real dt_max);

  int  NLIterationsTaken() {return nl_iterations_taken;}
  void SetNLIterationsTaken(int iters);
  int  MaxNLIterations() {return max_nl_iterations;}
  void ResetCounters();

  bool AdjustDt(Real      dt,
		NLSstatus nl_solver_status,
		Real&     new_dt);




  /*
    Time-step adjustment strategy (dt_new given dt)
    -----------------------------------------------

    Initialize dt_new to dt
    if first:
       num_consecutive_increases = num_consecutive_success = 0, prev_abs_err = -1, first = false

    if nl solved:
       last_chance = FALSE

       if nl_iterations_taken < min_nl_iterations_for_dt:
         num_consecutive_success++, num_consecutive_increases = 0

         if num_consecutive_success >= max_num_consecutive_success:
           fac = time_step_increase_factor

           if nl_iterations_taken < min_nl_iterations_for_dt_2:
             fac = time_step_increase_factor_2

           dt_new = dt * fac

       if nl_iterations_taken > max_nl_iterations_for_dt:
         (force J recalc next try)
         num_consecutive_increases++, num_consecutive_success=0

         if nl_iterations_taken > max_nl_iterations_for_dt:
           dt_new = dt * time_step_reduction_factor

       num_consecutive_failures_1 = num_consecutive_failures_2 = 0

    else // step was rejected

      if num_consecutive_failures_1 <= max_num_consecutive_failures_1:
        dt_new = dt * time_step_retry_factor

      else:
        num_consecutive_failures_2++

        if num_consecutive_failures_2 <= max_num_consecutive_failures_2:
          dt_new = dt * time_step_retry_factor_2

        else:

          if last_chance:
            return FALSE

          dt_new = dt * time_step_retry_factor_f

          last_chance = TRUE

      num_consecutive_success = 0

      (force J recalc next try)

    dt_new = std::min(max_time_step_size,dt_new)

    return TRUE
   */

  // Richard solver data
  RSdata* rs_data;

  // Control parameters (see strategy above)
  bool use_PETSc_snes;          // Use PETSc to solve the nonlinear problem
  bool use_fd_jac;              // Use matrix-free or finite-difference Jacobian
  bool use_dense_Jacobian;      // Use dense storage for Jacobian
  Real errfd;                   // Scale of perturbuation for fd Jacobian
  int  max_ls_iterations;       // Max line search attempts
  Real min_ls_factor;           // Minimum allowable line search scaling parameter
  Real ls_acceptance_factor;    // Maximun allowable reduction of residual in line search
  Real ls_reduction_factor;     // Scale to reduce line search scaling parameter for each try
  int  monitor_line_search;     // Write status information for line search
  int  maxit;                   // Maximum number of Newton iterations
  int  maxf;                    // Maximum number of function evaluations
  Real atol;                    // Absolute tolerance for norm(F)
  Real rtol;                    // Relative tolerance for norm(F)
  Real stol;                    // Absolute tolerance for norm(dp), dp is the Newton update
  bool scale_soln_before_solve; // Switch, scale P by 1/sigma before solve?
  bool centered_diff_J;         // Numerical J is computed with centered difference instead of forward
  bool ls_success;              // Pass-thru status from linesearch to convergence check
  std::string ls_reason;        // Status of linesearch to convergence success

  int max_nl_iterations;        // Max number of nonlinear (quasi-Newton iterations) to attempt
  Real max_nl_residual_norm;    // Largest workable residual value, abort solver if bigger
  int max_num_consecutive_success; // Number consecutive successful solves triggering increase dt
  int max_num_consecutive_failures_1; // Max allowable consecutive failures before retrying with dt reduced by time_step_retry_factor
  int max_num_consecutive_failures_2; // Max allowable consecutive failures before retrying with dt reduced by time_step_retry_factor_2
  int max_num_consecutive_increases; // Max allowable increases in residual beyond acceptance factor before retyring with dt reduced by consecutive_increase_reduction_factor
  Real consecutive_increase_reduction_factor;
  int min_nl_iterations_for_dt; // Maximum iterations to solve, below which next dt will be increased by time_step_increase_factor
  int min_nl_iterations_for_dt_2; // Maximum iterations to solve, below which next dt will be increased by time_step_increase_factor_2
  int max_nl_iterations_for_dt; // Maximum iterations to solve, above which next dt will not be modified
  Real max_time_step_size;      // Maximum time step to attempt
  Real time_step_increase_factor, time_step_increase_factor_2, time_step_reduction_factor;
  Real time_step_retry_factor, time_step_retry_factor_2, time_step_retry_factor_f;

private:

  // Status
  int nl_iterations_taken;
  Real nl_residual_norm;
  Real prev_abs_err;
  bool first;
  bool last_chance;
  int num_consecutive_success, num_consecutive_failures_1, num_consecutive_failures_2;
  int num_consecutive_increases;
};

struct RSdata
{
  /*
    Data structures required to evolve Richards equation
    Defines all the data that will be needed by the solver, but not who owns it or how to fill it.
    This separation will simplify the construction of a testing "harness" specific to the equation
    and data structures but one that does not rely on a fully constructed porous media integrator
    framework.
   */
  RSdata(int slev, int nlevs, Layout& layout, NLScontrol& nlsc);
  virtual ~RSdata();

  virtual void SetUpMemory(NLScontrol& nlsc);
  virtual void ResetRhoSat() = 0; // Called by the solver before each solve
  virtual void SetPressureBC(const BCRec& pbc) {pressure_bc = pbc;}
  virtual void SetInflowVelocity(PArray<MFTower>& velocity,
				 Real             time) = 0;
  virtual void FillStateBndry (Real time,
			       int  state_indx,
			       int  src_comp,
			       int  num_comp) = 0;
  virtual void calcInvPressure (MFTower&       N,
				const MFTower& P) const = 0;
  virtual void calcLambda (MFTower&       lbd,
			   const MFTower& N) = 0;
  virtual Array<int>& rinflowBCLo() = 0;
  virtual Array<int>& rinflowBCHi() = 0;


  virtual const MFTower* GetKappaCCdir(Real t) {return KappaCCdir;}
  virtual const MFTower* GetSource(Real time) {return Source;}


  // Algorithm switches for Richard Solver
  bool semi_analytic_J;                      // Explicitly use dPc/ds relationship in computing the Jacobian
  bool upwind_krel;                          // Upwind Krel based on sign of Darcy forcing
  bool subgrid_krel;                         // Where mixed cells, generate coarse krel by averaging refined krel
  Real variable_switch_saturation_threshold; // Newton iterate s if this is less than s, otherwise p

  int  pressure_maxorder;                    // Polynomial order of boundary condition extrapolation
  BCRec pressure_bc;                         // Boundary condition specifier for pressure
  Layout& layout;

  int nLevs;

  void SetCurrentTimestep(int step) {time_step = step;}
  int GetCurrentTimestep() const {return time_step;}


  // Cell-centered data
  Array<PArray<MultiFab> > ctmp;

  MFTower *InitialState;
  MFTower *Rhs;
  MFTower *RhoSatOld;
  MFTower *RhoSatNew;
  MFTower *Lambda;
  MFTower *Porosity;
  MFTower *SpecificStorage;
  MFTower *Pold;
  MFTower *Pnew;
  MFTower *KappaCCavg;
  MFTower *PCapParams;
  MFTower *KrParams;
  MFTower *Alpha;
  MFTower *CoeffCC;

  // Edge-centered data
  Array<PArray<MultiFab> > kappaEC;
  PArray<MFTower> KappaEC, DarcyVelocity, RichardCoefs;

  // Scalars
  Array<Real> gravity;
  Array<Real> density;
  Real time, dt;
  int time_step;
  int start_level, end_level;

  // Reuse Jacobian
  void SetMaxJacobianReuse(int max_num_reuse);
  void ResetJacobianCounter(int lev);
  void ResetJacobianCounter();
  bool UpdateJacobian(int lev);
  int max_num_Jacobian_reuses;
  Array<int> num_Jacobian_reuses_remaining;

protected:
  // These are "managed", access only via functions above
  MFTower *KappaCCdir;
  MFTower *Source;

  bool memory_is_set_up;
};

struct RSAMRdata
  : public RSdata
{
  /*
    An RSdata that manages memory and data evals via a PMAmr structure
   */
  RSAMRdata(int slev, int nlevs, Layout& layout, PMAmr* amrp, NLScontrol& nlsc);
  virtual ~RSAMRdata();
  virtual void SetUpMemory(NLScontrol& nlsc);
  virtual void ResetRhoSat();
  virtual void SetInflowVelocity(PArray<MFTower>& velocity,
				 Real             time);
  virtual void FillStateBndry (Real time,
			       int  state_indx,
			       int  src_comp,
			       int  num_comp);
  virtual void calcInvPressure (MFTower&       N,
				const MFTower& P) const;
  virtual void calcLambda (MFTower&       lbd,
			   const MFTower& N);
  virtual Array<int>& rinflowBCLo();
  virtual Array<int>& rinflowBCHi();

  enum RSdata_Property {
    RSdata_KappaCCdir = 0,
    RSdata_Num_Properties = 1
  };

  struct PropertyManagerCtx {
    MatFiller* matFiller;
    std::string property_name;
  };

  struct PropertyManager
  {
    PropertyManager(MFTower                           *pdat_ptr,
                    const std::set<PropertyManager*>& dependents,
                    PropertyManagerCtx                ctx)
      : property_dataPtr(pdat_ptr),
        property_dependents(dependents),
        property_ctx(ctx),
        eval_time(-1) {}
    virtual const MFTower* UpdateProperty(Real t);
    virtual void UpdateDependents(Real t)
      {
        for (std::set<PropertyManager*>::const_iterator it=property_dependents.begin(),
               End=property_dependents.end(); it!=End; ++it) {
          (*it)->UpdateProperty(t);
        }
      }
    MFTower* property_dataPtr;
    std::set<PropertyManager*> property_dependents;
    PropertyManagerCtx property_ctx;
    Real eval_time;
  };

  virtual const MFTower* GetKappaCCdir(Real t);
  virtual const MFTower* GetSource(Real t);

  PorousMedia& GetPMlevel(int level) {return pm[level];}
  PMAmr* pm_amr;
  PArray<PorousMedia> pm;
  std::map<RSdata_Property,PropertyManager*> managed_properties;
  Real source_eval_time;
};

class RichardSolver
{
public:

  RichardSolver(RSdata& rs_data, NLScontrol& nlsc);

  ~RichardSolver();

  void SetRecordFile(const std::string& file_name) {record_file = file_name;}

  const std::string& GetRecordFile() const {return record_file;}

  int Solve(Real        time,
            Real        dt,
	    int         timestep,
            NLScontrol& nlsc);

  void ComputeDarcyVelocity(MFTower& pressure,
                            Real     t);


  void ResetRhoSat() {rs_data.ResetRhoSat();}
  void SetCurrentTimestep(int step) {rs_data.SetCurrentTimestep(step);}
  int GetCurrentTimestep() const {return rs_data.GetCurrentTimestep();}

  Layout& GetLayout() {return rs_data.layout;}
  const Layout& GetLayout() const {return rs_data.layout;}
  void SetTime(Real time) {rs_data.time = time;}
  Real GetTime() const {return rs_data.time;}
  void SetDt(Real dt) {rs_data.dt = dt;}
  Real GetDt() const {return rs_data.dt;}
  Real GetNumLevels() const {return rs_data.nLevs;}
  MFTower& GetRhoSatN() {return *rs_data.RhoSatOld;}
  MFTower& GetRhoSatNp1() {return *rs_data.RhoSatNew;}
  MFTower& GetAlpha() {return *rs_data.Alpha;}
  MFTower& GetLambda() {return *rs_data.Lambda;}
  MFTower& GetKappaCCavg() {return *rs_data.KappaCCavg;}
  PArray<MFTower>& GetKappaEC() {return rs_data.KappaEC;}
  MFTower& GetPorosity() {return *rs_data.Porosity;}
  MFTower& GetSpecificStorage() {return *rs_data.SpecificStorage;}
  MFTower& GetPressureNp1() {return *rs_data.Pnew;}
  MFTower& GetPressureN() {return *rs_data.Pold;}
  MFTower& GetResidual() {return *rs_data.Rhs;}
  MFTower& GetPCapParams() {return *rs_data.PCapParams;}
  MFTower& GetKrParams() {return *rs_data.KrParams;}
  MFTower& GetCoeffCC() {return *rs_data.CoeffCC;}

  // Time-dependent properties
  const MFTower* GetKappaCCdir(Real t) {return rs_data.GetKappaCCdir(t);}

  PArray<MFTower>& GetDarcyVelocity() {return rs_data.DarcyVelocity;}
  PArray<MFTower>& GetRichardCoefs() {return rs_data.RichardCoefs;}
  const Array<Real>& GetDensity() {return rs_data.density;}
  const Array<Real>& GetGravity() {return rs_data.gravity;}
  const Array<BoxArray>& GridArray() const {return rs_data.layout.GridArray();}
  const Array<Geometry>& GeomArray() const {return rs_data.layout.GeomArray();}
  const Array<IntVect>& RefRatio() const {return rs_data.layout.RefRatio();}
  RSdata& GetRSdata() {return rs_data;}
  NLScontrol& GetNLScontrol() {return nlsc;}


  Vec& GetPressureV() {return SolnV;}
  Vec& GetResidualV() {return RhsV;}
  Vec& GetTrialResV() {return GV;}
  Vec& GetSolnTypV() {return SolnTypV;}
  Vec& GetSolnTypInvV() {return SolnTypInvV;}
  Vec& GetAlphaV() {return AlphaV;}
  SNES& GetSNES(){return snes;}
  MatFDColoring& GetMatFDColoring() {return matfdcoloring;}

  // Helper or internal functions
  void DivRhoU(MFTower& DivRhoU,
               MFTower& pressure,
               Real     t);

  void CalcResidual(MFTower& residual,
		    MFTower& pressure,
		    Real     t,
		    Real     dt);

  void CreateJac(Mat& J, 
		 MFTower& pressure,
		 Real dt);

  void ComputeRichardAlpha(Vec& Alpha,const Vec& Pressure);

  static void SetTheRichardSolver(RichardSolver* ptr);

  void ResetRemainingJacobianReuses();
  void UnsetRemainingJacobianReuses();
  bool ReusePreviousJacobian();

  void BuildMLPropEval();

  void BuildOpSkel(Mat& J);

  void YpambX(MFTower&           Y,
	      const MFTower&     X,
	      const Array<Real>& a,
	      const Array<Real>& b,
	      int                sComp,
	      int                dComp,
	      int                nComp,
	      int                nGrow) const;
  
  void FillPatch(MFTower& mft,
		 int sComp,
		 int nComp,
		 bool do_piecewise_constant=false);
  
  void CenterToEdgeUpwind(PArray<MFTower>&       mfte,
			  MFTower&               mftc,
			  const PArray<MFTower>& sgn,
			  int                    nComp,
                          const BCRec&           bc) const;
  
  void XmultYZ(MFTower&           X,
	       const MFTower&     Y,
	       const MFTower&     Z,
	       int                sCompY=0,
	       int                sCompZ=0,
	       int                dComp=0,
	       int                nComp=1,
	       int                nGrow=0);
  
  void CCtoECgradAdd(PArray<MFTower>&            mfte,
		     const MFTower&              mftc,
		     const Array<Array<Real> >&  a,
		     int                         sComp=0,
		     int                         dComp=0,
		     int                         nComp=1) const;
protected:
  
  MFTFillPatch* mftfp;
  RSdata& rs_data;
  NLScontrol& nlsc;

  SNES           snes;
  ISColoring     iscoloring;
  MatFDColoring  matfdcoloring;

  Vec RhsV, SolnV, SolnTypV, SolnTypInvV, GV, AlphaV;
  Mat Jac;

  int num_remaining_Jacobian_reuses;
  int current_timestep;
  std::string record_file;

  // Data for "subgrid" material model, off by default
  Array<BoxArray> state_to_fill;
  Array<BoxArray> derive_to_fill;
  PArray<MultiFab> phif, pcPf, kf, krf, pf, lf, lc;
};

#endif
