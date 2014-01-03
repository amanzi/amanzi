#ifndef _RichardSolver_H_
#define _RichardSolver_H_

#include <PorousMedia.H>
#include <PMAmr.H>
#include <MFTower.H>



class NLSreason
{
public:
  enum NLSreasonVal {
    NLS_SUCCESS,
    NLS_NONLINEAR_FAIL,
    NLS_LINEAR_FAIL,
    NLS_CATASTROPHIC_FAIL,
    NLS_INVALID };

  NLSreason(NLSreasonVal _val = NLS_INVALID)
    : val(_val)
  {
    BL_ASSERT(val>=0 && val<=4);
  }
  bool operator==(const NLSreason& rhs) {return val == rhs.val;}
  bool operator<(const NLSreason& rhs) {return val < rhs.val;}
  bool operator!=(const NLSreason& rhs) {return val != rhs.val;}
  NLSreasonVal val;
};



struct NLScontrol
{
  /*
    Control structure for generic time-stepper algorithm with Newton-based nonlinear solver
   */
  NLScontrol();

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
  void SetMaxDt(Real dt_max);

  bool UpdateJacobian(int lev);
  void ResetJacobianCounter(int lev);
  void ResetJacobianCounter();
  int  NLIterationsTaken() {return nl_iterations_taken;}
  void SetNLIterationsTaken(int iters);
  int  MaxNLIterations() {return max_nl_iterations;}
  void ResetCounters();

  bool AdjustDt(Real      dt, 
		NLSreason nl_solver_status, 
		Real&     new_dt);

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
  int max_num_Jacobian_reuses;

private:
  int start_level, end_level;
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
  Real max_time_step_size;
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
  virtual ~RSdata() {}

  virtual void InitBuild(NLScontrol& nlsc) = 0; // Called by the constructor to set up memory
  virtual void ResetRhoSat() = 0; // Called by the solver before each solve
  virtual void SetPressureBC(const BCRec& pbc) {pressure_bc = pbc;}
  virtual void SetInflowVelocity(PArray<MFTower>& velocity,
				 Real             time) = 0;
  virtual void FillStateBndry (Real time,
			       int  state_indx,
			       int  src_comp, 
			       int  num_comp,
			       int  amr_level) = 0;
  virtual void calcInvPressure (MultiFab&       N,
				const MultiFab& P,
				int             amr_level) const = 0;
  virtual void calcLambda (MultiFab*       lbd,
			   const MultiFab& N,
			   int             amr_level) = 0;
  virtual void getSource (MultiFab& source,
			  int       nGrow,
			  int       strt_comp,
			  int       num_comp,
			  Real      time,
			  int       amr_level,
			  bool      do_rho_scale = false) = 0;
  virtual Array<int>& rinflowBCLo() = 0;
  virtual Array<int>& rinflowBCHi() = 0;


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
  PArray <MultiFab>  initialState;

  MFTower *Rhs;
  MFTower *RhoSatOld;
  MFTower *RhoSatNew;
  MFTower *Lambda;
  MFTower *Porosity;
  MFTower *SpecificStorage;
  MFTower *Pold;
  MFTower *Pnew;
  MFTower *KappaCCavg;
  MFTower *KappaCCdir;
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
};

struct RSAMRdata
  : public RSdata
{
  /*
    An RSdata that manages memory and data evals via a PMAmr structure
   */
  RSAMRdata(int slev, int nlevs, Layout& layout, PMAmr* amrp, NLScontrol& nlsc);
  virtual ~RSAMRdata();
  virtual void InitBuild(NLScontrol& nlsc);
  virtual void ResetRhoSat();
  virtual void SetInflowVelocity(PArray<MFTower>& velocity,
				 Real             time);
  virtual void FillStateBndry (Real time,
			       int  state_indx,
			       int  src_comp, 
			       int  num_comp,
			       int  amr_level);
  virtual void calcInvPressure (MultiFab&       N,
				const MultiFab& P,
				int             amr_level) const;
  virtual void calcLambda (MultiFab*       lbd,
			   const MultiFab& N,
			   int             amr_level);
  virtual void getSource (MultiFab& source,
			  int       nGrow,
			  int       strt_comp,
			  int       num_comp,
			  Real      time,
			  int       amr_level,
			  bool      do_rho_scale = false);
  virtual Array<int>& rinflowBCLo();
  virtual Array<int>& rinflowBCHi();

  PorousMedia& GetPMlevel(int level) {return pm[level];}
  PMAmr* pm_amr;
  PArray<PorousMedia> pm;
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

  void ComputeDarcyVelocity(PArray<MFTower>&       darcy_vel,
			    MFTower&               pressure,
			    MFTower&               rhoSat,
			    MFTower&               lambda,
			    const PArray<MFTower>& kappa,
			    const Array<Real>&     density,
			    const Array<Real>&     gravity,
			    Real                   t);


  void UpdateDarcyVelocity(MFTower& pressure,
			   Real     t); // Set internal data


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
  MFTower& GetKappaCCdir() {return *rs_data.KappaCCdir;}
  PArray<MFTower>& GetKappaEC() {return rs_data.KappaEC;}
  MFTower& GetPorosity() {return *rs_data.Porosity;}
  MFTower& GetSpecificStorage() {return *rs_data.SpecificStorage;}
  MFTower& GetPressureNp1() {return *rs_data.Pnew;}
  MFTower& GetPressureN() {return *rs_data.Pold;}
  MFTower& GetResidual() {return *rs_data.Rhs;}
  MFTower& GetPCapParams() {return *rs_data.PCapParams;}
  MFTower& GetKrParams() {return *rs_data.KrParams;}
  MFTower& GetCoeffCC() {return *rs_data.CoeffCC;}

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

  void SetPermeability(Real cur_time);

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
  
  void CCtoECgradAdd(PArray<MFTower>& mfte,
		     const MFTower&   mftc,
		     const FArrayBox& a,
		     int              sComp=0,
		     int              dComp=0,
		     int              nComp=1) const;
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
