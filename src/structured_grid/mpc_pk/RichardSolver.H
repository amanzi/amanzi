#ifndef _RichardSolver_H_
#define _RichardSolver_H_

#include <PorousMedia.H>
#include <PMAmr.H>
#include <MFTower.H>

class RichardNLSdata;

class RichardSolver
{
public:
  struct RSParams
  {
    RSParams();
    bool use_fd_jac; // Use matrix-free or finite-difference Jacobian
    bool use_dense_Jacobian; // Use dense storage for Jacobian
    int upwind_krel; // Upwind Krel based on sign of Darcy forcing
    int pressure_maxorder; // Polynomial order of boundary condition extrapolation
    Real errfd; // Scale of perturbuation for fd Jacobian
    int max_ls_iterations; // Max line search attempts
    Real min_ls_factor; // Minimum allowable line search scaling parameter
    Real ls_acceptance_factor; // Maximun allowable reduction of residual in line search
    Real ls_reduction_factor; // Scale to reduce line search scaling parameter for each try
    int monitor_line_search; // Write status information for line search
    int  maxit; // Maximum number of Newton iterations
    int maxf; // Maximum number of function evaluations
    Real atol; // Absolute tolerance for norm(F)
    Real rtol; // Relative tolerance for norm(F)
    Real stol; // Absolute tolerance for norm(dp), dp is the Newton update
    bool scale_soln_before_solve; // Switch, scale P by 1/sigma before solve?
    bool semi_analytic_J; // Explicitly use dPc/ds relationship in computing the Jacobian
    bool centered_diff_J; // Numerical J is computed with centered difference instead of forward
  };

  RichardSolver(PMAmr&          pm_amr,
		const RSParams& params,
		Layout&         layout);
  ~RichardSolver();
  Layout& GetLayout() {return layout;}
  void ResetRhoSat();
  void SetTime(Real time);
  Real GetTime() const;
  void SetDt(Real dt);
  Real GetDt() const;
  Real GetNumLevels() const {return nLevs;}
  MFTower& GetRhoSatN() {return *RhoSatOld;}
  MFTower& GetRhoSatNp1() {return *RhoSatNew;}
  MFTower& GetAlpha() {return *Alpha;}
  Vec& GetAlphaV() {return AlphaV;}
  MFTower& GetLambda() {return *Lambda;}
  MFTower& GetKappaCC() {return *KappaCC;}
  PArray<MFTower>& GetKappa() {return Kappa;}
  MFTower& GetPorosity() {return *Porosity;}
  MFTower& GetPressure() {return *Pnew;}
  Vec& GetPressureV() {return SolnV;}
  MFTower& GetResidual() {return *Rhs;}
  Vec& GetResidualV() {return RhsV;}
  Vec& GetTrialResV() {return GV;}
  MFTower& GetPCapParams() {return *PCapParams;}
  Vec& GetSolnTypV() {return SolnTypV;}
  Vec& GetSolnTypInvV() {return SolnTypInvV;}
  PArray<MFTower>& GetDarcyVelocity() {return DarcyVelocity;}
  PArray<MFTower>& GetRichardCoefs() {return RichardCoefs;}
  const Array<Real>& GetDensity() {return density;}
  const Array<Real>& GetGravity() {return gravity;}
  SNES& GetSNES(){return snes;}
  PMAmr& GetPMAMR() {return pm_amr;}
  MatFDColoring& GetMatFDColoring() {return matfdcoloring;}
  PorousMedia& GetPMlevel(int level) {return pm[level];}
  const Array<BoxArray>& GridArray() const {return layout.GridArray();}
  const Array<Geometry>& GeomArray() const {return layout.GeomArray();}
  const Array<IntVect>& RefRatio() const {return layout.RefRatio();}
  const RichardSolver::RSParams& Parameters() {return params;}

  void UpdateDarcyVelocity(MFTower& pressure,
			   Real     t); // Set internal data

  void UpdateDarcyVelocity(PArray<MultiFab>& pressure,
			   Real              t); // Set internal data

  void ComputeDarcyVelocity(PArray<MFTower>&       darcy_vel,
			    MFTower&               pressure,
			    MFTower&               rhoSat,
			    MFTower&               lambda,
			    const PArray<MFTower>& kappa,
			    const Array<Real>&     density,
			    const Array<Real>&     gravity,
			    Real                   t);

  void DivRhoU(MFTower& DivRhoU,
               MFTower& pressure,
               Real     t);

  void DpDtResidual(MFTower& residual,
		    MFTower& pressure,
		    Real     t,
		    Real     dt);

  void CreateJac(Mat& J, 
		 MFTower& pressure,
		 Real dt);

  void ComputeRichardAlpha(Vec& Alpha,const Vec& Pressure);

  static void SetTheRichardSolver(RichardSolver* ptr);

  int Solve(Real            time,
            Real            dt,
	    int             timestep,
            RichardNLSdata& nl_data);
    

protected:
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
  
  void CCtoECgradAdd(PArray<MFTower>& mfte,
		     const MFTower&   mftc,
		     const FArrayBox& a,
		     int              sComp=0,
		     int              dComp=0,
		     int              nComp=1) const;
  
  void SetInflowVelocity(PArray<MFTower>& velocity,
			 Real             time);
  
  PMAmr& pm_amr;
  RSParams params;
  PArray<PorousMedia> pm;
  MFTFillPatch* mftfp;
  Layout& layout;

  SNES           snes;
  ISColoring     iscoloring;
  MatFDColoring  matfdcoloring;

  int nLevs;
  Array<PArray<MultiFab> > ktmp, ctmp; // Currently maintains its own copy of these, should be fixed
  MFTower *Rhs, *RhoSatOld, *RhoSatNew, *Lambda, *Porosity, *Pnew, *KappaCC, *PCapParams, *Alpha;
  PArray<MFTower> Kappa, DarcyVelocity, RichardCoefs;

  Vec RhsV, SolnV, SolnTypV, SolnTypInvV, GV, AlphaV;
  Mat Jac;

  Array<Real> gravity;
  Array<Real> density;
  Real mytime, mydt;
};

#endif
