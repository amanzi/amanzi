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
  };

  RichardSolver(PMAmr&          pm_amr,
		const RSParams& params);
  void ResetRhoSat();
  void SetTime(Real time);
  Real GetTime() const;
  void SetDt(Real dt);
  Real GetDt() const;
  MFTower& GetRhoSatN() {return *RhoSatOld;}
  MFTower& GetRhoSatNp1() {return *RhoSatNew;}
  MFTower& GetLambda() {return *Lambda;}
  MFTower& GetKappaCC() {return *KappaCC;}
  PArray<MFTower>& GetKappa() {return Kappa;}
  MFTower& GetPorosity() {return *Porosity;}
  MFTower& GetPressure() {return *Pnew;}
  MFTower& GetPStar() {return *PStar;}
  MFTower& GetRhoSatStar() {return *RhoSatStar;}
  MFTower& GetAlphaStar() {return *AlphaStar;}
  Vec& GetPressureV() {return SolnV;}
  MFTower& GetResidual() {return *Rhs;}
  Vec& GetResidualV() {return RhsV;}
  Vec& GetTrialResV() {return GV;}
  Vec& GetPCapParamsV() {return PCapParamsV;}
  PArray<MFTower>& GetDarcyVelocity() {return DarcyVelocity;}
  PArray<MFTower>& GetRichardCoefs() {return RichardCoefs;}
  const Array<Real>& GetDensity() {return density;}
  const Array<Real>& GetGravity() {return gravity;}
  SNES& GetSNES(){return snes;}
  PMAmr& GetPMAMR() {return pm_amr;}
  MatFDColoring& GetMatFDColoring() {return matfdcoloring;}
  PorousMedia& GetPMlevel(int level) {return pm[level];}
  const Array<BoxArray>& GridArray() const {return PMAmr::GetLayout().GridArray();}
  const Array<Geometry>& GeomArray() const {return PMAmr::GetLayout().GeomArray();}
  const Array<IntVect>& RefRatio() const {return PMAmr::GetLayout().RefRatio();}
  const RichardSolver::RSParams& Parameters() {return params;}

  void UpdateDarcyVelocity(MFTower& pressure,
			   Real     time); // Set internal data

  void UpdateDarcyVelocity(PArray<MultiFab>& pressure,
			   Real              time); // Set internal data

  void ComputeDarcyVelocity(PArray<MFTower>&       darcy_vel,
			    MFTower&               pressure,
			    MFTower&               rhoSat,
			    MFTower&               lambda,
			    const PArray<MFTower>& kappa,
			    const Array<Real>&     density,
			    const Array<Real>&     gravity,
			    Real                   time);

  void DpDtResidual(MFTower& residual,
		    MFTower& pressure,
		    Real     time,
		    Real     dt);

  void DpDtResidual_Alt(MFTower& residual,
                        MFTower& pressure,
                        Real     time,
                        Real     dt);

  void CreateJac(Mat& J, 
		 MFTower& pressure,
		 Real dt);

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

  SNES           snes;
  ISColoring     iscoloring;
  MatFDColoring  matfdcoloring;

  int nLevs;
  Array<PArray<MultiFab> > ktmp, ctmp; // Currently maintains its own copy of these, should be fixed
  MFTower *Rhs, *RhoSatOld, *RhoSatNew, *Lambda, *Porosity, *Pnew, *PStar, *RhoSatStar, *AlphaStar, *KappaCC;
  PArray<MFTower> Kappa, DarcyVelocity, RichardCoefs;

  Vec RhsV, SolnV, PCapParamsV, GV;
  Mat Jac;

  Array<Real> gravity;
  Array<Real> density;
  Real mytime, mydt;
};

#endif
