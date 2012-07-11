#ifndef _RichardSolver_H_
#define _RichardSolver_H_

#include <PorousMedia.H>
#include <PMAmr.H>
#include <MFTower.H>

class RichardSolver
{
public:
  RichardSolver(PMAmr& pm_amr);
  void ResetRhoSat();
  void SetTime(Real time);
  Real Time() const;
  void SetDt(Real dt);
  Real Dt() const;
  MFTower& GetRhoSatN() {return *RhoSatOld;}
  MFTower& GetRhoSatNp1() {return *RhoSatNew;}
  MFTower& GetLambda() {return *Lambda;}
  PArray<MFTower>& GetKappa() {return Kappa;}
  MFTower& GetPorosity() {return *Porosity;}
  MFTower& GetPressure() {return *Pnew;}
  Vec& GetPressureV() {return SolnV;}
  MFTower& GetResidual() {return *Rhs;}
  Vec& GetResidualV() {return RhsV;}
  PArray<MFTower>& GetDarcyVelocity() {return DarcyVelocity;}
  PArray<MFTower>& GetRichardCoefs() {return RichardCoefs;}
  const Array<Real>& GetDensity() {return density;}
  const Array<Real>& GetGravity() {return gravity;}
  SNES& GetSNES(){return snes;}
  PMAmr& GetPMAMR() {return pm_amr;}
  const Array<BoxArray>& GridArray() {return PMAmr::GetLayout().GridArray();}
  const Array<Geometry>& GeomArray() {return PMAmr::GetLayout().GeomArray();}
  const Array<IntVect>& RefRatio() {return PMAmr::GetLayout().RefRatio();}

  void DivU(MFTower& divu,
	    MFTower& rhoSat,
	    MFTower& pressure,
	    Real     time);
  
  void DpDtResidual(MFTower& residual,
		    MFTower& pressure,
		    Real     time,
		    Real     dt);

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
		 bool do_piecewise_constant=false);
  
  void CenterToEdgeUpwind(PArray<MFTower>&       mfte,
			  MFTower&               mftc,
			  const PArray<MFTower>& sgn,
			  int                    nComp=1) const;
  
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
  PArray<PorousMedia> pm;
  MFTFillPatch* mftfp;

  SNES           snes;
  ISColoring     iscoloring;
  MatFDColoring  matfdcoloring;

  int nLevs;
  Array<PArray<MultiFab> > ktmp, ctmp; // Currently maintains its own copy of these, should be fixed
  MFTower *Rhs, *RhoSatOld, *RhoSatNew, *Lambda, *Porosity, *Pnew, *Pinit;
  PArray<MFTower> Kappa, DarcyVelocity, RichardCoefs;

  Vec RhsV, SolnV;

  Array<Real> gravity;
  Array<Real> density;
  Real mytime, mydt;
};

#endif
