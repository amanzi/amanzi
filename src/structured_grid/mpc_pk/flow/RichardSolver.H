#ifndef _RichardSolver_H_
#define _RichardSolver_H_

#include <petsc.h>
#include <MFTower.H>
#include <NLScontrol.H>
#include <RSdata.H>

class RichardSolver
{
public:

  RichardSolver(RSdata& rs_data, NLScontrol& nlsc);

  ~RichardSolver();

  void SetRecordFile(const std::string& file_name) {record_file = file_name;}

  const std::string& GetRecordFile() const {return record_file;}

  int Solve(Real        old_time,
            Real        new_time,
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
  MFTower& GetPorosity() {return *rs_data.Porosity;}
  MFTower& GetSpecificStorage() {return *rs_data.SpecificStorage;}
  MFTower& GetPressureNp1() {return *rs_data.Pnew;}
  MFTower& GetPressureN() {return *rs_data.Pold;}
  MFTower& GetResidual() {return *rs_data.Rhs;}
  MFTower& GetPCapParams() {return *rs_data.PCapParams;}
  MFTower& GetKrParams() {return *rs_data.KrParams;}
  MFTower& GetCoeffCC() {return *rs_data.CoeffCC;}

  const iMultiFab& GetMaterialID(int lev) const {return materialID[lev];}

  // Time-dependent properties
  const MFTower* GetKappaCCdir(Real t) {return rs_data.GetKappaCCdir(t);}
  const PArray<MFTower>& GetKappaEC(Real t) {return rs_data.GetKappaEC(t);}

  PArray<MFTower>& GetDarcyVelocity() {return rs_data.DarcyVelocity;}
  PArray<MFTower>& GetRichardCoefs() {return rs_data.RichardCoefs;}
  const Array<Real>& GetDensity() {return rs_data.GetDensity();}
  const Array<Real>& GetGravity() {return rs_data.GetGravity();}
  const Array<Real>& GetViscosity() {return rs_data.GetViscosity();}
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
                 Real t,
		 Real dt);

  void ComputeRichardAlpha(Vec&       Alpha,
                           const Vec& Pressure,
                           Real       t);

  static void SetTheRichardSolver(RichardSolver* ptr);

  void ResetRemainingJacobianReuses();
  void UnsetRemainingJacobianReuses();
  bool ReusePreviousJacobian();

  void BuildOpSkel(Mat& J, bool calcSpace);

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

  PArray<iMultiFab> materialID;
};

#endif
