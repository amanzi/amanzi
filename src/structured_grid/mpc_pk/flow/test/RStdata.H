#ifndef _RStdata_H_
#define _RStdata_H_

#include <RichardSolver.H>

struct RStstruct
{
  Real kappa, phi, alpha, Sr, m, specific_storage, g, inflow_velocity, Pwt;
  Array<Real> rho, mu;
  bool saturated;
};

struct RStdata
  : public RSdata
{

  RStdata(int slev, int nlevs, Layout& layout, NLScontrol& nlsc, RStstruct& inputs);
  virtual ~RStdata();
  virtual void SetUpMemory(NLScontrol& nlsc);
  virtual void ResetRhoSat();
  virtual void SetInflowVelocity(PArray<MFTower>& velocity,
				 Real             time);
  virtual void FillStateBndry (MFTower& pressure,
                               Real time);
  virtual void calcInvPressure (MFTower&       N,
				const MFTower& P) const;
  virtual void calcLambda (MFTower&       lbd,
			   const MFTower& N);
  virtual void calcRichardAlpha(MFTower&       alpha,
                                const MFTower& rhoSat,
                                Real           time);
  virtual Array<int>& rinflowBCLo();
  virtual Array<int>& rinflowBCHi();
  virtual void SetDensity() {};
  virtual void SetGravity() {};
  virtual void SetViscosity() {};
  virtual void SetIsSaturated() {};
  virtual const MFTower* GetKappaCCdir(Real t);
  virtual const PArray<MFTower>& GetKappaEC(Real t);
  virtual const MFTower* GetSource(Real t);

  void SetConstantValue(MFTower& mft,
                        Real     val,
                        int      dComp,
                        int      nComp);
  void SetPCapParams();
  void SetPorosity();
  void SetSpecificStorage();

  bool IsNewTime(Real time) const;
  Real old_time, new_time;
  Real eval_time_for_KappaCCdir, eval_time_for_KappaEC, eval_time_for_source;

  static Real Pa_per_ATM;
  RStstruct inputs;
};

#endif
