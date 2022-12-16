/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#ifndef _RStdata_H_
#define _RStdata_H_

#include <RichardSolver.H>

struct RStstruct
{
  Real g, inflow_velocity, Pwt, Kr_smoothing_max_pcap;
  Array<Real> rho, mu;
  bool saturated;
};

struct RStdata
  : public RSdata
{

  RStdata(int slev, int nlevs, Layout& layout, NLScontrol& nlsc, RStstruct& inputs, const RockManager* rm);
  virtual ~RStdata();
  virtual void SetUpMemory(NLScontrol& nlsc);
  virtual void ResetRhoSat();
  virtual void SetInflowVelocity(PArray<MFTower>& velocity,
				 Real             time);
  virtual void FillStateBndry (MFTower& pressure,
                               Real time);
  virtual void calcInvPressure (MFTower&       N,
                                const MFTower& P,
                                Real           time,
                                int            sComp,
                                int            dComp,
                                int            nGrow) const;
  virtual void calcLambda (MFTower&       Lambda,
                           const MFTower& N,
                           Real           time,
                           int            sComp,
                           int            dComp,
                           int            nGrow) const;

  virtual void calcRichardAlpha (MFTower&       Alpha,
                                 const MFTower& N,
                                 Real           time,
                                 int            sComp,
                                 int            dComp,
                                 int            nGrow) const;

  virtual Array<int>& rinflowBCLo();
  virtual Array<int>& rinflowBCHi();
  virtual void SetDensity() {};
  virtual void SetGravity() {};
  virtual void SetViscosity() {};
  virtual void SetIsSaturated() {};
  virtual const MFTower* GetKappaCCdir(Real t);
  virtual const PArray<MFTower>& GetKappaEC(Real t);
  virtual const MFTower* GetSource(Real t);

  void SetPCapParams(Real t);
  void SetPorosity(Real t);
  void SetSpecificStorage(Real t);

  bool IsNewTime(Real time) const;
  Real old_time, new_time;
  Real eval_time_for_KappaCCdir, eval_time_for_KappaEC, eval_time_for_source;

  static Real Pa_per_ATM;
  RStstruct inputs;
  PArray<iMultiFab> materialID;
};

#endif
