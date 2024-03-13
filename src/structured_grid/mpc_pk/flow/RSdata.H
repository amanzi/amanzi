/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#ifndef _RSdata_H_
#define _RSdata_H_

#include <Layout.H>
#include <RockManager.H>

struct NLScontrol; // forward ref

struct RSdata
{
  /*
    Data structures required to evolve Richards equation
    Defines all the data that will be needed by the solver, but not who owns it or how to fill it.
    This separation will simplify the construction of a testing "harness" specific to the equation
    and data structures but one that does not rely on a fully constructed porous media integrator
    framework.
   */
  RSdata(int slev, int nlevs, Layout& layout, NLScontrol& nlsc, const RockManager* rm);
  virtual ~RSdata();

  virtual void SetUpMemory(NLScontrol& nlsc);
  virtual void ResetRhoSat() = 0; // Called by the solver before each solve
  virtual void SetPressureBC(const BCRec& pbc) {pressure_bc = pbc;}
  virtual void SetInflowVelocity(PArray<MFTower>& velocity,
				 Real             time) = 0;
  virtual void FillStateBndry (MFTower& press,
                               Real time) = 0;
  virtual void calcInvPressure (MFTower&       N,
                                const MFTower& P,
                                Real           time,
                                int            sComp,
                                int            dComp,
                                int            nGrow) const = 0;

  virtual void calcLambda (MFTower&       Lambda,
                           const MFTower& N,
                           Real           time,
                           int            sComp,
                           int            dComp,
                           int            nGrow) const = 0;

  virtual void calcRichardAlpha (MFTower&       Alpha,
                                 const MFTower& N,
                                 Real           time,
                                 int            sComp,
                                 int            dComp,
                                 int            nGrow) const = 0;

  virtual Array<int>& rinflowBCLo() = 0;
  virtual Array<int>& rinflowBCHi() = 0;

  virtual void SetDensity() = 0;
  virtual void SetGravity();
  virtual void SetViscosity() = 0;
  virtual void SetIsSaturated() = 0;

  virtual const PArray<MFTower>& GetKappaEC(Real t) {return KappaEC;}
  virtual const MFTower* GetKappaCCdir(Real t) {return KappaCCdir;}
  virtual const MFTower* GetSource(Real t) {return Source;}

  // Algorithm switches for Richard Solver
  bool semi_analytic_J;                      // Explicitly use dPc/ds relationship in computing the Jacobian
  std::string rel_perm_method;               // Upwind Krel based on sign of Darcy forcing, or average (harm/arith)
  Real variable_switch_saturation_threshold; // Newton iterate s if this is less than s, otherwise p

  int  pressure_maxorder;                    // Polynomial order of boundary condition extrapolation
  BCRec pressure_bc;                         // Boundary condition specifier for pressure
  Layout& layout;

  int nLevs;

  void SetCurrentTimestep(int step) {time_step = step;}
  int GetCurrentTimestep() const {return time_step;}

  const RockManager* GetRockManager() const {return rock_manager;}

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
  PArray<MFTower> KappaEC, DarcyVelocity, RichardCoefs;

  const Array<Real>& GetDensity() const {return density;}
  const Array<Real>& GetGravity();
  const Array<Real>& GetViscosity() const {return viscosity;}
  bool IsSaturated() const {return is_saturated;}
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
  Array<PArray<MultiFab> > kappaEC;
  bool memory_is_set_up;

  Array<Real> gravity;
  Array<Real> density;
  Array<Real> viscosity;
  bool is_saturated;

  const RockManager* rock_manager;
};

#endif
