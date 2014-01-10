#include <RSdata.H>
#include <NLScontrol.H>

// Richards Solver algorithm defaults
static bool semi_analytic_J_DEF                       = false;
static bool upwind_krel_DEF                           = false;
static bool subgrid_krel_DEF                          = false; // Not working yet
static Real variable_switch_saturation_threshold_DEF  = -0.9999;
static int pressure_maxorder_DEF                      = 3;
static int max_num_Jacobian_reuses_DEF                = 0; // This just doesnt seem to work very well....

RSdata::RSdata(int slev, int nlevs, Layout& _layout, NLScontrol& nlsc)
  : layout(_layout), nLevs(nlevs), start_level(slev), end_level(slev+nlevs-1),
    num_Jacobian_reuses_remaining(nlevs,0)
{
  semi_analytic_J = semi_analytic_J_DEF;
  variable_switch_saturation_threshold = variable_switch_saturation_threshold_DEF;
  upwind_krel = upwind_krel_DEF;
  subgrid_krel = subgrid_krel_DEF;
  pressure_maxorder = pressure_maxorder_DEF;
  max_num_Jacobian_reuses = max_num_Jacobian_reuses_DEF;

  // FIXME: Should we set default values for g and rho?

  nlsc.rs_data = this;
  ResetJacobianCounter();

  Porosity = 0;
  KappaCCavg = 0;
  KrParams = 0;
  Lambda = 0;
  PCapParams = 0;
  SpecificStorage = 0;
  Source = 0;
  KappaCCdir = 0;
  CoeffCC = 0;

  InitialState = 0;
  RhoSatOld = 0;
  RhoSatNew = 0;
  Pold = 0;
  Pnew = 0;
  Rhs = 0;
  Alpha = 0;

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
}

void RSdata::SetMaxJacobianReuse(int max_num_reuse) {max_num_Jacobian_reuses=max_num_reuse;}
void RSdata::ResetJacobianCounter(int lev) {num_Jacobian_reuses_remaining[lev]=max_num_Jacobian_reuses;}

void
RSdata::ResetJacobianCounter()
{
  int nlevs = end_level - start_level +1;
  for (int lev=0; lev<nlevs; ++lev) {
    ResetJacobianCounter(lev);
  }
}

bool
RSdata::UpdateJacobian(int lev)
{
  bool do_Jacobian_eval = false;
  num_Jacobian_reuses_remaining[lev]--;
  if (num_Jacobian_reuses_remaining[lev] <= 0) {
    do_Jacobian_eval = true;
    num_Jacobian_reuses_remaining[lev] = max_num_Jacobian_reuses;
  }
  return do_Jacobian_eval;
}

void
RSdata::SetUpMemory(NLScontrol& nlsc)
{
  InitialState = new MFTower(layout,IndexType(IntVect::TheZeroVector()),1,1,nLevs);
  SetDensity();
  SetGravity();
  SetViscosity();
  SetIsSaturated();
}

RSdata::~RSdata()
{
  delete InitialState;
}


