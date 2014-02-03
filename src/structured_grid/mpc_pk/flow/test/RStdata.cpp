#include <RStdata.H>

Real RStdata::Pa_per_ATM = 101325;
static bool saturated = false;

RStdata::RStdata(int slev, int nlevs, Layout& layout, NLScontrol& nlsc, RStstruct& _inputs)
  : RSdata(slev,nlevs,layout,nlsc), inputs(_inputs)
{
  Porosity = 0;
  KappaCCavg = 0;
  Lambda = 0;
  PCapParams = 0;
  SpecificStorage = 0;
  Source = 0;
  KappaCCdir = 0;
  CoeffCC = 0;

  Rhs = 0;
  Alpha = 0;

  eval_time_for_KappaCCdir = -1;
  eval_time_for_KappaEC = -1;
  eval_time_for_source = -1;

  // Set soil/fluid properties "known" by the base class
  gravity.resize(BL_SPACEDIM,0); gravity[BL_SPACEDIM-1] = inputs.g;
  density = inputs.rho;
  viscosity = inputs.mu;
  is_saturated = inputs.saturated;
}

RStdata::~RStdata()
{
  delete Porosity;
  delete KappaCCavg;
  delete Lambda;
  delete PCapParams;
  delete SpecificStorage;
  delete Source;
  delete KappaCCdir;
  delete CoeffCC;

  delete Rhs;
  delete Alpha;
  
  RichardCoefs.clear();
  DarcyVelocity.clear();
  KappaEC.clear();
}

void
RStdata::SetConstantValue(MFTower& mft,
                          Real     val,
                          int      dComp,
                          int      nComp)
{
  for (int lev=0; lev<nLevs; ++lev) {
    int nGrow = mft[lev].nGrow();
    mft[lev].setVal(val,dComp,nComp,nGrow);
  }
}

void
RStdata::SetPCapParams()
{
  Real aval = inputs.alpha * Pa_per_ATM;
  SetConstantValue(*PCapParams,3,0,1);
  SetConstantValue(*PCapParams,inputs.m,1,1);
  SetConstantValue(*PCapParams,aval,2,1);
  SetConstantValue(*PCapParams,inputs.Sr,3,1);
}

void
RStdata::SetPorosity()
{
  SetConstantValue(*Porosity,inputs.phi,0,1);
}

void
RStdata::SetSpecificStorage()
{
  SetConstantValue(*SpecificStorage,inputs.specific_storage,0,1);
}

void
RStdata::SetUpMemory(NLScontrol& nlsc)
{
  RSdata::SetUpMemory(nlsc);

  RhoSatOld = new MFTower(layout,IndexType(IntVect::TheZeroVector()),1,1,nLevs);
  RhoSatNew = new MFTower(layout,IndexType(IntVect::TheZeroVector()),1,1,nLevs);
  Pnew = new MFTower(layout,IndexType(IntVect::TheZeroVector()),1,1,nLevs);
  Pold = new MFTower(layout,IndexType(IntVect::TheZeroVector()),1,1,nLevs);

  Porosity = new MFTower(layout,IndexType(IntVect::TheZeroVector()),1,1,nLevs); SetPorosity();
  KappaCCavg = new MFTower(layout,IndexType(IntVect::TheZeroVector()),1,1,nLevs);
  Lambda = new MFTower(layout,IndexType(IntVect::TheZeroVector()),1,1,nLevs);
  PCapParams = new MFTower(layout,IndexType(IntVect::TheZeroVector()),4,1,nLevs); SetPCapParams();
  SpecificStorage = new MFTower(layout,IndexType(IntVect::TheZeroVector()),1,1,nLevs); SetSpecificStorage();
  Source = new MFTower(layout,IndexType(IntVect::TheZeroVector()),1,1,nLevs);
  KappaCCdir = new MFTower(layout,IndexType(IntVect::TheZeroVector()),BL_SPACEDIM,1,nLevs);
  CoeffCC    = new MFTower(layout,IndexType(IntVect::TheZeroVector()),BL_SPACEDIM,1,nLevs);

  Rhs = new MFTower(layout,IndexType(IntVect::TheZeroVector()),1,1,nLevs);
  Alpha = new MFTower(layout,IndexType(IntVect::TheZeroVector()),1,1,nLevs);
  
  RichardCoefs.resize(BL_SPACEDIM,PArrayManage);
  DarcyVelocity.resize(BL_SPACEDIM,PArrayManage);
  KappaEC.resize(BL_SPACEDIM,PArrayManage);
  for (int d=0; d<BL_SPACEDIM; ++d) {
    DarcyVelocity.set(d, new MFTower(layout,IndexType(BoxLib::BASISV(d)),1,0,nLevs));
    RichardCoefs.set(d, new MFTower(layout,IndexType(BoxLib::BASISV(d)),1,0,nLevs));
    KappaEC.set(d, new MFTower(layout,IndexType(BoxLib::BASISV(d)),1,0,nLevs));
  }
}

void
RStdata::ResetRhoSat()
{
  // No-op here since the data is owned by this clases
}

void
RStdata::SetInflowVelocity(PArray<MFTower>& velocity,
                           Real             time)
{
  int dir = 1;
  for (int lev=0; lev<nLevs; ++lev) {
    const Box bdryHi = BoxLib::bdryHi(layout.GeomArray()[lev].Domain(),dir,1);
    for (MFIter mfi(velocity[dir][lev]); mfi.isValid(); ++mfi) {
      Box ovlp = velocity[dir][lev][mfi].box() & bdryHi;
      if (ovlp.ok()) {
        velocity[dir][lev][mfi].setVal(inputs.inflow_velocity,ovlp,0,1);
      }
    }
  }
}

bool
RStdata::IsNewTime(Real time) const
{
  BL_ASSERT(old_time < new_time);
  Real teps = 1.e-6*(new_time - old_time);
  bool is_new_time = false;
  if (std::abs(time - new_time) <= teps) {
    is_new_time = true;
  }
  else if (std::abs(time - old_time) > teps)
  {
    BoxLib::Abort("Invalid time");
  }
  return is_new_time;
}

#  if defined(BL_FORT_USE_UPPERCASE)
#    define FORT_FILCC       FILCC
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_FILCC       filcc
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_FILCC       filcc_
#  endif

#include <ArrayLim.H>

extern "C" {
    void FORT_FILCC (const Real * q, ARLIM_P(q_lo), ARLIM_P(q_hi),
                     const int * domlo, const int * domhi,
                     const Real * dx_crse, const Real * xlo, 
                     const int * bc);
}
    
void
RStdata::FillStateBndry (MFTower& press,
                         Real time)
{
  MFTower& P = (IsNewTime(time) ? (*Pnew) : (*Pold) );
  int src_comp = 0;
  int num_comp = 1;
  int dir = BL_SPACEDIM-1;
  for (int lev=0; lev<nLevs; ++lev) {
    const Geometry& geom = layout.GeomArray()[lev];
    const Box& domain = geom.Domain();
    const int* domlo = domain.loVect();
    const int* domhi = domain.hiVect();
    const Real* dx = geom.CellSize();
    MultiFab& mf = P[lev];
    const BoxArray& ba = mf.boxArray();
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
      int         i       = mfi.index();
      RealBox     gridloc = RealBox(mfi.validbox(),geom.CellSize(),geom.ProbLo());
      FArrayBox&  fab     = mf[mfi];
      const int*  flo     = fab.loVect();
      const int*  fhi     = fab.hiVect();
      const Real* xlo     = gridloc.lo();
      FORT_FILCC(fab.dataPtr(), ARLIM(flo), ARLIM(fhi),
                 domlo, domhi, dx, xlo,
                 pressure_bc.vect());
    }
    mf.FillBoundary(src_comp,num_comp);
    geom.FillPeriodicBoundary(mf, src_comp, num_comp);

    Box gbox = domain;
    for (int d=0; d<BL_SPACEDIM; ++d) {
      if (d!=dir) gbox.grow(d,P.NGrow());
    }
    const Box adjCellLo = BoxLib::adjCellLo(gbox,dir,1);
    for (MFIter mfi(P[lev]); mfi.isValid(); ++mfi) {
      Box ovlp = P[lev][mfi].box() & adjCellLo;
      if (ovlp.ok()) {
        P[lev][mfi].setVal(inputs.Pwt,ovlp,0,1);
      }
    }
  }
}

void
RStdata::calcInvPressure (MFTower&       N,
                          const MFTower& P) const
{
  Real rho_loc = GetDensity()[0];
  for (int lev=0; lev<nLevs; ++lev) {
    int nGrow = 1;
    BL_ASSERT(N[lev].nGrow()>=nGrow);
    BL_ASSERT(P[lev].nGrow()>=nGrow);
    for (MFIter mfi(N[lev]); mfi.isValid(); ++mfi) {
      FArrayBox& rhosat = N[lev][mfi];
      const FArrayBox& pressure = P[lev][mfi];
      const FArrayBox& pc_params = (*PCapParams)[lev][mfi];
      const Box& box = Box(mfi.validbox()).grow(nGrow);
      BL_ASSERT(rhosat.box().contains(box));
      BL_ASSERT(pressure.box().contains(box));
      for (IntVect iv=box.smallEnd(), End=box.bigEnd(); iv<=End; box.next(iv)) {
        Real p = pressure(iv,0);
        Real seff = 1;
        Real m_loc = pc_params(iv,1);
        Real a_loc = pc_params(iv,2);
        Real Sr_loc = pc_params(iv,3);
        if (p < 0) {
          seff = std::pow(1 + std::pow(-p*a_loc, 1./(1-m_loc)), -m_loc);
        }
        rhosat(iv,0) = rho_loc*(seff*(1-Sr_loc)+Sr_loc);
      }
    }
  }
}

void
RStdata::calcLambda (MFTower&       lbd,
                     const MFTower& N)
{
  Real mu_loc = GetViscosity()[0];
  Real rho_loc = GetDensity()[0];
  for (int lev=0; lev<nLevs; ++lev) {
    int nGrow = 1;
    BL_ASSERT(N[lev].nGrow()>=nGrow);
    BL_ASSERT(lbd[lev].nGrow()>=nGrow);
    for (MFIter mfi(N[lev]); mfi.isValid(); ++mfi) {
      const FArrayBox& rhosat = N[lev][mfi];
      FArrayBox& lambda = lbd[lev][mfi];
      const FArrayBox& pc_params = (*PCapParams)[lev][mfi];
      const Box& box = Box(mfi.validbox()).grow(nGrow);
      BL_ASSERT(rhosat.box().contains(box));
      BL_ASSERT(lambda.box().contains(box));
      for (IntVect iv=box.smallEnd(), End=box.bigEnd(); iv<=End; box.next(iv)) {
        Real m_loc = pc_params(iv,1);
        Real a_loc = pc_params(iv,2);
        Real Sr_loc = pc_params(iv,3);
        Real seff = (rhosat(iv,0)/rho_loc - Sr_loc)/(1 - Sr_loc);
        lambda(iv,0) = std::sqrt(seff) * std::pow(1-std::pow(1-std::pow(seff,1/m_loc),m_loc),2) / mu_loc;
      }
    }
  }
}

void
RStdata::calcRichardAlpha(MFTower&       alpha,
                          const MFTower& N,
                          Real           t)
{
  Real mu_loc = GetViscosity()[0];
  Real rho_loc = GetDensity()[0];
  for (int lev=0; lev<nLevs; ++lev) {
    int nGrow = 1;
    BL_ASSERT(N[lev].nGrow()>=nGrow);
    BL_ASSERT(alpha[lev].nGrow()>=nGrow);
    for (MFIter mfi(N[lev]); mfi.isValid(); ++mfi) {
      FArrayBox& afab = alpha[lev][mfi];
      const FArrayBox& Nfab = N[lev][mfi];
      const FArrayBox& phifab = (*Porosity)[lev][mfi];
      const FArrayBox& pc_params = (*PCapParams)[lev][mfi];
      const Box& box = Box(mfi.validbox()).grow(nGrow);
      BL_ASSERT(Nfab.box().contains(box));
      BL_ASSERT(phifab.box().contains(box));
      for (IntVect iv=box.smallEnd(), End=box.bigEnd(); iv<=End; box.next(iv)) {
        Real m_loc = pc_params(iv,1);
        Real a_loc = pc_params(iv,2);
        Real Sr_loc = pc_params(iv,3);
        Real seff = (Nfab(iv,0)/rho_loc - Sr_loc)/(1 - Sr_loc);

        seff = std::max(1.e-3,seff);

        Real n_loc = -1/m_loc;
        Real phi_loc = phifab(iv,0);

        Real sn = std::pow(seff,n_loc);
        Real dpc_dN = n_loc * (1-m_loc) * std::pow(sn-1,-m_loc) * sn / (seff * a_loc);

        afab(iv,0) = - phi_loc / dpc_dN * rho_loc * (1-inputs.Sr);
      }
    }
  }
}

Array<int>&
RStdata::rinflowBCLo()
{}

Array<int>&
RStdata::rinflowBCHi()
{}

const MFTower*
RStdata::GetKappaCCdir(Real t)
{
  if (t!=eval_time_for_KappaCCdir) {
    Real kappaval = inputs.kappa * Pa_per_ATM;
    for (int lev=0; lev<nLevs; ++lev) {
      int nGrow = (*KappaCCdir)[lev].nGrow();
      (*KappaCCdir)[lev].setVal(kappaval,0,BL_SPACEDIM,nGrow);
    }
    eval_time_for_KappaCCdir = t;
  }
  return KappaCCdir;
}

const MFTower*
RStdata::GetSource(Real t)
{
  if (t!=eval_time_for_source) {
    for (int lev=0; lev<nLevs; ++lev) {
      int nGrow = (*Source)[lev].nGrow();
      (*Source)[lev].setVal(0,0,1,nGrow);
    }
    eval_time_for_source = t;
  }
  return Source;
}

const PArray<MFTower>&
RStdata::GetKappaEC(Real t)
{
  static bool do_harm = true;
  if (t!=eval_time_for_KappaEC) {
    GetKappaCCdir(t);
    MFTower::CCtoECavg(KappaEC,*KappaCCdir,1.0,0,0,-1,do_harm,nLevs);
    eval_time_for_KappaEC = t;
  }
  return KappaEC;
}

