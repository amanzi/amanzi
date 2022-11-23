#include <RStdata.H>

Real RStdata::Pa_per_ATM = 101325;
static bool saturated = false;

RStdata::RStdata(int slev, int nlevs, Layout& layout, NLScontrol& nlsc, RStstruct& _inputs, const RockManager* rm)
  : RSdata(slev,nlevs,layout,nlsc,rm), inputs(_inputs)
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

  bool ignore_mixed = true;
  int nLevs = layout.NumLevels();
  materialID.resize(nLevs,PArrayManage);
  for (int lev=0; lev<nLevs; ++lev) {
    const BoxArray& ba = layout.GridArray()[lev];
    materialID.set(lev,new iMultiFab(ba,1,3));
    rock_manager->GetMaterialID(lev,materialID[lev],materialID[lev].nGrow(),ignore_mixed);
  }
}

RStdata::~RStdata()
{
  delete RhoSatOld;
  delete RhoSatNew;
  delete Pnew;
  delete Pold;

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

static std::string CP_model_None = "None";
static std::string CP_model_vG = "VanGenuchten";
static std::string CP_model_BC = "BrooksCorey";

// FIXME: Pasted in from RockManager.cpp
static int CPL_MODEL_ID = 0;

static int VG_M                   = 1;
static int VG_ALPHA               = 2;
static int VG_SR                  = 3;
static int VG_ELL                 = 4;
static int VG_KR_MODEL_ID         = 5;
static int VG_KR_SMOOTHING_MAX_PC = 6;

static int BC_LAMBDA              = 1;
static int BC_ALPHA               = 2;
static int BC_SR                  = 3;
static int BC_ELL                 = 4;
static int BC_KR_MODEL_ID         = 5;
static int BC_KR_SMOOTHING_MAX_PC = 6;

void
RStdata::SetPCapParams(Real t)
{
  if (rock_manager->CanDerive("capillary_pressure")) {
    int nComp = rock_manager->NComp("capillary_pressure");
    for (int lev=0; lev<nLevs; ++lev) {
      int nGrow = (*PCapParams)[lev].nGrow();
      MultiFab pcp((*PCapParams)[lev].boxArray(),nComp,nGrow);
      bool ret = rock_manager->GetProperty(t,lev,pcp,"capillary_pressure",0,nGrow);
      if (!ret) BoxLib::Abort("Failed to build capillary_pressure");
      for (MFIter mfi(pcp); mfi.isValid(); ++mfi) {
        FArrayBox& res = (*PCapParams)[lev][mfi];
        const FArrayBox& pc_params = pcp[mfi];
        const Box& box = pc_params.box();
        for (IntVect iv(box.smallEnd()), End=box.bigEnd(); iv <= End; box.next(iv)) {
          bool is_vG = rock_manager->Is_CP_model_XX(pc_params(iv,CPL_MODEL_ID),CP_model_vG);
          bool is_BC = rock_manager->Is_CP_model_XX(pc_params(iv,CPL_MODEL_ID),CP_model_BC);
          if (is_vG) {
            res(iv,0) = 0;
            res(iv,1) = pc_params(iv,VG_M);
            res(iv,2) = pc_params(iv,VG_ALPHA) / Pa_per_ATM;
            res(iv,3) = pc_params(iv,VG_SR);
            res(iv,4) = pc_params(iv,VG_KR_SMOOTHING_MAX_PC);
          }
          else if (is_BC) {
            res(iv,0) = 0;
            res(iv,1) = pc_params(iv,BC_LAMBDA);
            res(iv,2) = pc_params(iv,BC_ALPHA) / Pa_per_ATM;
            res(iv,3) = pc_params(iv,BC_SR);
            res(iv,4) = pc_params(iv,BC_KR_SMOOTHING_MAX_PC);
          }
          else {
            BoxLib::Abort("Unknown PC model");
          }
        }
      }
    }
  }
}

void
RStdata::SetPorosity(Real t)
{
  for (int lev=0; lev<nLevs; ++lev) {
    int nGrow = (*Porosity)[lev].nGrow();
    bool ret = rock_manager->GetProperty(t,lev,(*Porosity)[lev],"porosity",0,nGrow);
    if (!ret) BoxLib::Abort("Failed to build porosity");
  }
}

void
RStdata::SetSpecificStorage(Real t)
{
  for (int lev=0; lev<nLevs; ++lev) {
    int nGrow = (*SpecificStorage)[lev].nGrow();
    bool ret = rock_manager->GetProperty(t,lev,(*SpecificStorage)[lev],"specific_storage",0,nGrow);
    if (!ret) {
      (*SpecificStorage)[lev].setVal(0);
    }
  }
}

void
RStdata::SetUpMemory(NLScontrol& nlsc)
{
  RSdata::SetUpMemory(nlsc);

  Real t = 0;

  RhoSatOld = new MFTower(layout,IndexType(IntVect::TheZeroVector()),1,1,nLevs);
  RhoSatNew = new MFTower(layout,IndexType(IntVect::TheZeroVector()),1,1,nLevs);
  Pnew = new MFTower(layout,IndexType(IntVect::TheZeroVector()),1,1,nLevs);
  Pold = new MFTower(layout,IndexType(IntVect::TheZeroVector()),1,1,nLevs);

  Porosity = new MFTower(layout,IndexType(IntVect::TheZeroVector()),1,1,nLevs); SetPorosity(t);
  KappaCCavg = new MFTower(layout,IndexType(IntVect::TheZeroVector()),1,1,nLevs);
  if (!is_saturated) {
    Lambda = new MFTower(layout,IndexType(IntVect::TheZeroVector()),1,1,nLevs);
    PCapParams = new MFTower(layout,IndexType(IntVect::TheZeroVector()),5,1,nLevs); SetPCapParams(t);
  }
  SpecificStorage = new MFTower(layout,IndexType(IntVect::TheZeroVector()),1,1,nLevs); SetSpecificStorage(t);
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
  int dir = BL_SPACEDIM-1;
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
  int src_comp = 0;
  int num_comp = 1;
  int dir = BL_SPACEDIM-1;
  for (int lev=0; lev<nLevs; ++lev) {
    const Geometry& geom = layout.GeomArray()[lev];
    const Box& domain = geom.Domain();
    const int* domlo = domain.loVect();
    const int* domhi = domain.hiVect();
    const Real* dx = geom.CellSize();
    MultiFab& mf = press[lev];
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

    Box gbox = domain;
    for (int d=0; d<BL_SPACEDIM; ++d) {
      if (d!=dir) gbox.grow(d,press.NGrow());
    }
    const Box adjCellLo = BoxLib::adjCellLo(gbox,dir,1);
    for (MFIter mfi(press[lev]); mfi.isValid(); ++mfi) {
      Box ovlp = press[lev][mfi].box() & adjCellLo;
      if (ovlp.ok()) {
        press[lev][mfi].setVal(inputs.Pwt,ovlp,0,1);
      }
    }
  }
}

void
RStdata::calcInvPressure (MFTower&       N,
                          const MFTower& P,
                          Real           time,
                          int            sComp,
                          int            dComp,
                          int            nGrow) const
{
  //
  // Pcap = Pgas - Pwater, then get N=s.rho from Pcap(s)^{-1}
  //
  BL_ASSERT(N.NGrow() >= nGrow  && P.NGrow() >= nGrow);
  BL_ASSERT(N.NComp() >= dComp && P.NComp() >= sComp);
  Real rho_loc = GetDensity()[0];

  for (int lev=0; lev<nLevs; ++lev) {
    BL_ASSERT(N[lev].boxArray() == P[lev].boxArray());
    MultiFab pc(P[lev].boxArray(),1,nGrow);
    pc.setVal(inputs.Pwt,0,1,nGrow);
    MultiFab::Subtract(pc,P[lev],0,0,1,nGrow);
    rock_manager->InverseCapillaryPressure(pc,materialID[lev],time,N[lev],0,dComp,nGrow);    
    N[lev].mult(rho_loc,dComp,1,nGrow);
    if (nGrow > 0) {
      N[lev].FillBoundary(dComp,1);
    }
  }
}

static Real vgKr(Real seff, Real m, Real ell) {
  return std::pow(seff, ell) * std::pow(1-std::pow(1-std::pow(seff,1/m),m),2);
}

void
RStdata::calcLambda (MFTower&       Lambda,
                     const MFTower& N,
                     Real           time,
                     int            sComp,
                     int            dComp,
                     int            nGrow) const
{
  Real mu_loc = GetViscosity()[0];
  Real rho_loc = GetDensity()[0];

  for (int lev=0; lev<nLevs; ++lev) {
    BL_ASSERT(N[lev].boxArray() == Lambda[lev].boxArray());
    MultiFab sat(N[lev].boxArray(),1,nGrow);
    MultiFab::Copy(sat,N[lev],0,0,1,nGrow);
    sat.mult(1/rho_loc,0,1,nGrow);
    rock_manager->RelativePermeability(sat,materialID[lev],time,Lambda[lev],0,dComp,nGrow);
    Lambda[lev].mult(1/mu_loc,dComp,1,nGrow);
    if (nGrow > 0) {
      Lambda[lev].FillBoundary(dComp,1);
    }
  }
}

void
RStdata::calcRichardAlpha (MFTower&       Alpha,
                           const MFTower& N,
                           Real           time,
                           int            sComp,
                           int            dComp,
                           int            nGrow) const
{
  BL_ASSERT(N.NGrow() >= nGrow); // Assumes that boundary cells have been properly filled
  BL_ASSERT(Alpha.NGrow() >= nGrow); // Fill boundary cells (in F)

  Real rho_loc = GetDensity()[0];
  for (int lev=0; lev<nLevs; ++lev) {
    const BoxArray& grids = N[lev].boxArray();
    MultiFab sat(grids,1,nGrow);
    MultiFab::Copy(sat,N[lev],0,0,1,nGrow);
    sat.mult(1/rho_loc,0,1,nGrow);
    rock_manager->DInverseCapillaryPressure(sat,materialID[lev],time,Alpha[lev],sComp,dComp,nGrow);
    Alpha[lev].mult(-rho_loc,dComp,1,nGrow);
    MultiFab::Multiply(Alpha[lev],(*Porosity)[lev],0,dComp,1,nGrow);
    if (nGrow > 0) {
      Alpha[lev].FillBoundary(dComp,1);
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
    for (int lev=0; lev<nLevs; ++lev) {
      int nGrow = (*KappaCCdir)[lev].nGrow();
      bool ret = rock_manager->GetProperty(t,lev,(*KappaCCdir)[lev],"permeability",0,nGrow);
      if (!ret) BoxLib::Abort("Failed to build permeability");
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

