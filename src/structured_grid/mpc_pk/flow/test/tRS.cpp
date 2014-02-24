#include <RichardSolver.H>
#include <ParmParse.H>
#include <VisMF.H>
#include <Utility.H>

static Real inflow_velocity = -1.1091e-10;
static Real Pwt = 0;
static Real dt = 1; // Initial
static int Niter = 60;
static std::string plt_name = "";
static bool verbose = false;

// Soil
static Real kappa = 2.87e-13;
static Real phi = 0.38;
static Real alpha = 3.02e-4;
static Real Sr = 0.354;
static Real m = 0.291;
static Real specific_storage = 0;

// Water
static Real rho = 998.2;
static Real mu = 0.001005;
static Real g = 9.81117;

static bool saturated = false;

static Real H = 107.52;
static Real W = H/16;
// static int Nx = 16;
// static int Ny = Nx*16;
// static int Nz = 64;
static int Nx = 2;
static int Ny = Nx*16;
static int Nz = 64;

static Real Pa_per_ATM = 101325;

struct RStdata
  : public RSdata
{
  RStdata(int slev, int nlevs, Layout& layout, NLScontrol& nlsc);
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
  virtual void SetDensity();
  virtual void SetGravity();
  virtual void SetViscosity();
  virtual void SetIsSaturated();
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
};

void
RStdata::SetDensity()
{
  density.resize(1,rho);
}

void
RStdata::SetViscosity()
{
  viscosity.resize(1,mu);
}

void
RStdata::SetGravity()
{
  gravity.resize(BL_SPACEDIM,0);
  gravity[BL_SPACEDIM-1] = g / Pa_per_ATM;
}

void
RStdata::SetIsSaturated()
{
  is_saturated = saturated;
}

RStdata::RStdata(int slev, int nlevs, Layout& layout, NLScontrol& nlsc)
  : RSdata(slev,nlevs,layout,nlsc)
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

  delete RhoSatOld;
  delete RhoSatNew;
  delete Pnew;
  delete Pold;
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
  Real aval = alpha * Pa_per_ATM;
  SetConstantValue(*PCapParams,   3,0,1);
  SetConstantValue(*PCapParams,   m,1,1);
  SetConstantValue(*PCapParams,aval,2,1);
  SetConstantValue(*PCapParams,  Sr,3,1);
}

void
RStdata::SetPorosity()
{
  SetConstantValue(*Porosity,phi,0,1);
}

void
RStdata::SetSpecificStorage()
{
  SetConstantValue(*SpecificStorage,specific_storage,0,1);
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
        velocity[dir][lev][mfi].setVal(inflow_velocity,ovlp,0,1);
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
      
#ifdef _OPENMP
#pragma omp parallel for
#endif
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
        P[lev][mfi].setVal(Pwt,ovlp,0,1);
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

        afab(iv,0) = - phi_loc / dpc_dN * rho_loc * (1-Sr);
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
    Real kappaval = kappa * Pa_per_ATM;
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


//
// Components are  Interior, Inflow, Outflow, Symmetry, SlipWall, NoSlipWall.
//
static int press_bc[] =
{
  INT_DIR, FOEXTRAP, EXT_DIR, REFLECT_EVEN, FOEXTRAP, FOEXTRAP
};

static
void
set_pressure_bc (BCRec&       bc,
                 const BCRec& phys_bc)
{
  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();
  for (int i = 0; i < BL_SPACEDIM; i++) {
    bc.setLo(i,press_bc[lo_bc[i]]);
    bc.setHi(i,press_bc[hi_bc[i]]);
  }
}

static void
GradFill (MFTower&               mft,
          int                    nLevs,
          const Array<Real>&     grad,
          const Array<Geometry>& geoms)
{
  BL_ASSERT(nLevs<=mft.NumLevels());
  for (int lev=0; lev<nLevs; ++lev) {
    int nGrow = mft.NGrow();
    BL_ASSERT(mft[lev].nGrow()>=nGrow);
    for (MFIter mfi(mft[lev]); mfi.isValid(); ++mfi) {
      FArrayBox& fab = mft[lev][mfi];
      const Box& box = Box(mfi.validbox()).grow(nGrow);
      for (IntVect iv=box.smallEnd(), End=box.bigEnd(); iv<=End; box.next(iv)) {
        Real val = 0;
        for (int d=0; d<BL_SPACEDIM; ++d) {
          val += grad[d]*(iv[d]+0.5)*geoms[lev].CellSize()[d];
        }
        fab(iv,0) = val;
      }
    }
  }
}

static void
MultMFT (MFTower&           mft,
         Real               a,
         int                sComp = 0,
         int                nComp = -1,
         int                nLevs = -1)
{
  nLevs = (nLevs==-1 ? mft.NumLevels() : std::max(nLevs,mft.NumLevels()));
  nComp = (nComp==-1 ? mft.NComp() : std::max(nComp,mft.NComp()));
  for (int lev=0; lev<nLevs; ++lev) {
    int nGrow = mft.NGrow();
    BL_ASSERT(mft[lev].nGrow()>=nGrow);
    mft[lev].mult(a,sComp,nComp);
  }
}

static void
PlusMFT (MFTower&           mft,
         Real               a,
         int                sComp = 0,
         int                nComp = -1,
         int                nLevs = -1)
{
  nLevs = (nLevs==-1 ? mft.NumLevels() : std::max(nLevs,mft.NumLevels()));
  nComp = (nComp==-1 ? mft.NComp() : std::max(nComp,mft.NComp()));
  for (int lev=0; lev<nLevs; ++lev) {
    int nGrow = mft.NGrow();
    BL_ASSERT(mft[lev].nGrow()>=nGrow);
    mft[lev].plus(a,sComp,nComp);
  }
}

static void
CopyMFT (MFTower&           dMFT,
         const MFTower&     sMFT,
         int                sComp = 0,
         int                dComp = 0,
         int                nComp = -1,
         int                nLevs = -1)
{
  nLevs = (nLevs==-1 ? sMFT.NumLevels() : std::max(nLevs,sMFT.NumLevels()));
  BL_ASSERT(nLevs <= dMFT.NumLevels());
  nComp = (nComp==-1 ? sMFT.NComp() : std::max(nComp,sMFT.NComp()));
  BL_ASSERT(nComp <= dMFT.NComp());
  for (int lev=0; lev<nLevs; ++lev) {
    BL_ASSERT(sMFT[lev].boxArray()==dMFT[lev].boxArray());
    int nGrow = sMFT.NGrow();
    BL_ASSERT(dMFT[lev].nGrow()>=nGrow);
    MultiFab::Copy(dMFT[lev],sMFT[lev],sComp,dComp,nComp,nGrow);
  }
}

int
main (int   argc,
      char* argv[])
{
  BoxLib::Initialize(argc,argv);

  ParmParse pp;

#ifdef BL_USE_PETSC
  std::string petsc_help = "Amanzi-S passthrough access to PETSc help option\n";
  std::string petsc_file_str = "Petsc_Options_File";
  std::string petsc_options_file = ".petsc";
  if (pp.countval(petsc_file_str.c_str())) {
    pp.get(petsc_file_str.c_str(),petsc_options_file);
  }
  PetscInitialize(&argc,&argv,petsc_options_file.c_str(),petsc_help.c_str());
#endif

  // Build simple single-level Layout
  int nLevs = 1;
  Array<IntVect>  refRatio_array;
  Array<BoxArray> grid_array(nLevs);
  Array<Geometry> geom_array(nLevs);

  Box domain(IntVect(D_DECL(0,0,0)),IntVect(D_DECL(Nx-1,Ny-1,Nz-1)));
  grid_array[0] = BoxArray(domain);

  Array<Real> problo(BL_SPACEDIM,0);
  Array<Real> probhi(BL_SPACEDIM);
  for (int d=0; d<BL_SPACEDIM; ++d) {
    probhi[d] = (d==BL_SPACEDIM-1 ? H : W);
  }
  RealBox rb(problo.dataPtr(),probhi.dataPtr());
  int coord = 0; // Cartesian
  Array<int> is_per(BL_SPACEDIM,0); // Not periodic
  geom_array[0] = Geometry(domain,&rb,coord,is_per.dataPtr());
  
  Layout layout(refRatio_array,grid_array,geom_array,nLevs);

  NLScontrol nlsc;

  nlsc.max_nl_iterations=15;
  nlsc.time_step_reduction_factor=0.8;
  nlsc.time_step_retry_factor=0.5;
  //nlsc.max_num_consecutive_success=3;
  nlsc.max_num_consecutive_success=1;
  //nlsc.min_nl_iterations_for_dt=10;
  nlsc.min_nl_iterations_for_dt=13;
  nlsc.time_step_increase_factor=1.4;
  nlsc.ls_acceptance_factor=5;
  nlsc.monitor_line_search=0;
  RStdata rs_data(0,nLevs,layout,nlsc);
  rs_data.upwind_krel=1;
  rs_data.semi_analytic_J=true;

  pp.query("inflow_velocity",inflow_velocity);
  pp.query("Niter",Niter);
  pp.query("dt",dt);
  pp.query("plt_name",plt_name);
  pp.query("verbose",verbose);

  pp.query("do_upwind",rs_data.upwind_krel);

  // Boundary conditions:
  // Components are  0:Interior, 1:Inflow, 2:Outflow, 3:Symmetry, 4:SlipWall, 5:NoSlipWall.
  //                  INT_DIR,   FOEXTRAP,  EXT_DIR, REFLECT_EVEN, FOEXTRAP,    FOEXTRAP
  static Array<int> lo_bc(BL_SPACEDIM,3), hi_bc(BL_SPACEDIM,3);
  lo_bc[BL_SPACEDIM-1] = 2;
  hi_bc[BL_SPACEDIM-1] = 1;

  BCRec phys_bc;
  for (int i = 0; i < BL_SPACEDIM; i++) {
    phys_bc.setLo(i,lo_bc[i]);
    phys_bc.setHi(i,hi_bc[i]);
  }

  BCRec pressure_bc; set_pressure_bc(pressure_bc,phys_bc);
  rs_data.SetPressureBC(pressure_bc);
  rs_data.SetUpMemory(nlsc);

  //RichardSolver rs(rs_data,nlsc);
  RichardSolver* rs = new RichardSolver(rs_data,nlsc);

  MFTower& Pold = *(rs_data.Pold);
  MFTower& Pnew = *(rs_data.Pnew);
  MFTower& RSold = *(rs_data.RhoSatOld);
  MFTower& RSnew = *(rs_data.RhoSatNew);

  Array<Real> grad(BL_SPACEDIM,0); grad[BL_SPACEDIM-1] = - rho*g / Pa_per_ATM;

  rs_data.old_time = 0;
  rs_data.new_time = dt;
  int step = 0;
  rs_data.SetCurrentTimestep(step);

  GradFill(Pold,nLevs,grad,geom_array);
  GradFill(Pnew,nLevs,grad,geom_array);

  rs_data.FillStateBndry(Pnew,rs_data.new_time);
  rs_data.calcInvPressure(RSnew,Pnew);

  Array<MFTower*> output_set;
  Array<std::string> output_names;
  output_set.push_back(&Pnew); output_names.push_back("Pnew");
  output_set.push_back(&RSnew); output_names.push_back("RSnew");

  // Write initial data
  if (plt_name!="") {
    std::string outfile=BoxLib::Concatenate(plt_name,step,3);
    MultMFT(RSnew,1/rs_data.GetDensity()[0]);
    MultMFT(Pnew,Pa_per_ATM);
    PlusMFT(Pnew,Pa_per_ATM);
    MFTower::WriteSet(outfile,output_set,output_names,rs_data.new_time);
    PlusMFT(Pnew,-Pa_per_ATM);
    MultMFT(Pnew,1/Pa_per_ATM);
    MultMFT(RSnew,rs_data.GetDensity()[0]);
  }

  Real dt_new;
  NLSstatus ret;
  for ( ; step < Niter; ++step) {
    bool cont = true;
    int retCode = -1;
    while (cont && retCode < 0) {
      rs_data.SetCurrentTimestep(step+1);
      rs_data.new_time = rs_data.old_time+dt;
      if (verbose && ParallelDescriptor::IOProcessor()) {
        std::cout << "................ attempting dt = " << dt << " step " << step << std::endl;
      }
      retCode = rs->Solve(rs_data.old_time,rs_data.new_time,step,nlsc);
      if (retCode > 0) {
        if (verbose && ParallelDescriptor::IOProcessor()) {
          std::cout << "................ SUCCEEDED dt = " << dt << std::endl;
        }
        ret = NLSstatus::NLS_SUCCESS;
      } 
      else {
        if (verbose && ParallelDescriptor::IOProcessor()) {
          std::cout << "................ FAILED dt = " << dt << std::endl;
        }
        if (retCode == -3 || retCode == 0) {
          ret = NLSstatus::NLS_LINEAR_FAIL;
        }
        else {
          ret = NLSstatus::NLS_NONLINEAR_FAIL;
        }
        CopyMFT(Pnew,Pold);
      }
      cont = nlsc.AdjustDt(dt,ret,dt_new);
      dt = dt_new;
    }

    rs_data.FillStateBndry(Pnew,rs_data.new_time);
    rs_data.calcInvPressure(RSnew,Pnew);

    if (plt_name!="") {
      std::string outfile_loc=BoxLib::Concatenate(plt_name,step,3);
      MultMFT(RSnew,1/rs_data.GetDensity()[0]);
      MultMFT(Pnew,Pa_per_ATM);
      PlusMFT(Pnew,Pa_per_ATM);
      MFTower::WriteSet(outfile_loc,output_set,output_names,rs_data.new_time);
      PlusMFT(Pnew,-Pa_per_ATM);
      MultMFT(Pnew,1/Pa_per_ATM);
      MultMFT(RSnew,rs_data.GetDensity()[0]);
    }

    CopyMFT(Pold,Pnew);
  } 

  delete rs;

  // Must clear prior to PetscFinalize to explicitly delete petsc data structures, or 
  // the PetscFinalize will do it, and the Layout destructor will then generate a seg fault.
  layout.Clear();

#ifdef BL_USE_PETSC
    PetscFinalize();
#endif

  ParallelDescriptor::EndParallel();
  return 0;
}
