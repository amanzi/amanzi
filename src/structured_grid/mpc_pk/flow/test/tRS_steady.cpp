// This is needed to define some external symbols for linking. 
// It's kind of silly to do things this way. -JNJ
#include <VerboseObject_objs.hh> 

#include <RichardSolver.H>
#include <RStdata.H>

#include <ParmParse.H>
#include <VisMF.H>
#include <Utility.H>

// static Real H = 107.52;
// static Real W = H/16;
// static int NyFAC = 16;

static Real dt_init = 1;
static int Niter = 60;
static std::string plt_name = "";
static bool verbose = false;

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
      const Box box = Box(mfi.validbox()).grow(nGrow);
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

  int nLevs = 1;
  Array<IntVect>  refRatio_array;
  Array<BoxArray> grid_array(nLevs);
  Array<Geometry> geom_array(nLevs);

  NLScontrol nlsc;

  nlsc.max_nl_iterations=15;
  nlsc.time_step_reduction_factor=0.8;
  nlsc.time_step_retry_factor=0.5;
  nlsc.max_num_consecutive_success=1;
  nlsc.min_nl_iterations_for_dt=13;
  nlsc.time_step_increase_factor=1.4;
  nlsc.ls_acceptance_factor=5;
  nlsc.monitor_line_search=0;

  RStstruct inputs;
  inputs.rho.resize(1,998.2);
  inputs.mu.resize(1,0.001005);
  inputs.g = 9.81117 / RStdata::Pa_per_ATM;
  inputs.saturated = false;
  inputs.inflow_velocity = -1.1091e-10;
  inputs.Pwt = 1;

  pp.query("inflow_velocity",inputs.inflow_velocity);
  pp.query("Niter",Niter);
  pp.query("plt_name",plt_name);
  pp.query("verbose",verbose);
  pp.query("dt",dt_init);

  bool abort_on_nl_fail = false; pp.query("abort_on_nl_fail",abort_on_nl_fail);

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
  Array<Real> grad(BL_SPACEDIM,0); grad[BL_SPACEDIM-1] = - inputs.rho[0]*inputs.g;

  int Nlev = 4; pp.query("Nlev",Nlev);
  if (Nlev<3) {
    BoxLib::Abort("Must have Nlev>=3 in order to compute convergence rate");
  }
  int Ny0 = 32; pp.query("Ny0",Ny0);
  Array<int> Ny(Nlev,Ny0);

  PArray<FArrayBox> Psave(Nlev);
  for (int lev=1; lev<Nlev; ++lev) {
    Ny[lev] = Ny[lev-1]*2;
  }

  std::string krel_upwind_method = "upwind-darcy_velocity"; pp.query("krel_upwind_method",krel_upwind_method);

  for (int lev = 0; lev < Nlev; ++lev) {
    Real dt = dt_init;
    int Nx = 2;
    Box domain(IntVect(D_DECL(0,0,0)),IntVect(D_DECL(Nx-1,Ny[lev]-1,1)));
    grid_array[0] = BoxArray(domain);
    geom_array[0] = Geometry(domain);

    Array<Real> plo(BL_SPACEDIM), phi(BL_SPACEDIM);
    for (int i = 0; i < BL_SPACEDIM; ++i) {
      plo[i] = (Geometry::ProbLo())[i];
      phi[i] = (Geometry::ProbLo())[i];
    }
    RegionManager rm(plo,phi);
  
    Layout layout(refRatio_array,grid_array,geom_array,nLevs);

    int nGrow = 0;
    RockManager rockManager(&rm);
    rockManager.FinalizeBuild(geom_array,refRatio_array,nGrow,false);

    RStdata rs_data(0,nLevs,layout,nlsc,inputs,&rockManager);
    rs_data.rel_perm_method = krel_upwind_method;
    rs_data.semi_analytic_J=true; pp.query("semi_analytic",rs_data.semi_analytic_J);

    rs_data.SetPressureBC(pressure_bc);
    rs_data.SetUpMemory(nlsc);

    RichardSolver* rs = new RichardSolver(rs_data,nlsc);
    if (pp.countval("record_file") ){
      std::string rf; pp.get("record_file",rf);
      rs->SetRecordFile(rf);
    }

    MFTower& Pold = *(rs_data.Pold);
    MFTower& Pnew = *(rs_data.Pnew);
    MFTower& RSold = *(rs_data.RhoSatOld);
    MFTower& RSnew = *(rs_data.RhoSatNew);

    rs_data.old_time = 0;
    rs_data.new_time = dt;
    int step = 0;
    rs_data.SetCurrentTimestep(step);

    GradFill(Pold,nLevs,grad,geom_array);
    GradFill(Pnew,nLevs,grad,geom_array);

    PlusMFT(Pnew,inputs.Pwt);
    PlusMFT(Pold,inputs.Pwt);
    rs_data.calcInvPressure(RSnew,Pnew,rs_data.new_time,0,0,0);

    Array<MFTower*> output_set;
    Array<std::string> output_names;
    output_set.push_back(&Pnew); output_names.push_back("Pnew");
    output_set.push_back(&RSnew); output_names.push_back("Snew");

    // Write initial data
    if (plt_name!="") {
      std::string outfile=BoxLib::Concatenate(plt_name,step,3);
      MultMFT(RSnew,1/rs_data.GetDensity()[0]); // S = N / rho
      MFTower::WriteSet(outfile,output_set,output_names,rs_data.new_time);
      MultMFT(RSnew,rs_data.GetDensity()[0]); // N = S.rho
    }

    Real dt_new;
    NLSstatus ret;
    for (step++ ; step <= Niter; ++step) {
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
	    std::cout << "................ FAILED dt = " << dt
		      << " return code = " << retCode << std::endl;
	  }
	  if (retCode == -3 || retCode == 0) {
	    ret = NLSstatus::NLS_LINEAR_FAIL;
	  }
	  else {
          ret = NLSstatus::NLS_NONLINEAR_FAIL;
	  }
	  if (abort_on_nl_fail) {
	    BoxLib::Abort("Nonlinear solver failed");
	  }
	  CopyMFT(Pnew,Pold);
	}
	cont = nlsc.AdjustDt(dt,ret,dt_new);
	dt = dt_new;
      }

      rs_data.FillStateBndry(Pnew,rs_data.new_time);
      rs_data.calcInvPressure(RSnew,Pnew,rs_data.new_time,0,0,0);

      if (plt_name!="") {
	std::string outfile_loc=BoxLib::Concatenate(plt_name,step,3);
	MultMFT(RSnew,1/rs_data.GetDensity()[0]); // S = N/rho
	MFTower::WriteSet(outfile_loc,output_set,output_names,rs_data.new_time);
	MultMFT(RSnew,rs_data.GetDensity()[0]); // N = S.rho
      }

      CopyMFT(Pold,Pnew);
    }

    Psave.set(lev, new FArrayBox(Pnew[0].boxArray()[0],1));
    Psave[lev].copy(Pnew[0][0]);

#if 0
    if (lev==Nlev-1) {
      Box box=Pnew[0].boxArray()[0];
      IntVect be = box.bigEnd(); be[0]=box.smallEnd()[0];
      box.setBig(be);
      Box ebox = box.surroundingNodes(1);
      FArrayBox u(ebox,1);
      u.copy(rs_data.DarcyVelocity[1][0][0]);
      std::cout << u << std::endl;
    }
#endif

    delete rs;

    // Must clear prior to PetscFinalize to explicitly delete petsc data structures, or 
    // the PetscFinalize will do it, and the Layout destructor will then generate a seg fault.
    layout.Clear();
  }

  Array<Array<Real> > norms(Nlev-1,Array<Real>(3));
  IntVect ref(IntVect(D_DECL(1,2,1)));
  for (int lev=0; lev<Nlev-1; ++lev) {
    const FArrayBox& ffab = Psave[lev+1];
    const FArrayBox& cfab = Psave[lev];
    const Box& fbox = ffab.box();
    const Box& cbox = cfab.box();
    FArrayBox cf(cbox,1); cf.setVal(0);
    if (Box(cbox).refine(ref) != fbox) {
      BoxLib::Abort("Boxes not correct");
    }
    for (IntVect ivc = cbox.smallEnd(); ivc <= cbox.bigEnd(); cbox.next(ivc)) {
      Box ivfbx = Box(ivc,ivc).refine(ref);
      for (IntVect ivf = ivfbx.smallEnd(); ivf <= ivfbx.bigEnd(); ivfbx.next(ivf)) {
	cf(ivc,0) += ffab(ivf,0);
      }
      cf(ivc,0) *= 1.0/ivfbx.numPts();
    }
    FArrayBox err(cbox,1);
    err.copy(cf);
    err.minus(cfab);
    norms[lev][0] = err.norm(0); // infinity norm (max norm).
    norms[lev][1] = err.norm(1) / err.box().numPts(); // 1-norm 
    // norms[lev][2] = err.norm(2) / err.box().numPts(); // 2-norm
    norms[lev][2] = norms[lev][1]; // missing implementation of 2-norm (use dot ???)
  }

  int retVal = 0;
  Array<Array<Real> > rates(Nlev-2,Array<Real>(3));
  for (int lev=0; lev<Nlev-2; ++lev) {
    for (int j=0; j<3; ++j) {
      rates[lev][j] = std::log(norms[lev][j] / norms[lev+1][j]) / std::log(2);
    }
    if (verbose) {
      std::cout << "Level " << lev << ": (r0,r1,r2) = " << rates[lev][0] << " " << rates[lev][1] << " " << rates[lev][2] << std::endl;
    }
    bool do_upwind = krel_upwind_method == "upwind-darcy_velocity";
    if (!do_upwind && (rates[lev][0]< 1.9 || rates[lev][1]< 1.9 || rates[lev][2]< 1.9)) retVal = 1; // If harmonic averaging used, we expect second order
    if (do_upwind &&  (rates[lev][0]< 0.9 || rates[lev][1]< 0.9 || rates[lev][2]< 0.9)) retVal = 1; // If upwind averaging used, we expect first order
  }

#ifdef BL_USE_PETSC
    PetscFinalize();
#endif

  if (retVal != 0) {
    BoxLib::Abort("Steady Richard differencing fails accuracy test");
  }
  ParallelDescriptor::EndParallel();

  return 0;
}
