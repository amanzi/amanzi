/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

// This is needed to define some external symbols for linking.
// It's kind of silly to do things this way. -JNJ
#include <VerboseObject_objs.hh>

#include <RichardSolver.H>
#include <RStdata.H>

#include <ParmParse.H>
#include <VisMF.H>
#include <Utility.H>

#if BL_SPACEDIM==2
static Real H = 107.52;
static Real W = H/16;
static int Nx = 2;
static int Ny = Nx*16;
#else
static Real H = 107.52;
static Real W = H/16;
static int Nx = 2;
static int Ny = 2;
static int Nz = Nx*16;
#endif

static Real dt = 1; // Initial
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

  nlsc.verbosity = verbose ? 5 : 0;
  nlsc.max_nl_iterations=15;
  nlsc.time_step_reduction_factor=0.8;
  nlsc.time_step_retry_factor=0.2;
  nlsc.max_num_consecutive_success=0;
  nlsc.min_nl_iterations_for_dt=10;
  nlsc.time_step_increase_factor=1.6;
  nlsc.ls_acceptance_factor=1.4;
  nlsc.monitor_line_search=0;
  nlsc.scale_soln_before_solve = 1;

  RStstruct inputs;
  inputs.rho.resize(1,998.2);
  inputs.mu.resize(1,0.001005);
  inputs.g = 9.807 / RStdata::Pa_per_ATM;
  inputs.saturated = false;
  inputs.inflow_velocity = -1.109096373e-10;
  inputs.Pwt = 1;

  pp.query("inflow_velocity",inputs.inflow_velocity);

  RegionManager rm(problo,probhi);

#if 0
  if (ParallelDescriptor::IOProcessor()) {
    std::cout << "The Regions: " << std::endl;
    const Array<const Region*> regions = rm.RegionPtrArray();
    for (int i=0; i<regions.size(); ++i) {
      std::cout << *(regions[i]) << std::endl;
    }
  }
#endif

  int nGrow = 0;
  RockManager rockManager(&rm);
  rockManager.FinalizeBuild(geom_array,refRatio_array,nGrow,false);


  RStdata rs_data(0,nLevs,layout,nlsc,inputs,&rockManager);
  rs_data.rel_perm_method = "upwind-darcy_velocity";
  rs_data.semi_analytic_J=true;

  pp.query("Niter",Niter);
  pp.query("dt",dt);
  pp.query("plt_name",plt_name);
  pp.query("verbose",verbose);
  pp.query("semi_analytic",rs_data.semi_analytic_J);
  bool abort_on_nl_fail = false; pp.query("abort_on_nl_fail",abort_on_nl_fail);

  pp.query("rel_perm_method",rs_data.rel_perm_method);
  pp.query("pressure_maxorder",rs_data.pressure_maxorder);

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

  RichardSolver* rs = new RichardSolver(rs_data,nlsc);

  if (pp.countval("record_file") ){
    std::string rf; pp.get("record_file",rf);
    rs->SetRecordFile(rf);
  }

  MFTower& Pold = *(rs_data.Pold);
  MFTower& Pnew = *(rs_data.Pnew);
  MFTower& RSold = *(rs_data.RhoSatOld);
  MFTower& RSnew = *(rs_data.RhoSatNew);

  Array<Real> grad(BL_SPACEDIM,0); grad[BL_SPACEDIM-1] = - inputs.rho[0]*inputs.g;

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
  int retCode;
  for (step++ ; step <= Niter; ++step) {
    bool cont = true;
    retCode = -1;
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

  delete rs;

  // Must clear prior to PetscFinalize to explicitly delete petsc data structures, or
  // the PetscFinalize will do it, and the Layout destructor will then generate a seg fault.
  layout.Clear();

#ifdef BL_USE_PETSC
  PetscFinalize();
#endif

  ParallelDescriptor::EndParallel();
  return (retCode < 0 ? 1 : 0);
}
