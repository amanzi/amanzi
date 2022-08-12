#include <TensorDiffusion_PK.H>
#include <tTensorDiffusion_PK_F.H>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

// closing DSO objects
#include "VerboseObject_objs.hh"

#include <VisMF.H>
#include <Utility.H>
#include <MCCGSolver.H>
#include <MCMultiGrid.H>

enum bc_t {Periodic = 0,
	   Dirichlet = LO_DIRICHLET, 
	   Neumann = LO_NEUMANN};

BCRec defaultBC()
{
  return BCRec(D_DECL(EXT_DIR,EXT_DIR,EXT_DIR),
	       D_DECL(EXT_DIR,EXT_DIR,EXT_DIR));
}

const Real visc_tol = 1.e-8;
bool use_cg_solve = false;
bool use_mg_precond_flag = false;

enum SolveMode {PREDICTOR, CORRECTOR, ONEPASS};

Real get_scaled_abs_tol (const MultiFab& rhs,
			 Real            reduction)
{
  Real norm_est = 0;
  for (MFIter Rhsmfi(rhs); Rhsmfi.isValid(); ++Rhsmfi)
    norm_est = std::max(norm_est, rhs[Rhsmfi].norm(0));
  ParallelDescriptor::ReduceRealMax(norm_est);
  return norm_est * reduction;
}

void
loadBndryData (MCInterpBndryData&  bd,
	       MultiFab&           S_fine, int sComp_S_fine,
	       MultiFab*           S_crse, int sComp_S_crse,
	       const Array<BCRec>& bc,
	       const Geometry&     geom,
	       int                 ratio,
	       int                 nComp)
{
  if (S_crse == 0) {
    bd.setBndryValues(S_fine,sComp_S_fine,0,nComp,bc);
  } else {
    BoxArray cgrids = BoxArray(S_fine.boxArray()).coarsen(ratio);
    BndryRegister crse_br(cgrids,0,1,2,nComp);
    crse_br.copyFrom(*S_crse,S_crse->nGrow(),sComp_S_crse,0,nComp);
    bd.setBndryValues(crse_br,0,S_fine,sComp_S_fine,0,nComp,ratio,bc);
  }
}

using Amanzi::AmanziTransport::getOp;

void
diffuse_tracer(Real                   t_old,
	       Real                   t_new,
	       Real                   be_cn_theta,
	       const MultiFab&        S_old,     int sComp_S_old,
	       MultiFab&              S_new,     int sComp_S_new,
	       MultiFab*              W_old,     int sComp_W_old,
	       MultiFab*              W_new,     int sComp_W_new,
	       MultiFab*              W_half,    int sComp_W_half,
	       int                    W_flag,
	       const MCInterpBndryData& bd_old,  int sComp_bd_old,
	       const MCInterpBndryData& bd_new,  int sComp_bd_new,
	       MultiFab* const*       fluxn,
	       MultiFab* const*       fluxnp1,   int dComp_flux,
	       MultiFab*              delta_rhs, int sComp_rhs,
	       const MultiFab*        alpha_in,  int sComp_alpha_in,
	       const MultiFab* const* betan,     int sComp_betan,
	       const MultiFab* const* betanp1,   int sComp_betanp1,
	       const MultiFab* const* beta1n,    int sComp_beta1n,
	       const MultiFab* const* beta1np1,  int sComp_beta1np1,
	       int                    nComp,
	       const SolveMode&       solve_mode,
	       int                    max_order,
	       bool                   add_old_time_divFlux,
	       int                    mg_verbose,
	       int                    mg_usecg)
{
  const BoxArray& grids = S_old.boxArray();
  BL_ASSERT(solve_mode==ONEPASS || (delta_rhs && delta_rhs->boxArray()==grids));
  if (alpha_in) {
    BL_ASSERT(alpha_in->nComp() >= nComp + sComp_alpha_in);
  }

  Real dt = t_new - t_old;
  BL_ASSERT(dt > 0);
  const Geometry& geom = bd_old.getGeom();
  const Real* dx = geom.CellSize();
  //
  // Set up Rhs.
  //
  MultiFab Rhs(grids,nComp,0), Soln(grids,nComp,1);
  MultiFab volume; geom.GetVolume(volume,grids,1);
  MultiFab area[BL_SPACEDIM];
  for (int d = 0; d < BL_SPACEDIM; d++) {
    geom.GetFaceArea(area[d],grids,d,0);
  }

  if (add_old_time_divFlux) {
    Real a = 0.0;
    Real b = -(1.0-be_cn_theta)*dt;
    TensorOp* op_old = getOp(a,b,bd_old,sComp_bd_old,1,W_old,sComp_W_old,1,W_half,sComp_W_half,W_flag,
			     betan,sComp_betan,1,beta1n,sComp_beta1n,1,volume,area,alpha_in,sComp_alpha_in);
    op_old->maxOrder(max_order);
    
    MultiFab::Copy(Soln,S_old,sComp_S_old,0,nComp,0);

    op_old->apply(Rhs,Soln);
    op_old->compFlux(D_DECL(*fluxn[0],*fluxn[1],*fluxn[2]),Soln,MCInhomogeneous_BC,
		     0,dComp_flux,nComp,sComp_bd_old);
    delete op_old;
    
    for (int d = 0; d < BL_SPACEDIM; ++d) {
      for (int n=0; n<nComp; ++n) {
	(*fluxn[d]).mult(-b/(dt*dx[d]),dComp_flux+n,1,0);
      }
    }
  } else {
    Rhs.setVal(0);
  }

  //
  // Add body sources
  //
  if (delta_rhs != 0) {
    FArrayBox tmpfab;
    for (MFIter mfi(*delta_rhs); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      tmpfab.resize(box,nComp);
      tmpfab.copy((*delta_rhs)[mfi],box,sComp_rhs,box,0,nComp);
      tmpfab.mult(dt,box,0,nComp);
      for (int n=0; n<nComp; ++n) {
	tmpfab.mult(volume[mfi],box,0,n,1);
      }
      Rhs[mfi].plus(tmpfab,box,0,0,nComp);
    }
  }

  //
  // Increment Rhs with S_old*V (or S_old*V*rho_half if rho_flag==1
  //                             or S_old*V*rho_old  if rho_flag==3)
  //
  for (MFIter mfi(Soln); mfi.isValid(); ++mfi) {
    const Box& box = mfi.validbox();

    for (int n=0; n<nComp; ++n) {

      Soln[mfi].copy(S_old[mfi],sComp_S_new+n,n,1);

      if (W_flag == 2) {
	Soln[mfi].mult((*W_old)[mfi],box,sComp_W_old,n,1);
      }
      else if (W_flag == 1) {
	Soln[mfi].mult((*W_half)[mfi],box,sComp_W_half,n,1);
      }

      Soln[mfi].mult(volume[mfi],box,0,n,1);

      if (alpha_in!=0) {
	Soln[mfi].mult((*alpha_in)[mfi],box,sComp_alpha_in+n,n,1);
      }

      Rhs[mfi].plus(Soln[mfi],box,n,n,1);
    }
  }

  MultiFab::Copy(Soln,S_new,sComp_S_new,0,nComp,0);

  //
  // Construct viscous operator with bndry data at time N+1.
  //
  Real a = 1.0;
  Real b = be_cn_theta*dt;
  TensorOp* op_new = getOp(a,b,bd_new,sComp_bd_new,1,W_new,sComp_W_new,1,W_half,sComp_W_half,W_flag,
			   betanp1,sComp_betanp1,1,beta1np1,sComp_beta1np1,1,volume,area,alpha_in,sComp_alpha_in);
  op_new->maxOrder(max_order);

  Real new_scale = 1.0/op_new->aCoefficients().max(0);
  op_new->setScalars(a*new_scale, b*new_scale);
  Rhs.mult(new_scale,0,nComp);

#if 0
  // test: build relaxation coeff on all cells
  MultiFab Ax(Soln.boxArray(),1,0);
  MultiFab Axp(Soln.boxArray(),1,0);
  MultiFab coef(Soln.boxArray(),1,0);
  op_new->apply(Ax,Soln);
  for (MFIter mfi(Soln); mfi.isValid(); ++mfi) {
    const Box& box = mfi.validbox();
    FArrayBox& s=Soln[mfi];
    for (IntVect iv=box.smallEnd(); iv<=box.bigEnd(); box.next(iv)) {
      Real saveVal = s(iv,0);
      Real old = s(iv,0);
      Real eps = .0002;
      s(iv,0) += eps;
      op_new->apply(Axp,Soln);
      s(iv,0) = old;
      coef[mfi](iv,0) = (Axp[mfi](iv,0) - Ax[mfi](iv,0))/eps;
    }
  }
  VisMF::Write(coef,"COEF");
  BoxLib::Abort();
#endif

  //
  // Construct solver and call it.
  //
  const Real S_tol     = visc_tol;
  const Real S_tol_abs = get_scaled_abs_tol(Rhs, visc_tol);

  if (use_cg_solve) {
    MCCGSolver cg(*op_new,use_mg_precond_flag);
    cg.solve(Soln,Rhs,S_tol,S_tol_abs);
  }
  else
    {
      MCMultiGrid mg(*op_new);
      mg.setUseCG(mg_usecg);
      mg.setVerbose(mg_verbose);
      mg.solve(Soln,Rhs,S_tol,S_tol_abs);
    }

  Rhs.clear();
  //
  // Get extensivefluxes from new-time op
  //
  op_new->compFlux(D_DECL(*fluxnp1[0],*fluxnp1[1],*fluxnp1[2]),Soln,MCInhomogeneous_BC,
		   0,dComp_flux,nComp,sComp_bd_new);
  for (int i = 0; i < BL_SPACEDIM; ++i) {
    for (int n=0; n<nComp; ++n) {
      (*fluxnp1[i]).mult(b/(dt*dx[i]),dComp_flux+n,1,0);
    }
  }
  delete op_new;
  //
  // Copy into state variable at new time, without bc's
  //
  MultiFab::Copy(S_new,Soln,0,sComp_S_new,nComp,0);
}

void loadSoln(FArrayBox& fab);

void WritePlotfile(const std::string         &pfversion,
                   const PArray<MultiFab>    &data,
                   const Real                 time,
                   const Real*                probLo,
                   const Real*                probHi,
                   const Array<int>          &refRatio,
                   const Array<Box>          &probDomain,
                   const Array<Array<Real> > &dxLevel,
                   const int                  coordSys,
                   const std::string         &oFile,
                   const Array<std::string>  &names,
                   const bool                 verbose,
		   const bool                 isCartGrid,
		   const Real                *vfeps,
		   const int                 *levelSteps)
{
    if(ParallelDescriptor::IOProcessor()) {
      if( ! BoxLib::UtilCreateDirectory(oFile,0755)) {
         BoxLib::CreateDirectoryFailed(oFile);
      }
    }
    //
    // Force other processors to wait till directory is built.
    //
    ParallelDescriptor::Barrier();
    
    std::string oFileHeader(oFile);
    oFileHeader += "/Header";
    
    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
    
    std::ofstream os;
    
    //os.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
    
    if(verbose && ParallelDescriptor::IOProcessor()) {
      std::cout << "Opening file = " << oFileHeader << '\n';
    }
    
    os.open(oFileHeader.c_str(), std::ios::out|std::ios::binary);
    
    if(os.fail()) {
      BoxLib::FileOpenFailed(oFileHeader);
    }
    //
    // Start writing plotfile.
    //
    os << pfversion << '\n';
    int n_var = data[0].nComp();
    os << n_var << '\n';
    for (int n = 0; n < n_var; n++) os << names[n] << '\n';
    os << BL_SPACEDIM << '\n';
    os << std::setprecision(30) << time << '\n';
    const int finestLevel = data.size() - 1;
    os << finestLevel << '\n';
    for (int i = 0; i < BL_SPACEDIM; i++) os << probLo[i] << ' ';
    os << '\n';
    for (int i = 0; i < BL_SPACEDIM; i++) os << probHi[i] << ' ';
    os << '\n';
    for (int i = 0; i < finestLevel; i++) os << refRatio[i] << ' ';
    os << '\n';
    for (int i = 0; i <= finestLevel; i++) os << probDomain[i] << ' ';
    os << '\n';
    if(levelSteps != 0) {
      for (int i = 0; i <= finestLevel; i++) os << levelSteps[i] << ' ';
    } else {
      for (int i = 0; i <= finestLevel; i++) os << 0 << ' ';
    }
    os << '\n';
    for(int i = 0; i <= finestLevel; i++) {
      for(int k = 0; k < BL_SPACEDIM; k++) {
            os << dxLevel[i][k] << ' ';
      }
      os << '\n';
    }
    if(isCartGrid) {
      for(int i(0); i <= finestLevel; i++) {
        os << vfeps[i] << ' ';
      }
      os << '\n';
    }
    os << coordSys << '\n';
    os << 0 << '\n';                  // --------------- The bndry data width.
    //
    // Write out level by level.
    //
    for(int iLevel(0); iLevel <= finestLevel; ++iLevel) {
        //
        // Write state data.
        //
        const BoxArray &ba = data[iLevel].boxArray();
        int nGrids = ba.size();
        char buf[64];
        sprintf(buf, "Level_%d", iLevel);
        
        if(ParallelDescriptor::IOProcessor()) {
            os << iLevel << ' ' << nGrids << ' ' << time << '\n';
            if(levelSteps != 0) {
              os << levelSteps[iLevel] << '\n';
	    } else {
              os << 0 << '\n';
	    }
            
            for(int i(0); i < nGrids; ++i) {
              const Box &b = ba[i];
              for(int n(0); n < BL_SPACEDIM; ++n) {
                Real glo = b.smallEnd()[n] * dxLevel[iLevel][n];
                Real ghi = (b.bigEnd()[n]+1) * dxLevel[iLevel][n];
                os << glo << ' ' << ghi << '\n';
              }
            }
            //
            // Build the directory to hold the MultiFabs at this level.
            //
            std::string Level(oFile);
            Level += '/';
            Level += buf;
            
            if( ! BoxLib::UtilCreateDirectory(Level, 0755)) {
              BoxLib::CreateDirectoryFailed(Level);
	    }
        }
        //
        // Force other processors to wait till directory is built.
        //
        ParallelDescriptor::Barrier();
        //
        // Now build the full relative pathname of the MultiFab.
        //
        static const std::string MultiFabBaseName("MultiFab");
        
        std::string PathName(oFile);
        PathName += '/';
        PathName += buf;
        PathName += '/';
        PathName += MultiFabBaseName;
        
        if(ParallelDescriptor::IOProcessor()) {
            //
            // The full name relative to the Header file.
            //
            std::string RelativePathName(buf);
            RelativePathName += '/';
            RelativePathName += MultiFabBaseName;
            os << RelativePathName << '\n';
        }
        VisMF::Write(data[iLevel], PathName);
    }
    
    os.close();
}

int
main (int   argc,
      char* argv[])
{
  BoxLib::Initialize(argc,argv);

  std::cout << std::setprecision(10);

  // Set up a default set of boxes
  ParmParse pp;
  int ratio=2; pp.query("ratio", ratio);
  Array<Box> domains(2);
  domains[0] = Box(IntVect(D_DECL(0,0,0)),IntVect(D_DECL(31,31,31)));
  domains[1] = BoxLib::refine(domains[0],ratio);

  Box fine_box = BoxLib::refine(BoxLib::grow(domains[0],-6),ratio);
  BoxArray bs(fine_box);

  // allocate/init soln and rhs
  int Ncomp=1;
  int Nghost=1;
  MultiFab soln(bs, Ncomp, Nghost, Fab_allocate); soln.setVal(0.0);
  MultiFab out(bs, Ncomp, Nghost, Fab_allocate); 

  if (ParallelDescriptor::IOProcessor())
    std::cout << "grids:\n" << bs << std::endl;

  // This says we are using Cartesian coordinates
  int coord = 0;
    
  // This defines the physical size of the box.  Right now the box is [0,1] in each direction.
  RealBox real_box;
  for (int n = 0; n < BL_SPACEDIM; n++) {
    real_box.setLo(n, 0.0);
    real_box.setHi(n, 1.0);
  }

  // This sets the boundary conditions to be periodic or not
  int* is_per = new int[BL_SPACEDIM];
  bc_t bc_type = Dirichlet;
    
  std::string bc_type_s = "Dirichlet";
  pp.query("bc_type",bc_type_s);
  if (bc_type_s == "Dirichlet") {
    bc_type = Dirichlet;
  }
  else if (bc_type_s == "Neumann") {
    bc_type = Neumann;
  }
  else if (bc_type_s == "Periodic") {
    bc_type = Periodic;
  }
  else {
    if (ParallelDescriptor::IOProcessor()) {
      std::cout << "Don't know this boundary type: " << bc_type << std::endl;
    }
    BoxLib::Error("");
  }

  if (bc_type == Dirichlet || bc_type == Neumann) {
    for (int n = 0; n < BL_SPACEDIM; n++) is_per[n] = 0;
  } 
  else {
    for (int n = 0; n < BL_SPACEDIM; n++) is_per[n] = 1;
  }
    
  // This defines a Geometry object which is useful for writing the plotfiles
  Geometry geom(domains[0],&real_box,coord,is_per);
  Array<Array<Real> > dx(2,Array<Real>(BL_SPACEDIM));
  for ( int n=0; n<BL_SPACEDIM; n++ ) {
    dx[1][n] = ( geom.ProbHi(n) - geom.ProbLo(n) )/domains[1].length(n);
    dx[0][n] = ratio * dx[1][n];
  }

  // Create the BCRec's interpreted by ViscBndry objects
  // Create the boundary object
  int nBndComp = MCLinOp::bcComponentsNeeded(Ncomp);
  Array<BCRec> pbcarray(nBndComp,defaultBC());
  if (bc_type == Neumann) {
    for (int d=0; d<BL_SPACEDIM; d++ ) {
      for (int n=0; n<Ncomp; ++n) {
        pbcarray[n].setLo(d,REFLECT_EVEN);
        pbcarray[n].setHi(d,REFLECT_EVEN);
      }
    }
  }

  // Create "background coarse data"
  int nGrowC = 3; pp.query("nGrowC",nGrowC);
  BoxArray cba(domains[0]);
  cba.maxSize(32);

  int crse_flag = 0; pp.query("crse_flag",crse_flag);
  MultiFab crse(cba, Ncomp, 0);
  for (MFIter mfi(crse); mfi.isValid(); ++mfi)
  {
    FORT_FILLCRSE(crse[mfi].dataPtr(),
                  ARLIM(crse[mfi].loVect()),ARLIM(crse[mfi].hiVect()),
                  dx[0].dataPtr(),&Ncomp,domains[0].loVect(),domains[0].hiVect(),geom.ProbLo(),&crse_flag);
  }

  MultiFab fine(bs, Ncomp, 2);
  for (MFIter mfi(fine); mfi.isValid(); ++mfi)
  {
    const Box& box = fine[mfi].box();
    FORT_FILLCRSE(fine[mfi].dataPtr(),
                  ARLIM(fine[mfi].loVect()),ARLIM(fine[mfi].hiVect()),
                  dx[1].dataPtr(),&Ncomp,box.loVect(),box.hiVect(),geom.ProbLo(),&crse_flag);
  }

  MultiFab rhs(bs, Ncomp, 0, Fab_allocate);
  for(MFIter rhsmfi(rhs); rhsmfi.isValid(); ++rhsmfi)
  {
    FORT_FILLRHS(rhs[rhsmfi].dataPtr(),
                 ARLIM(rhs[rhsmfi].loVect()),ARLIM(rhs[rhsmfi].hiVect()),
                 dx[1].dataPtr(),&Ncomp,domains[1].loVect(),domains[1].hiVect(),geom.ProbLo());
  }

  int fine_flag = 0; pp.query("fine_flag",fine_flag);
  Nghost = 1; // need space for bc info
  MultiFab fineg(BoxArray(bs).grow(Nghost+2),Ncomp,0);

  fine.FillBoundary();
  TensorDiffusionBndry bd(bs,Ncomp,geom);
  loadBndryData(bd,fine,0,&crse,0,pbcarray,geom,ratio,Ncomp);

  Real t_start = 0; pp.query("t_start",t_start);
  Real delta_t_init = .01; pp.query("delta_t_init",delta_t_init);  
  Real t_end = 1.e15; pp.query("t_end",t_end);
  Real theta = 0.5; pp.query("theta",theta);
  int max_order = 3; pp.query("max_order",max_order);

  Nghost = 1;
  MultiFab soln_new(bs, Ncomp, Nghost);
  MultiFab W_new(bs, Ncomp, Nghost); W_new.setVal(2);
  MultiFab W_old(bs, Ncomp, Nghost); W_old.setVal(2);
  MultiFab dRhs(bs, Ncomp, 0);  dRhs.setVal(0);
  Real W_flag = 2;
  MultiFab* W_half = 0;
  MultiFab *fluxn[BL_SPACEDIM], *fluxnp1[BL_SPACEDIM];
  MultiFab *beta[BL_SPACEDIM], *beta1[BL_SPACEDIM];
  for (int d=0; d<BL_SPACEDIM; d++) {
    BoxArray fbox(bs); fbox.surroundingNodes(d);
    fluxn[d] = new MultiFab(fbox,Ncomp,Nghost);
    fluxnp1[d] = new MultiFab(fbox,Ncomp,Nghost);
    beta[d] = new MultiFab(fbox,Ncomp,Nghost);
    beta1[d] = new MultiFab(fbox,Ncomp,Nghost);
  }
  MultiFab* alpha = 0;
  MultiFab::Copy(soln,fine,0,0,Ncomp,1);

  Real aT = 0; pp.query("aT",aT);
  Real aL = 1; pp.query("aL",aL);
  Real u = -1; pp.query("u",u);
  Real v = 1; pp.query("v",v);

  for (int d=0; d<BL_SPACEDIM; ++d)
  {
    for(MFIter bmfi(*beta[d]); bmfi.isValid(); ++bmfi) {
      FORT_MAKEMU((*beta[d])[bmfi].dataPtr(),
                  ARLIM((*beta[d])[bmfi].loVect()),ARLIM((*beta[d])[bmfi].hiVect()),dx[1].dataPtr(),d,
		  &aT, &aL, &u, &v);
    }

    for(MFIter b1mfi(*beta1[d]); b1mfi.isValid(); ++b1mfi) {
      FORT_MAKEMU1((*beta1[d])[b1mfi].dataPtr(),
                   ARLIM((*beta1[d])[b1mfi].loVect()),ARLIM((*beta1[d])[b1mfi].hiVect()),dx[1].dataPtr(),d,
		   &aT, &aL, &u, &v);
    }
  }
  
  MultiFab::Copy(soln_new,soln,0,0,Ncomp,soln.nGrow());

  Real old_time = t_start;
  Real delta_t = delta_t_init;
  Real new_time = std::min(t_end,t_start+delta_t);
  int step = 0;
  int max_step = 10; pp.query("max_step",max_step);

  bool do_output=false; pp.query("do_output",do_output);

  std::string pfversion = "PorousMedia-V1.1";
  PArray<MultiFab> data(2,PArrayNoManage);
  data.set(0, &crse);
  data.set(1, &soln_new);
  Array<std::string> varnames(1,"var");
  int plv=0;
  std::string pltname = BoxLib::Concatenate("plt",step,3);
  int isCartGrid = 0;
  Real vfeps=0;
  Array<int> levelSteps(2,0); levelSteps[1]=step;
  if (do_output) {
    WritePlotfile(pfversion,data,new_time,geom.ProbLo(),geom.ProbHi(),Array<int>(1,ratio),
		  domains,dx,coord,pltname,varnames,plv,isCartGrid,&vfeps,levelSteps.dataPtr());
  }
  int mg_verbose=0;pp.query("mg_verbose",mg_verbose);
  int mg_usecg=0;pp.query("mg_usecg",mg_usecg);

  step++;
  while (old_time < t_end && step <= max_step) {

    std::cout << "Taking step " << step << " from " << old_time << " to " << new_time << std::endl;

    bool add_old_time_divFlux = true;
    data[1].setVal(0);
    diffuse_tracer(old_time, new_time, theta, soln, 0, soln_new, 0, &W_old, 0, &W_new, 0,
                   W_half, 0, W_flag, bd, 0, bd, 0, fluxn, fluxnp1, 0, &dRhs, 0, 
                   alpha, 0, beta, 0, beta, 0, beta1, 0, beta1, 0, Ncomp, 
                   ONEPASS, max_order, add_old_time_divFlux,mg_verbose,mg_usecg);

    if (do_output) {
      levelSteps[1]=step;
      pltname = BoxLib::Concatenate("plt",step,3);
      WritePlotfile(pfversion,data,new_time,geom.ProbLo(),geom.ProbHi(),Array<int>(1,ratio),
		    domains,dx,coord,pltname,varnames,plv,isCartGrid,&vfeps,levelSteps.dataPtr());
    }
    old_time = new_time;
    delta_t *= 1.5;
    new_time = std::min(t_end,old_time+delta_t);
    step++;
    MultiFab::Copy(soln,soln_new,0,0,Ncomp,0);
  }

  bool dump_current_soln = false; pp.query("dump_current_soln",dump_current_soln);
  if (dump_current_soln) {
    VisMF::Write(soln,"SOLN");
  }
  bool doCheck=true; pp.query("doCheck",doCheck);
  if (doCheck>0) {

    // Try to make sure we are doing the correct problem!
    BL_ASSERT(pp.countval("boxes")==0);
    BL_ASSERT(bc_type_s == "Dirichlet");
    BL_ASSERT(ratio == 2);
    BL_ASSERT(t_start == 0);
    BL_ASSERT(delta_t_init == .01);
    BL_ASSERT(t_end > 0.1133300781);
    BL_ASSERT(theta == 0.5);
    BL_ASSERT(max_order == 3);
    BL_ASSERT(max_step == 10);
    BL_ASSERT(crse_flag == 0);
    BL_ASSERT(aT == 0);
    BL_ASSERT(aL == 1);
    BL_ASSERT(u == -1);
    BL_ASSERT(v == 1);
    
    MultiFab mf;
    VisMF::Read(mf,"tTensorDiffusion_TestRes");
    MultiFab::Subtract(mf,soln_new,0,0,1,0);
    Real normDiff = mf.max(0);
    std::cout << "Max of difference: " << normDiff << std::endl;
    
    if (normDiff > 1.e-2) {// How close is close?
      std::cerr << "TensorDiffusion hosed!" << std::endl;
      throw std::exception();
    }
  }
  
  for (int d=0; d<BL_SPACEDIM; d++) {
    delete fluxn[d];
    delete fluxnp1[d];
    delete beta[d];
    delete beta1[d];
  }
  ParallelDescriptor::EndParallel();
  }
