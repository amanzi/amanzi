#include <TensorDiffusion_PK.H>
#include <tTensorDiffusion_PK_F.H>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

#include <VisMF.H>
#include <Utility.H>

enum bc_t {Periodic = 0,
	   Dirichlet = LO_DIRICHLET, 
	   Neumann = LO_NEUMANN};

BCRec defaultBC()
{
  return BCRec(D_DECL(EXT_DIR,EXT_DIR,EXT_DIR),
	       D_DECL(EXT_DIR,EXT_DIR,EXT_DIR));
}

using Amanzi::AmanziTransport::loadBndryData;
using Amanzi::AmanziTransport::diffuse_tracer;
using Amanzi::AmanziTransport::ONEPASS;

int
main (int   argc,
      char* argv[])
{
  BoxLib::Initialize(argc,argv);

  std::cout << std::setprecision(10);

  if (argc < 2)
  {
    std::cerr << "usage:  " << argv[0] << " inputsfile [options]" << '\n';
    exit(-1);
  }

  ParmParse pp;
    
#if BL_SPACEDIM == 2
  Box domain(IntVect(0,0),IntVect(11,11));
  std::string boxfile("grids/gr2D") ;
#elif BL_SPACEDIM == 3
  Box domain(IntVect(0,0,0),IntVect(11,11,11));
  std::string boxfile("grids/gr.3_2x3x4") ;
#endif
  pp.query("boxes", boxfile);

  std::ifstream ifs(boxfile.c_str(), std::ios::in);

  if (!ifs)
  {
    std::string msg = "problem opening grids file: ";
    msg += boxfile.c_str();
    BoxLib::Abort(msg.c_str());
  }

  ifs >> domain;

  if (ParallelDescriptor::IOProcessor())
    std::cout << "domain: " << domain << std::endl;

  BoxArray bs;
  bs.readFrom(ifs);

  // allocate/init soln and rhs
  int Ncomp=1;
  int Nghost=1;
  MultiFab soln(bs, Ncomp, Nghost, Fab_allocate); soln.setVal(0.0);
  MultiFab out(bs, Ncomp, Nghost, Fab_allocate); 

  Nghost = 0;

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
    
  std::string bc_type_s;
  pp.get("bc_type",bc_type_s);
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
  Geometry geom(domain,&real_box,coord,is_per);
  Real dx[BL_SPACEDIM];
  for ( int n=0; n<BL_SPACEDIM; n++ ) {
    dx[n] = ( geom.ProbHi(n) - geom.ProbLo(n) )/domain.length(n);
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
  int ratio=2; pp.query("ratio", ratio);

  MultiFab rhs(bs, Ncomp, 0, Fab_allocate);
  for(MFIter rhsmfi(rhs); rhsmfi.isValid(); ++rhsmfi)
  {
    FORT_FILLRHS(rhs[rhsmfi].dataPtr(),
                 ARLIM(rhs[rhsmfi].loVect()),ARLIM(rhs[rhsmfi].hiVect()),
                 dx,&Ncomp,domain.loVect(),domain.hiVect(),geom.ProbLo());
  }
    
  Nghost = 1; // need space for bc info
  MultiFab fine(bs,Ncomp,Nghost,Fab_allocate);
  for(MFIter finemfi(fine); finemfi.isValid(); ++finemfi)
  {
    FORT_FILLFINE(fine[finemfi].dataPtr(),
                  ARLIM(fine[finemfi].loVect()),ARLIM(fine[finemfi].hiVect()),
                  dx,&Ncomp,domain.loVect(),domain.hiVect(),geom.ProbLo());
  }

  // Create "background coarse data"
  //Box crse_bx = Box(domain).coarsen(ratio).grow(1);
  Box crse_bx = Box(domain).coarsen(ratio);
  BoxArray cba(crse_bx);
  cba.maxSize(32);
  Real h_crse[BL_SPACEDIM];
  for (int d=0; d<BL_SPACEDIM; d++) h_crse[d] = dx[d]*ratio;

  MultiFab crse(cba, Ncomp, 0);
  for (MFIter mfi(crse); mfi.isValid(); ++mfi)
  {
    FORT_FILLCRSE(crse[mfi].dataPtr(),
                  ARLIM(crse[mfi].loVect()),ARLIM(crse[mfi].hiVect()),
                  h_crse,&Ncomp,crse_bx.loVect(),crse_bx.hiVect(),geom.ProbLo());
  }

  TensorDiffusionBndry bd(bs,Ncomp,geom);
  bool single_level=true; pp.query("single_level",single_level);
  MultiFab* c_ptr = single_level ? 0 : &crse;
  loadBndryData(bd,fine,0,c_ptr,0,pbcarray,geom,ratio,Ncomp);

  Real t_start = 0; pp.query("t_start",t_start);
  Real delta_t = 1; pp.query("delta_t",delta_t);  
  Real t_end = 10*delta_t; pp.query("t_end",t_end);
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
  MultiFab::Multiply(soln,W_old,0,0,Ncomp,0);

  soln.setVal(0,0,Ncomp,0);

  for (int d=0; d<BL_SPACEDIM; ++d)
  {
    for(MFIter bmfi(*beta[d]); bmfi.isValid(); ++bmfi) {
      FORT_MAKEMU((*beta[d])[bmfi].dataPtr(),
                  ARLIM((*beta[d])[bmfi].loVect()),ARLIM((*beta[d])[bmfi].hiVect()),dx,d);
    }

    for(MFIter b1mfi(*beta1[d]); b1mfi.isValid(); ++b1mfi) {
      FORT_MAKEMU1((*beta1[d])[b1mfi].dataPtr(),
                   ARLIM((*beta1[d])[b1mfi].loVect()),ARLIM((*beta1[d])[b1mfi].hiVect()),dx,d);
    }
  }
  
  MultiFab::Copy(soln_new,soln,0,0,Ncomp,soln.nGrow());

  Real old_time = t_start;
  Real new_time = std::min(t_end,t_start+delta_t);
  int step = 0;
  int max_step = 10; pp.query("max_step",max_step);

  int dump_soln_new=0; pp.query("dump_soln_new",dump_soln_new);

  if (dump_soln_new) {
    VisMF::Write(soln_new,BoxLib::Concatenate("soln_new_",step,3));
  }
  step++;
  while (old_time < t_end && step <= max_step) {

    std::cout << "Taking step " << step << " from " << old_time << " to " << new_time << std::endl;

    bool add_old_time_divFlux = true;
    diffuse_tracer(old_time, new_time, theta, soln, 0, soln_new, 0, &W_old, 0, &W_new, 0,
                   W_half, 0, W_flag, bd, 0, bd, 0, fluxn, fluxnp1, 0, &dRhs, 0, 
                   alpha, 0, beta, 0, beta, 0, beta1, 0, beta1, 0, Ncomp, 
                   ONEPASS, max_order, add_old_time_divFlux);
    
    if (dump_soln_new) {
      VisMF::Write(soln_new,BoxLib::Concatenate("soln_new_",step,3));
    }

    old_time = new_time;
    delta_t *= 1.5;
    new_time = std::min(t_end,old_time+delta_t);
    step++;
    MultiFab::Copy(soln,soln_new,0,0,Ncomp,0);
  }

  for (int d=0; d<BL_SPACEDIM; d++) {
    delete fluxn[d];
    delete fluxnp1[d];
    delete beta[d];
    delete beta1[d];
  }
  ParallelDescriptor::EndParallel();
}
