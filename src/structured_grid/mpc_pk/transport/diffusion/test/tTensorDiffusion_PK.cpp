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

const Real visc_tol = 1.e-12;
bool use_cg_solve = false;
bool use_mg_precond_flag = false;

namespace Amanzi {
  namespace AmanziTransport {

    enum SolveMode {PREDICTOR, CORRECTOR, ONEPASS};


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
                   bool                   add_old_time_divFlux)
    {
#if 0
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
      MultiFab Rhs(grids,nComp,0), Soln(grids,nComp,MCLinOp_grow);
      MultiFab volume; geom.GetVolume(volume,grids,Geom_Grow);
      PArray<MultiFab> area(BL_SPACEDIM);
      for (int d = 0; d < BL_SPACEDIM; d++) {
        area.set(d, new MultiFab());
        geom.GetFaceArea(area[d],grids,d,0);
      }

      if (add_old_time_divFlux)
      {
        Real a = 0.0;
        Real b = -(1.0-be_cn_theta)*dt;
        Real scale_old = 0;
        TensorOp* op_old = getOp(a,b,bd_old,sComp_bd_old,W_old,sComp_W_old,W_half,sComp_W_half,W_flag,scale_old,
                                 betan,sComp_betan,beta1n,sComp_beta1n,volume,area,alpha_in,sComp_alpha_in,nComp);
        op_old->maxOrder(max_order);

        MultiFab::Copy(Soln,S_old,sComp_S_old,0,nComp,0);
        if (W_flag == 2) {
          for (MFIter Smfi(Soln); Smfi.isValid(); ++Smfi) {
            for (int n=0; n<nComp; ++n) {
              Soln[Smfi].divide((*W_old)[Smfi],Smfi.validbox(),sComp_W_old,n,1);
            }
          }
        }

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
      // If this is a predictor step, put "explicit" updates passed via S_new
      // into delta_rhs after scaling by rho_half if reqd, so they dont get lost,
      // pull it off S_new to avoid double counting
      //   (for rho_flag == 1:
      //       S_new = S_old - dt.(U.Grad(phi)); want Rhs -= rho_half.(U.Grad(phi)),
      //    else
      //       S_new = S_old - dt.Div(U.Phi),   want Rhs -= Div(U.Phi) )
      //
      if (solve_mode == PREDICTOR) {
        FArrayBox tmpfab;
        for (MFIter Smfi(S_new); Smfi.isValid(); ++Smfi) {
          const Box& box = Smfi.validbox();
          tmpfab.resize(box,nComp);
          tmpfab.copy(S_new[Smfi],box,sComp_S_new,box,0,nComp);
          tmpfab.minus(S_old[Smfi],box,sComp_S_old,0,nComp);
          S_new[Smfi].minus(tmpfab,box,0,sComp_S_new,nComp); // Remove this term from S_new
          tmpfab.mult(1.0/dt,box,0,nComp);
          if (W_flag == 1) {
            for (int n=0; n<nComp; ++n) {
              tmpfab.mult((*W_half)[Smfi],box,sComp_W_half,n,1);
            }
          }
          if (alpha_in!=0) {
            tmpfab.mult((*alpha_in)[Smfi],box,sComp_alpha_in,0,nComp);
          }
          (*delta_rhs)[Smfi].plus(tmpfab,box,0,sComp_rhs,nComp);
        }
      }

      //
      // Add body sources
      //
      if (delta_rhs != 0) {
        FArrayBox tmpfab;
        for (MFIter mfi(*delta_rhs); mfi.isValid(); ++mfi)
        {
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
      //  (Note: here S_new holds S_old, but also maybe an explicit increment
      //         from advection if solve_mode != PREDICTOR)
      //
      MultiFab::Copy(Soln,S_new,sComp_S_new,0,nComp,0);

      for (MFIter mfi(Soln); mfi.isValid(); ++mfi) {
        const Box& box = mfi.validbox();
        for (int n=0; n<nComp; ++n) {
          Soln[mfi].mult(volume[mfi],box,0,n,1);
        }
        if (W_flag == 1) {
          for (int n=0; n<nComp; ++n) {
            Soln[mfi].mult((*W_half)[mfi],box,sComp_W_half,n,1);
          }
        } else  if (W_flag == 3) {
          for (int n=0; n<nComp; ++n) {
            Soln[mfi].mult((*W_old)[mfi],box,sComp_W_old,n,1);
          }
        }
        if (alpha_in!=0) {
          for (int n=0; n<nComp; ++n) {
            Soln[mfi].mult((*alpha_in)[mfi],box,sComp_alpha_in+n,n,1);
          }
          Rhs[mfi].plus(Soln[mfi],box,0,0,nComp);
        }
      }

      MultiFab::Copy(Soln,S_new,sComp_S_new,0,nComp,0);
      if (W_flag == 2) {
        for (MFIter Smfi(Soln); Smfi.isValid(); ++Smfi) {
          for (int n=0; n<nComp; ++n) {
            Soln[Smfi].divide((*W_new)[Smfi],Smfi.validbox(),sComp_W_new,n,1);
          }
        }
      }

      //
      // Construct viscous operator with bndry data at time N+1.
      //
      Real a = 1.0;
      Real b = be_cn_theta*dt;
      Real scale_new = 1;
      TensorOp* op_new = getOp(a,b,bd_new,sComp_bd_new,W_new,sComp_W_new,W_half,sComp_W_half,W_flag,scale_new,
                               betanp1,sComp_betanp1,beta1np1,sComp_beta1np1,volume,area,alpha_in,sComp_alpha_in,nComp);
      op_new->maxOrder(max_order);
      Rhs.mult(scale_new,0,nComp);

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
    
      if (W_flag == 2) {
        for (MFIter Smfi(S_new); Smfi.isValid(); ++Smfi) {
          for (int n=0; n<nComp; ++n) {
            S_new[Smfi].mult((*W_new)[Smfi],Smfi.validbox(),sComp_W_new,sComp_S_new+n,1);
          }
        }
      }
#endif
    }

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
                   const InterpBndryData& bd_old,    int sComp_bd_old,
                   const InterpBndryData& bd_new,    int sComp_bd_new,
                   MultiFab* const*       fluxn,
                   MultiFab* const*       fluxnp1,   int dComp_flux,
                   MultiFab*              delta_rhs, int sComp_rhs,
                   const MultiFab*        alpha_in,  int sComp_alpha_in,
                   const MultiFab* const* betan,     int sComp_betan,
                   const MultiFab* const* betanp1,   int sComp_betanp1,
                   int                    nComp,
                   const SolveMode&       solve_mode,
                   int                    max_order,
                   bool                   add_old_time_divFlux)
    {
#if 0
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
      MultiFab Rhs(grids,nComp,0), Soln(grids,nComp,MCLinOp_grow);
      MultiFab volume; geom.GetVolume(volume,grids,Geom_Grow);
      PArray<MultiFab> area(BL_SPACEDIM,PArrayManage);
      for (int d = 0; d < BL_SPACEDIM; d++) {
        area.set(d, new MultiFab());
        geom.GetFaceArea(area[d],grids,d,0);
      }

      if (add_old_time_divFlux)
      {
        Real a = 0.0;
        Real b = -(1.0-be_cn_theta)*dt;
        Real scale_old = 0;
        ABecHelper* op_old = getOp(a,b,bd_old,sComp_bd_old,W_old,sComp_W_old,1,W_half,sComp_W_half,W_flag,
                                   betan,sComp_betan,nComp,volume,area,alpha_in,sComp_alpha_in,nComp);
        op_old->maxOrder(max_order);

        MultiFab::Copy(Soln,S_old,sComp_S_old,0,nComp,0);
        if (W_flag == 2) {
          for (MFIter Smfi(Soln); Smfi.isValid(); ++Smfi) {
            for (int n=0; n<nComp; ++n) {
              Soln[Smfi].divide((*W_old)[Smfi],Smfi.validbox(),sComp_W_old,n,1);
            }
          }
        }

        op_old->apply(Rhs,Soln);
        op_old->compFlux(D_DECL(*fluxn[0],*fluxn[1],*fluxn[2]),Soln,LinOp::Inhomogeneous_BC,
                         0,dComp_flux,nComp,sComp_bd_old);
        delete op_old;

        for (int i = 0; i < BL_SPACEDIM; ++i) {
          for (int n=0; n<nComp; ++n) {
            (*fluxn[i]).mult(-b/(dt*dx[i]),dComp_flux+n,1,0);
          }
        }
      } else {
        Rhs.setVal(0);
      }

      //
      // If this is a predictor step, put "explicit" updates passed via S_new
      // into delta_rhs after scaling by rho_half if reqd, so they dont get lost,
      // pull it off S_new to avoid double counting
      //   (for rho_flag == 1:
      //       S_new = S_old - dt.(U.Grad(phi)); want Rhs -= rho_half.(U.Grad(phi)),
      //    else
      //       S_new = S_old - dt.Div(U.Phi),   want Rhs -= Div(U.Phi) )
      //
      if (solve_mode == PREDICTOR) {
        FArrayBox tmpfab;
        for (MFIter Smfi(S_new); Smfi.isValid(); ++Smfi) {
          const Box& box = Smfi.validbox();
          tmpfab.resize(box,nComp);
          tmpfab.copy(S_new[Smfi],box,sComp_S_new,box,0,nComp);
          tmpfab.minus(S_old[Smfi],box,sComp_S_old,0,nComp);
          S_new[Smfi].minus(tmpfab,box,0,sComp_S_new,nComp); // Remove this term from S_new
          tmpfab.mult(1.0/dt,box,0,nComp);
          if (W_flag == 1) {
            for (int n=0; n<nComp; ++n) {
              tmpfab.mult((*W_half)[Smfi],box,sComp_W_half,n,1);
            }
          }
          if (alpha_in!=0) {
            tmpfab.mult((*alpha_in)[Smfi],box,sComp_alpha_in,0,nComp);
          }
          (*delta_rhs)[Smfi].plus(tmpfab,box,0,sComp_rhs,nComp);
        }
      }

      //
      // Add body sources
      //
      if (delta_rhs != 0) {
        FArrayBox tmpfab;
        for (MFIter mfi(*delta_rhs); mfi.isValid(); ++mfi)
        {
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
      //  (Note: here S_new holds S_old, but also maybe an explicit increment
      //         from advection if solve_mode != PREDICTOR)
      //
      MultiFab::Copy(Soln,S_new,sComp_S_new,0,nComp,0);

      for (MFIter mfi(Soln); mfi.isValid(); ++mfi) {
        const Box& box = mfi.validbox();
        for (int n=0; n<nComp; ++n) {
          Soln[mfi].mult(volume[mfi],box,0,n,1);
        }
        if (W_flag == 1) {
          for (int n=0; n<nComp; ++n) {
            Soln[mfi].mult((*W_half)[mfi],box,sComp_W_half,n,1);
          }
        } else  if (W_flag == 3) {
          for (int n=0; n<nComp; ++n) {
            Soln[mfi].mult((*W_old)[mfi],box,sComp_W_old,n,1);
          }
        }
        if (alpha_in!=0) {
          for (int n=0; n<nComp; ++n) {
            Soln[mfi].mult((*alpha_in)[mfi],box,sComp_alpha_in+n,n,1);
          }
          Rhs[mfi].plus(Soln[mfi],box,0,0,nComp);
        }
      }

      MultiFab::Copy(Soln,S_new,sComp_S_new,0,nComp,0);
      if (W_flag == 2) {
        for (MFIter Smfi(Soln); Smfi.isValid(); ++Smfi) {
          for (int n=0; n<nComp; ++n) {
            Soln[Smfi].divide((*W_new)[Smfi],Smfi.validbox(),sComp_W_new,n,1);
          }
        }
      }

      //
      // Construct viscous operator with bndry data at time N+1.
      //
      Real a = 1.0;
      Real b = be_cn_theta*dt;
      Real scale_new = 1;
      ABecLaplacian* op_new = getOp(a,b,bd_new,sComp_bd_new,W_new,sComp_W_new,W_half,sComp_W_half,W_flag,scale_new,
                                    betanp1,sComp_betanp1,volume,area,alpha_in,sComp_alpha_in);
      op_new->maxOrder(max_order);
      Rhs.mult(scale_new,0,nComp);

      //
      // Construct solver and call it.
      //
      const Real S_tol     = visc_tol;
      const Real S_tol_abs = get_scaled_abs_tol(Rhs, visc_tol);

      if (use_cg_solve) {
        CGSolver cg(*op_new,use_mg_precond_flag);
        cg.solve(Soln,Rhs,S_tol,S_tol_abs);
      }
      else
      {
        MultiGrid mg(*op_new);
        mg.solve(Soln,Rhs,S_tol,S_tol_abs);
      }

      Rhs.clear();
      //
      // Get extensivefluxes from new-time op
      //
      op_new->compFlux(D_DECL(*fluxnp1[0],*fluxnp1[1],*fluxnp1[2]),Soln,LinOp::Inhomogeneous_BC,
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
    
      if (W_flag == 2) {
        for (MFIter Smfi(S_new); Smfi.isValid(); ++Smfi) {
          for (int n=0; n<nComp; ++n) {
            S_new[Smfi].mult((*W_new)[Smfi],Smfi.validbox(),sComp_W_new,sComp_S_new+n,1);
          }
        }
      }
#endif
    }
  } /* AmanziTransport */
} /* Amanzi */

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
