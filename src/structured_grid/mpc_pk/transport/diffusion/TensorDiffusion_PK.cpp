
#include <TensorDiffusion_PK.H>

#include <MCMultiGrid.H>
#include <MultiGrid.H>
#include <MCCGSolver.H>
#include <ParallelDescriptor.H>


namespace Amanzi {
  namespace AmanziTransport {

    const Real visc_tol = 1.e-12;
    bool use_cg_solve = false;
    bool use_mg_precond_flag = false;


    Real get_scaled_abs_tol (const MultiFab& rhs,
                             Real            reduction)
    {
      Real norm_est = 0;
      for (MFIter Rhsmfi(rhs); Rhsmfi.isValid(); ++Rhsmfi)
        norm_est = std::max(norm_est, rhs[Rhsmfi].norm(0));
      ParallelDescriptor::ReduceRealMax(norm_est);
      return norm_est * reduction;
    }

    TensorOp*
    getOp(Real                     a,
          Real                     b,
          const MCInterpBndryData& bd,       int sComp_bd,
          MultiFab*                W,        int sComp_W,
          MultiFab*                W_half,   int sComp_W_half,
          int                      W_flag,
          Real                     scale,
          const MultiFab* const*   beta,     int sComp_beta,
          const MultiFab* const*   beta1,    int sComp_beta1,
          const MultiFab&          volume,
          const PArray<MultiFab>&  area,
          const MultiFab*          alpha_in, int sComp_alpha_in,
          int                      nComp);

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
        crse_br.copyFrom(S_fine,S_fine.nGrow(),sComp_S_fine,0,nComp);
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
    }

    TensorOp*
    getOp(Real                     a,
          Real                     b,
          const MCInterpBndryData& bd,       int sComp_bd,
          MultiFab*                W,        int sComp_W,
          MultiFab*                W_half,   int sComp_W_half,
          int                      W_flag,
          Real                     scale,
          const MultiFab* const*   beta,     int sComp_beta,
          const MultiFab* const*   beta1,    int sComp_beta1,
          const MultiFab&          volume,
          const PArray<MultiFab>&  area,
          const MultiFab*          alpha_in, int sComp_alpha_in,
          int                      nComp)
    {
      const Real* dx = bd.getGeom().CellSize();
      TensorOp* op = new TensorOp(bd,dx,nComp);
      int nComp_alpha = (W_flag > 1 || alpha_in!=0 ?  nComp  :  1); 
      const BoxArray& grids = volume.boxArray();
      MultiFab alpha(grids,nComp_alpha,Geom_Grow);

      for (int n=0; n<nComp_alpha; ++n) {
        MultiFab::Copy(alpha,volume,0,n,1,Geom_Grow);
      }
      if (W_flag == 1) {
        for (int n=0; n<nComp_alpha; ++n) {
          MultiFab::Multiply(alpha,*W_half,sComp_W_half,n,1,alpha.nGrow());
        }
      }
      else if (W_flag == 2 || W_flag == 3) {
        for (int n=0; n<nComp; ++n) {
          MultiFab::Multiply(alpha,*W,sComp_W,n,1,alpha.nGrow());
        }
      }
      if (alpha_in) {
        MultiFab::Multiply(alpha,*alpha_in,sComp_alpha_in,0,nComp,alpha.nGrow());
      }
  
      scale = (scale != 0 ? 1.0/alpha.max(0) : 1.0);
      op->setScalars(a*scale, b*scale);
      op->aCoefficients(alpha);
      alpha.clear();

      for (int d = 0; d < BL_SPACEDIM; d++) {
        MultiFab bcoeffs(area[d].boxArray(),nComp,0);
        for (int n=0; n<nComp; ++n) {
          MultiFab::Copy(bcoeffs,area[d],0,n,1,0);
        }
        for (MFIter mfi(bcoeffs); mfi.isValid(); ++mfi) {
          bcoeffs[mfi].mult((*beta[d])[mfi],sComp_beta,0,nComp);
          bcoeffs[mfi].mult(dx[d]);
        }
        op->bCoefficients(bcoeffs,d);
    
        for (int n=0; n<nComp; ++n) {
          MultiFab::Copy(bcoeffs,area[d],0,n,1,0);
        }
        for (MFIter mfi(bcoeffs); mfi.isValid(); ++mfi) {
          bcoeffs[mfi].mult((*beta1[d])[mfi],sComp_beta1,0,nComp);
          bcoeffs[mfi].mult(dx[d]);
        }
        op->b1Coefficients(bcoeffs,d);
      }
      return op;
    }

  } /* AmanziTransport */
} /* Amanzi */
