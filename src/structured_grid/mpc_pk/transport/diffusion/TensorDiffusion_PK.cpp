/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

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

    TensorOp*
    getOp(Real                     a,
          Real                     b,
          const MCInterpBndryData& bd,       int sComp_bd, int nComp_bd,
          MultiFab*                W,        int sComp_W,  int nComp_W,
          MultiFab*                W_half,   int sComp_W_half,
          int                      W_flag,
          const MultiFab* const    beta[BL_SPACEDIM],     int sComp_beta, int nComp_beta,
          const MultiFab* const    beta1[BL_SPACEDIM],    int sComp_beta1, int nComp_beta1,
          const MultiFab&          volume,
          MultiFab*   area,
          const MultiFab*          alpha_in, int sComp_alpha_in)
    {
      int nComp = 1;
      const Real* dx = bd.getGeom().CellSize();
      TensorOp* op = new TensorOp(bd,dx,nComp);
      int nComp_alpha = (W_flag > 1 || alpha_in!=0 ?  nComp_W  :  1);
      const BoxArray& grids = volume.boxArray();
      if (a!=0) {
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
          for (int n=0; n<nComp_alpha; ++n) {
            MultiFab::Multiply(alpha,*W,sComp_W,n,1,alpha.nGrow());
          }
        }
        if (alpha_in) {
          for (int n=0; n<nComp_alpha; ++n) {
            MultiFab::Multiply(alpha,*alpha_in,sComp_alpha_in,n,1,alpha.nGrow());
          }
        }
        op->aCoefficients(alpha);
        alpha.clear();
      }
      op->setScalars(a, b);
      for (int d = 0; d < BL_SPACEDIM; d++) {
        MultiFab bcoeffs(area[d].boxArray(),nComp_beta,0);
        for (int n=0; n<nComp_beta; ++n) {
          MultiFab::Copy(bcoeffs,area[d],0,n,1,0);
        }
        for (MFIter mfi(bcoeffs); mfi.isValid(); ++mfi) {
          bcoeffs[mfi].mult((*beta[d])[mfi],sComp_beta,0,nComp_beta);
          bcoeffs[mfi].mult(dx[d]);
        }
        op->bCoefficients(bcoeffs,d);

        for (int n=0; n<nComp_beta1; ++n) {
          MultiFab::Copy(bcoeffs,area[d],0,n,1,0);
        }
        for (MFIter mfi(bcoeffs); mfi.isValid(); ++mfi) {
          bcoeffs[mfi].mult((*beta1[d])[mfi],sComp_beta1,0,nComp_beta1);
          bcoeffs[mfi].mult(dx[d]);
        }
        op->b1Coefficients(bcoeffs,d);
      }
      return op;
    }

    ABecHelper*
    getOp(Real                     a,
          Real                     b,
          const ViscBndry&         bd,       int sComp_bd,
          MultiFab*                W,        int sComp_W, int nComp_W,
          MultiFab*                W_half,   int sComp_W_half,
          int                      W_flag,
          const MultiFab* const  beta[BL_SPACEDIM],     int sComp_beta, int nComp_beta,
          const MultiFab&          volume,
          MultiFab*  area,
          const MultiFab*          alpha_in, int sComp_alpha_in,
          int                      nComp)
    {
      const Real* dx = bd.getGeom().CellSize();
      ABecHelper* op = new ABecHelper(bd,dx);
      int nComp_alpha = (W_flag > 1 || alpha_in!=0 ?  nComp_W  :  1);
      const BoxArray& grids = volume.boxArray();
      if (a!=0) {
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
          for (int n=0; n<nComp_alpha; ++n) {
            MultiFab::Multiply(alpha,*W,sComp_W,n,1,alpha.nGrow());
          }
        }
        if (alpha_in) {
          for (int n=0; n<nComp_alpha; ++n) {
            MultiFab::Multiply(alpha,*alpha_in,sComp_alpha_in,n,1,alpha.nGrow());
          }
        }
        op->aCoefficients(alpha);
        alpha.clear();
      }
      op->setScalars(a, b);

      for (int d = 0; d < BL_SPACEDIM; d++) {
        MultiFab bcoeffs(area[d].boxArray(),nComp_beta,0);
        for (int n=0; n<nComp_beta; ++n) {
          MultiFab::Copy(bcoeffs,area[d],0,n,1,0);
        }
        for (MFIter mfi(bcoeffs); mfi.isValid(); ++mfi) {
          bcoeffs[mfi].mult((*beta[d])[mfi],sComp_beta,0,nComp_beta);
          bcoeffs[mfi].mult(dx[d]);
        }
        op->bCoefficients(bcoeffs,d);
      }
      return op;
    }

  } /* AmanziTransport */
} /* Amanzi */
