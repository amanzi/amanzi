/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#ifndef _ABecHelper_H_
#define _ABecHelper_H_

#include <ABecLaplacian.H>

class ABecHelper
  : public ABecLaplacian
{
public:
    ABecHelper (const BndryData& _bd,
                const Real*      _h)
      : ABecLaplacian(_bd,_h), default_bndryComp(0), default_alphaComp(0), default_betaComp(0) {}
  void set_bndryComp(int comp) {default_bndryComp = comp;}
  void set_alphaComp(int comp) {default_alphaComp = comp;}
  void set_betaComp(int comp) {default_betaComp = comp;}
  int bndryComp() const {return default_bndryComp;}
  int alphaComp() const {return default_alphaComp;}
  int betaComp() const {return default_betaComp;}
  virtual ~ABecHelper () {}

  virtual void compFlux (D_DECL(MultiFab &xflux, MultiFab &yflux, MultiFab &zflux),
                         MultiFab& in, const BC_Mode& bc_mode=LinOp::Inhomogeneous_BC,
                         int src_comp=0, int dst_comp=0, int num_comp=1, int bnd_comp=-1);

  void compFlux (D_DECL(MultiFab &xflux, MultiFab &yflux, MultiFab &zflux),
                 MultiFab& in, bool do_ApplyBC, const BC_Mode& bc_mode=LinOp::Inhomogeneous_BC,
                 int src_comp=0, int dst_comp=0, int num_comp=1, int bnd_comp=-1);


  //
  // Applies level LinOp to "in", returns "out", uses BC_mode flag.
  //
  virtual void apply (MultiFab&      out,
                      MultiFab&      in,
                      int            level   = 0,
                      LinOp::BC_Mode bc_mode = LinOp::Inhomogeneous_BC,
                      bool           local   = false,
                      int            src_comp = 0,
                      int            dst_comp = 0,
                      int            num_comp = 1,
                      int            bndry_comp = -1);
  //
  // Fills level boundary cells using BC_mode flag, int. BC data if reqd.
  //
  virtual void applyBC (MultiFab&      inout,
                        int            src_comp,
                        int            num_comp,
                        int            level   = 0,
                        LinOp::BC_Mode bc_mode = LinOp::Inhomogeneous_BC,
                        bool           local   = false,
                        int            bndryComp = -1);
  //
  // Compute the level residual = rhsL - L(solnL).
  //
  virtual void residual (MultiFab&       residL,
                         const MultiFab& rhsL,
                         MultiFab&       solnL,
                         int             level   = 0,
                         LinOp::BC_Mode  bc_mode = LinOp::Inhomogeneous_BC,
                         bool            local   = false);
  //
  // Smooth the level system L(solnL)=rhsL.
  //
  void smooth (MultiFab&       solnL,
               const MultiFab& rhsL,
               int             level   = 0,
               LinOp::BC_Mode  bc_mode = LinOp::Inhomogeneous_BC);

  void jacobi_smooth (MultiFab&       solnL,
                      const MultiFab& rhsL,
                      int             level   = 0,
                      LinOp::BC_Mode  bc_mode = LinOp::Inhomogeneous_BC);
  //
  // Estimate the norm of the operator.
  //
  virtual Real norm (int nm = 0, int level = 0, const bool local = false);

protected:
  //
  // compute out=L(in) at level=level
  //
  virtual void Fapply (MultiFab&       out,
                       const MultiFab& in,
                       int             level);
  virtual void Fapply (MultiFab&       out,
                       int             dst_comp,
                       const MultiFab& in,
                       int             src_comp,
                       int             num_comp,
                       int             level);
  //
  // apply GSRB smoother to improve residual to L(solnL)=rhsL
  //
  virtual void Fsmooth (MultiFab&       solnL,
                        const MultiFab& rhsL,
                        int             level,
                        int             rgbflag);
  //
  // apply Jacobi smoother to improve residual to L(solnL)=rhsL
  //
  virtual void Fsmooth_jacobi (MultiFab&       solnL,
                               const MultiFab& rhsL,
                               int             level);

protected:
  int default_bndryComp;
  int default_alphaComp;
  int default_betaComp;
};

#endif
