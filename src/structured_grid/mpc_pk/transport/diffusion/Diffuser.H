/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#ifndef _Diffuser_H_
#define _Diffuser_H_

//namespace Amanzi {
//  namespace AmanziTransport {

    template <class Vector, class Op>
    class DiffuserOp
    {
    public:
      DiffuserOp(Op* _op, Vector* _W) : op(_op), W(_W) {}

      Real ApplyOperator(Vector& Ax, Vector& x) {op->apply(Ax,x); return 0.0;}
      Real ScaleOperator() {
        if (op==0) return 1;
        Real a = op->get_alpha();
        Real b = op->get_beta();
        Real new_scale = 1.0/op->aCoefficients().max(0);
        op->setScalars(a*new_scale, b*new_scale);
        return new_scale;
      }
      void UnscaleOperator(Real scale) {
        if (op!=0) {
          Real a = op->get_alpha();
          Real b = op->get_beta();
          op->setScalars(a/scale, b/scale);
        }
      }
      void UnweightVector(const Vector& xw, Vector& xunw) {
	if (&xunw != &xw) {
	  xunw.COPY(xw);
	}
	if (W != 0) {
	  xunw.MULTAYI(*W);
	}
      }
      void WeightVector(const Vector& xunw, Vector& xw) {
	if (&xunw != &xw) {
	  xw.COPY(xunw);
	}
	if (W != 0) {
	  xw.MULTAY(*W);
	}
      }
      bool isValid() const {return op != 0;}
      Op& LinOp() {BL_ASSERT(op!=0); return *op;}
    protected:
      Op* op;
      Vector* W;
    };

    template <class Vector, class DiffOp>
    class LinSolver
    {
    public:
      LinSolver(DiffOp& _dop)
        : dop(_dop) {}
      void Solve(Vector& X, const Vector& Rhs, Real abs_tol, Real rel_tol) {BoxLib::Abort("No solver implemented for this DiffOp");}
      DiffOp& DiffuseOp() {return dop;}
      bool isValid() const {return dop.isValid();}
    protected:
      DiffOp& dop;
    };

    template <class Vector, class Op>
    class Diffuser
    {
    public:
      Diffuser(Op* _oldOp, Op* _newOp, const Vector& volume, Vector* _Wold, Vector* _Wnew)
        : opOld(new DiffuserOp<Vector,Op>(_oldOp,_Wold)),
          opNew(new DiffuserOp<Vector,Op>(_newOp,_Wnew)),
          Volume(volume), X(volume), Y(volume) {}
      ~Diffuser() {delete opOld; delete opNew;}
      void Diffuse(const Vector* Sold, Vector& Snew, const Vector& dRhs, Real t_old, Real t_new, Real abs_tol, Real rel_tol) {
        Real dt = t_new - t_old;
        if (Sold != 0) {
          X.COPY(*Sold);                     // S (=X) = Cn
          if (opOld->isValid()) {
            opOld->ApplyOperator(Y,X);       // Rhs (=Y) = Ln(Cn)
            opOld->WeightVector(X,X);        // S (=X) = phi.s.Cn
            X.AXPBY(dRhs,1,dt);              // Rhs (=X) = phi.s.Cn + dt * dRhs
          }
          X.MULTAY(Volume,1);                // Rhs (=X) = ( phi.s.Cn + dt * dRhs ) * Vol
          Y.AXPBY(X,1,1);                    // Rhs (=Y) = Ln(Cn) + ( phi.s.Cn + dt * dRhs ) * Vol
        } else {
          Y.COPY(dRhs);
          Y.MULTAY(Volume,dt);
        }
        Real scale = opNew->ScaleOperator(); // a,b <- a,b/scale
        Y.SCALE(scale);                      // Rhs <- Rhs.scale
        LinSolver<Vector, DiffuserOp<Vector, Op> > solver(*opNew);
        solver.Solve(Snew,Y,abs_tol,rel_tol);// Cnp1: phi.s.Cnp1.Vol + Lnp1(Cnp1) =  Ln(Cn) + (phi.s.Cn + dt * dRhs ) * Vol
        opNew->UnscaleOperator(scale);       // a,b <- a,b.scale
      }
    protected:
      DiffuserOp<Vector,Op> *opOld, *opNew;
      const Vector& Volume;
      Vector X, Y;
    };

//  } /* AmanziTransport */
//} /* Amanzi */

#endif /* _Diffuser_H_ */
