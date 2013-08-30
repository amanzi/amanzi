#ifndef _Diffuser_H_
#define _Diffuser_H_

//namespace Amanzi {
//  namespace AmanziTransport {

    template <class Vector, class Op>
    class DiffuserOp
    {
    public:
      DiffuserOp(Op& _op, Vector& _W) : op(_op), W(_W) {}
      Real ApplyOperator(Vector& Ax, Vector& x) {op.apply(Ax,x);}
      Real ScaleOperator() {
        Real a = op.get_alpha();
        Real b = op.get_beta();
        Real new_scale = 1.0/op.aCoefficients().max(0);
        op.setScalars(a*new_scale, b*new_scale);
        return new_scale;
      }
      void UnscaleOperator(Real scale) {
        Real a = op.get_alpha();
        Real b = op.get_beta();
        op.setScalars(a/scale, b/scale);
      }
      void UnweightVector(const Vector& xw, Vector& xunw) {
        xunw.COPY(xw);
        xunw.MULTAYI(W);
      }
      void WeightVector(const Vector& xunw, Vector& xw) {
        xw.COPY(xunw);
        xw.MULTAY(W);
      }
      Op& LinOp() {return op;}
    protected:
      Op& op;
      Vector& W;
    };

    template <class Vector, class DiffOp>
    class LinSolver
    {
    public:
      LinSolver(DiffOp& _dop)
        : dop(_dop) {}
      void solve(Vector& X, const Vector& Rhs) {}
    protected:
      DiffOp& dop;
    };

#include <VisMF.H>

    template <class Vector, class Op>
    class Diffuser
    {
    public:
      Diffuser(Op& _oldOp, Op& _newOp, const Vector& volume, Vector& _Wold, Vector& _Wnew)
        : opOld(new DiffuserOp<Vector,Op>(_oldOp,_Wold)),
          opNew(new DiffuserOp<Vector,Op>(_newOp,_Wnew)),
          Volume(volume), X(volume.Clone()), Y(volume.Clone()) {}
      ~Diffuser() {delete X; delete Y;}
      void Diffuse(const Vector& Sold, Vector& S_new, const Vector& dRhs, Real t_old, Real t_new) {
        Real dt = t_new - t_old;
        opOld->UnweightVector(Sold,*X);      // S (=X) = phi.s.Cn/phi.s = Cn
        opOld->ApplyOperator(*Y,*X);         // Rhs (=Y) = Ln(Cn)
        opOld->WeightVector(*X,*X);          // S (=X) = phi.s.Cn
        X->AXPBY(dRhs,1,dt);                 // Rhs (=X) = phi.s.Cn + dt * dRhs
        X->MULTAY(Volume,1);                 // Rhs (=X) = ( phi.s.Cn + dt * dRhs ) * Vol
        Y->AXPBY(*X,1,1);                    // Rhs (=Y) = Ln(Cn) + ( phi.s.Cn + dt * dRhs ) * Vol
        opNew->UnweightVector(S_new,*Y);     // Sguess = phi.s.Cnp1/phi.s = Cnp1
        Real scale = opNew->ScaleOperator(); //   a <- a/scale
        X->SCALE(scale);                     // Rhs <- Rhs.scale 
        LinSolver<Vector, DiffuserOp<Vector, Op> > solver(*opNew);
        solver.solve(*Y,*X);                 // Cnp1: Lnp1(Cnp1) =  Ln(Cn) + (phi.s.Cn + dt * dRhs ) * Vol
        opNew->UnscaleOperator(scale);       // a <- a.scale
        opNew->WeightVector(*Y,S_new);       // Snew = phi.s.Cnp1
      }
      void getForcing(const Vector& S, Real t_eval) {}
    protected:
      DiffuserOp<Vector,Op> *opOld, *opNew;
      const Vector& Volume;
      Vector *X, *Y;
    };

//  } /* AmanziTransport */
//} /* Amanzi */

#endif /* _Diffuser_H_ */
