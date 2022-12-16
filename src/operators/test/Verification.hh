/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  This is the operators component of the Amanzi code.

*/

#ifndef AMANZI_OPERATOR_VERIFICATION_HH_
#define AMANZI_OPERATOR_VERIFICATION_HH_

#include "Mesh.hh"

#include "CompositeVector.hh"
#include "DenseMatrix.hh"
#include "Operator.hh"
#include "TreeOperator.hh"
#include "TreeVector.hh"

template <class Vector, class Operator>
class Verification {
 public:
  Verification(Teuchos::RCP<const Operator> op) : op_(op){};
  ~Verification(){};

  void
  CheckMatrixSPD(bool symmetry = true, bool pos_def = true, int nloops = 2, bool assembled = false)
  {
    Vector a(op_->DomainMap()), ha(a), b(a), hb(a);

    for (int n = 0; n < nloops; n++) {
      a.Random();
      b.Random();
      if (assembled) {
        op_->ApplyAssembled(a, ha);
        op_->ApplyAssembled(b, hb);
      } else {
        op_->Apply(a, ha);
        op_->Apply(b, hb);
      }

      double ahb, bha, aha, bhb;
      a.Dot(hb, &ahb);
      b.Dot(ha, &bha);
      a.Dot(ha, &aha);
      b.Dot(hb, &bhb);

      if (a.Comm()->MyPID() == 0) {
        if (n == 0) std::cout << "Matrix (assembled=" << assembled << "):\n";
        if (symmetry) printf("  Symmetry test: %21.14e = %21.14e\n", ahb, bha);
        if (pos_def) std::cout << "  Positivity test: " << aha << " " << bhb << std::endl;
      }
      if (symmetry) CHECK_CLOSE(ahb, bha, 1e-12 * fabs(ahb));
      if (pos_def) {
        CHECK(aha > 0.0);
        CHECK(bhb > 0.0);
      }
    }
  }

  void CheckPreconditionerSPD(double rtol = 1e-12, bool symmetry = true, bool pos_def = true)
  {
    Vector a(op_->DomainMap()), ha(a), b(a), hb(a);

    a.Random();
    b.Random();
    op_->ApplyInverse(a, ha);
    op_->ApplyInverse(b, hb);

    double ahb, bha, aha, bhb;
    a.Dot(hb, &ahb);
    b.Dot(ha, &bha);
    a.Dot(ha, &aha);
    b.Dot(hb, &bhb);

    if (a.Comm()->MyPID() == 0) {
      int size = (op_->A() != Teuchos::null) ? op_->A()->NumGlobalRows() : -1;
      std::cout << "Preconditioner: size=" << size << "\n";
      if (symmetry) printf("  Symmetry test: %21.14e = %21.14e\n", ahb, bha);
      if (pos_def) std::cout << "  Positivity test: " << aha << " " << bhb << std::endl;
    }
    if (symmetry) CHECK_CLOSE(ahb, bha, rtol * fabs(ahb));
    if (pos_def) {
      CHECK(aha > 0.0);
      CHECK(bhb > 0.0);
    }
  }

  void CheckResidual(const Vector x, const Vector b, double tol)
  {
    Vector r(b);

    op_->ApplyAssembled(x, r);
    r.Update(1.0, b, -1.0);

    double tmp, xnorm;
    r.Dot(r, &tmp);
    x.Dot(x, &xnorm);
    CHECK_CLOSE(0.0, tmp, tol * tol * xnorm * xnorm);
  }

  void CheckResidual(const Vector x, double tol)
  {
    Vector b(*op_->rhs());
    CheckResidual(x, b, tol);
  }


  void CheckSpectralBounds(int kernel_dim = 1)
  {
    double emin = 1e+99, emax = -1e+99;
    double cndmin = 1e+99, cndmax = 1.0, cndavg = 0.0;

    for (auto it = op_->begin(); it != op_->end(); ++it) {
      auto& matrices = (*it)->matrices;
      for (int i = 0; i < matrices.size(); i++) {
        Amanzi::WhetStone::DenseMatrix Acell(matrices[i]);
        int n = Acell.NumRows();

        int info, ldv(1), lwork = 4 * n;
        double VL, VR, dmem2[n], dwork[lwork];
        Amanzi::WhetStone::DenseVector dmem1(n);

        Amanzi::WhetStone::DGEEV_F77("N",
                                     "N",
                                     &n,
                                     Acell.Values(),
                                     &n,
                                     dmem1.Values(),
                                     dmem2,
                                     &VL,
                                     &ldv,
                                     &VR,
                                     &ldv,
                                     dwork,
                                     &lwork,
                                     &info);

        OrderByIncrease_(n, dmem1.Values());

        int m = std::min(kernel_dim, dmem1.NumRows() - 1);
        double e, a = dmem1(m), b = dmem1(m); // skip the kernel
        for (int k = m + 1; k < n; k++) {
          e = dmem1(k);
          a = std::min(a, e);
          b = std::max(b, e);
        }

        emin = std::min(emin, a);
        emax = std::max(emax, b);

        double cnd = fabs(b) / (fabs(a) + 1e-16);
        cndmin = std::min(cndmin, cnd);
        cndmax = std::max(cndmax, cnd);
        cndavg += cnd;
      }
      cndavg /= matrices.size();

      int MyPID = 0;
      if (MyPID == 0) {
        printf("  Matrices:\n");
        printf("    eigenvalues (min,max) = %8.3g %8.3g\n", emin, emax);
        printf("    conditioning (min,max,avg) = %8.2g %8.2g %8.2g\n", cndmin, cndmax, cndavg);
      }
    }
  }

 private:
  void OrderByIncrease_(int n, double* mem)
  {
    for (int i = 0; i < n; i++) {
      for (int j = 1; j < n - i; j++) {
        if (mem[j - 1] > mem[j]) {
          double tmp = mem[j];
          mem[j] = mem[j - 1];
          mem[j - 1] = tmp;
        }
      }
    }
  }

 private:
  Teuchos::RCP<const Operator> op_;
};

typedef Verification<Amanzi::CompositeVector, Amanzi::Operators::Operator> VerificationCV;
typedef Verification<Amanzi::TreeVector, Amanzi::Operators::TreeOperator> VerificationTV;

#endif
