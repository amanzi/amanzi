/*
  This is the operators component of the Amanzi code.

  License: BSD
  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_OPERATOR_VERIFICATION_HH_
#define AMANZI_OPERATOR_VERIFICATION_HH_

#include "Mesh.hh"

#include "CompositeVector.hh"
#include "Operator.hh"
#include "TreeOperator.hh"
#include "TreeVector.hh"

template <class Vector, class Operator>
class Verification {
 public:
  Verification(Teuchos::RCP<const Operator> op) : op_(op) {};
  ~Verification() {};

  void CheckMatrixSPD(bool symmetry = true, bool pos_def = true) {
    Vector a(op_->DomainMap()), ha(a), b(a), hb(a);

    a.Random();
    b.Random();
    op_->Apply(a, ha);
    op_->Apply(b, hb);

    double ahb, bha, aha, bhb;
    a.Dot(hb, &ahb);
    b.Dot(ha, &bha);
    a.Dot(ha, &aha);
    b.Dot(hb, &bhb);

    if (a.Comm().MyPID() == 0) {
      std::cout << "Matrix:\n";
      if (symmetry)
          std::cout << "  Symmetry test: " << ahb << " = " << bha << std::endl;
      if (pos_def)
          std::cout << "  Positivity test: " << aha << " " << bhb << std::endl;
    } 
    if (symmetry) CHECK_CLOSE(ahb, bha, 1e-12 * fabs(ahb));
    if (pos_def) {
      CHECK(aha > 0.0);
      CHECK(bhb > 0.0);
    }
  }

  void CheckPreconditionerSPD(bool symmetry = true, bool pos_def = true) {
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

    if (a.Comm().MyPID() == 0) {
      int size = (op_->A() != Teuchos::null) ? op_->A()->NumGlobalRows() : -1;
      std::cout << "Preconditioner: size=" << size << "\n";
      if (symmetry)
          std::cout << "  Symmetry test: " << ahb << " = " << bha << std::endl;
      if (pos_def)
          std::cout << "  Positivity test: " << aha << " " << bhb << std::endl;
    } 
    if (symmetry) CHECK_CLOSE(ahb, bha, 1e-12 * fabs(ahb));
    if (pos_def) {
      CHECK(aha > 0.0);
      CHECK(bhb > 0.0);
    }
  }

  void CheckResidual(const Vector x, const Vector b, double tol) {
    Vector r(b);

    op_->ApplyAssembled(x, r);
    r.Update(1.0, b, -1.0);

    double tmp;
    r.Dot(r, &tmp);
    CHECK_CLOSE(0.0, tmp, tol * tol);
  }

  void CheckResidual(const Vector x, double tol) {
    Vector b(*op_->rhs());
    CheckResidual(x, b, tol);
  }

 private:
  Teuchos::RCP<const Operator> op_;
};

typedef Verification<Amanzi::CompositeVector, Amanzi::Operators::Operator> VerificationCV;
typedef Verification<Amanzi::TreeVector, Amanzi::Operators::TreeOperator> VerificationTV;

#endif

