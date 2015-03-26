/*
  This is the operators component of the Amanzi code.

  License: BSD
  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_OPERATOR_VERIFICATION_HH_
#define AMANZI_OPERATOR_VERIFICATION_HH_

#include "Mesh.hh"

class Verification {
 public:
  Verification(Teuchos::RCP<const Amanzi::Operators::Operator> op) : op_(op) {
    MyPID = op_->A()->DomainMap().Comm().MyPID();
  };
  ~Verification() {};

  void CheckMatrixSPD(bool symmetry = true, bool pos_def = true) {
    const Amanzi::CompositeVectorSpace& cvs = op_->DomainMap();
    Amanzi::CompositeVector a(cvs), ha(cvs), b(cvs), hb(cvs);

    a.Random();
    b.Random();
    op_->Apply(a, ha);
    op_->Apply(b, hb);

    double ahb, bha, aha, bhb;
    a.Dot(hb, &ahb);
    b.Dot(ha, &bha);
    a.Dot(ha, &aha);
    b.Dot(hb, &bhb);

    if (MyPID == 0) {
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
    const Amanzi::CompositeVectorSpace& cvs = op_->DomainMap();
    Amanzi::CompositeVector a(cvs), ha(cvs), b(cvs), hb(cvs);

    a.Random();
    b.Random();
    op_->ApplyInverse(a, ha);
    op_->ApplyInverse(b, hb);

    double ahb, bha, aha, bhb;
    a.Dot(hb, &ahb);
    b.Dot(ha, &bha);
    a.Dot(ha, &aha);
    b.Dot(hb, &bhb);

    if (MyPID == 0) {
      std::cout << "Preconditioner: size=" << op_->A()->NumGlobalRows() << "\n";
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

 private:
  Teuchos::RCP<const Amanzi::Operators::Operator> op_;
  int MyPID;
};

#endif

