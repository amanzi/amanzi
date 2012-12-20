/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  A high level test class for the MatrixMFD operator.

  License: BSD
  Authors: Ethan Coon (ATS version) (ecoon@lanl.gov)
*/

#ifndef PK_TESTS_DIVGRAD_TEST_HH_
#define PK_TESTS_DIVGRAD_TEST_HH_

#include "boundary_function.hh"
#include "matrix_mfd.hh"

#include "pk_factory.hh"
#include "pk_physical_base.hh"

namespace Amanzi {

// forward declarations


namespace TestPKs {

class DivGradTest : public PKPhysicalBase {

public:
  DivGradTest(Teuchos::ParameterList& plist,
           const Teuchos::RCP<TreeVector>& solution) :
      PKDefaultBase(plist,solution),
      PKPhysicalBase(plist, solution) {
    // set a few parameters before setup
    plist_.set("solution key", "solution");
  }

  // Virtual destructor
  virtual ~DivGradTest() {}

  // main methods
  // -- Setup data.
  virtual void setup(const Teuchos::Ptr<State>& S);

  // -- Initialize owned (dependent) variables.
  virtual void initialize(const Teuchos::Ptr<State>& S);

  // -- Commit any secondary (dependent) variables.
  virtual void commit_state(double dt, const Teuchos::RCP<State>& S) {}

  // -- Update diagnostics for vis.
  virtual void calculate_diagnostics(const Teuchos::RCP<State>& S) {}

  virtual bool advance(double dt) {}

  virtual double get_dt() { return 1.e99; }

protected:
  // boundary condition members
  virtual void UpdateBoundaryConditions_();
  virtual void ApplyBoundaryConditions_(const Teuchos::RCP<CompositeVector>& pres);
  virtual bool TestRegularFaceValues_(const Teuchos::RCP<CompositeVector>& pres);

protected:
  // mathematical operators
  Teuchos::RCP<Operators::MatrixMFD> matrix_;

  // boundary condition data
  Teuchos::RCP<Functions::BoundaryFunction> bc_dirichlet_;
  Teuchos::RCP<Functions::BoundaryFunction> bc_neumann_;
  std::vector<Operators::Matrix_bc> bc_markers_;
  std::vector<double> bc_values_;

 private:
  // factory registration
  static RegisteredPKFactory<DivGradTest> reg_;
};

}  // namespace
}  // namespace Amanzi

#endif
