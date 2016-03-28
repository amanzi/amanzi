/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Interface for the Constant Temp PK.  This PK simply provides a constant
temperature, and is provided for testing with other PKs that depend upon an
energy equation.  This could easily be provided by the state as an independent
variable, but this is nice for testing the full hierarchy with a simple PK.

Example usage:

  <ParameterList name="energy">
    <Parameter name="PK model" type="string" value="Constant Temperature"/>
    <Parameter name="Constant Temperature" type="double" value="290.0"/>
  </ParameterList>

------------------------------------------------------------------------- */

#ifndef PKS_ENERGY_CONSTANT_TEMPERATURE_HH_
#define PKS_ENERGY_CONSTANT_TEMPERATURE_HH_

// #include "pk_factory_ats.hh"
// #include "pk_physical_bdf_base.hh"
#include "PK_Factory.hh"
#include "PK_PhysicalBDF_ATS.hh"

namespace Amanzi {
namespace Energy {

class ConstantTemperature : public PK_PhysicalBDF_ATS {

public:

  ConstantTemperature(Teuchos::ParameterList& FElist,
                      const Teuchos::RCP<Teuchos::ParameterList>& plist,
                      const Teuchos::RCP<State>& S,
                      const Teuchos::RCP<TreeVector>& solution) :
    PK_Default(plist, FElist, solution),
    PK_PhysicalBDF_ATS(FElist, plist, S, solution) {
    plist_->set("solution key", "temperature");
  }

  // Virtual destructor
  virtual ~ConstantTemperature() {}

  // ConstantTemperature is a PK
  // -- Setup data
  virtual void Setup(const Teuchos::Ptr<State>& S);

  // -- Initialize owned (dependent) variables.
  virtual void Initialize(const Teuchos::Ptr<State>& S);

  // -- Commit any secondary (dependent) variables.
  virtual void CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S){};

  // -- Update diagnostics for vis.
  virtual void CalculateDiagnostics(const Teuchos::RCP<State>& S){};

  // -- advance via one of a few methods
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit);

  // ConstantTemperature is a BDFFnBase
  // computes the non-linear functional f = f(t,u,udot)
  virtual void Functional(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                   Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> f);

  // applies preconditioner to u and returns the result in Pu
  virtual int ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);

  // updates the preconditioner
  virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h);

private:
  // A few options for advance
  bool advance_analytic_(double dt);
  bool advance_bdf_(double dt);

  // initial temperature
  Teuchos::RCP<CompositeVector> temp0_;

  // misc setup information
  Teuchos::ParameterList energy_plist_;

  // time integration
  double atol_;
  double rtol_;

  // factory registration
  static RegisteredPKFactory<ConstantTemperature> reg_;
};

} // namespace
} // namespace

#endif
