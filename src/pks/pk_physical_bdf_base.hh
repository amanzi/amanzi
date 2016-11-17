/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
//! Standard base for most implemented PKs, this combines both domains/meshes of PKPhysicalBase and BDF methods of PKBDFBase.

/*
  ATS is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/


/*!

A base class for all PKs that are both physical, in the sense that they
implement an equation and are not couplers, and support the implicit
integration interface.  This largely just supplies a default error norm based
on error in conservation relative to the extent of the conserved quantity.

* `"absolute error tolerance`" [double] **1.0**

  Absolute tolerance, :math:`a_tol` in the equation below.

* `"relative error tolerance`" [double] **1.0**

  Relative tolerance, :math:`r_tol` in the equation below.

By default, the error norm used by solvers is given by:

:math:`ENORM(u, du) = |du| / ( a_tol + r_tol * |u| )`

The defaults here are typically good, or else good defaults are set in the
code, so these need not be supplied.


NOTE: ``PKPhysicalBDFBase -->`` PKBDFBase_
      ``PKPhysicalBDFBase -->`` PKPhysicalBase_
      ``PKPhysicalBDFBase (v)-->`` PKDefaultBase_

*/


#ifndef AMANZI_PK_PHYSICAL_BDF_BASE_HH_
#define AMANZI_PK_PHYSICAL_BDF_BASE_HH_

#include "errors.hh"
#include "pk_default_base.hh"
#include "pk_bdf_base.hh"
#include "pk_physical_base.hh"

#include "Operator.hh"

namespace Amanzi {

class PKPhysicalBDFBase : public PKBDFBase, public PKPhysicalBase {

 public:
  PKPhysicalBDFBase(const Teuchos::RCP<Teuchos::ParameterList>& plist,
                    Teuchos::ParameterList& FElist,
                    const Teuchos::RCP<TreeVector>& solution) :
      PKDefaultBase(plist, FElist, solution),
      PKPhysicalBase(plist, FElist, solution),
      PKBDFBase(plist, FElist, solution) {}

  virtual void setup(const Teuchos::Ptr<State>& S);

  // initialize.  Note both BDFBase and PhysicalBase have initialize()
  // methods, so we need a unique overrider.
  virtual void initialize(const Teuchos::Ptr<State>& S);

  // Default preconditioner is Picard
  virtual int ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu) {
    *Pu = *u;
    return 0;
  }

  // updates the preconditioner, default does nothing
  virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h) {}

  // Default implementations of BDFFnBase methods.
  // -- Compute a norm on u-du and return the result.
  virtual double ErrorNorm(Teuchos::RCP<const TreeVector> u,
                       Teuchos::RCP<const TreeVector> du);

  // -- Experimental approach -- calling this indicates that the time
  //    integration scheme is changing the value of the solution in
  //    state.
  virtual void ChangedSolution();

  virtual double BoundaryValue(const Teuchos::RCP<const Amanzi::CompositeVector>& solution, int face_id);
  virtual int BoundaryDirection(int face_id);
  virtual void ApplyBoundaryConditions_(const Teuchos::Ptr<CompositeVector>& u);
  
  // PC operator access
  Teuchos::RCP<Operators::Operator> preconditioner() { return preconditioner_; }

  // BC access
  std::vector<int>& bc_markers() { return bc_markers_; }
  std::vector<double>& bc_values() { return bc_values_; }
  Teuchos::RCP<Operators::BCs> BCs() { return bc_; }

 protected:
  // PC
  Teuchos::RCP<Operators::Operator> preconditioner_;

  // BCs
  std::vector<int> bc_markers_;
  std::vector<double> bc_values_;
  Teuchos::RCP<Operators::BCs> bc_;

  // error criteria
  Key conserved_key_;
  Key cell_vol_key_;
  double atol_, rtol_, fluxtol_;

};


} // namespace

#endif
