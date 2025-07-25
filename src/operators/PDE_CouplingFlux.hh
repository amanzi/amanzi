/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

  A coupling operator discretizes the two-point flux between
  physics fields. At the moment, we test the flux between two PKs.
*/

#ifndef AMANZI_OPERATOR_PDE_COUPLING_FLUX_HH_
#define AMANZI_OPERATOR_PDE_COUPLING_FLUX_HH_

#include <memory>
#include <string>
#include <vector>

// Amanzi
#include "BilinearForm.hh"
#include "Polynomial.hh"
#include "Tensor.hh"
#include "WhetStoneDefs.hh"

// Operators
#include "PDE_HelperDiscretization.hh"

namespace Amanzi {
namespace Operators {

class PDE_CouplingFlux : public PDE_HelperDiscretization {
 public:
  PDE_CouplingFlux(Teuchos::ParameterList& plist,
                   const Teuchos::RCP<const CompositeVectorSpace>& cvs_row,
                   const Teuchos::RCP<const CompositeVectorSpace>& cvs_col,
                   std::shared_ptr<const std::vector<std::vector<int>>> row_inds,
                   std::shared_ptr<const std::vector<std::vector<int>>> col_inds,
                   const Teuchos::RCP<Operator> global_op = Teuchos::null)
  {
    global_op_ = global_op;
    Init_(plist, cvs_row, cvs_col, row_inds, col_inds);
  }
  ~PDE_CouplingFlux() {};

  // main members
  // -- required by the interface
  using PDE_HelperDiscretization::UpdateMatrices;
  virtual void UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& u,
                              const Teuchos::Ptr<const CompositeVector>& p) override;

  // -- setup
  void Setup(std::shared_ptr<const std::vector<double>> K, double factor)
  {
    K_ = K;
    factor_ = factor;
  }

  // optional calculation of flux from potential p
  virtual void UpdateFlux(const Teuchos::Ptr<const CompositeVector>& p,
                          const Teuchos::Ptr<CompositeVector>& u) override {};

  // not implemented
  virtual void ApplyBCs(bool primary, bool eliminate, bool essential_eqn) override
  {
    Errors::Message msg("Coupling operator does not support boundary conditions.");
    Exceptions::amanzi_throw(msg);
  }

 private:
  void Init_(Teuchos::ParameterList& plist,
             const Teuchos::RCP<const CompositeVectorSpace>& cvs_row,
             const Teuchos::RCP<const CompositeVectorSpace>& cvs_col,
             std::shared_ptr<const std::vector<std::vector<int>>>& row_inds,
             std::shared_ptr<const std::vector<std::vector<int>>>& col_inds);

 protected:
  std::shared_ptr<const std::vector<double>> K_;

 private:
  double factor_;
};

} // namespace Operators
} // namespace Amanzi

#endif
