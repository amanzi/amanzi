/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  A coupling operator discretizes the two-point flux between 
  physics fields. At the moment, we tes the flux between two PKs.
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
                   std::shared_ptr<const std::vector<std::vector<int> > > row_inds,
                   std::shared_ptr<const std::vector<std::vector<int> > > col_inds,
                   const Teuchos::RCP<Operator> global_op = Teuchos::null) :
      PDE_HelperDiscretization(global_op)
  {
    Init_(plist, cvs_row, cvs_col, row_inds, col_inds);
  }
  ~PDE_CouplingFlux() {};

  // main members 
  // -- required by the interface
  virtual void UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& u,
          const Teuchos::Ptr<const CompositeVector>& p);
  
  // -- setup
  void Setup(std::shared_ptr<const std::vector<double> > K, double factor) {
    K_ = K;
    factor_ = factor;
  }

 private:
  void Init_(Teuchos::ParameterList& plist,
             const Teuchos::RCP<const CompositeVectorSpace>& cvs_row,
             const Teuchos::RCP<const CompositeVectorSpace>& cvs_col,
             std::shared_ptr<const std::vector<std::vector<int> > >& row_inds,
             std::shared_ptr<const std::vector<std::vector<int> > >& col_inds);

 protected:
  std::shared_ptr<const std::vector<double> > K_;

 private:
  double factor_;
};

}  // namespace Operators
}  // namespace Amanzi

#endif

