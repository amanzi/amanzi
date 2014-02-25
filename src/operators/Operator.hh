/*
  This is the Operator component of the Amanzi code.

  License: BSD
  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)

  The operators can be initialized from other operators.
  Since data are never copied by default, we have to track 
  down the ownership of data.
*/

#ifndef AMANZI_OPERATOR_HH_
#define AMANZI_OPERATOR_HH_

#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_FECrsMatrix.h"

#include "Mesh.hh"
#include "DenseMatrix.hh"
#include "CompositeVector.hh"
#include "CompositeVectorSpace.hh"
#include "Preconditioner.hh"
#include "OperatorTypeDefs.hh"

namespace Amanzi {
namespace Operators {

class Operator {
 public:
  Operator() {};
  Operator(Teuchos::RCP<const CompositeVectorSpace> cvs, int dummy);
  Operator(const Operator& op);
  ~Operator() {};

  // main members
  void Init();
  int Apply(const CompositeVector& X, CompositeVector& Y) const;
  virtual int ApplyInverse(const CompositeVector& X, CompositeVector& Y) const { Y = X; }

  void SymbolicAssemble();
  void SymbolicAssembleNodes();
  void SymbolicAssembleFaces();

  void AssembleStencilMFD_Nodes();
  void AssembleStencilMFD_Faces();

  void ApplyBCs(std::vector<int>& bc_model, std::vector<double>& bc_values);

  const CompositeVectorSpace& DomainMap() const { return *cvs_; }
  const CompositeVectorSpace& RangeMap() const { return *cvs_; }

  // preconditioners
  void InitPreconditioner(const std::string& prec_name, const Teuchos::ParameterList& plist);

  // access
  Teuchos::RCP<CompositeVector>& rhs() { return rhs_; }

 public:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::RCP<const CompositeVectorSpace> cvs_;
  mutable bool data_validity_;

  std::vector<Teuchos::RCP<std::vector<WhetStone::DenseMatrix> > > matrix_blocks_;
  std::vector<int> matrix_blocks_type_;
  Teuchos::RCP<CompositeVector> diagonal_block_;

  Teuchos::RCP<CompositeVector> rhs_;

 public:
  int ncells_owned, nfaces_owned, nnodes_owned;
  int ncells_wghost, nfaces_wghost, nnodes_wghost;
 
  Teuchos::RCP<Epetra_FECrsMatrix> A_;
  Teuchos::RCP<AmanziPreconditioners::Preconditioner> preconditioner_;
};

}  // namespace Operators
}  // namespace Amanzi


#endif


