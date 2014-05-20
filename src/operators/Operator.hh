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

/* ******************************************************************
1. Operator is a linear operator acting from linear space X to linear
space Y. These spaces are described by CompositeVectors (CV). A few
maps X->Y is supported. 

At the moment X = Y. Extension to TreeVectors should not be done in 
this class.

2. Operator is an un-ordered additive collection of lower-rank (or 
equal) simple operators. During its construction, an operator can 
only grow by assimulating more operators. 

At the moment, an operator cannot be split into two operators, but
there are no desing restriction for doing it in the future.

3. A simple operator (a set of 1 operators) is defined by triple:
scheme + elemental matrices + diagonal matrix. The schema specifies
structure of elemental matrices, e.g. cell-based matrices 
representing interation between face-based unknowns.

4. Operator can be converted to Epetra_FECrsMatrix matrix to generate
a preconditioner. This operation cannot be applied to a subset of
defining operators. 
 
****************************************************************** */ 

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
  void Clone(const Operator& op);
  int Apply(const CompositeVector& X, CompositeVector& Y) const;
  virtual int ApplyInverse(const CompositeVector& X, CompositeVector& Y) const;

  void SymbolicAssembleMatrix(int schema);
  virtual void AssembleMatrix(int schema);

  virtual void ApplyBCs(std::vector<int>& bc_model, std::vector<double>& bc_values);

  const CompositeVectorSpace& DomainMap() const { return *cvs_; }
  const CompositeVectorSpace& RangeMap() const { return *cvs_; }

  void CreateCheckPoint();
  void RestoreCheckPoint();

  void AddAccumulationTerm(const CompositeVector& u0, const CompositeVector& ss, double dT);
  void AddAccumulationTerm(const CompositeVector& u0, const CompositeVector& ss);

  // preconditioners
  virtual void InitPreconditioner(const std::string& prec_name, const Teuchos::ParameterList& plist,
                                  std::vector<int>& bc_model, std::vector<double>& bc_values);

  // access
  Teuchos::RCP<CompositeVector>& rhs() { return rhs_; }
  bool data_validity() { return data_validity_; }

 public:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::RCP<const CompositeVectorSpace> cvs_;
  mutable bool data_validity_;

  std::vector<Teuchos::RCP<std::vector<WhetStone::DenseMatrix> > > blocks_;
  std::vector<Teuchos::RCP<std::vector<WhetStone::DenseMatrix> > > blocks_shadow_;
  std::vector<int> blocks_schema_;
  Teuchos::RCP<CompositeVector> diagonal_, diagonal_checkpoint_;
  Teuchos::RCP<CompositeVector> rhs_, rhs_checkpoint_;

 public:
  int ncells_owned, nfaces_owned, nnodes_owned;
  int ncells_wghost, nfaces_wghost, nnodes_wghost;
 
  Teuchos::RCP<Epetra_FECrsMatrix> A_;
  Teuchos::RCP<AmanziPreconditioners::Preconditioner> preconditioner_;
  int offset_global_[3], offset_my_[3];
};

}  // namespace Operators
}  // namespace Amanzi


#endif


