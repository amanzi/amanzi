/*
  Multiphase PK

  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_FLATTENED_TREE_OPERATOR_HH_
#define AMANZI_FLATTENED_TREE_OPERATOR_HH_

#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "TreeVector.hh"
#include "TreeVectorSpace.hh"
#include "TreeOperator.hh"
#include "VerboseObject.hh"

#include "Operator.hh"

/* ******************************************************************
  TreeOperator is the block analogue of Operators -- it provides a
  linear operator acting on a TreeVectorSpace (TVS). This class 
  flattens a TVS by creting as many blocks as there are copies of
  a CompositeVectorSpace in the TVS.
****************************************************************** */ 

namespace Amanzi {
namespace Operators {

class SuperMap;

class FlattenedTreeOperator : public Operators::TreeOperator {
 public:
  using Vector_t = TreeVector;
  using VectorSpace_t = TreeVector::VectorSpace_t;

  FlattenedTreeOperator() {};
  FlattenedTreeOperator(Teuchos::RCP<const TreeVectorSpace> tvs);
  ~FlattenedTreeOperator() {};

  // modified algorithms that use two supermaps
  virtual void SymbolicAssembleMatrix();
  virtual void AssembleMatrix();

  // only assembled matrix is allowed
  virtual int Apply(const TreeVector& X, TreeVector& Y) const;

 private:
  Teuchos::RCP<TreeVectorSpace> tvs_flat_;
  Teuchos::RCP<SuperMap> smap_flat_;
};

}  // namespace Operators
}  // namespace Amanzi


#endif


