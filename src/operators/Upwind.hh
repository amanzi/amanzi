/*
  This is the operators component of the Amanzi code. 

  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_UPWIND_HH_
#define AMANZI_UPWIND_HH_

#include <string>
#include <vector>

#include "Epetra_IntVector.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "CompositeVector.hh"
#include "Mesh.hh"
#include "VerboseObject.hh"

#include "OperatorDefs.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* The base class for all upwind methods. 
****************************************************************** */ 

template<class Model>
class Upwind {
 public:
  Upwind() {};
  Upwind(Teuchos::RCP<const AmanziMesh::Mesh> mesh,
         Teuchos::RCP<const Model> model)
      : mesh_(mesh), model_(model) {};
  ~Upwind() {};

  // main methods
  virtual void Init(Teuchos::ParameterList& plist) = 0;

  virtual void Compute(const CompositeVector& flux,
                       const std::vector<int>& bc_model, const std::vector<double>& bc_value,
                       const CompositeVector& field, CompositeVector& field_upwind,
                       double (Model::*Value)(int, double) const) = 0;

 protected:
  Teuchos::RCP<VerboseObject> vo_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::RCP<const Model> model_;
};

}  // namespace Operators
}  // namespace Amanzi

#endif

