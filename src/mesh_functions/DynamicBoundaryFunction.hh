/*
  Mesh Functions

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: (v1) Neil Carlson
           (v2) Ethan Coon

  Function applied to a mesh component with at most one function 
  application per entity.
 
  Aamanzi is no longer using this function!
*/

#ifndef AMANZI_DYNAMICBOUNDARY_FUNCTION_HH_
#define AMANZI_DYNAMICBOUNDARY_FUNCTION_HH_

#include <vector>
#include <map>
#include <string>

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "MultiFunction.hh"
#include "UniqueMeshFunction.hh"
#include "BoundaryFunction.hh"

namespace Amanzi {
namespace Functions {

class DynamicBoundaryFunction : public BoundaryFunction {

public:
  DynamicBoundaryFunction(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) :
    BoundaryFunction(mesh) {};

  void AddFunction(const Teuchos::RCP<BoundaryFunction>& f);
  
  Teuchos::RCP<BoundaryFunction> GetFunction(int id){return func_[id];}

  int Func_ID(double time);

  void Compute(double time);

protected:
  //std::map<int,double> value_;
  //bool finalized_;
  std::vector< Teuchos::RCP<BoundaryFunction> > func_;
};

} // namespace
} // namespace


#endif
