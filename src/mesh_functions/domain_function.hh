/*
This is the flow component of the Amanzi code. 

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided Reconstruction.cppin the top-level COPYRIGHT file.

Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_DOMAIN_FUNCTION_HH_
#define AMANZI_DOMAIN_FUNCTION_HH_

#include <vector>
#include <string>

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "function.hh"
#include "mesh_function.hh"

namespace Amanzi {

class DomainFunction : public MeshFunction {
 public:
  DomainFunction(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) { mesh_ = mesh; }

  void Define(const std::vector<std::string>& regions, const Teuchos::RCP<const Function>& f);
  void Compute(double T);
  void ComputeDistribute(double T);
  void ComputeDistribute(double T, double* weight);
};

}  // namespace Amanzi

#endif  // AMANZI_DOMAIN_FUNCTION_HH_
