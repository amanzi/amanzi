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
#include "Function.hh"
#include "mesh_function.hh"

namespace Amanzi {

const int DOMAIN_FUNCTION_ACTION_NONE = 0;
const int DOMAIN_FUNCTION_ACTION_DISTRIBUTE_VOLUME = 1;
const int DOMAIN_FUNCTION_ACTION_DISTRIBUTE_PERMEABILITY = 2;

class DomainFunction : public MeshFunction {
 public:
  DomainFunction(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) { mesh_ = mesh; }

  void Define(const std::vector<std::string>& regions, 
              const Teuchos::RCP<const Function>& f, int action);
  void Compute(double T);
  void ComputeDistribute(double T);
  void ComputeDistribute(double T, double* weight);

  // experimental support of multi-valued functions
  void DefineMultiValue(const std::vector<std::string>& regions,
                        const Teuchos::RCP<const Function>& f, int action, 
                        const std::string& name);

  void ComputeMultiValue(double T, const std::string& name);
  void ComputeDistributeMultiValue(double T, const std::string& name);
  void ComputeDistributeMultiValue(double T, const std::string& name, double* weight);

  // extract internal information
  int CollectActionsList();

 private:
  std::vector<int> actions_;
  std::vector<std::string> names_;
};

}  // namespace Amanzi

#endif  // AMANZI_DOMAIN_FUNCTION_HH_
