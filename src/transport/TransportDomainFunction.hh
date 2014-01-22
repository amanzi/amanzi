/*
  This is the flow component of the Amanzi code. 

  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided Reconstruction.cppin the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_TRANSPORT_DOMAIN_FUNCTION_HH_
#define AMANZI_TRANSPORT_DOMAIN_FUNCTION_HH_

#include <vector>
#include <string>

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "unique_mesh_function.hh"

namespace Amanzi {
namespace AmanziTransport {

namespace TransportActions {
const int DOMAIN_FUNCTION_ACTION_NONE = 0;
const int DOMAIN_FUNCTION_ACTION_DISTRIBUTE_VOLUME = 1;
const int DOMAIN_FUNCTION_ACTION_DISTRIBUTE_PERMEABILITY = 2;
}

class TransportDomainFunction : public Functions::UniqueMeshFunction {
 public:
  explicit TransportDomainFunction(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) : 
      UniqueMeshFunction(mesh),
      finalized_(false) {}
  
  void Define(const std::vector<std::string>& regions,
              const Teuchos::RCP<const MultiFunction>& f,
              int action, 
              const std::string& name);

  void Define(const std::string region,
              const Teuchos::RCP<const MultiFunction>& f,
              int action,
              const std::string& name);
  
  void Compute(double time);
  void ComputeDistribute(double T);
  void ComputeDistribute(double T, double* weight);
  void ComputeDistributeMultiValue(double T, const std::string& name);
  void ComputeDistributeMultiValue(double T, const std::string& name, double* weight);
  
  void Finalize() {}
 
  // extract internal information
  int CollectActionsList();
 
  // iterator methods
  typedef std::map<int,double>::const_iterator Iterator;
  Iterator begin() const { return value_.begin(); }
  Iterator end() const  { return value_.end(); }
  Iterator find(const int j) const { return value_.find(j); }

 private:
  std::vector<int> actions_;
  std::vector<std::string> names_;

 protected:
  std::map<int,double> value_;
  bool finalized_;
};


}  // namespace AmanziTransport
}  // namespace Amanzi

#endif  // AMANZI_TRANSPORT_DOMAIN_FUNCTION_HH_
