/*
  This is the flow component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_FLOW_DOMAIN_FUNCTION_HH_
#define AMANZI_FLOW_DOMAIN_FUNCTION_HH_

#include <string>
#include <vector>

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "unique_mesh_function.hh"

namespace Amanzi {
namespace Functions {

const int DOMAIN_FUNCTION_ACTION_NONE = 0;
const int DOMAIN_FUNCTION_ACTION_DISTRIBUTE_VOLUME = 1;
const int DOMAIN_FUNCTION_ACTION_DISTRIBUTE_PERMEABILITY = 2;

const int DOMAIN_FUNCTION_SUBMODEL_RATE = 0;
const int DOMAIN_FUNCTION_SUBMODEL_INTEGRAND = 1;

class FlowDomainFunction : public UniqueMeshFunction {
 public:
  explicit FlowDomainFunction(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) : 
      UniqueMeshFunction(mesh),
      finalized_(false) {}

  void Define(const std::vector<std::string>& regions,
              const Teuchos::RCP<const MultiFunction>& f,
              int action, int submodel);

  void Define(std::string& region,
              const Teuchos::RCP<const MultiFunction>& f,
              int action, int submodel);

  // source term on time interval (t0, t1]
  void Compute(double t0, double t1);
  void ComputeDistribute(double t0, double t1);
  void ComputeDistribute(double t0, double t1, double* weight);
  
  void Finalize() {};
 
  // iterator methods
  typedef std::map<int,double>::const_iterator Iterator;
  Iterator begin() const { return value_.begin(); }
  Iterator end() const  { return value_.end(); }
  Iterator find(const int j) const { return value_.find(j); }
  std::map<int,double>::size_type size() { return value_.size(); }

  // extract internal information
  int CollectActionsList();

 private:
  std::vector<int> actions_;
  std::vector<int> submodel_;

 protected:
  std::map<int,double> value_;
  bool finalized_;
};

}  // namespace Functions
}  // namespace Amanzi

#endif

