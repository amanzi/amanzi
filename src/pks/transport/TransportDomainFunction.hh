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

#include "CommonDefs.hh"
#include "Mesh.hh"
#include "unique_mesh_function.hh"

namespace Amanzi {
namespace Transport {

class TransportDomainFunction : public Functions::UniqueMeshFunction {
 public:
  explicit TransportDomainFunction(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) : 
      UniqueMeshFunction(mesh),
      finalized_(false) {}
  
  void Define(const std::vector<std::string>& regions,
              const Teuchos::RCP<const MultiFunction>& f,
              int action, int submodel,
              const std::string& name);

  void Define(const std::string region,
              const Teuchos::RCP<const MultiFunction>& f,
              int action, int submodel,
              const std::string& name);
  
  void Compute(double t0, double t1);
  void ComputeDistribute(double t0, double t1);
  void ComputeDistribute(double t0, double t1, double* weight);
  void ComputeDistributeMultiValue(double t0, double t1);
  void ComputeDistributeMultiValue(double t0, double t1, double* weight);
  
  void Finalize() {}
 
  // extract internal information
  int CollectActionsList();
 
  // iterator methods
  typedef std::map<int,double>::const_iterator Iterator;
  Iterator begin() const { return value_.begin(); }
  Iterator end() const  { return value_.end(); }
  Iterator find(const int j) const { return value_.find(j); }

  // access
  const std::string& tcc_name() { return tcc_name_; }
  int tcc_index() { return tcc_index_; }
  void set_tcc_index(int i) { tcc_index_ = i; }

 private:
  std::vector<int> actions_;
  std::vector<int> submodel_;
  std::string tcc_name_;
  int tcc_index_; // index the global list of components

 protected:
  std::map<int,double> value_;
  bool finalized_;
};


}  // namespace Transport
}  // namespace Amanzi

#endif  // AMANZI_TRANSPORT_DOMAIN_FUNCTION_HH_
