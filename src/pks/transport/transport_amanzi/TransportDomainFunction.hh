/*
  Transport PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_TRANSPORT_DOMAIN_FUNCTION_HH_
#define AMANZI_TRANSPORT_DOMAIN_FUNCTION_HH_

#include <vector>
#include <string>

#include "Teuchos_RCP.hpp"

#include "CommonDefs.hh"
#include "Mesh.hh"
#include "PK_DomainFunction.hh"

namespace Amanzi {
namespace Transport {

class TransportDomainFunction : public PK_DomainFunction {
 public:
  TransportDomainFunction(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) : 
      PK_DomainFunction(mesh) {};
  
  void Define(const std::vector<std::string>& regions,
              const Teuchos::RCP<const MultiFunction>& f,
              int action, int submodel,
              const std::string& name);

  void Define(const std::string& region,
              const Teuchos::RCP<const MultiFunction>& f,
              int action, int submodel,
              const std::string& name);
  
  // access
  const std::string& tcc_name() { return tcc_name_; }
  int tcc_index() { return tcc_index_; }
  void set_tcc_index(int i) { tcc_index_ = i; }

 private:
  std::string tcc_name_;
  int tcc_index_; // index the global list of components
};


}  // namespace Transport
}  // namespace Amanzi

#endif  // AMANZI_TRANSPORT_DOMAIN_FUNCTION_HH_
