/*
  Transport PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Daniil Svyatsky (dasvyat@lanl.gov)
*/

#ifndef AMANZI_TRANSPORT_SOURCE_FUNCTION_CONCENTRATION_HH_
#define AMANZI_TRANSPORT_SOURCE_FUNCTION_CONCENTRATION_HH_

#include <vector>
#include <map>
#include <string>

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "FunctionTabularString.hh"
#include "TransportDomainFunction.hh"

namespace Amanzi {
namespace Transport {

class TransportSourceFunction_Concentration : public TransportDomainFunction {
 public:
  TransportSourceFunction_Concentration(const Teuchos::ParameterList& plist,
                                        const Teuchos::RCP<const AmanziMesh::Mesh>& mesh);
  ~TransportSourceFunction_Concentration();
  
  void Compute(double t_old, double t_new);

  // require by the case class
  virtual std::string name() const { return "concentration source"; } 

 private:
  void Init_(const std::vector<std::string> &regions);

 private:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  std::vector<std::string> components_;

};

}  // namespace Transport
}  // namespace Amanzi

#endif


#endif
