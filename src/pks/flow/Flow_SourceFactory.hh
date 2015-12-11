/*
  Flow PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_FLOW_SOURCE_FACTORY_HH_
#define AMANZI_FLOW_SOURCE_FACTORY_HH_

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "Point.hh"

#include "FlowDomainFunction.hh"

namespace Amanzi {
namespace Flow {

class FlowSourceFactory {
 public:
  FlowSourceFactory(const Teuchos::RCP<const AmanziMesh::Mesh> mesh,
                    const Teuchos::RCP<Teuchos::ParameterList> plist)
     : mesh_(mesh), plist_(plist) {};
  ~FlowSourceFactory() {};
  
  void Create(std::vector<FlowDomainFunction*>& srcs) const;

 private:
  void ProcessSourceSpec(Teuchos::ParameterList& list, FlowDomainFunction* src) const;
  void ProcessStringActions(const std::string& name, int* method) const;
     
 private:
  const Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  const Teuchos::RCP<Teuchos::ParameterList> plist_;
};

}  // namespace Flow
}  // namespace Amanzi

#endif

