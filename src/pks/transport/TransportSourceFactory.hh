/*
  Transport PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_TRANSPORT_SOURCE_FACTORY_HH_
#define AMANZI_TRANSPORT_SOURCE_FACTORY_HH_

#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Point.hh"
#include "Mesh.hh"
#include "TransportSourceFunction.hh"

namespace Amanzi {
namespace Transport {

class TransportSourceFactory {
 public:
  TransportSourceFactory(const Teuchos::RCP<const AmanziMesh::Mesh> mesh,
                         const Teuchos::RCP<Teuchos::ParameterList> plist)
     : mesh_(mesh), plist_(plist) {};
  ~TransportSourceFactory() {};
  
  void Create(std::vector<TransportSourceFunction*>& srcs);

 private:
  void ProcessSourceSpec(Teuchos::ParameterList& list, 
                         const std::string& name, 
                         TransportSourceFunction* src) const;
  void ProcessStringActions(const std::string& name, int* method) const;
     
 private:
  const Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  const Teuchos::RCP<Teuchos::ParameterList> plist_;
};

}  // namespace Transport
}  // namespace Amanzi

#endif // TRANSPORT_SOURCE_FACTORY_HH_
