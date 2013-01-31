#ifndef __FLOW_SOURCE_FACTORY_HH__
#define __FLOW_SOURCE_FACTORY_HH__

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Point.hh"
#include "Mesh.hh"
#include "domain_function.hh"


namespace Amanzi {
namespace AmanziFlow {

class FlowSourceFactory {
 public:
  FlowSourceFactory(const Teuchos::RCP<const AmanziMesh::Mesh> mesh,
                    const Teuchos::RCP<Teuchos::ParameterList> params)
     : mesh_(mesh), params_(params) {};
  ~FlowSourceFactory() {};
  
  DomainFunction* createSource() const;

 private:
  void ProcessSourceSpec(Teuchos::ParameterList& list, DomainFunction* src) const;
  void ProcessStringActions(const std::string& name, int* method) const;
     
 private:
  const Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  const Teuchos::RCP<Teuchos::ParameterList> params_;
};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif // AMANZI_FLOW_SOURCE_FACTORY_HH_
