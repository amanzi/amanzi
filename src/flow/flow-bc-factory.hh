#ifndef AMANZI_FLOW_BC_FACTORY_HH_
#define AMANZI_FLOW_BC_FACTORY_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Mesh.hh"

namespace Amanzi {

class BoundaryFunction; // forward declaration

class FlowBCFactory {
 public:
  FlowBCFactory(const Teuchos::RCP<const AmanziMesh::Mesh> &mesh,
                const Teuchos::RCP<Teuchos::ParameterList> &params)
     : mesh_(mesh), params_(params) {}
  ~FlowBCFactory() {};
  
  BoundaryFunction* CreatePressure() const;
  BoundaryFunction* CreateMassFlux() const;
  BoundaryFunction* CreateStaticHead(double, double, double) const;

 private: // auxillary functions
  void process_pressure_list(Teuchos::ParameterList&, BoundaryFunction*) const;
  void process_pressure_spec(Teuchos::ParameterList&, BoundaryFunction*) const;
  void process_mass_flux_list(Teuchos::ParameterList&, BoundaryFunction*) const;
  void process_mass_flux_spec(Teuchos::ParameterList&, BoundaryFunction*) const;
  void process_static_head_list(double, double, double, Teuchos::ParameterList&,
      BoundaryFunction*) const;
  void process_static_head_spec(double, double, double,
      Teuchos::ParameterList&, BoundaryFunction*) const;
     
 private: // data
  const Teuchos::RCP<const AmanziMesh::Mesh> &mesh_;
  const Teuchos::RCP<Teuchos::ParameterList> &params_;
};

} // namespace Amanzi

#endif // AMANZI_FLOW_BC_FACTORY_HH_
