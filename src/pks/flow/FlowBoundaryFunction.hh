/*
  Flow PK
 
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Neil Carlson (version 1)
           Ethan Coon (version 2)
           Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)

  Function applied to a mesh component with at most one function 
  application per entity.
*/

#ifndef AMANZI_FLOW_BOUNDARY_FUNCTION_HH_
#define AMANZI_FLOW_BOUNDARY_FUNCTION_HH_

#include <string>
#include <vector>

#include "Epetra_Vector.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "PK_DomainFunction.hh"

#include "FlowDefs.hh"


namespace Amanzi {
namespace Flow {

class FlowBoundaryFunction : public PK_DomainFunction {
 public:
  FlowBoundaryFunction()
    : rainfall_(false),
      relative_to_top_(false),
      relative_to_bottom_(false),
      no_flow_above_water_table_(false),
      seepage_flux_threshold_(0.0),
      bc_name_("underfined"),
      seepage_model_("")
  {}

  FlowBoundaryFunction(const Teuchos::ParameterList& plist);
  void ComputeSubmodel(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh);

  // modifiers and access
  void set_bc_name(const std::string& name) { bc_name_ = name; }
  std::string get_bc_name() { return bc_name_; }

  bool no_flow_above_water_table() const { return no_flow_above_water_table_; }
  std::string seepage_model() const { return seepage_model_; }
  double ref_pressure() const { return ref_pressure_; }
  double seepage_flux_threshold() const { return seepage_flux_threshold_; }

 private:
  void CalculateShiftWaterTable_(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                                 const std::string& region);

  void set_intersection_(const std::vector<AmanziMesh::Entity_ID>& v1,
                         const std::vector<AmanziMesh::Entity_ID>& v2,
                         std::vector<AmanziMesh::Entity_ID>* vv);

 private:
  bool rainfall_;
  bool relative_to_top_, relative_to_bottom_;
  bool no_flow_above_water_table_;

  double rho_, g_;
  double ref_pressure_, seepage_flux_threshold_;
  std::string bc_name_, seepage_model_;

  std::vector<std::string> regions_;

  Teuchos::RCP<Epetra_Vector> shift_water_table_;
  int nedges_;
};

} // namespace Flow
} // namespace Amanzi

#endif
