/*
  Transport PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lbl.gov)
*/

#ifndef AMANZI_TRANSPORT_BOUNDARY_FUNCTION_CHEMISTRY_HH_
#define AMANZI_TRANSPORT_BOUNDARY_FUNCTION_CHEMISTRY_HH_

#include <vector>
#include <map>
#include <string>

#include "Teuchos_RCP.hpp"

#include "Beaker.hh"
#include "Mesh.hh"
#include "TransportDomainFunction.hh"

#include "Amanzi_PK.hh"

namespace Amanzi {
namespace Transport {

class TransportBoundaryFunction_Chemistry : public TransportDomainFunction {
 public:
  TransportBoundaryFunction_Chemistry(){};
  TransportBoundaryFunction_Chemistry(const Teuchos::ParameterList& plist)
    : TransportDomainFunction(plist)
  {
    amanzi_pk_ = Teuchos::rcp_dynamic_cast<AmanziChemistry::Amanzi_PK>(
      plist.get<Teuchos::RCP<AmanziChemistry::Chemistry_PK>>("chemical pk"));
    chem_engine_ = amanzi_pk_->get_engine();

    constraints_ = plist.get<Teuchos::Array<std::string>>("names").toVector();
  }

  virtual void ComputeSubmodel(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                               Teuchos::RCP<CompositeVector> tcc)
  {
    AmanziMesh::Entity_ID_List cells;
    auto beaker_state = amanzi_pk_->beaker_state();
    const auto& beaker_parameters = amanzi_pk_->beaker_parameters();

    for (auto it = begin(); it != end(); ++it) {
      int f = it->first;
      mesh->face_get_cells(f, AmanziMesh::Parallel_type::OWNED, &cells);
      int c = cells[0];

      amanzi_pk_->CopyCellStateToBeakerState(c, tcc->ViewComponent("cell", true));

      auto& values = it->second;
      chem_engine_->EnforceConstraint(&beaker_state, beaker_parameters, constraints_, values);

      for (int i = 0; i < values.size(); i++) { values[i] = beaker_state.total.at(i); }
    }
  }

 private:
  Teuchos::RCP<AmanziChemistry::Amanzi_PK> amanzi_pk_;
  std::shared_ptr<AmanziChemistry::Beaker> chem_engine_;
  std::vector<std::string> constraints_;
};

} // namespace Transport
} // namespace Amanzi

#endif
