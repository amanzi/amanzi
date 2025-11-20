/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Transport PK

*/

#ifndef AMANZI_TRANSPORT_SOURCE_FUNCTION_ALQUIMIA_HH_
#define AMANZI_TRANSPORT_SOURCE_FUNCTION_ALQUIMIA_HH_

#include <vector>
#include <map>
#include <string>

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "FunctionTabularString.hh"
#include "TransportDomainFunction.hh"

#ifdef ALQUIMIA_ENABLED
#include "Alquimia_PK.hh"
#include "ChemistryEngine.hh"

namespace Amanzi {
namespace Transport {

class TransportSourceFunction_Alquimia : public TransportDomainFunction {
 public:
  TransportSourceFunction_Alquimia(const Teuchos::ParameterList& plist,
                                   const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                                   Teuchos::RCP<AmanziChemistry::Alquimia_PK> alquimia_pk,
                                   Teuchos::RCP<AmanziChemistry::ChemistryEngine> chem_engine);
  ~TransportSourceFunction_Alquimia();

  virtual void Compute(double t_old, double t_new) override;

  DomainFunction_kind getType() const override { return DomainFunction_kind::ALQUIMIA; }

 protected:
  void Init_(const std::vector<std::string>& regions);

 protected:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;

  // string function of geochemical conditions
  Teuchos::RCP<FunctionTabularString> f_;

  // Chemistry state and engine.
  Teuchos::RCP<AmanziChemistry::Alquimia_PK> alquimia_pk_;
  Teuchos::RCP<AmanziChemistry::ChemistryEngine> chem_engine_;

  // Containers for interacting with the chemistry engine.
  AmanziChemistry::AlquimiaBeaker beaker_;
};

} // namespace Transport
} // namespace Amanzi

#endif


#endif
