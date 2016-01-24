/*
  Transport PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_TRANSPORT_BC_FACTORY_HH_
#define AMANZI_TRANSPORT_BC_FACTORY_HH_

#include <vector>

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "TransportBoundaryFunction.hh"
#include "TransportBoundaryFunction_Tracer.hh"

#ifdef ALQUIMIA_ENABLED
#include "TransportBoundaryFunction_Alquimia.hh"
#include "Alquimia_PK.hh"
#include "ChemistryEngine.hh"
#endif

namespace Amanzi {
namespace Transport {

class TransportBCFactory {
 public:
#ifdef ALQUIMIA_ENABLED
  // Alquimia-enabled constructors.
  TransportBCFactory(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                     const Teuchos::RCP<Teuchos::ParameterList>& list) :
      mesh_(mesh),
      list_(list),
      chem_pk_(Teuchos::null),
      chem_engine_(Teuchos::null) {};

  TransportBCFactory(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                     const Teuchos::RCP<Teuchos::ParameterList>& list,
                     const Teuchos::RCP<AmanziChemistry::Alquimia_PK>& chem_pk,
                     const Teuchos::RCP<AmanziChemistry::ChemistryEngine>& chem_engine) :
      mesh_(mesh),
      list_(list), 
      chem_pk_(chem_pk),
      chem_engine_(chem_engine) {};
#else
  TransportBCFactory(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                     const Teuchos::RCP<Teuchos::ParameterList>& list) :
      mesh_(mesh),
      list_(list) {};
#endif

  ~TransportBCFactory() {};
  
  void Create(std::vector<TransportBoundaryFunction*>& bcs) const;

  // non-reactive components
  void ProcessTracerList(std::vector<TransportBoundaryFunction*>& bcs) const;
  void ProcessTracerSpec(
      Teuchos::ParameterList& spec, TransportBoundaryFunction_Tracer* bc) const;

  // reactive components
  void ProcessGeochemicalConditionList(std::vector<TransportBoundaryFunction*>& bcs) const;

 private:
  const Teuchos::RCP<const AmanziMesh::Mesh>& mesh_;
  const Teuchos::RCP<Teuchos::ParameterList>& list_;
#ifdef ALQUIMIA_ENABLED
  const Teuchos::RCP<AmanziChemistry::Alquimia_PK> chem_pk_;
  const Teuchos::RCP<AmanziChemistry::ChemistryEngine> chem_engine_;
#endif
};

}  // namespace Transport
}  // namespace Amanzi

#endif
