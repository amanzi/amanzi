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

// TPLs
#include "Teuchos_RCP.hpp"

// amanzi
#ifdef ALQUIMIA_ENABLED
#include "Alquimia_PK.hh"
#include "ChemistryEngine.hh"
#endif
#include "State.hh"
#include "Mesh.hh"

// Transport
#include "TransportBoundaryFunction.hh"
#include "TransportBoundaryFunction_Tracer.hh"
#include "TransportBoundaryFunction_Coupler.hh"
#ifdef ALQUIMIA_ENABLED
#include "TransportBoundaryFunction_Alquimia.hh"
#endif

namespace Amanzi {
namespace Transport {

class TransportBCFactory {
 public:
#ifdef ALQUIMIA_ENABLED
  // Alquimia-enabled constructors.
  TransportBCFactory(const Teuchos::RCP<const State>& S,
                     const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                     const Teuchos::RCP<Teuchos::ParameterList>& list,
                     const Teuchos::RCP<AmanziChemistry::Alquimia_PK>& chem_pk,
                     const Teuchos::RCP<AmanziChemistry::ChemistryEngine>& chem_engine) :
    S_(S),
      mesh_(mesh),
      list_(list), 
      chem_pk_(chem_pk),
      chem_engine_(chem_engine) {};
#else
  TransportBCFactory(const Teuchos::RCP<const State>& S,
                     const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                     const Teuchos::RCP<Teuchos::ParameterList>& list) :
    S_(S),
      mesh_(mesh),
      list_(list) {};
#endif

  ~TransportBCFactory() {};
  
  void Create(std::vector<TransportBoundaryFunction*>& bcs) const;

 private:
  // non-reactive components
  void ProcessTracerList_(std::vector<TransportBoundaryFunction*>& bcs) const;
  void ProcessTracerSpec_(
      Teuchos::ParameterList& spec, TransportBoundaryFunction_Tracer* bc) const;

  // reactive components
  void ProcessGeochemicalConditionList_(std::vector<TransportBoundaryFunction*>& bcs) const;

  // domain coupling
  void ProcessCouplingConditionList_(std::vector<TransportBoundaryFunction*>& bcs) const;
  void ProcessCouplerSpec_(Teuchos::ParameterList& spec, TransportBoundaryFunction_Coupler* bc) const;

 private:
  const Teuchos::RCP<const State>& S_;
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
