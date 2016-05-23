/*
  Transport PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Jeffrey Johnson (jnjohnson@lbl.gov)
*/

#ifndef AMANZI_TRANSPORT_BOUNDARY_FUNCTION_ALQUIMIA_HH_
#define AMANZI_TRANSPORT_BOUNDARY_FUNCTION_ALQUIMIA_HH_

#include <vector>
#include <map>
#include <string>

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "MultiFunction.hh"
#include "TransportBoundaryFunction.hh"

#ifdef ALQUIMIA_ENABLED
#include "Alquimia_PK.hh"
#include "ChemistryEngine.hh"

namespace Amanzi {
namespace Transport {

class TransportBoundaryFunction_Alquimia : public TransportBoundaryFunction {
 public:
  TransportBoundaryFunction_Alquimia(const std::vector<double>& times,
                                     const std::vector<std::string>& cond_names, 
                                     const Teuchos::RCP<const State>& S,
                                     const Teuchos::RCP<const AmanziMesh::Mesh> &mesh,
                                     Teuchos::RCP<AmanziChemistry::Alquimia_PK> chem_pk,
                                     Teuchos::RCP<AmanziChemistry::ChemistryEngine> chem_engine);
  ~TransportBoundaryFunction_Alquimia();
  
  void Compute(double time);

  void Define(const std::vector<std::string> &regions);
  void Define(std::string region);

 private:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  const Teuchos::RCP<const State>& S_;

  // The geochemical conditions we are enforcing, and the times we are enforcing them at.
  std::vector<double> times_;
  std::vector<std::string> cond_names_;

  // Chemistry state and engine.
  Teuchos::RCP<AmanziChemistry::Alquimia_PK> chem_pk_;
  Teuchos::RCP<AmanziChemistry::ChemistryEngine> chem_engine_;

  // Containers for interacting with the chemistry engine.
  AlquimiaState alq_state_;
  AlquimiaMaterialProperties alq_mat_props_;
  AlquimiaAuxiliaryData alq_aux_data_;
  AlquimiaAuxiliaryOutputData alq_aux_output_;

  // A mapping of boundary face indices to interior cells.
  std::map<int, int> cell_for_face_;
};

}  // namespace Transport
}  // namespace Amanzi

#endif


#endif
