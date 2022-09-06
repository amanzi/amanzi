/*
  Process Kernels

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Miscalleneous collection of simple non-member functions.
*/

#ifndef AMANZI_PK_UTILS_HH_
#define AMANZI_PK_UTILS_HH_

#include <map>
#include <string>
#include <vector>

#include "CompositeVector.hh"
#include "State.hh"
#include "Tag.hh"
#include "VerboseObject.hh"

namespace Amanzi {

class StateArchive {
 public:
  StateArchive() = delete;
  StateArchive(Teuchos::RCP<State>& S, Teuchos::RCP<VerboseObject>& vo) : S_(S), vo_(vo) {};

  void Add(std::vector<std::string> fields,
           std::vector<std::string> evals,
           std::vector<std::string> primary,
           const Tag& tag = Tags::DEFAULT);

  void Restore(const std::string& passwd);

  // access
  const CompositeVector get(const std::string& name);

 private:
  Teuchos::RCP<State> S_;
  Teuchos::RCP<VerboseObject> vo_;

  Tag tag_;
  std::vector<std::string> primary_;
  std::map<std::string, CompositeVector> fields_, evals_;
};


// Miscalleneous functions.
// Average permeability tensor in horizontal direction.
void PKUtils_CalculatePermeabilityFactorInWell(
    const Teuchos::Ptr<State>& S, Teuchos::RCP<Epetra_Vector>& Kxy);

AmanziGeometry::Point PKUtils_EntityCoordinates(
    int id, AmanziMesh::Entity_ID kind, const AmanziMesh::Mesh& mesh);

}  // namespace Amanzi

#endif
