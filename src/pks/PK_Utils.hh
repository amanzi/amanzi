/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Process Kernels

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

/*
  Keeps copies of fields which could be restored. e.g. when a time integration step
  fails. The list is populated by

    - input fields which are primary evaluators
    - input fields which have no evaluators
    - fields with previous prev_ when they are overwritten
*/

class StateArchive {
 public:
  StateArchive() = delete;
  StateArchive(Teuchos::RCP<State>& S, Teuchos::RCP<VerboseObject>& vo) : S_(S), vo_(vo){};

  void Add(std::vector<std::string>& fields, const Tag& tag);

  void Restore(const std::string& passwd);

  void CopyFieldsToPrevFields(std::vector<std::string>& fields, const std::string& passwd);

  // access
  const CompositeVector& get(const std::string& name);

 private:
  Teuchos::RCP<State> S_;
  Teuchos::RCP<VerboseObject> vo_;

  Tag tag_;
  std::map<std::string, CompositeVector> fields_;
};


// Miscalleneous functions.
// Average permeability tensor in horizontal direction.
void
PKUtils_CalculatePermeabilityFactorInWell(const Teuchos::Ptr<State>& S,
                                          Teuchos::RCP<Epetra_MultiVector>& Kxy);

AmanziGeometry::Point
PKUtils_EntityCoordinates(int id, AmanziMesh::Entity_ID kind, const AmanziMesh::Mesh& mesh);

} // namespace Amanzi

#endif
