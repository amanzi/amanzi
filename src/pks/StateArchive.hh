/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Process Kernels

  Keeps copies of fields which could be restored. e.g. when a time
  integration step fails.
*/

#ifndef AMANZI_PK_STATE_ARCHIVE_HH_
#define AMANZI_PK_STATE_ARCHIVE_HH_

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
  StateArchive(Teuchos::RCP<State>& S, Teuchos::RCP<VerboseObject>& vo)
    : S_(S), vo_(vo) {};

  void Add(const std::vector<std::string>& fields, const Tag& tag);

  void Restore(const std::string& passwd);

  void CopyFieldsToPrevFields(std::vector<std::string>& fields,
                              const std::string& passwd,
                              bool add);

  // access
  const CompositeVector& get(const std::string& name);

 private:
  Teuchos::RCP<State> S_;
  Teuchos::RCP<VerboseObject> vo_;

  Tag tag_;
  std::map<std::string, CompositeVector> fields_;
};

} // namespace Amanzi

#endif
