/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
    Ethan Coon (coonet@ornl.gov)
*/

//!

#pragma once

#include <string>

#include "Teuchos_RCP.hpp"
#include "Mesh.hh"

namespace Amanzi {
namespace Functions {

class BCsFunction {
 public:
  BCsFunction(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
      : mesh_(mesh) {}

  void Define(const std::vector<std::string>& regions,
              const Teuchos::RCP<const MultiFunction>& f);

  void Define(std::string& region, const Teuchos::RCP<const MultiFunction>& f);

  void Compute(double time);

  void Finalize();

  // iterator methods
  typedef std::map<int, double>::const_iterator Iterator;
  Iterator begin() const { return value_.begin(); }
  Iterator end() const { return value_.end(); }
  Iterator find(const int j) const { return value_.find(j); }
  std::map<int, double>::size_type size() { return value_.size(); }

 protected:
  std::map<int, double> value_;
  bool finalized_;
};

} // namespace Functions
} // namespace Amanzi


