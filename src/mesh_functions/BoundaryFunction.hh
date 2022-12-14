/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: (v1) Neil Carlson
      (v2) Ethan Coon
*/

/*
  Mesh Functions

  Function applied to a mesh component with at most one function
  application per entity.

  Amanzi is no longer using this function!
*/

#ifndef AMANZI_BOUNDARY_FUNCTION_HH_
#define AMANZI_BOUNDARY_FUNCTION_HH_

#include <vector>
#include <map>
#include <string>

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "MultiFunction.hh"
#include "UniqueMeshFunction.hh"

namespace Amanzi {
namespace Functions {

class BoundaryFunction : public UniqueMeshFunction {
 public:
  BoundaryFunction(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
    : UniqueMeshFunction(mesh), finalized_(false){};

  void Define(const std::vector<std::string>& regions, const Teuchos::RCP<const MultiFunction>& f);

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


#endif
