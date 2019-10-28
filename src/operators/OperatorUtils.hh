/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon (coonet@ornl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#ifndef AMANZI_OPERATORS_UTILS_HH_
#define AMANZI_OPERATORS_UTILS_HH_

#include "Teuchos_RCP.hpp"
#include "Schema.hh"

namespace Amanzi {
namespace Operators {

class Schema;

// Estimate the max number of unknowns per row. Note this can be an
// overestimate, but shouldn't be an underestimate.
unsigned int
MaxRowSize(const AmanziMesh::Mesh& mesh, int schema, unsigned int n_dofs = 1);
unsigned int
MaxRowSize(const AmanziMesh::Mesh& mesh, Schema& schema);


} // namespace Operators
} // namespace Amanzi


#endif
