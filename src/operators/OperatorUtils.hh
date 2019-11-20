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


// Do a Kokkos local mat-mult.
//
// This likely needs to be a better implementation...  FIXME --etc
//
// Calculates:
//   Y[Y_inds] += mat[elem] * X[X_inds]
template<class MatView_type, typename GO, class VecView_type>
void
LocalMatMultAdd(const MatView_type& mat,
                const GO elem,
                const typename VecView_type::const_type& X,
                const AmanziMesh::Entity_ID_View& X_ids,
                const VecView_type& Y,
                const AmanziMesh::Entity_ID_View& Y_ids)
{
  auto num_cols = X_ids.extent(0);
  auto num_rows = Y_ids.extent(0);
  Kokkos::parallel_for(
      "local_mat_mult_row_loop",
      num_rows,
      KOKKOS_LAMBDA(const int& i) {
        double result = 0.0;
        Kokkos::parallel_reduce("local_mat_mult_inner_reduction",
                X_ids.extent(0),
                [&](const int& j, double& lsum) {
                  lsum += mat(elem, j+i*num_cols) * X(X_ids(j), 0);
                }, result);
        Kokkos::atomic_add(Y(Y_ids(i), 0), result);
      });
}

} // namespace Operators
} // namespace Amanzi


#endif
