/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon (coonet@ornl.gov)
*/

//! Container for local matrices.

/*

  Op classes are a small container that generically stores local matrices and
  provides metadata about the structure of those local matrices via visitor
  pattern.

*/

#include "Op.hh"

#if 0 
namespace Amanzi {
namespace Operators {

void
Op::Zero()
{
  Kokkos::MDRangePolicy<Kokkos::Rank<2>> policy({0,0},
          {data.extent(0), data.extent(1)});
  Kokkos::parallel_for(policy,
                       KOKKOS_LAMBDA(const int i, const int j) {
                         data(i,j) = 0.0; });
}



// Restore pristine value of the matrices, i.e. before BCs.
int
Op::CopyShadowToMaster()
{
  if (shadow.extent(0) == data.extent(0)) {
    Kokkos::deep_copy(data, shadow);
    
  } else {
    // sparse copy
    Kokkos::MDRangePolicy<Kokkos::Rank<2>> policy({0,0},
            {shadow.extent(0), shadow.extent(1)});
    Kokkos::parallel_for(policy,
                         KOKKOS_LAMBDA(const int i, const int j) {
                           data(shadow_indices(i), j) = shadow(i,j); });
  }
  return 0;
}
    

// Matching rules for schemas.
bool
Op::Matches(int match_schema, int matching_rule)
{
  if (matching_rule == OPERATOR_SCHEMA_RULE_EXACT) {
    if ((match_schema & schema_old) == schema_old) return true;
  } else if (matching_rule == OPERATOR_SCHEMA_RULE_SUBSET) {
    if (match_schema & schema_old) return true;
  }
  return false;
}


// -- rescale local matrices in the container using a double
void
Op::Rescale(double scaling)
{
  Kokkos::MDRangePolicy<Kokkos::Rank<2>> policy({0,0},
          {data.extent(0), data.extent(1)});
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const int i,const int j) {
        data(i,j) *= scaling; });
}

} // namespace Operators
} // namespace Amanzi
#endif 

