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

namespace Amanzi {
namespace Operators {

void
Op::Zero()
{
  if (matrices.size()) {
    Kokkos::parallel_for(matrices.size(),
                         KOKKOS_LAMBDA(const int i) {
                           Zero(i);
                         });
  }
  if (diag.get()) diag->putScalar(0.);
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
void Op::Rescale(double scaling)
{
  if (matrices.size())
    for (int i = 0; i != matrices.size(); ++i) matrices[i] *= scaling;
  if (diag.get()) diag->scale(scaling);
}

} // namespace Operators
} // namespace Amanzi



