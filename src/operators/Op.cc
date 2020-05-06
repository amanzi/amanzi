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
  if (A.size()) {
    Kokkos::parallel_for(
      "Op::Zero",
      A.size(),
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
  AMANZI_ASSERT(0);
  if (A.size()) {
    A.update_entries_host();
    for (int i = 0; i != A.size(); ++i){
      WhetStone::DenseMatrix<DefaultHostMemorySpace> lm(
        A.at_host(i),A.size_host(i,0),A.size_host(i,1)); 
      lm *= scaling; 
    }
    A.update_entries_device();
  }
  
  if (diag.get()) diag->scale(scaling);
}


// allocate work vector space
void Op::PreallocateWorkVectors() const
{
  int nfaces_owned = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  AMANZI_ASSERT(A.size() == nfaces_owned);
  v = CSR_Vector(A.size());
  Av = CSR_Vector(A.size());

  for (int i=0; i!=A.size(); ++i) {
    v.sizes_.view_host()(i,0) = A.size_host(i,0);
    v.row_map_.view_host()(i) = A.size_host(i,0);
    Av.sizes_.view_host()(i,0) = A.size_host(i,1);
    Av.row_map_.view_host()(i) = A.size_host(i,1);
  }
  
  v.prefix_sum(); 
  Av.prefix_sum();
}    


} // namespace Operators
} // namespace Amanzi



