/*
  License: BSD
  Ethan Coon (ecoon@lanl.gov)

  MatrixMFD_Surf provides for solving the p-lambda system where a subset of
  the lambdas are more densely coupled for overland flow.

*/

#ifndef OPERATORS_MATRIX_MFD_SURF_HH_
#define OPERATORS_MATRIX_MFD_SURF_HH_

#include "matrix_mfd_tpfa.hh"
#include "matrix_mfd.hh"

namespace Amanzi {
namespace Operators {


class MatrixMFD_Surf : public MatrixMFD {

 public:
  MatrixMFD_Surf(Teuchos::ParameterList& plist,
                 MFD_method method,
                 const Teuchos::RCP<const AmanziMesh::Mesh> mesh,
                 const Teuchos::RCP<const AmanziMesh::Mesh> surface_mesh);

  virtual void SymbolicAssembleGlobalMatrices();

  void AssembleGlobalMatricesWithSurface(const MatrixMFD_TPFA& surface_A);
  virtual void ComputeSchurComplement(const std::vector<Matrix_bc>& bc_markers,
          const std::vector<double>& bc_values);

 protected:
  Teuchos::RCP<const AmaniMesh::Mesh> surface_mesh_;
};


} //namespace
} //namespace
