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
                 const Teuchos::RCP<const AmanziMesh::Mesh> mesh,
                 const Teuchos::RCP<const AmanziMesh::Mesh> surface_mesh);

  // NOTE this is not a copy constructor!
  MatrixMFD_Surf(const MatrixMFD& other,
                 const Teuchos::RCP<const AmanziMesh::Mesh> surface_mesh);

  virtual void SymbolicAssembleGlobalMatrices();

  virtual void AssembleGlobalMatrices();

  virtual void ApplyBoundaryConditions(const std::vector<Matrix_bc>& subsurface_markers,
          const std::vector<double>& subsurface_values);

  virtual void ComputeSchurComplement(const std::vector<Matrix_bc>& bc_markers,
          const std::vector<double>& bc_values);

  void set_surface_A(const Teuchos::RCP<MatrixMFD_TPFA>& surface_A) {
    surface_A_ = surface_A; }

 protected:
  Teuchos::RCP<const AmanziMesh::Mesh> surface_mesh_;
  Teuchos::RCP<MatrixMFD_TPFA> surface_A_;
};


} //namespace
} //namespace


#endif
