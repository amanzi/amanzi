/*
  License: BSD
  Ethan Coon (ecoon@lanl.gov)

  MatrixMFD_Surf provides for solving the p-lambda system where a subset of
  the lambdas are more densely coupled for overland flow.

*/

#ifndef OPERATORS_MATRIX_MFD_SURF_HH_
#define OPERATORS_MATRIX_MFD_SURF_HH_

#include "MatrixMFD_TPFA.hh"
#include "MatrixMFD.hh"

namespace Amanzi {
namespace Operators {

class MatrixMFD_Surf : virtual public MatrixMFD {

 public:
  MatrixMFD_Surf(Teuchos::ParameterList& plist,
                 const Teuchos::RCP<const AmanziMesh::Mesh>& mesh);

  virtual void AssembleGlobalMatrices();

  virtual void ApplyBoundaryConditions(const std::vector<MatrixBC>& subsurface_markers,
          const std::vector<double>& subsurface_values);

  virtual void ComputeSchurComplement(const std::vector<MatrixBC>& bc_markers,
          const std::vector<double>& bc_values);

  virtual void SetSurfaceOperator(const Teuchos::RCP<MatrixMFD_TPFA>& surface_A) {
    surface_mesh_ = surface_A->Mesh();
    surface_A_ = surface_A; }
  virtual void GetSurfaceOperator(Teuchos::RCP<MatrixMFD_TPFA>& surface_A) {
    surface_A = surface_A_; }

 protected:
  virtual void FillMatrixGraphs_(const Teuchos::Ptr<Epetra_CrsGraph> cf_graph,
          const Teuchos::Ptr<Epetra_FECrsGraph> ff_graph);

 protected:
  Teuchos::RCP<const AmanziMesh::Mesh> surface_mesh_;
  Teuchos::RCP<MatrixMFD_TPFA> surface_A_;

};


} //namespace
} //namespace


#endif
