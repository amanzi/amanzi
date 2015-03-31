/*
  License: BSD
  Ethan Coon (ecoon@lanl.gov)
  Daniil Svyatskiy (dasvyat@lanl.gov)

  Matrix_TPFA_Surf provides for solving the cell-centered system with lambdas
  on the boundary system where a subset of
  the lambdas are more densely coupled for overland flow.

*/

#ifndef OPERATORS_MATRIX_TPFA_SURF_HH_
#define OPERATORS_MATRIX_TPFA_SURF_HH_

#include "Matrix_TPFA.hh"
#include "MatrixMFD_TPFA.hh"
#include "MatrixMFD.hh"

namespace Amanzi {
namespace Operators {

class Matrix_TPFA_Surf : public Matrix_TPFA {

 public:
  Matrix_TPFA_Surf(Teuchos::ParameterList& plist,
		   const Teuchos::RCP<const AmanziMesh::Mesh>& mesh);

  virtual void AssembleGlobalMatrices(){};
  virtual void ApplyBoundaryConditions(const std::vector<MatrixBC>& subsurface_markers,
          const std::vector<double>& subsurface_values, bool ADD_BC_FLUX=true);
  virtual void ComputeSchurComplement(const std::vector<MatrixBC>& bc_markers,
          const std::vector<double>& bc_values);

  virtual void SetSurfaceOperator(const Teuchos::RCP<MatrixMFD_TPFA>& surface_A) {
    surface_mesh_ = surface_A->Mesh();
    surface_A_ = surface_A; }
  virtual Teuchos::RCP<MatrixMFD_TPFA> GetSurfaceOperator() {
    return surface_A_; }

  virtual void SymbolicAssembleGlobalMatrices() {
    // This must be protected from being called too early.
    if (surface_mesh_ != Teuchos::null) {
      Matrix_TPFA::SymbolicAssembleGlobalMatrices();
    }
  }
  // virtual int Apply(const CompositeVector& X,
  //                    CompositeVector& Y); 
  // virtual int ApplyInverse(const CompositeVector& X,
  //                           CompositeVector& Y);

  virtual void AssembleRHS_() const;
  virtual void AssembleSchur_() const;


  virtual void ComputeNegativeResidual(const CompositeVector& solution,
				       const Teuchos::Ptr<CompositeVector>& residual) const;

  
 protected:
  virtual void FillMatrixGraphs_(const Teuchos::Ptr<Epetra_CrsGraph> cf_graph,
          const Teuchos::Ptr<Epetra_FECrsGraph> ff_graph);

 protected:
  Teuchos::RCP<const AmanziMesh::Mesh> surface_mesh_;
  Teuchos::RCP<MatrixMFD_TPFA> surface_A_;

  bool fill_graph;

  friend class MatrixMFD_Coupled_Surf;
};


} //namespace
} //namespace


#endif
