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

  virtual void ApplyBoundaryConditions(const std::vector<MatrixBC>& subsurface_markers,
				       const std::vector<double>& subsurface_values, 
				       bool ADD_BC_FLUX=true);

  virtual void SetSurfaceOperator(const Teuchos::RCP<MatrixMFD_TPFA>& surface_A);
  virtual Teuchos::RCP<MatrixMFD_TPFA> GetSurfaceOperator() {
    return surface_A_; }

  virtual void SymbolicAssembleGlobalMatrices() {
    // This must be protected from being called too early.
    if (surface_mesh_ != Teuchos::null) {
      MatrixMFD::SymbolicAssembleGlobalMatrices();
    }
  }

  virtual int Apply(const CompositeVector& X,
                     CompositeVector& Y) const;

 protected:
  virtual void FillMatrixGraphs_(const Teuchos::Ptr<Epetra_CrsGraph> cf_graph,
          const Teuchos::Ptr<Epetra_FECrsGraph> ff_graph);

  virtual void AssembleAff_() const;
  virtual void AssembleSchur_() const;
  virtual void AssembleRHS_() const;

 protected:
  Teuchos::RCP<const AmanziMesh::Mesh> surface_mesh_;
  Teuchos::RCP<MatrixMFD_TPFA> surface_A_;
  bool dump_schur_;

  // TRILINOS FAIL
  //  Teuchos::RCP<const Epetra_Import> surf_importer_;
  Teuchos::RCP<const Epetra_Map> surf_map_in_subsurf_;
  
  friend class MatrixMFD_Coupled_Surf;
};


} //namespace
} //namespace


#endif
