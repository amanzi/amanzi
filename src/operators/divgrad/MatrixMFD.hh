/*
  License: BSD
  Authors: Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
           Ethan Coon (ecoon@lanl.gov) (ATS version)
  MatrixMFD provides a mimetic discretization for the elliptic operator div K grad u.

*/

#ifndef OPERATORS_MATRIX_MFD_HH_
#define OPERATORS_MATRIX_MFD_HH_

#include <strings.h>

#include "Teuchos_RCP.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"


#include "Epetra_Map.h"
#include "Epetra_Operator.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_SerialDenseVector.h"

#include "Mesh.hh"
#include "Point.hh"
#include "CompositeVectorSpace.hh"
#include "CompositeVector.hh"
#include "CompositeMatrix.hh"
#include "EpetraMatrix.hh"
#include "DenseMatrix.hh"
#include "mfd3d_diffusion.hh"
#include "Preconditioner.hh"
#include "PreconditionerFactory.hh"
#include "EpetraMatrixDefault.hh"
#include "LinearOperator.hh"
#include "LinearOperatorFactory.hh"

#include "MatrixMFD_Defs.hh"

namespace Amanzi {
namespace Operators {

class MatrixMFD : public CompositeMatrix {
 public:

  MatrixMFD(Teuchos::ParameterList& plist,
            const Teuchos::RCP<const AmanziMesh::Mesh>& mesh);

  MatrixMFD(const MatrixMFD& other);

  MatrixMFD& operator=(const MatrixMFD& other);
  
  // Virtual destructor
  virtual ~MatrixMFD() {};

  // CompositeMatrix methods
  virtual const CompositeVectorSpace& DomainMap() const {
    return *space_; }

  virtual const CompositeVectorSpace& RangeMap() const {
    return *space_; }

  virtual Teuchos::RCP<CompositeMatrix> Clone() const {
    return Teuchos::rcp(new MatrixMFD(*this)); }

  virtual int Apply(const CompositeVector& X,
                     CompositeVector& Y) const;

  virtual int ApplyInverse(const CompositeVector& X,
                            CompositeVector& Y) const;


  // Access to local matrices for external tweaking.
  std::vector<double>& Acc_cells() {
    MarkLocalMatricesAsChanged_();
    return Acc_cells_;
  }
  std::vector<Teuchos::SerialDenseMatrix<int, double> >& Aff_cells() {
    MarkLocalMatricesAsChanged_();
    return Aff_cells_;
  }
  std::vector<Epetra_SerialDenseVector>& Acf_cells() {
    MarkLocalMatricesAsChanged_();
    return Acf_cells_;
  }
  std::vector<Epetra_SerialDenseVector>& Afc_cells() {
    MarkLocalMatricesAsChanged_();
    return Afc_cells_;
  }
  const std::vector<double>& Acc_cells() const { return Acc_cells_; }
  const std::vector<Teuchos::SerialDenseMatrix<int, double> >& Aff_cells() const { return Aff_cells_; }
  const std::vector<Epetra_SerialDenseVector>& Acf_cells() const { return Acf_cells_; }
  const std::vector<Epetra_SerialDenseVector>& Afc_cells() const { return Afc_cells_; }

  // Access to local rhs
  std::vector<double>& Fc_cells() {
    assembled_rhs_ = false;
    return Fc_cells_;
  }
  std::vector<Epetra_SerialDenseVector>& Ff_cells() {
    assembled_rhs_ = false;
    return Ff_cells_;
  }

  // Const access to assembled matrices.
  Teuchos::RCP<const Epetra_Vector> Acc() {
    return Acc_;
  }
  Teuchos::RCP<const Epetra_FECrsMatrix> Aff() {
    if (!assembled_operator_) AssembleAff_();
    return Aff_;
  }
  Teuchos::RCP<const Epetra_FECrsMatrix> Schur() {
    if (!assembled_schur_) AssembleSchur_();
    return Sff_;
  }
  Teuchos::RCP<const CompositeVector> rhs() {
    if (!assembled_rhs_) AssembleRHS_();
    return rhs_;
  }

  // Other accessors/mutators.
  MFDMethod method(){return method_;}
  bool symmetric() { return flag_symmetry_; }
  void set_symmetric(bool flag_symmetry) { flag_symmetry_ = flag_symmetry; }
  const Epetra_Comm& Comm() const { return *(mesh_->get_comm()); }
  Teuchos::RCP<const AmanziMesh::Mesh> Mesh() const { return mesh_; }

  // Main computational methods
  // -- local matrices
  virtual void CreateMFDmassMatrices(
      const Teuchos::Ptr<std::vector<WhetStone::Tensor> >& K);
  virtual void CreateMFDstiffnessMatrices(
      const Teuchos::Ptr<const CompositeVector>& Krel);
  virtual void CreateMFDrhsVectors();
  virtual void Add2MFDstiffnessMatrices(std::vector<double>* Acc_ptr,
					std::vector<Teuchos::SerialDenseMatrix<int, double> >* Aff_ptr,
					std::vector<Epetra_SerialDenseVector>* Acf_ptr,
					std::vector<Epetra_SerialDenseVector>* Afc_ptr);

  virtual void ApplyBoundaryConditions(const std::vector<MatrixBC>& bc_markers,
				       const std::vector<double>& bc_values,
				       bool ADD_BC_FLUX=true);

  // -- global matrices
  virtual void SymbolicAssembleGlobalMatrices();

  // Operator methods.
  virtual void ComputeResidual(const CompositeVector& X,
                       const Teuchos::Ptr<CompositeVector>& F) const;
  virtual void ComputeNegativeResidual(const CompositeVector& X,
          const Teuchos::Ptr<CompositeVector>& F) const;

  // Solver methods.
  virtual void InitPreconditioner();

  // First derivative quantities.
  virtual void DeriveFlux(const CompositeVector& solution,
                          const Teuchos::Ptr<CompositeVector>& flux) const;
  // void DeriveFluxFV(const CompositeVector& solution,
  // 		    const Teuchos::Ptr<CompositeVector>& flux,
  // 		    const Teuchos::Ptr<const CompositeVector>& rel_perm){};
  virtual void DeriveCellVelocity(const CompositeVector& flux,
          const Teuchos::Ptr<CompositeVector>& velocity) const;

  // Consistency methods
  virtual void UpdateConsistentFaceConstraints(const Teuchos::Ptr<CompositeVector>& u);
  virtual void UpdateConsistentFaceCorrection(const CompositeVector& u,
          const Teuchos::Ptr<CompositeVector>& Pu);
  virtual void UpdateConsistentCellCorrection(const CompositeVector& u,
          const Teuchos::Ptr<CompositeVector>& Pu);

  // note X == Y is valid
  int ApplyAfc(const CompositeVector& X, CompositeVector& Y, double scalar) const;
  int ApplyAcf(const CompositeVector& X, CompositeVector& Y, double scalar) const;

 protected:
  // Assertions of assembly process
  void AssertAssembledOperator_or_die_() const;
  void AssertAssembledSchur_or_die_() const;
  void AssertAssembledRHS_or_die_() const;
  virtual void MarkLocalMatricesAsChanged_() {
    assembled_operator_ = false;
    assembled_schur_ = false;
    assembled_rhs_ = false;
  }

  int ApplyAfc_(const Epetra_MultiVector& X, CompositeVector& Y, double scalar) const;
  int ApplyAcf_(const CompositeVector& X, Epetra_MultiVector& Y, double scalar) const;

  void InitializeFromPList_();
  virtual void UpdatePreconditioner_() const;

  virtual void FillMatrixGraphs_(const Teuchos::Ptr<Epetra_CrsGraph> cf_graph,
          const Teuchos::Ptr<Epetra_FECrsGraph> ff_graph);
  virtual void CreateMatrices_(const Epetra_CrsGraph& cf_graph,
          const Epetra_FECrsGraph& ff_graph);

  virtual void AssembleAff_() const;
  virtual void AssembleRHS_() const;
  virtual void AssembleSchur_() const;
  
 private:
  // These are dangerous -- they require ghosted vectors, and a
  // communication either before or after application.
  int ApplyAfc_(const Epetra_MultiVector& X, Epetra_MultiVector& Y, double scalar) const;
  int ApplyAcf_(const Epetra_MultiVector& X, Epetra_MultiVector& Y, double scalar) const;



 protected:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::ParameterList plist_;
  bool flag_symmetry_;

  // lazy assembly
  mutable bool assembled_operator_;
  mutable bool assembled_schur_;
  mutable bool assembled_rhs_;

  MFDMethod method_;

  // local matrices
  std::vector<WhetStone::DenseMatrix > Mff_cells_;
  std::vector<Teuchos::SerialDenseMatrix<int, double> > Aff_cells_;
  std::vector<Epetra_SerialDenseVector> Acf_cells_, Afc_cells_;
  std::vector<double> Acc_cells_;  // duplication may be useful later
  std::vector<Epetra_SerialDenseVector> Ff_cells_;
  std::vector<double> Fc_cells_;

  // boundary condition flags
  std::vector<MatrixBC> bc_markers_;

  // global matrices
  Teuchos::RCP<Epetra_Vector> Acc_;
  mutable Teuchos::RCP<Epetra_FECrsMatrix> Aff_;
  mutable Teuchos::RCP<Epetra_FECrsMatrix> Sff_;  // Schur complement

  // global rhs
  mutable Teuchos::RCP<CompositeVector> rhs_;

  // diagnostics
  int nokay_;
  int npassed_; // performance of algorithms generating mass matrices

  // VectorSpace describing both domain and range
  Teuchos::RCP<CompositeVectorSpace> space_;

  // preconditioner for Schur complement
  mutable Teuchos::RCP<AmanziPreconditioners::Preconditioner> S_pc_;
  mutable Teuchos::RCP<AmanziPreconditioners::Preconditioner> Aff_pc_;

  // LinearOperator and Preconditioner for solving face system
  // Aff * x_f = r_Aff c - Afc * x_c for x_f
  Teuchos::RCP<EpetraMatrixDefault<Epetra_FECrsMatrix> > Aff_op_;
  Teuchos::RCP<EpetraMatrix> Aff_solver_;

  // verbose object
  Teuchos::RCP<VerboseObject> vo_;

  friend class MatrixMFD_Coupled;
  friend class MatrixMFD_Coupled_Surf;

};


}  // namespace AmanziFlow
}  // namespace Amanzi

#endif
