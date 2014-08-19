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

  // Constructor
  MatrixMFD(Teuchos::ParameterList& plist,
            const Teuchos::RCP<const AmanziMesh::Mesh>& mesh);

  MatrixMFD(const MatrixMFD& other);


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
    return Acc_cells_;
  }
  std::vector<Teuchos::SerialDenseMatrix<int, double> >& Aff_cells() {
    return Aff_cells_;
  }
  std::vector<Epetra_SerialDenseVector>& Acf_cells() {
    return Acf_cells_;
  }
  std::vector<Epetra_SerialDenseVector>& Afc_cells() {
    return Afc_cells_;
  }

  // Access to local rhs
  std::vector<double>& Fc_cells() {
    return Fc_cells_;
  }
  std::vector<Epetra_SerialDenseVector>& Ff_cells() {
    return Ff_cells_;
  }

  // Const access to assembled matrices.
  Teuchos::RCP<const Epetra_Vector> Acc() {
    AssertAssembledOperator_or_die_();
    return Acc_;
  }
  Teuchos::RCP<const Epetra_FECrsMatrix> Aff() {
    AssertAssembledOperator_or_die_();
    return Aff_;
  }
  Teuchos::RCP<const Epetra_CrsMatrix> Afc() {
    AssertAssembledOperator_or_die_();
    return Afc_;
  }
  Teuchos::RCP<const Epetra_CrsMatrix> Acf() {
    AssertAssembledOperator_or_die_();
    return Acf_;
  }
  Teuchos::RCP<const Epetra_FECrsMatrix> Schur() {
    AssertAssembledSchur_or_die_();
    return Sff_;
  }
  Teuchos::RCP<const CompositeVector> rhs() {
    AssertAssembledOperator_or_die_();
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
  virtual void AssembleGlobalMatrices();
  virtual void ComputeSchurComplement(const std::vector<MatrixBC>& bc_markers,
          const std::vector<double>& bc_values);

  // Operator methods.
  virtual void ComputeResidual(const CompositeVector& X,
                       const Teuchos::Ptr<CompositeVector>& F) const;
  virtual void ComputeNegativeResidual(const CompositeVector& X,
          const Teuchos::Ptr<CompositeVector>& F) const;

  // Solver methods.
  virtual void InitPreconditioner();
  virtual void UpdatePreconditioner();

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


 protected:
  void InitializeFromPList_();

  // Assertions of assembly process
  void AssertAssembledOperator_or_die_() const;
  void AssertAssembledSchur_or_die_() const;

  virtual void FillMatrixGraphs_(const Teuchos::Ptr<Epetra_CrsGraph> cf_graph,
          const Teuchos::Ptr<Epetra_FECrsGraph> ff_graph);
  virtual void CreateMatrices_(const Epetra_CrsGraph& cf_graph,
          const Epetra_FECrsGraph& ff_graph);


 protected:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::ParameterList plist_;
  bool flag_symmetry_;

  // lazy assembly
  bool assembled_operator_;
  bool assembled_schur_;

  MFDMethod method_;

  // local matrices
  std::vector<WhetStone::DenseMatrix > Mff_cells_;
  std::vector<Teuchos::SerialDenseMatrix<int, double> > Aff_cells_;
  std::vector<Epetra_SerialDenseVector> Acf_cells_, Afc_cells_;
  std::vector<double> Acc_cells_;  // duplication may be useful later
  std::vector<Epetra_SerialDenseVector> Ff_cells_;
  std::vector<double> Fc_cells_;

  // global matrices
  Teuchos::RCP<Epetra_Vector> Acc_;
  Teuchos::RCP<Epetra_CrsMatrix> Acf_;
  Teuchos::RCP<Epetra_CrsMatrix> Afc_;  // We generate transpose of this matrix block.
  Teuchos::RCP<Epetra_FECrsMatrix> Aff_;
  Teuchos::RCP<Epetra_FECrsMatrix> Sff_;  // Schur complement

  // global rhs
  Teuchos::RCP<CompositeVector> rhs_;

  // diagnostics
  int nokay_;
  int npassed_; // performance of algorithms generating mass matrices

  // VectorSpace describing both domain and range
  Teuchos::RCP<CompositeVectorSpace> space_;

  // preconditioner for Schur complement
  Teuchos::RCP<AmanziPreconditioners::Preconditioner> S_pc_;

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
