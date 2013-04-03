/*
  License: BSD
  Authors: Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
  Ethan Coon (ecoon@lanl.gov)
  MatrixMFD provides a mimetic discretization for the elliptic operator div K grad u.

*/

#ifndef OPERATORS_MATRIX_MFD_HH_
#define OPERATORS_MATRIX_MFD_HH_

#include <strings.h>

#include "Epetra_Map.h"
#include "Epetra_Operator.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_FECrsMatrix.h"
#include "ml_MultiLevelPreconditioner.h"

#include "Ifpack.h"
#include "Ifpack_ILU.h"
#include "Ifpack_AdditiveSchwarz.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_LAPACK.hpp"

#include "Mesh.hh"
#include "Point.hh"
#include "tree_vector.hh"
#include "composite_vector.hh"
#include "mfd3d.hpp"

#include "matrix.hh"


namespace Amanzi {
namespace Operators {

enum MFD_method {
  MFD_NULL = 0,
  MFD_POLYHEDRA = 1,
  MFD_POLYHEDRA_MONOTONE = 2,  // under development
  MFD_HEXAHEDRA_MONOTONE = 3,
  MFD_TWO_POINT_FLUX = 4, // without consistency
  MFD_SUPPORT_OPERATOR = 5, // rc1 compatibility
  MFD_OPTIMIZED = 6
};

const int MFD_HEX_FACES = 6;  // Hexahedron is the common element
const int MFD_HEX_NODES = 8;
const int MFD_HEX_EDGES = 12;

const int MFD_QUAD_FACES = 4;  // Quadrilateral is the common element
const int MFD_QUAD_NODES = 4;
const int MFD_QUAD_EDGES = 4;

const int MFD_MAX_FACES = 14;  // Kelvin's tetrakaidecahedron
const int MFD_MAX_NODES = 47;  // These polyhedron parameters must
const int MFD_MAX_EDGES = 60;  // be calculated in Init().



class MatrixMFD : public Matrix {

 public:
  MatrixMFD(Teuchos::ParameterList& plist,
            const Teuchos::RCP<const AmanziMesh::Mesh> mesh);

  MatrixMFD(const MatrixMFD& other);

  virtual ~MatrixMFD() {};

  void InitializeFromPList_();

  // access to data for updating manually
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
  Teuchos::RCP<Epetra_FECrsMatrix> Schur() {
    return Sff_;
  }


  std::vector<double>& Fc_cells() {
    return Fc_cells_;
  }
  std::vector<Epetra_SerialDenseVector>& Ff_cells() {
    return Ff_cells_;
  }
  Teuchos::RCP<const Epetra_Vector> Acc() {
    return Acc_;
  }
  Teuchos::RCP<const Epetra_FECrsMatrix> Aff() {
    return Aff_;
  }

  // performance of algorithms generating mass matrices
  int nokay() {
    return nokay_;
  }
  int npassed() {
    return npassed_;
  }

  // main computational methods
  void SetSymmetryProperty(bool flag_symmetry) {
    flag_symmetry_ = flag_symmetry;
  }

  void CreateMFDmassMatrices(const Teuchos::Ptr<std::vector<WhetStone::Tensor> >& K);
  virtual void CreateMFDstiffnessMatrices(const Teuchos::Ptr<const CompositeVector>& Krel);
  void RescaleMFDstiffnessMatrices(const Epetra_Vector& old_scale,
          const Epetra_Vector& new_scale);
  void CreateMFDrhsVectors();

  Teuchos::RCP<CompositeVector>& rhs() {
    return rhs_;
  }
  void InitializeSuperVecs(const CompositeVector& sample);

  virtual void ApplyBoundaryConditions(const std::vector<Matrix_bc>& bc_markers,
          const std::vector<double>& bc_values);

  virtual void SymbolicAssembleGlobalMatrices();
  virtual void AssembleGlobalMatrices();
  virtual void ComputeSchurComplement(const std::vector<Matrix_bc>& bc_markers,
          const std::vector<double>& bc_values);

  // operator methods: MatrixMFD is an Epetra_Operator
  // int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;
  // int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;
  bool UseTranspose() const {
    return false;
  }
  int SetUseTranspose(bool) {
    return 1;
  }

  const Epetra_Comm& Comm() const {
    return *(mesh_->get_comm());
  }
  const Epetra_Map& OperatorDomainMap() const {
    return *supermap_;
  }
  const Epetra_Map& OperatorRangeMap() const {
    return *supermap_;
  }

  const char* Label() const {
    return strdup("Matrix MFD");
  }
  double NormInf() const {
    return 0.0;
  }
  bool HasNormInf() const {
    return false;
  }

  // operator methods for CompositeVectors
  virtual void Apply(const CompositeVector& X,
                     const Teuchos::Ptr<CompositeVector>& Y) const;
  virtual void ApplyInverse(const CompositeVector& X,
                            const Teuchos::Ptr<CompositeVector>& Y) const;

  virtual void Apply(const TreeVector& X,
                     const Teuchos::Ptr<TreeVector>& Y) const {
    return Apply(*X.data(), Y->data().ptr());
  }

  virtual void ApplyInverse(const TreeVector& X,
                            const Teuchos::Ptr<TreeVector>& Y) const {
    return ApplyInverse(*X.data(), Y->data().ptr());
  }

  void ComputeResidual(const CompositeVector& X,
                       const Teuchos::Ptr<CompositeVector>& F) const;
  void ComputeNegativeResidual(const CompositeVector& X,
          const Teuchos::Ptr<CompositeVector>& F) const;

  // extra methods for preconditioning
  virtual void InitPreconditioner();
  virtual void UpdatePreconditioner();

  // extra methods for convenience
  void DeriveFlux(const CompositeVector& solution,
                  const Teuchos::Ptr<CompositeVector>& flux) const;
  void DeriveCellVelocity(const CompositeVector& flux,
                          const Teuchos::Ptr<CompositeVector>& velocity) const;

  // development methods
  virtual void UpdateConsistentFaceConstraints(const Teuchos::Ptr<CompositeVector>& u);
  virtual void UpdateConsistentFaceCorrection(const CompositeVector& u,
          const Teuchos::Ptr<CompositeVector>& Pu);

 protected:
  virtual void FillMatrixGraphs_(const Teuchos::Ptr<Epetra_CrsGraph> cf_graph,
          const Teuchos::Ptr<Epetra_FECrsGraph> ff_graph);
  virtual void CreateMatrices_(const Epetra_CrsGraph& cf_graph,
          const Epetra_FECrsGraph& ff_graph);

 protected:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::ParameterList plist_;
  bool flag_symmetry_;
  MFD_method method_;

  std::vector<Teuchos::SerialDenseMatrix<int, double> > Mff_cells_;
  std::vector<Teuchos::SerialDenseMatrix<int, double> > Aff_cells_;
  std::vector<Epetra_SerialDenseVector> Acf_cells_, Afc_cells_;
  std::vector<double> Acc_cells_;  // duplication may be useful later

  std::vector<Epetra_SerialDenseVector> Ff_cells_;
  std::vector<double> Fc_cells_;

  Teuchos::RCP<Epetra_Vector> Acc_;
  Teuchos::RCP<Epetra_CrsMatrix> Acf_;
  Teuchos::RCP<Epetra_CrsMatrix> Afc_;  // We generate transpose of this matrix block.
  Teuchos::RCP<Epetra_FECrsMatrix> Aff_;
  Teuchos::RCP<Epetra_FECrsMatrix> Sff_;  // Schur complement

  Teuchos::RCP<CompositeVector> rhs_;

  int nokay_;
  int npassed_; // performance of algorithms generating mass matrices

  enum PrecMethod { PREC_METHOD_NULL,
                    TRILINOS_ML,
                    TRILINOS_ILU,
                    TRILINOS_BLOCK_ILU,
                    HYPRE_AMG,
                    HYPRE_EUCLID,
                    HYPRE_PARASAILS };
  PrecMethod prec_method_;

  Teuchos::RCP<ML_Epetra::MultiLevelPreconditioner> ml_prec_;
  Teuchos::ParameterList ml_plist_;

#ifdef HAVE_HYPRE
  Teuchos::RCP<Ifpack_Hypre> IfpHypre_Sff_;
  Teuchos::ParameterList hypre_plist_;
  int hypre_ncycles_, hypre_nsmooth_;
  double hypre_tol_, hypre_strong_threshold_;
#endif

  Teuchos::RCP<Ifpack_ILU> ilu_prec_;
  Teuchos::ParameterList ilu_plist_;

  Teuchos::RCP<Ifpack_Preconditioner> ifp_prec_;
  Teuchos::ParameterList ifp_plist_;

  Teuchos::RCP<const Epetra_Map> supermap_;
  Teuchos::RCP<CompositeVector> vector_x_; // work vectors for AztecOO
  Teuchos::RCP<CompositeVector> vector_y_;

  friend class MatrixCoupledMFD;
};


}  // namespace AmanziFlow
}  // namespace Amanzi

#endif
