/*
  License: BSD
  Authors: Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)

  MatrixMFD provides a mimetic discretization for the elliptic operators div K grad u.

*/

#ifndef OPERATORS_MATRIX_MFD_HH_
#define OPERATORS_MATRIX_MFD_HH_

#include <strings.h>

#include "Epetra_Map.h"
#include "Epetra_Operator.h"
#include "Epetra_MultiVector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_FECrsMatrix.h"
#include "ml_MultiLevelPreconditioner.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_LAPACK.hpp"

#include "Mesh.hh"
#include "Point.hh"
#include "boundary-function.hh"
#include "mfd3d.hpp"

namespace Amanzi {
namespace Operator {

enum MFD_method {
  MFD_NULL = 0,
  MFD_POLYHEDRA = 1,
  MFD_POLYHEDRA_MONOTONE = 2,  // under development
  MFD_HEXAHEDRA_MONOTONE = 3
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

enum Matrix_bc {
  MATRIX_BC_NULL = 0,
  MATRIX_BC_DIRICHLET,
  MATRIX_BC_FLUX
};

class MatrixMFD {

public:
  Matrix_MFD(Teuchos::ParameterList& plist, Teuchos::RCP<AmanziMesh::Mesh> mesh) :
    plist_(plist), mesh_(mesh) {}

  // main methods
  void SetSymmetryProperty(bool flag_symmetry) { flag_symmetry_ = flag_symmetry; }

  void CreateMFDmassMatrices(MFD_method method, std::vector<WhetStone::Tensor>& K);
  void CreateMFDrhsVectors();
  void CreateMFDstiffnessMatrices(MFD_method method,
          const std::vector<WhetStone::Tensor>& K, const CompositeVector& K_faces);
  void RescaleMFDstiffnessMatrices(const Epetra_Vector& old_scale,
          const Epetra_Vector& new_scale);
  void ApplyBoundaryConditions(const std::vector<Matrix_bc>& bc_markers,
          const std::vector<double>& bc_values);

  void SymbolicAssembleGlobalMatrices();
  void AssembleGlobalMatrices();
  void ComputeSchurComplement(const std::vector<Matrix_bc>& bc_markers,
          const std::vector<double>& bc_values);

  void Apply(const CompositeVector& X,
             const Teuchos::RCP<CompositeVector>& Y) const;
  void ApplyInverse(const CompositeVector& X,
                    const Teuchos::RCP<CompositeVector>& Y) const;
  void ComputeResidual(const CompositeVector& X,
                       const Teuchos::RCP<CompositeVector>& F) const;
  void ComputeNegativeResidual(const CompositeVector& X,
                               const Teuchos::RCP<CompositeVector>& F) const;

  void DeriveFlux(const CompositeVector& solution,
                  const Teuchos::RCP<CompositeVector>& flux) const;
  void DeriveCellVelocity(const CompositeVector& flux,
                          const Teuchos::RCP<CompositeVector>& velocity) const;

  void InitMLPreconditioner(Teuchos::ParameterList& ml_plist_);
  void UpdateMLPreconditioner();

  bool UseTranspose() const { return false; }
  int SetUseTranspose(bool) { return 1; }

  // access methods
  Teuchos::RCP<Epetra_Vector>& rhs() { return rhs_; }

private:
  Teuchos::RCP<AmanziMesh::Mesh> mesh_;
  bool flag_symmetry_;

  std::vector<Teuchos::SerialDenseMatrix<int, double> > Aff_cells_;
  std::vector<Epetra_SerialDenseVector> Acf_cells_, Afc_cells_;
  std::vector<double> Acc_cells;  // duplication may be useful later

  std::vector<Epetra_SerialDenseVector> Ff_cells_;
  std::vector<double> Fc_cells_;

  Teuchos::RCP<Epetra_Vector> Acc_;
  Teuchos::RCP<Epetra_CrsMatrix> Acf_;
  Teuchos::RCP<Epetra_CrsMatrix> Afc_;  // We generate transpose of this matrix block.
  Teuchos::RCP<Epetra_FECrsMatrix> Aff_;
  Teuchos::RCP<Epetra_FECrsMatrix> Sff_;  // Schur complement

  Teuchos::RCP<CompositeVector> rhs_;

  Teuchos::RCP<ML_Epetra::MultiLevelPreconditioner> ml_prec_;
  Teuchos::ParameterList ml_plist_;
};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif
