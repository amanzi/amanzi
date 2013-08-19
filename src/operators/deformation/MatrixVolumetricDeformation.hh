/*
  License:
  Authors: Ethan Coon (ecoon@lanl.gov) (ATS version)

  Deformation optimization matrix
*/

#ifndef OPERATORS_MATRIX_VOLUMETRIC_DEFORMATION_HH_
#define OPERATORS_MATRIX_VOLUMETRIC_DEFORMATION_HH_


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

#include "Mesh.hh"
#include "composite_vector.hh"
#include "composite_matrix.hh"


namespace Amanzi {
namespace Operators {

class MatrixVolumetricDeformation : public CompositeMatrix {

 public:

  MatrixVolumetricDeformation(Teuchos::ParameterList& plist,
          const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
          const Teuchos::RCP<const AmanziMesh::Entity_ID_List>& fixed_nodes);

  MatrixVolumetricDeformation(const MatrixVolumetricDeformation& other);


  // Vector space of the Matrix's domain.
  virtual Teuchos::RCP<const CompositeVectorFactory> domain() const {
    return domain_; }

  // Vector space of the Matrix's range.
  virtual Teuchos::RCP<const CompositeVectorFactory> range() const {
    return range_; }

  // Virtual copy constructor.
  virtual Teuchos::RCP<CompositeMatrix> Clone() const {
    return Teuchos::rcp(new MatrixVolumetricDeformation(*this)); }

  // Apply matrix, b <-- Ax
  virtual void Apply(const CompositeVector& x,
                     const Teuchos::Ptr<CompositeVector>& b) const;

  // Apply the inverse, x <-- A^-1 b
  virtual void ApplyInverse(const CompositeVector& b,
                            const Teuchos::Ptr<CompositeVector>& x) const;

  // This is a Normal equation, so we need to apply N^T to the rhs
  void ApplyRHS(const CompositeVector& x_cell,
                const Teuchos::Ptr<CompositeVector>& x_node) const;


 protected:
  void InitializeFromOptions_();
  void Assemble_();
  void UpdateInverse_();

 protected:
  // solver methods
  Teuchos::RCP<Ifpack_Hypre> IfpHypre_;
  Teuchos::ParameterList hypre_plist_;
  int hypre_ncycles_, hypre_nsmooth_;
  double hypre_tol_, hypre_strong_threshold_;
  int hypre_coarsen_type_, hypre_relax_type_;
  int hypre_verbose_, hypre_cycle_type_;

  // local data
  Teuchos::RCP<CompositeVectorFactory> range_;
  Teuchos::RCP<CompositeVectorFactory> domain_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::ParameterList plist_;
  Teuchos::RCP<Epetra_CrsMatrix> operator_;
  Teuchos::RCP<Epetra_CrsMatrix> dVdz_;
  Teuchos::RCP<const AmanziMesh::Entity_ID_List> fixed_nodes_;

  // parameters to control optimization
  double diagonal_shift_;
  double smoothing_;
};


}
}



#endif
