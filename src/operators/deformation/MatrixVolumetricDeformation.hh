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
          const Teuchos::RCP<std::vector<std::string> >& bottom_region_list);

  MatrixVolumetricDeformation(const MatrixVolumetricDeformation& other);


  // Vector space of the Matrix's domain.
  virtual Teuchos::RCP<const CompositeVectorFactory> domain() {
    return domain_; }

  // Vector space of the Matrix's range.
  virtual Teuchos::RCP<const CompositeVectorFactory> range() {
    return range_; }

  // Virtual copy constructor.
  virtual Teuchos::RCP<CompositeMatrix> Clone() const {
    return Teuchos::rcp(new MatrixVolumetricDeformation(*this)); }

  // Apply matrix, b <-- Ax
  virtual void Apply(const CompositeVector& x,
                     const Teuchos::Ptr<CompositeVector>& b);

  // Apply the inverse, x <-- A^-1 b
  virtual void ApplyInverse(const CompositeVector& b,
                            const Teuchos::Ptr<CompositeVector>& x);

 protected:
  void InitializeFromOptions_();
  void Assemble_();
  void UpdateInverse_();

 protected:
  enum PrecMethod { PREC_METHOD_NULL,
                    TRILINOS_ML,
                    TRILINOS_ILU,
                    TRILINOS_BLOCK_ILU,
                    HYPRE_AMG,
                    HYPRE_EUCLID,
                    HYPRE_PARASAILS };
  PrecMethod prec_method_;

  // solver methods
  Teuchos::RCP<ML_Epetra::MultiLevelPreconditioner> ml_prec_;
  Teuchos::ParameterList ml_plist_;
#ifdef HAVE_HYPRE
  Teuchos::RCP<Ifpack_Hypre> IfpHypre_;
  Teuchos::ParameterList hypre_plist_;
  int hypre_ncycles_, hypre_nsmooth_;
  double hypre_tol_, hypre_strong_threshold_;
#endif
  Teuchos::RCP<Ifpack_ILU> ilu_prec_;
  Teuchos::ParameterList ilu_plist_;
  Teuchos::RCP<Ifpack_Preconditioner> ifp_prec_;
  Teuchos::ParameterList ifp_plist_;

  // local data
  Teuchos::RCP<CompositeVectorFactory> range_;
  Teuchos::RCP<CompositeVectorFactory> domain_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::ParameterList plist_;
  Teuchos::RCP<Epetra_CrsMatrix> operator_;
  Teuchos::RCP<Epetra_CrsMatrix> dVdz_;
  Teuchos::RCP<std::vector<std::string> > bottom_region_list_;

  // parameters to control optimization
  double diagonal_shift_;
  double smoothing_;
};


}
}



#endif
