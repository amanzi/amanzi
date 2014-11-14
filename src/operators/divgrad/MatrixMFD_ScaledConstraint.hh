/*
  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)

  MatrixMFD_ScaledConstraint provides an MFD implementation where the
  constraint equation can be rescaled (by, for instance, the faces rel perm)
  which allows zero diffusion coefficient.
*/

#ifndef OPERATORS_MATRIX_MFD_SCALED_CONSTRAINT_HH_
#define OPERATORS_MATRIX_MFD_SCALED_CONSTRAINT_HH_

#include <strings.h>

#include "MatrixMFD.hh"

namespace Amanzi {
namespace Operators {

class MatrixMFD_ScaledConstraint : virtual public MatrixMFD {
 public:

  // Constructor
  MatrixMFD_ScaledConstraint(Teuchos::ParameterList& plist,
                             const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) :
      MatrixMFD(plist, mesh) {
    flag_symmetry_ = false;
  }

  MatrixMFD_ScaledConstraint(const MatrixMFD_ScaledConstraint& other) :
      MatrixMFD(other) {
    if (other.Krel_ != Teuchos::null) {
      Krel_ = Teuchos::rcp(new Epetra_Vector(*other.Krel_));
    }
    flag_symmetry_ = false;
  }

  // By definition, this is NOT symmetric
  void set_symmetric(bool flag_symmetry);

  // Main computational methods
  // -- local matrices
  virtual void CreateMFDstiffnessMatrices(
      const Teuchos::Ptr<const CompositeVector>& Krel);

  virtual void ApplyBoundaryConditions(const std::vector<MatrixBC>& bc_markers,
				       const std::vector<double>& bc_values, 
				       bool ADD_BC_FLUX=true);

  // First derivative quantities.
  virtual void DeriveFlux(const CompositeVector& solution,
                          const Teuchos::Ptr<CompositeVector>& flux) const;

 protected:
  virtual void CreateMatrices_(const Epetra_CrsGraph& cf_graph,
          const Epetra_FECrsGraph& ff_graph);

 protected:
  Teuchos::RCP<Epetra_Vector> Krel_;

  friend class MatrixMFD_Coupled;
  friend class MatrixMFD_Permafrost;

};


}  // namespace AmanziFlow
}  // namespace Amanzi

#endif
