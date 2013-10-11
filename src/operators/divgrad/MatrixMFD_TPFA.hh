/*
  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
  Daniil Svyatskiy (dasvyat@lanl.gov)
         Ethan Coon (ATS version) (ecoon@lanl.gov)

  The class provides a different implementation of solvers than in
  the base class. In particular, Lagrange multipliers are elliminated
  from the DAE system and short vectors are used in the nonlinear solver.
*/

#ifndef OPERATORS_MATRIX_MFD_TPFA_HH__
#define OPERATORS_MATRIX_MFD_TPFA_HH__

#include <strings.h>

#include "Teuchos_RCP.hpp"
#include "Ifpack.h"

#include "MatrixMFD.hh"
#include "MatrixMFD_Factory.hh"

namespace Amanzi {
namespace Operators {

class MatrixMFD_TPFA : virtual public MatrixMFD {
 public:
  MatrixMFD_TPFA(Teuchos::ParameterList& plist,
                 const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) :
      MatrixMFD(plist,mesh) {}

  // override main methods of the base class
  virtual void CreateMFDstiffnessMatrices(const Teuchos::Ptr<const CompositeVector>& Krel);
  virtual void SymbolicAssembleGlobalMatrices();
  virtual void AssembleGlobalMatrices();
  virtual void ComputeSchurComplement(const std::vector<MatrixBC>& bc_markers,
          const std::vector<double>& bc_values) { Sff_ = Spp_; }

  virtual int Apply(const CompositeVector& X,
                     CompositeVector& Y) const;
  virtual int ApplyInverse(const CompositeVector& X,
                            CompositeVector& Y) const;

  void AnalyticJacobian(const CompositeVector& height,
                        const CompositeVector& potential,
                        const CompositeVector& Krel,
                        const CompositeVector& dKrel_dp,
                        const CompositeVector& Krel_cell,
                        const CompositeVector& dKrel_cell_dp);

  Teuchos::RCP<Epetra_FECrsMatrix> TPFA() {
    AssertAssembledOperator_or_die_();
    return Spp_;
  }

  virtual void UpdateConsistentFaceConstraints(const Teuchos::Ptr<CompositeVector>& u);
  virtual void UpdateConsistentFaceCorrection(const CompositeVector& u,
          const Teuchos::Ptr<CompositeVector>& Pu);

 protected:
  void ComputeJacobianLocal(int mcells,
                            int face_id,
                            double dist,
                            double *height,
                            double *potential,
                            double *k_rel,
                            double *dk_rel_dp,
                            double *k_rel_cell,
                            double *dk_rel_cell_dp,
                            Teuchos::SerialDenseMatrix<int, double>& Jpp);

 protected:
  Teuchos::RCP<CompositeVector> Dff_;
  Teuchos::RCP<Epetra_FECrsMatrix> Spp_;  // Explicit Schur complement

 private:
  MatrixMFD_TPFA(const MatrixMFD& other);
  void operator=(const MatrixMFD_TPFA& matrix);

  friend class MatrixMFD_Surf;

private:
  // factory registration
  static RegisteredMatrixMFD_Factory<MatrixMFD_TPFA> reg_;

};

}  // namespace Operators
}  // namespace Amanzi

#endif
