/*
  This is the flow component of the Amanzi code.

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

#include "matrix_mfd.hh"


namespace Amanzi {
namespace Operators {

class MatrixMFD_Surf;


class MatrixMFD_TPFA : public MatrixMFD {
 public:
  MatrixMFD_TPFA(Teuchos::ParameterList& plist,
                 const Teuchos::RCP<const AmanziMesh::Mesh> mesh) :
      MatrixMFD(plist,mesh) {}

  MatrixMFD_TPFA(const MatrixMFD& other) :
    MatrixMFD(other) {}

  virtual ~MatrixMFD_TPFA() {};

  // override main methods of the base class
  virtual void CreateMFDstiffnessMatrices(const Teuchos::Ptr<const CompositeVector>& Krel);
  virtual void SymbolicAssembleGlobalMatrices();
  virtual void AssembleGlobalMatrices();
  virtual void ComputeSchurComplement(const std::vector<Matrix_bc>& bc_markers,
          const std::vector<double>& bc_values) {}


  virtual void Apply(const CompositeVector& X,
             const Teuchos::Ptr<CompositeVector>& Y) const;
  virtual void ApplyInverse(const CompositeVector& X,
                    const Teuchos::Ptr<CompositeVector>& Y) const;

  virtual void InitPreconditioner(Teuchos::ParameterList& prec_list);
  virtual void UpdatePreconditioner();

  const char* Label() const {
    return strdup("Matrix MFD_TPFA");
  }

  Teuchos::RCP<Epetra_FECrsMatrix> TPFA() {
    return Spp_;
  }

 protected:

  Teuchos::RCP<CompositeVector> Dff_;
  Teuchos::RCP<Epetra_FECrsMatrix> Spp_;  // Explicit Schur complement
  Teuchos::RCP<Epetra_FECrsMatrix> NumJac_;  // Numerical Jacobian

#ifdef HAVE_HYPRE
  Teuchos::RCP<Ifpack_Hypre> IfpHypre_Spp_;
#endif

 private:
  void operator=(const MatrixMFD_TPFA& matrix);

  friend class MatrixMFD_Surf;
};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif
