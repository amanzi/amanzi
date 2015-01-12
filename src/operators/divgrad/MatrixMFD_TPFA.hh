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

#ifndef OPERATORS_MATRIX_MFD_TPFA_HH_
#define OPERATORS_MATRIX_MFD_TPFA_HH_

#include <strings.h>

#include "Teuchos_RCP.hpp"
#include "MatrixMFD.hh"
#include "upwinding.hh"

namespace Amanzi {
namespace Operators {

class MatrixMFD_TPFA : virtual public MatrixMFD {
 public:
  MatrixMFD_TPFA(Teuchos::ParameterList& plist,
                 const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) :
      MatrixMFD(plist,mesh),
      assembled_app_(false) {
    cells_only_ = plist.get<bool>("TPFA use cells only", false);
    if (cells_only_) {
      space_ = Teuchos::rcp(new CompositeVectorSpace());
      space_->SetMesh(mesh_)->SetComponent("cell", AmanziMesh::CELL, 1);
    }
  }

  // override main methods of the base class
  virtual void CreateMFDstiffnessMatrices(const Teuchos::Ptr<const CompositeVector>& Krel);
  virtual void SymbolicAssembleGlobalMatrices();

  virtual int Apply(const CompositeVector& X,
                     CompositeVector& Y) const;
  virtual int ApplyInverse(const CompositeVector& X,
                            CompositeVector& Y) const;

  int Apply(const Epetra_MultiVector& X,
            Epetra_MultiVector& Y) const;
  
  void AnalyticJacobian(const Upwinding& upwinding,
			const Teuchos::Ptr<State>& S,
			std::string potential_key,
			const CompositeVector& dconductivity,
			const std::vector<MatrixBC>& bc_markers,
			const std::vector<double>& bc_values);

  Teuchos::RCP<Epetra_FECrsMatrix> TPFA() {
    if (!assembled_app_) AssembleApp_();
    return App_;
  }

  virtual void UpdateConsistentFaceConstraints(const Teuchos::Ptr<CompositeVector>& u);
  virtual void UpdateConsistentFaceCorrection(const CompositeVector& u,
          const Teuchos::Ptr<CompositeVector>& Pu);

  Teuchos::RCP<const CompositeVector> Dff() const {
    return Dff_;
  }

 protected:
  virtual void MarkLocalMatricesAsChanged_() {
    assembled_operator_ = false;
    assembled_schur_ = false;
    assembled_rhs_ = false;
    assembled_dff_ = false;
    assembled_app_ = false;
  }

  virtual void AssembleApp_() const;
  virtual void AssembleRHS_() const;
  virtual void AssembleSchur_() const;
  virtual void AssembleDff_() const;

 protected:
  mutable Teuchos::RCP<CompositeVector> Dff_;
  mutable Teuchos::RCP<Epetra_FECrsMatrix> App_;  // Explicit Schur complement
  bool cells_only_;
  mutable bool assembled_app_;
  mutable bool assembled_dff_;

 private:
  MatrixMFD_TPFA(const MatrixMFD& other);
  void operator=(const MatrixMFD_TPFA& matrix);

};

}  // namespace Operators
}  // namespace Amanzi

#endif
