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

#ifndef OPERATORS_MATRIX_TPFA_HH_
#define OPERATORS_MATRIX_TPFA_HH_

#include <strings.h>

#include "Teuchos_RCP.hpp"
#include "MatrixMFD.hh"
//#include "BlockMatrix.hh"
#include "upwinding.hh"

namespace Amanzi {
namespace Operators {

class Matrix_TPFA :  public MatrixMFD {
 public:
  Matrix_TPFA(Teuchos::ParameterList& plist,
	      const Teuchos::RCP<const AmanziMesh::Mesh>& mesh);

  // override main methods of the base class
  virtual void CreateMFDmassMatrices(const Teuchos::Ptr<std::vector<WhetStone::Tensor> >& K);
  virtual void CreateMFDstiffnessMatrices(const Teuchos::Ptr<const CompositeVector>& Krel);
  // virtual void CreateMFDrhsVectors();
  virtual void SymbolicAssembleGlobalMatrices();
  virtual void AssembleGlobalMatrices(){};
  virtual void ComputeSchurComplement(const std::vector<MatrixBC>& bc_markers,
				      const std::vector<double>& bc_values);
  virtual void ApplyBoundaryConditions(const std::vector<MatrixBC>& bc_markers,
				       const std::vector<double>& bc_values, bool ADD_BC_FLUX = true);

  virtual void ComputeNegativeResidual(const CompositeVector& solution,
				       const Teuchos::Ptr<CompositeVector>& residual) const;

  virtual void FillMatrixGraphs_(const Teuchos::Ptr<Epetra_CrsGraph> cf_graph,
				 const Teuchos::Ptr<Epetra_FECrsGraph> ff_graph);
  virtual int Apply(const CompositeVector& X,
                     CompositeVector& Y) const;
  virtual int ApplyInverse(const CompositeVector& X,
                            CompositeVector& Y) const;

  // void DeriveFluxFV(const CompositeVector& solution,
  // 		    const Teuchos::Ptr<CompositeVector>& flux,
  // 		    const Teuchos::Ptr<const CompositeVector>& rel_perm);

  virtual void DeriveFlux(const CompositeVector& solution,
			    const Teuchos::Ptr<CompositeVector>& flux) const;

  virtual void UpdatePreconditioner_() const;


  void ComputeTransmissibilities_(const Teuchos::Ptr<std::vector<WhetStone::Tensor> >& K);
  Teuchos::RCP<Epetra_Vector>  gravity_terms() {return gravity_term_;}


  int Apply(const Epetra_MultiVector& X,
            Epetra_MultiVector& Y) const;
  
  void AnalyticJacobian(const Upwinding& upwinding,
			const Teuchos::Ptr<State>& S,
			std::string potential_key,
			const CompositeVector& dconductivity,
			const std::vector<MatrixBC>& bc_markers,
			const std::vector<double>& bc_values);

  Teuchos::RCP<Epetra_FECrsMatrix> TPFA() {
    AssertAssembledOperator_or_die_();
    return Spp_;
  }


  //  virtual void AssembleAff_() const;
  virtual void AssembleRHS_() const;
  virtual void AssembleSchur_() const;

  virtual void UpdateConsistentFaceConstraints(const Teuchos::Ptr<CompositeVector>& u);
  virtual void UpdateConsistentFaceCorrection(const CompositeVector& u,
          const Teuchos::Ptr<CompositeVector>& Pu);

  Teuchos::RCP<const CompositeVector> Dff() const {
    return Dff_;
  }

  double BoundaryValue(const Amanzi::CompositeVector& solution, int face_id) const;
  void SetBoundaryValue(Amanzi::CompositeVector& solution, int face_id, double value);

 protected:
  mutable Teuchos::RCP<CompositeVector> Dff_;
  mutable Teuchos::RCP<Epetra_FECrsMatrix> Spp_;  // Explicit Schur complement
  mutable Teuchos::RCP<Epetra_CrsMatrix> Afc_;
  mutable Teuchos::RCP<Epetra_CrsMatrix> Acf_;
  //mutable Teuchos::RCP<AmanziPreconditioners::Preconditioner> Aff_pc_;
  bool cells_only_;

  Teuchos::RCP<Epetra_Vector> transmissibility_;
  Teuchos::RCP<Epetra_Vector> rel_perm_transmissibility_;
  Teuchos::RCP<Epetra_Vector> gravity_term_;
  std::vector<int> face_flag_;  

 private:
  Matrix_TPFA(const MatrixMFD& other);
  void operator=(const Matrix_TPFA& matrix);

};

}  // namespace Operators
}  // namespace Amanzi

#endif
