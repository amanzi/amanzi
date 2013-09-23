/*
This is the flow component of the Amanzi code.  

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided in the top-level COPYRIGHT file.

Authors: Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
         Daniil Svyatskiy (dasvyat@lanl.gov)

The class provides a different implementation of solvers than in 
the base class. In particular, Lagrange multipliers are elliminated
from the DAE system and short vectors are used in the nonlinear solver.
*/

#ifndef AMANZI_MATRIX_TPFA_HH__
#define AMANZI_MATRIX_TPFA_HH__

#include <strings.h>

#include "Teuchos_RCP.hpp"
#include "Ifpack.h" 

#include "Matrix_MFD.hh"


namespace Amanzi {
namespace AmanziFlow {

class Matrix_MFD_TPFA : public Matrix_MFD {
 public:
  Matrix_MFD_TPFA(Teuchos::RCP<Flow_State> FS, Teuchos::RCP<const Epetra_Map> map);
  Matrix_MFD_TPFA(Teuchos::RCP<Flow_State> FS, 
		  Teuchos::RCP<const Epetra_Map> map,
		  Teuchos::RCP<Epetra_Vector> Krel_faces,
		  Teuchos::RCP<Epetra_Vector> Trans_faces,
		  Teuchos::RCP<Epetra_Vector> Grav_faces);
  ~Matrix_MFD_TPFA() {};

  void Set_Krel_faces (Teuchos::RCP<Epetra_Vector> Krel_faces) { Krel_faces_ = Krel_faces;}
  void Set_Trans_faces(Teuchos::RCP<Epetra_Vector> Trans_faces) { trans_on_faces_ = Trans_faces;}
  void Set_Grav_faces (Teuchos::RCP<Epetra_Vector> Grav_faces) { grav_on_faces_ = Grav_faces;}

  // override main methods of the base class
  virtual void CreateMFDstiffnessMatrices(RelativePermeability& rel_perm);
  virtual void SymbolicAssembleGlobalMatrices(const Epetra_Map& super_map);
  virtual void AssembleGlobalMatrices();
  virtual void AssembleSchurComplement(std::vector<int>& bc_model, std::vector<bc_tuple>& bc_values);
  // void AssembleGlobalMatrices(const Epetra_Vector& Krel_faces, const Epetra_Vector& Trans_faces);
  // void AssembleSchurComplement(const Epetra_Vector& Krel_faces, const Epetra_Vector& Trans_faces);
  
  double ComputeNegativeResidual(const Epetra_Vector& solution,  
				 Epetra_Vector& residual);

  void UpdatePreconditioner() { preconditioner_->Update(Spp_); } 


  void AnalyticJacobian(const Epetra_Vector& solution, int dim,
                        std::vector<int>& bc_markers, std::vector<bc_tuple>& bc_values,
			const Epetra_Vector& trans_faces,
			const Epetra_Vector& grav_term_faces,
                        RelativePermeability& rel_perm); 




  void ApplyBoundaryConditions(std::vector<int>& bc_model,
			       std::vector<bc_tuple>& bc_values); 

  void DeriveDarcyMassFlux(const Epetra_Vector& solution,
			   std::vector<int>& bc_model, 
			   std::vector<bc_tuple>& bc_values,
			   Epetra_Vector& darcy_mass_flux);



  int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;
  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;



 private:
  void ComputeJacobianLocal(int mcells,
                            int face_id,
                            int Krel_method,
                            std::vector<int>& bc_markers,
                            std::vector<bc_tuple>& bc_values,
                            double *pres,
                            double *dk_dp_cell,
			    const Epetra_Vector& trans_faces,
			    const Epetra_Vector& grav_term_faces,
                            Teuchos::SerialDenseMatrix<int, double>& Jpp);
         
  Teuchos::RCP<Epetra_Vector> Dff_;
  Teuchos::RCP<Epetra_FECrsMatrix> Spp_;  // Explicit Schur complement

  Teuchos::RCP<Epetra_Vector> Krel_faces_;
  Teuchos::RCP<Epetra_Vector> trans_on_faces_;
  Teuchos::RCP<Epetra_Vector> grav_on_faces_;


 private:
  void operator=(const Matrix_MFD_TPFA& matrix);
};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif
