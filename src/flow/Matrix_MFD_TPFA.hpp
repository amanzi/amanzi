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

#ifndef __MATRIX_MFD_TPFA_HPP__
#define __MATRIX_MFD_TPFA_HPP__

#include <strings.h>

#include "Teuchos_RCP.hpp"
#include "Ifpack.h" 

#include "Matrix_MFD.hpp"

namespace Amanzi {
namespace AmanziFlow {

class Matrix_MFD_TPFA : public Matrix_MFD {
 public:
  Matrix_MFD_TPFA(Teuchos::RCP<Flow_State> FS_, const Epetra_Map& map_) : Matrix_MFD(FS_, map_) {};
  ~Matrix_MFD_TPFA() {};

  // override main methods of the base class
  void CreateMFDstiffnessMatrices(Epetra_Vector& Krel_cells, Epetra_Vector& Krel_faces);
  void SymbolicAssembleGlobalMatrices(const Epetra_Map& super_map);
  void AssembleGlobalMatrices();
  void ComputeSchurComplement(std::vector<int>& bc_model, std::vector<bc_tuple>& bc_values) {};
  
   void AnalyticJacobian(const Epetra_Vector& solution,
                        int dim,
                        int Krel_method,
                        std::vector<int>& bc_markers, 
                        std::vector<bc_tuple>& bc_values,
                        Epetra_Vector& Krel_cells, 
                        Epetra_Vector& dK_dP_cells,
                        Epetra_Vector& Krel_faces, 
                        Epetra_Vector& dK_dP_faces);

  int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;
  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  void InitPreconditioner(int method, Teuchos::ParameterList& prec_list);
  void UpdatePreconditioner();
  
  void AddCol2NumJacob(int irow, Epetra_Vector& r);
  void CompareJacobians();

  const char* Label() const { return strdup("Matrix MFD_TPFA"); }
  bool IsPureNewton() const { return true; } 

 private:
  void ComputeJacobianLocal(int mcells,
                            int face_id,
                            int dim,
                            int Krel_method,
                            std::vector<int>& bc_markers,
                            std::vector<bc_tuple>& bc_values,
                            double dist,
                            double *pres,
                            double *perm_abs_vert,
                            double *perm_abs_horz,
                            double *k_rel,
                            double *dk_dp_cell,
                            AmanziGeometry::Point& normal,
                            Teuchos::SerialDenseMatrix<int, double>& Jpp);
         
  Teuchos::RCP<Epetra_Vector> Dff_;
  Teuchos::RCP<Epetra_FECrsMatrix> Spp_;  // Explicit Schur complement
  Teuchos::RCP<Epetra_FECrsMatrix> NumJac_;  // Numerical Jacobian

#ifdef HAVE_HYPRE
  Teuchos::RCP<Ifpack_Hypre> IfpHypre_Spp_;
#endif

 private:
  void operator=(const Matrix_MFD_TPFA& matrix);
};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif
