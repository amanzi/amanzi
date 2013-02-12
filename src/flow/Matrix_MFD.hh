/*
This is the flow component of the Amanzi code. 

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided in the top-level COPYRIGHT file.

Authors: Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
*/

#ifndef __MATRIX_MFD_HH__
#define __MATRIX_MFD_HH__

#include <strings.h>

#include "Epetra_Map.h"
#include "Epetra_Operator.h"
#include "Epetra_MultiVector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_FECrsMatrix.h"
#include "ml_MultiLevelPreconditioner.h"

#include "Ifpack.h" 
// note that if trilinos is compiled with hypre support, then
// including Ifpack.h results in the definition of HAVE_HYPRE

#include "Teuchos_RCP.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_LAPACK.hpp"

#include "Mesh.hh"
#include "Point.hh"
#include "BoundaryFunction.hh"
#include "mfd3d.hh"

#include "Flow_State.hh"
#include "Flow_typedefs.hh"


namespace Amanzi {
namespace AmanziFlow {

class Matrix_MFD : public Epetra_Operator {
 public:
   Matrix_MFD(Teuchos::RCP<Flow_State> FS_, const Epetra_Map& map_);
  ~Matrix_MFD();

  // main methods
  void SetSymmetryProperty(bool flag_symmetry) { flag_symmetry_ = flag_symmetry; }
  void AddActionProperty(int action) { actions_ |= action; } 

  void CreateMFDmassMatrices(int mfd3d_method, std::vector<WhetStone::Tensor>& K);
  void CreateMFDrhsVectors();
  void ApplyBoundaryConditions(std::vector<int>& bc_model, std::vector<bc_tuple>& bc_values);

  virtual void CreateMFDstiffnessMatrices(Epetra_Vector& Krel_cells, Epetra_Vector& Krel_faces, int method);
  virtual void SymbolicAssembleGlobalMatrices(const Epetra_Map& super_map);
  virtual void AssembleGlobalMatrices();
  virtual void AssembleSchurComplement(std::vector<int>& bc_model, std::vector<bc_tuple>& bc_values);

  double ComputeResidual(const Epetra_Vector& solution, Epetra_Vector& residual);
  double ComputeNegativeResidual(const Epetra_Vector& solution, Epetra_Vector& residual);

  int ReduceGlobalSystem2LambdaSystem(Epetra_Vector& u);

  void DeriveDarcyMassFlux(const Epetra_Vector& solution, 
                           const Epetra_Import& face_importer, 
                           Epetra_Vector& darcy_mass_flux);

  virtual void InitPreconditioner(int method, Teuchos::ParameterList& prec_list);
  virtual void UpdatePreconditioner();

  // required methods
  int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;
  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;
  bool UseTranspose() const { return false; }
  int SetUseTranspose(bool) { return 1; }

  const Epetra_Comm& Comm() const { return *(mesh_->get_comm()); }
  const Epetra_Map& OperatorDomainMap() const { return map_; }
  const Epetra_Map& OperatorRangeMap() const { return map_; }

  const char* Label() const { return strdup("Matrix MFD"); }
  double NormInf() const { return 0.0; }
  bool HasNormInf() const { return false; }

  // development methods
  void CreateMFDmassMatrices_ScaledStability(int method, double factor, std::vector<WhetStone::Tensor>& K);
  void RescaleMFDstiffnessMatrices(const Epetra_Vector& old_scale, const Epetra_Vector& new_scale);

  // access methods
  std::vector<Teuchos::SerialDenseMatrix<int, double> >& Aff_cells() { return Aff_cells_; }
  std::vector<Epetra_SerialDenseVector>& Acf_cells() { return Acf_cells_; }
  std::vector<Epetra_SerialDenseVector>& Afc_cells() { return Afc_cells_; }
  std::vector<double>& Acc_cells() { return Acc_cells_; }
  std::vector<Epetra_SerialDenseVector>& Ff_cells() { return Ff_cells_; }
  std::vector<double>& Fc_cells() { return Fc_cells_; }
  Teuchos::RCP<Epetra_Vector>& rhs() { return rhs_; }
  Teuchos::RCP<Epetra_Vector>& rhs_faces() { return rhs_faces_; }

  Teuchos::RCP<Epetra_FECrsMatrix>& Aff() { return Aff_; }
  Teuchos::RCP<Epetra_FECrsMatrix>& Sff() { return Sff_; }
  Teuchos::RCP<Epetra_Vector>& Acc() { return Acc_; }
  Teuchos::RCP<Epetra_CrsMatrix>& Acf() { return Acf_; }
  Teuchos::RCP<Epetra_CrsMatrix>& Afc() { return Afc_; }

#ifdef HAVE_HYPRE
  Teuchos::RCP<Ifpack_Hypre> IfpHypre_Sff() { return IfpHypre_Sff_; }
#endif

  int nokay() { return nokay_; }
  int npassed() { return npassed_; }

 protected:
  Teuchos::RCP<Flow_State> FS_;
  Teuchos::RCP<AmanziMesh::Mesh> mesh_;
  Epetra_Map map_;

  bool flag_symmetry_;
  int actions_;  // applly, apply inverse, or both

  std::vector<Teuchos::SerialDenseMatrix<int, double> > Mff_cells_;
  std::vector<Teuchos::SerialDenseMatrix<int, double> > Aff_cells_;
  std::vector<Epetra_SerialDenseVector> Acf_cells_, Afc_cells_;
  std::vector<double> Acc_cells_;  // duplication may be useful later

  std::vector<Epetra_SerialDenseVector> Ff_cells_;
  std::vector<double> Fc_cells_;

  Teuchos::RCP<Epetra_Vector> Acc_;
  Teuchos::RCP<Epetra_CrsMatrix> Acf_; 
  Teuchos::RCP<Epetra_CrsMatrix> Afc_;  // We generate transpose of this matrix block. 
  Teuchos::RCP<Epetra_FECrsMatrix> Aff_;
  Teuchos::RCP<Epetra_FECrsMatrix> Sff_;  // Schur complement

  Teuchos::RCP<Epetra_Vector> rhs_;
  Teuchos::RCP<Epetra_Vector> rhs_cells_;
  Teuchos::RCP<Epetra_Vector> rhs_faces_;

  int method_;  // Preconditioners
  Teuchos::RCP<ML_Epetra::MultiLevelPreconditioner> MLprec;
  Teuchos::ParameterList ML_list;
  
  Teuchos::RCP<Ifpack_Preconditioner> ifp_prec_;
  Teuchos::ParameterList ifp_plist_;

  int nokay_, npassed_;  // performance of algorithms generating mass matrices 

#ifdef HAVE_HYPRE
  Teuchos::RCP<Ifpack_Hypre> IfpHypre_Sff_;
  double hypre_tol, hypre_strong_threshold;
  int hypre_nsmooth, hypre_ncycles;
#endif

 private:
  void operator=(const Matrix_MFD& matrix);
};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif
