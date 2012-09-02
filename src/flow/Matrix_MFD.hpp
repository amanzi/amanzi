/*
This is the flow component of the Amanzi code. 
License: BSD
Authors: Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
*/

#ifndef __MATRIX_MFD_HPP__
#define __MATRIX_MFD_HPP__

#include <strings.h>

#include "Epetra_Map.h"
#include "Epetra_Operator.h"
#include "Epetra_MultiVector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_FECrsMatrix.h"
#include "ml_MultiLevelPreconditioner.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_LAPACK.hpp"

#include "Mesh.hh"
#include "Point.hh"
#include "boundary_function.hh"
#include "mfd3d.hpp"

#include "Flow_State.hpp"

#include "Ifpack.h"
#include "Ifpack_Hypre.h"

namespace Amanzi {
namespace AmanziFlow {

class Matrix_MFD : public Epetra_Operator {
 public:
   Matrix_MFD(Teuchos::RCP<Flow_State> FS_, const Epetra_Map& map_) : FS(FS_), map(map_) { mesh_ = FS->mesh(); }
  ~Matrix_MFD();

  // main methods
  void SetSymmetryProperty(bool flag_symmetry) { flag_symmetry_ = flag_symmetry; }
  void CreateMFDmassMatrices(int mfd3d_method, std::vector<WhetStone::Tensor>& K);
  void CreateMFDrhsVectors();
  void CreateMFDstiffnessMatrices(Epetra_Vector& Krel_cells, Epetra_Vector& Krel_faces);
  void RescaleMFDstiffnessMatrices(const Epetra_Vector& old_scale, const Epetra_Vector& new_scale);
  void ApplyBoundaryConditions(std::vector<int>& bc_markers, std::vector<double>& bc_values);

  void SymbolicAssembleGlobalMatrices(const Epetra_Map& super_map);
  void AssembleGlobalMatrices();
  void ComputeSchurComplement(std::vector<int>& bc_markers, std::vector<double>& bc_values);

  double ComputeResidual(const Epetra_Vector& solution, Epetra_Vector& residual);
  double ComputeNegativeResidual(const Epetra_Vector& solution, Epetra_Vector& residual);

  void DeriveDarcyMassFlux(const Epetra_Vector& solution, 
                           const Epetra_Import& face_importer, 
                           Epetra_Vector& darcy_mass_flux);

  void InitPreconditioner(int method, Teuchos::ParameterList& prec_list);
  void UpdatePreconditioner();

  // required methods
  int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;
  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;
  bool UseTranspose() const { return false; }
  int SetUseTranspose(bool) { return 1; }

  const Epetra_Comm& Comm() const { return *(mesh_->get_comm()); }
  const Epetra_Map& OperatorDomainMap() const { return map; }
  const Epetra_Map& OperatorRangeMap() const { return map; }

  const char* Label() const { return strdup("Matrix MFD"); }
  double NormInf() const { return 0.0; }
  bool HasNormInf() const { return false; }

  // access methods
  std::vector<Teuchos::SerialDenseMatrix<int, double> >& Aff_cells() { return Aff_cells_; }
  std::vector<double>& Acc_cells() { return Acc_cells_; }
  std::vector<Epetra_SerialDenseVector>& Ff_cells() { return Ff_cells_; }
  std::vector<double>& Fc_cells() { return Fc_cells_; }
  Teuchos::RCP<Epetra_Vector>& rhs() { return rhs_; }
  Teuchos::RCP<Epetra_Vector>& rhs_faces() { return rhs_faces_; }

  Teuchos::RCP<Epetra_FECrsMatrix>& Aff() { return Aff_; }
  Teuchos::RCP<Epetra_Vector>& Acc() { return Acc_; }
  Teuchos::RCP<Epetra_CrsMatrix>& Acf() { return Acf_; }
  Teuchos::RCP<Epetra_CrsMatrix>& Afc() { return Afc_; }

  int nokay() { return nokay_; }
  int npassed() { return npassed_; }

 private:
  Teuchos::RCP<Flow_State> FS;
  Teuchos::RCP<AmanziMesh::Mesh> mesh_;
  Epetra_Map map;

  bool flag_symmetry_;

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
  ML_Epetra::MultiLevelPreconditioner* MLprec;
  Teuchos::ParameterList ML_list;
  
  Teuchos::RCP<Ifpack_Preconditioner> ifp_prec_;
  Teuchos::ParameterList ifp_plist_;

#ifdef HAVE_HYPRE_API
  Teuchos::RCP<Ifpack_Hypre> IfpHypre_Sff_;
  double hypre_tol, hypre_strong_threshold;
  int hypre_nsmooth, hypre_ncycles;
#endif

  int nokay_, npassed_;  // performance of algorithms generating mass matrices 

 private:
  void operator=(const Matrix_MFD& matrix);
};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif
