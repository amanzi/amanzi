/*
This is the flow component of the Amanzi code. 

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided in the top-level COPYRIGHT file.

Authors: Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
*/

#ifndef AMANZI_STIFFNESS_MFD_HH_
#define AMANZI_STIFFNESS_MFD_HH_

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
#include "flow-boundary-function.hh"

#include "Flow_State.hh"
#include "Flow_typedefs.hh"


namespace Amanzi {
namespace AmanziFlow {

class Stiffness_MFD : public Epetra_Operator {
 public:
   Stiffness_MFD(Teuchos::RCP<Flow_State> FS_, const Epetra_Map& map);
  ~Stiffness_MFD() {};

  // main methods
  void CreateMFDstiffnessMatrices(std::vector<WhetStone::Tensor>& K, double factor);
  void CreateMFDrhsVectors();
  void ApplyBoundaryConditions(std::vector<int>& bc_model, std::vector<double>& bc_values);

  void SymbolicAssembleGlobalMatrices();
  void AssembleGlobalMatrices();

  void InitPreconditioner(int method, Teuchos::ParameterList& prec_list);
  void UpdatePreconditioner();

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

  // access methods
  std::vector<WhetStone::DenseMatrix>& Avv_cells() { return Avv_cells_; }
  std::vector<Epetra_SerialDenseVector>& Fv_cells() { return Fv_cells_; }

  Teuchos::RCP<Epetra_FECrsMatrix>& Avv() { return Avv_; }
  Teuchos::RCP<Epetra_Vector>& rhs() { return rhs_; }

#ifdef HAVE_HYPRE
  Teuchos::RCP<Ifpack_Hypre> IfpHypre_Sff() { return IfpHypre_Sff_; }
#endif

 protected:
  Teuchos::RCP<Flow_State> FS_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Epetra_Map map_;

  std::vector<WhetStone::DenseMatrix> Avv_cells_;
  std::vector<Epetra_SerialDenseVector> Fv_cells_;

  Teuchos::RCP<Epetra_FECrsMatrix> Avv_;
  Teuchos::RCP<Epetra_Vector> rhs_;

  int method_;  // Preconditioners
  Teuchos::RCP<ML_Epetra::MultiLevelPreconditioner> MLprec;
  Teuchos::ParameterList ML_list;
  
  Teuchos::RCP<Ifpack_Preconditioner> ifp_prec_;
  Teuchos::ParameterList ifp_plist_;

#ifdef HAVE_HYPRE
  Teuchos::RCP<Ifpack_Hypre> IfpHypre_Sff_;
  double hypre_tol, hypre_strong_threshold;
  int hypre_nsmooth, hypre_ncycles, hypre_verbosity;
#endif
};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif
