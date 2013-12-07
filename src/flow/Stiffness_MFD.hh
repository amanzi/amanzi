/*
  This is the flow component of the Amanzi code. 

  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_STIFFNESS_MFD_HH_
#define AMANZI_STIFFNESS_MFD_HH_

#include <strings.h>

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
#include "flow_boundary_function.hh"

#include "Flow_typedefs.hh"


namespace Amanzi {
namespace AmanziFlow {

class Stiffness_MFD {
 public:
  Stiffness_MFD() {};
  Stiffness_MFD(Teuchos::RCP<const AmanziMesh::Mesh> mesh) : mesh_(mesh) {};
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
  const Epetra_BlockMap& DomainMap() const {
    return *map_;
  }
  const Epetra_BlockMap& RangeMap() const {
    return *map_;
  }

  // access methods
  std::vector<WhetStone::DenseMatrix>& Avv_cells() { return Avv_cells_; }
  std::vector<Epetra_SerialDenseVector>& Fv_cells() { return Fv_cells_; }

  Teuchos::RCP<Epetra_FECrsMatrix>& Avv() { return Avv_; }
  Teuchos::RCP<Epetra_Vector>& rhs() { return rhs_; }

#ifdef HAVE_HYPRE
  Teuchos::RCP<Ifpack_Hypre> IfpHypre_Sff() { return IfpHypre_Sff_; }
#endif

 private:
  void CombineGhostNode2MasterNode_(Epetra_Vector& v, Epetra_CombineMode mode = Insert);

 protected:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::RCP<const Epetra_Map> map_;

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
