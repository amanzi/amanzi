/*
  This is the flow component of the Amanzi code. 

  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
*/

#ifndef AMANZI_MATRIX_MFD_HH_
#define AMANZI_MATRIX_MFD_HH_

#include <strings.h>

#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_FECrsMatrix.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_LAPACK.hpp"

#include "Preconditioner.hh"
#include "CompositeVector.hh"
#include "CompositeVectorSpace.hh"

#include "Mesh.hh"
#include "Point.hh"
#include "DenseMatrix.hh"
#include "boundary_function.hh"

#include "State.hh"
#include "Matrix.hh"
#include "RelativePermeability.hh"

#include "FlowTypeDefs.hh"

namespace Amanzi {
namespace AmanziFlow {

class Matrix_MFD : public Matrix<CompositeVector, CompositeVectorSpace> {
 public:
  Matrix_MFD() {};
  Matrix_MFD(Teuchos::RCP<State> S, Teuchos::RCP<RelativePermeability> rel_perm);
  ~Matrix_MFD();

  // main methods (required methods)
  void Init() {};
  void CreateMassMatrices(int mfd3d_method, std::vector<WhetStone::Tensor>& K);
  void CreateRHSVectors();

  void CreateStiffnessMatricesDarcy(int mfd3d_method, std::vector<WhetStone::Tensor>& K);
  void CreateStiffnessMatricesRichards();

  void SymbolicAssemble();
  void Assemble();
  void AssembleSchurComplement(std::vector<int>& bc_model, std::vector<bc_tuple>& bc_values);

  void AddGravityFluxesDarcy(double rho, const AmanziGeometry::Point& gravity,
                             std::vector<WhetStone::Tensor>& K);
  void AddGravityFluxesRichards(double rho, const AmanziGeometry::Point& gravity,
                                std::vector<WhetStone::Tensor>& K);

  void AddTimeDerivative(
      const Epetra_MultiVector& p, const Epetra_MultiVector& phi, double rho, double dT);
  void AddTimeDerivativeSpecificStorage(
      const Epetra_MultiVector& p, const Epetra_MultiVector& ss, double g, double dT);
  void AddTimeDerivativeSpecificYield(
      const Epetra_MultiVector& p, const Epetra_MultiVector& sy, double g, double dT);

  void ApplyBoundaryConditions(std::vector<int>& bc_model, std::vector<bc_tuple>& bc_values);

  void InitPreconditioner(const std::string& name, const Teuchos::ParameterList& plist);
  virtual void UpdatePreconditioner() { preconditioner_->Update(Sff_); } 
  void DestroyPreconditioner() { preconditioner_->Destroy(); }

  void DeriveMassFlux(const CompositeVector& solution, CompositeVector& darcy_mass_flux,
                      std::vector<int>& bc_model, std::vector<bc_tuple>& bc_values);

  virtual int Apply(const CompositeVector& X, CompositeVector& Y) const;
  virtual int ApplyInverse(const CompositeVector& X, CompositeVector& Y) const;
  const CompositeVectorSpace& DomainMap() const {
    return cvs_;
  }
  const CompositeVectorSpace& RangeMap() const {
    return cvs_;
  }

  // control methods
  bool CheckActionProperty(int action) { return actions_ & action; }

  // other main methods
  int ReduceGlobalSystem2LambdaSystem(CompositeVector& u);

  // development methods
  void CreateMassMatrices_ScaledStability(int method, double factor, std::vector<WhetStone::Tensor>& K);
  int PopulatePreconditioner(Matrix_MFD& matrix);
  void RescaleMFDstiffnessMatrices(const Epetra_Vector& old_scale, const Epetra_Vector& new_scale);

  // access methods
  Teuchos::RCP<Epetra_FECrsMatrix>& Aff() { return Aff_; }
  Teuchos::RCP<Epetra_FECrsMatrix>& Sff() { return Sff_; }
  Teuchos::RCP<Epetra_Vector>& Acc() { return Acc_; }
  Teuchos::RCP<Epetra_CrsMatrix>& Acf() { return Acf_; }
  Teuchos::RCP<Epetra_CrsMatrix>& Afc() { return Afc_; }

 protected:
  CompositeVectorSpace cvs_;

  std::vector<WhetStone::DenseMatrix> Mff_cells_;
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
  Teuchos::RCP<AmanziPreconditioners::Preconditioner> preconditioner_;

 private:
  void operator=(const Matrix_MFD& matrix);
};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif
