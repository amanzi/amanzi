/*
This is the transport component of the Amanzi code. 

Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided in the top-level COPYRIGHT file.

Author: Konstantin Lipnikov (lipnikov@lanl.gov)
Usage: 
*/

#ifndef __MATRIX_DISPERSION_HH__
#define __MATRIX_DISPERSION_HH__

#include "Teuchos_RCP.hpp"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_Vector.h"

#include "Preconditioner.hh"

namespace Amanzi {
namespace AmanziTransport {

class Dispersion_Specs {
 public:
  Dispersion_Specs() {
    model = TRANSPORT_DISPERSIVITY_MODEL_NULL;
    dispersivity_longitudinal = 0.0;
    dispersivity_transverse = 0.0;
    method = TRANSPORT_DISPERSION_METHOD_TPFA; 
    preconditioner = "identity";
  }
  ~Dispersion_Specs() {};

 public:
  int model, method;
  double dispersivity_longitudinal, dispersivity_transverse;
  string preconditioner;
};


class Matrix_Dispersion {
 public:
  Matrix_Dispersion() {};
  Matrix_Dispersion(Teuchos::RCP<const AmanziMesh::Mesh> mesh) : mesh_(mesh) {};
  ~Matrix_Dispersion() {};

  // primary members
  void Init(Dispersion_Specs& specs);
  void Apply(const Epetra_Vector& v,  Epetra_Vector& av) const;
  void ApplyInverse(const Epetra_Vector& v,  Epetra_Vector& hv) const;

  void CalculateDispersionTensor(const Epetra_Vector& darcy_flux,
                                 const Epetra_Vector& porosity,
                                 const Epetra_Vector& saturation);
  void SymbolicAssembleGlobalMatrix();
  void AssembleGlobalMatrix();
  void AddTimeDerivative(double dT, const Epetra_Vector& porosity, 
                         const Epetra_Vector& saturation);

  void UpdatePreconditioner() { preconditioner_->Update(App_); }
  void InitPreconditioner(const Teuchos::ParameterList& list) { preconditioner_->Init(list); }

 private:
  void PopulateHarmonicPoints();
  void ExtractBoundaryConditions(const int component,
                                 std::vector<int>& bc_face_id,
                                 std::vector<double>& bc_face_value);

 private:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  int dim;

  int ncells_owned, ncells_wghost;
  int nfaces_owned, nfaces_wghost;

  Dispersion_Specs* specs_;

  std::vector<AmanziGeometry::Point> hap_points_;
  std::vector<double> hap_weights_;

  std::vector<WhetStone::Tensor> D;

  Teuchos::RCP<Epetra_FECrsMatrix> App_;
  Teuchos::RCP<AmanziPreconditioners::Preconditioner> preconditioner_;
};

}  // namespace AmanziTransport
}  // namespace Amanzi

#endif

