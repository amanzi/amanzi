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
#include "Teuchos_ParameterList.hpp"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_Vector.h"

#include "Preconditioner.hh"
#include "Transport_State.hh"

namespace Amanzi {
namespace AmanziTransport {

class DispersionModel {
 public:
  DispersionModel() {
    model = TRANSPORT_DISPERSIVITY_MODEL_NULL;
    alphaL = 0.0;
    alphaT = 0.0;
    D = 0.0;
    tau = 0.0;
  }
  ~DispersionModel() {};

 public:
  int model;
  double alphaL, alphaT, D, tau;
  std::vector<std::string> regions;
};


class Matrix_Dispersion {
 public:
  Matrix_Dispersion() {};
  Matrix_Dispersion(Teuchos::RCP<const AmanziMesh::Mesh> mesh) : mesh_(mesh) {};
  ~Matrix_Dispersion() {};

  // primary members
  void Init(std::vector<Teuchos::RCP<DispersionModel> >& specs, 
            const std::string& preconditioner, 
            const Teuchos::ParameterList& prec_list);
  void Apply(const Epetra_Vector& v,  Epetra_Vector& av) const;
  void ApplyInverse(const Epetra_Vector& v,  Epetra_Vector& hv) const;

  void CalculateDispersionTensor(const Epetra_Vector& darcy_flux,
                                 const Epetra_Vector& porosity,
                                  const Epetra_Vector& saturation);
  void SymbolicAssembleGlobalMatrix();
  void AssembleGlobalMatrixTPFA(const Teuchos::RCP<Transport_State>& TS);
  void AssembleGlobalMatrixNLFV(const Teuchos::RCP<Transport_State>& TS);
  void AddTimeDerivative(double dT, const Epetra_Vector& porosity, 
                         const Epetra_Vector& saturation);

  void UpdatePreconditioner() { preconditioner_->Update(App_); }

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

  std::vector<Teuchos::RCP<DispersionModel> >* specs_;

  std::vector<AmanziGeometry::Point> hap_points_;
  std::vector<double> hap_weights_;

  std::vector<WhetStone::Tensor> D;

  Teuchos::RCP<Epetra_FECrsMatrix> App_;
  Teuchos::RCP<AmanziPreconditioners::Preconditioner> preconditioner_;
};

}  // namespace AmanziTransport
}  // namespace Amanzi

#endif

