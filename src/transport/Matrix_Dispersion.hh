/*
  This is the transport component of the Amanzi code. 

  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
  Usage: 
*/

#ifndef AMANZI_MATRIX_DISPERSION_HH_
#define AMANZI_MATRIX_DISPERSION_HH_

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


struct FaceStencil {
 public:
  FaceStencil() {};
  ~FaceStencil() {};

  void Init(int d) {
    p.init(d);
    weights.resize(2 * d);
    stencil.resize(2 * d);
    faces.resize(2 * d);
  }

 public:
  AmanziGeometry::Point p;  // harmonic averaging point
  double gamma;  // coefficient for cell with the lowest id 
  std::vector<double> weights;  // weights in positive decompositions
  std::vector<int> stencil;  // ids of cells in positive decompositions
  std::vector<int> faces;  // ids of interface faces
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
  void AddTimeDerivative(double dT, const Epetra_Vector& porosity, 
                         const Epetra_Vector& saturation);

  // TPFA members
  void AssembleGlobalMatrixTPFA(const Teuchos::RCP<Transport_State>& TS);

  // NLFV members
  void InitNLFV();  // additional initialization of nonlinear scheme
  void CreateFluxStencils();
  void AssembleGlobalMatrixNLFV(const Epetra_Vector& p);

  void UpdatePreconditioner() { preconditioner_->Update(App_); }

 private:
  void PopulateHarmonicPoints_();

 private:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  int dim;

  int ncells_owned, ncells_wghost;
  int nfaces_owned, nfaces_wghost;

  std::vector<Teuchos::RCP<DispersionModel> >* specs_;

  std::vector<WhetStone::Tensor> D;
  std::vector<FaceStencil> stencil_;

  Teuchos::RCP<Epetra_FECrsMatrix> App_;
  Teuchos::RCP<AmanziPreconditioners::Preconditioner> preconditioner_;
};

}  // namespace AmanziTransport
}  // namespace Amanzi

#endif

