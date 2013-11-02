/*
  This is the transport component of the Amanzi code. 

  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
  Usage: 
*/

#ifndef AMANZI_DISPERSION_NLFV_HH_
#define AMANZI_DISPERSION_NLFV_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_Vector.h"

#include "Preconditioner.hh"
#include "Dispersion.hh"

namespace Amanzi {
namespace AmanziTransport {

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


class Dispersion_NLFV : public Dispersion {
 public:
  Dispersion_NLFV() {};
  Dispersion_NLFV(std::vector<Teuchos::RCP<DispersionModel> >* specs,
                  Teuchos::RCP<const AmanziMesh::Mesh> mesh, Teuchos::RCP<State> S)
      : Dispersion(specs, mesh, S) {};
  ~Dispersion_NLFV() {};

  // primary members
  void Apply(const Epetra_Vector& v,  Epetra_Vector& av) const;
  int ApplyInverse(const Epetra_Vector& v,  Epetra_Vector& hv) const;

  void SymbolicAssembleMatrix();
  void ModifySymbolicAssemble();
  void AssembleMatrix(const Epetra_MultiVector& p);

  // additional members
  void InitNLFV();  // additional initialization of nonlinear scheme
  void CreateFluxStencils();

 private:
  void PopulateHarmonicPoints_();

 private:
  std::vector<FaceStencil> stencil_;
};

}  // namespace AmanziTransport
}  // namespace Amanzi

#endif

