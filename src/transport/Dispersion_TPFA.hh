/*
  This is the transport component of the Amanzi code. 

  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
  Usage: 
*/

#ifndef AMANZI_DISPERSION_TPFA_HH_
#define AMANZI_DISPERSION_TPFA_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_Vector.h"

#include "Preconditioner.hh"
#include "Transport_State.hh"
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


class Dispersion_TPFA : public Dispersion {
 public:
  Dispersion_TPFA() {};
  Dispersion_TPFA(std::vector<Teuchos::RCP<DispersionModel> >* specs,
                  Teuchos::RCP<const AmanziMesh::Mesh> mesh, Teuchos::RCP<Transport_State> TS)
      : Dispersion(specs, mesh, TS) {};
  ~Dispersion_TPFA() {};

  // primary members
  void Apply(const Epetra_Vector& v,  Epetra_Vector& av) const;
  int ApplyInverse(const Epetra_Vector& v,  Epetra_Vector& hv) const;

  void SymbolicAssembleGlobalMatrix();
  void AssembleGlobalMatrix(const Epetra_Vector& p) { AssembleGlobalMatrixTPFA(TS_); }
  void AssembleGlobalMatrixTPFA(const Teuchos::RCP<Transport_State>& TS);

  // NLFV members
  void InitNLFV();  // additional initialization of nonlinear scheme
  void CreateFluxStencils();
  void AssembleGlobalMatrixNLFV(const Epetra_Vector& p);

 private:
  void PopulateHarmonicPoints_();

 private:
  std::vector<FaceStencil> stencil_;
};

}  // namespace AmanziTransport
}  // namespace Amanzi

#endif

