/*
  This is the transport component of the Amanzi code. 

  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
  Usage: 
*/

#ifndef AMANZI_DISPERSION_HH_
#define AMANZI_DISPERSION_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_Vector.h"

#include "tensor.hh"
#include "State.hh"
#include "Preconditioner.hh"
#include "TransportDefs.hh"

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


class Dispersion {
 public:
  Dispersion() {};
  Dispersion(std::vector<Teuchos::RCP<DispersionModel> >* specs,
             Teuchos::RCP<const AmanziMesh::Mesh> mesh, Teuchos::RCP<State> S) 
      : specs_(specs), mesh_(mesh), S_(S) {};
  ~Dispersion() {};

  // required members
  virtual void SymbolicAssembleMatrix() {};  // It fixes a large stencil S.
  virtual void ModifySymbolicAssemble() {};  // It allows to tweak the stencil a little.
  virtual void AssembleMatrix(const Epetra_MultiVector& p) {};

  virtual int Apply(const Epetra_Vector& v,  Epetra_Vector& av) const { return 0; };
  virtual int ApplyInverse(const Epetra_Vector& v,  Epetra_Vector& hv) const { return 0; }

  // generic members
  void Init();
  void CalculateDispersionTensor(
      const Epetra_MultiVector& darcy_flux,
      const Epetra_MultiVector& porosity, const Epetra_MultiVector& saturation);

  void AddTimeDerivative(
      double dT, const Epetra_MultiVector& porosity, const Epetra_MultiVector& saturation);

  const Epetra_Map& RangeMap() const { return App_->RangeMap(); }
  const Epetra_Map& DomainMap() const { return App_->DomainMap(); }

  void InitPreconditioner(const std::string& prec_name, const Teuchos::ParameterList& prec_list);
  void UpdatePreconditioner() { preconditioner_->Update(App_); }

 public:
  int nprimary, nsecondary;  // number of matrices
  int num_simplex_itrs;

 protected:
  std::vector<Teuchos::RCP<DispersionModel> >* specs_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::RCP<State> S_;

  int dim;
  int ncells_owned, ncells_wghost;
  int nfaces_owned, nfaces_wghost;

  std::vector<WhetStone::Tensor> D;
  Teuchos::RCP<Epetra_FECrsMatrix> App_;
  Teuchos::RCP<AmanziPreconditioners::Preconditioner> preconditioner_;
};

}  // namespace AmanziTransport
}  // namespace Amanzi

#endif

