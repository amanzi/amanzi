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
#include "Dispersion.hh"

namespace Amanzi {
namespace AmanziTransport {

class Dispersion_TPFA : public Dispersion {
 public:
  Dispersion_TPFA() {};
  Dispersion_TPFA(std::vector<Teuchos::RCP<DispersionModel> >* specs,
                  Teuchos::RCP<const AmanziMesh::Mesh> mesh, Teuchos::RCP<State> S)
      : Dispersion(specs, mesh, S) {};
  ~Dispersion_TPFA() {};

  // primary members
  int Apply(const Epetra_Vector& v,  Epetra_Vector& av) const;
  int ApplyInverse(const Epetra_Vector& v,  Epetra_Vector& hv) const;

  void SymbolicAssembleMatrix();
  void AssembleMatrix(const Epetra_MultiVector& p);
};

}  // namespace AmanziTransport
}  // namespace Amanzi

#endif

