/*
  This is the transport component of the Amanzi code. 

  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
  Usage: 
*/

#ifndef AMANZI_DISPERSION_MFD_HH_
#define AMANZI_DISPERSION_MFD_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_Map.h"

#include "Preconditioner.hh"
#include "Dispersion.hh"

namespace Amanzi {
namespace AmanziTransport {

class Dispersion_MFD : public Dispersion {
 public:
  Dispersion_MFD() {};
  Dispersion_MFD(std::vector<Teuchos::RCP<DispersionModel> >* specs,
                Teuchos::RCP<const AmanziMesh::Mesh> mesh, Teuchos::RCP<State> S)
      : Dispersion(specs, mesh, S) {};
  ~Dispersion_MFD() {};

  // primary members
  void Apply(const Epetra_Vector& v,  Epetra_Vector& av) const;
  int ApplyInverse(const Epetra_Vector& v,  Epetra_Vector& hv) const;

  void SymbolicAssembleMatrix();
  void AssembleMatrix(const Epetra_Vector& p);

  // access (now only for unit tests)
  const Epetra_Map& super_map() { return *super_map_; }

 private:
  Teuchos::RCP<Epetra_Map> super_map_;
};

}  // namespace AmanziTransport
}  // namespace Amanzi

#endif

