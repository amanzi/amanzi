/*
  This is the transport component of the Amanzi code. 

  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
  Usage: 
*/

#ifndef AMANZI_DISPERSION_MODEL_HH_
#define AMANZI_DISPERSION_MODEL_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_Vector.h"

#include "tensor.hh"
#include "State.hh"
#include "Preconditioner.hh"
#include "TransportDefs.hh"

namespace Amanzi {
namespace Transport {

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

}  // namespace Transport
}  // namespace Amanzi

#endif

