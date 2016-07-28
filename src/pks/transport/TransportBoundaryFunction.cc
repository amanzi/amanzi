/*
  Transport PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <vector>
#include <string>

#include "Teuchos_RCP.hpp"

#include "TransportBoundaryFunction.hh"

namespace Amanzi {
namespace Transport {

/* ******************************************************************
* Constructor.
****************************************************************** */
TransportBoundaryFunction::TransportBoundaryFunction(const Teuchos::ParameterList& plist)
    : domain_volume_(-1.0)
{
  mass_ratio_ = false;
  if (plist.isSublist("boundary mass ratio")) {
    mass_ratio_ = true;
    molar_mass_ = plist.get<double>("molar mass");
  }
}


/* ******************************************************************
* Process additional parameters for BC submodels. 
****************************************************************** */
void TransportBoundaryFunction::ComputeSubmodel(
    const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
    Teuchos::RCP<const Epetra_MultiVector>& field)
{
  if (mass_ratio_) {
    double area(1.0);
    bool flag = (name() != "simple");

    for (Iterator it = begin(); it != end(); ++it) {
      int f = it->first;
      if (flag) area = mesh->face_area(f);
      it->second *= (area / molar_mass_) * fabs((*field)[0][f]);
    }
  }
}

}  // namespace Transport
}  // namespace Amanzi

