/*
  Transport PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <cmath>
#include <string>

// Amanzi
#include "Point.hh"
#include "Tensor.hh"

#include "MDM_Isotropic.hh"

namespace Amanzi {
namespace Transport {

/* ******************************************************************
* Setup fundamental parameters for this model.                                            
****************************************************************** */
MDM_Isotropic::MDM_Isotropic(Teuchos::ParameterList& plist)
{
  alpha_ = plist.get<double>("alpha", 0.0);
}


/* ******************************************************************
* Isotropic tensor of rank 1.
****************************************************************** */
WhetStone::Tensor MDM_Isotropic::mech_dispersion(
    const AmanziGeometry::Point& u, int axi_symmetric, double s, double phi) const
{
  WhetStone::Tensor D(dim_, 1);
  D(0, 0) = alpha_ * s * phi * norm(u);
  return D;
}

}  // namespace Transport
}  // namespace Amanzi

