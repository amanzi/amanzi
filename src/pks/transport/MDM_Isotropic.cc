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
#include "exceptions.hh"
#include "errors.hh"
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
  int count = 0;
  alpha_ = 0.;
  dispersivity_ = false;

  if (plist.isParameter("dispersivity [m]")) {
    alpha_ = plist.get<double>("dispersivity [m]");
    dispersivity_ = true;
    count++;
  }
  if (plist.isParameter("dispersion coefficient [m^2 s^-1]")) {
    alpha_ = plist.get<double>("dispersion coefficient [m^2 s^-1]");
    dispersivity_ = false;
    count++;
  }

  if (plist.isParameter("dispersion coefficient")) {
    alpha_ = plist.get<double>("dispersion coefficient");
    dispersivity_ = false;
    count++;
  }
  if (count == 0 || plist.isParameter("alpha")) {
    alpha_ = plist.get<double>("alpha", 0.0);
    dispersivity_ = true;
    count++;
  }

  if (count > 1) {
    Errors::Message msg("Only one of \"alpha\" and \"dispersion coefficient\" may be supplied.");
    Exceptions::amanzi_throw(msg);
  }
  if (count == 0) {
    Errors::Message msg("Either \"alpha\" or \"dispersion coefficient\" must be supplied.");
    Exceptions::amanzi_throw(msg);
  }
}


/* ******************************************************************
* Isotropic tensor of rank 1.
****************************************************************** */
WhetStone::Tensor
MDM_Isotropic::mech_dispersion(const AmanziGeometry::Point& u,
                               int axi_symmetric,
                               double s,
                               double phi) const
{
  WhetStone::Tensor D(dim_, 1);
  D(0, 0) = dispersivity_ ? alpha_ * s * phi * norm(u) : alpha_ * s * phi;
  return D;
}

} // namespace Transport
} // namespace Amanzi
