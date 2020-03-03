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

#include "MDM_BurnettFrind.hh"

namespace Amanzi {
namespace Transport {

/* ******************************************************************
* Setup fundamental parameters for this model.                                            
****************************************************************** */
MDM_BurnettFrind::MDM_BurnettFrind(Teuchos::ParameterList& plist)
{
  alphaL_ = plist.get<double>("alphaL", 0.0);
  alphaTV_ = plist.get<double>("alphaTV", 0.0);
  alphaTH_ = plist.get<double>("alphaTH", 0.0);
}


/* ******************************************************************
* Anisotropic tensor defined by three parameters.
****************************************************************** */
WhetStone::Tensor MDM_BurnettFrind::mech_dispersion(
    const AmanziGeometry::Point& u, int axi_symmetry, double s, double phi) const
{
  WhetStone::Tensor D(dim_, 2);
  D.PutScalar(0.0);

  // calculate pore velocity
  AmanziGeometry::Point v(u);
  v /= phi;
  double vnorm = norm(v);
  if (vnorm != 0.0) {
    double a1, a2, a3;
    a1 = alphaTV_ * vnorm;
    a2 = (alphaL_ - alphaTV_) / vnorm;
    a3 = (alphaTH_ - alphaTV_) / vnorm;

    for (int i = 0; i < dim_; i++) {
      D(i, i) = a1;
      for (int j = i; j < dim_; j++) {
        D(i, j) += a2 * v[i] * v[j];
        D(j, i) = D(i, j);
      }
    }
    D(0, 0) += a3 * v[1] * v[1];
    D(1, 1) += a3 * v[0] * v[0];

    D(0, 1) -= a3 * v[0] * v[1];
    D(1, 0) -= a3 * v[0] * v[1];

    D *= phi * s;
  }

  return D;
}

}  // namespace Transport
}  // namespace Amanzi

