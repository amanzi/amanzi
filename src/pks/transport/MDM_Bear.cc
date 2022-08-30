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

#include "MDM_Bear.hh"

namespace Amanzi {
namespace Transport {

/* ******************************************************************
* Setup fundamental parameters for this model.                                            
****************************************************************** */
MDM_Bear::MDM_Bear(Teuchos::ParameterList& plist)
{
  alphaL_ = plist.get<double>("alpha_l");
  alphaT_ = plist.get<double>("alpha_t");
}


/* ******************************************************************
* Anisotropic tensor defined by two parameters.
****************************************************************** */
WhetStone::Tensor MDM_Bear::mech_dispersion(
    const AmanziGeometry::Point& u, int axi_symmetry, double wc, double phi) const
{
  WhetStone::Tensor D(dim_, 2);
  D.PutScalar(0.0);

  // calculate pore velocity
  AmanziGeometry::Point v(u);
  v /= phi;
  double vnorm = norm(v);

  if (vnorm != 0.0) {
    double anisotropy = (alphaL_ - alphaT_) / vnorm;
    for (int i = 0; i < dim_; i++) {
      D(i, i) = alphaT_ * vnorm;
      for (int j = i; j < dim_; j++) {
        D(j, i) = D(i, j) += anisotropy * v[i] * v[j];
      }
    }

    D *= wc;
  }

  return D;
}

}  // namespace Transport
}  // namespace Amanzi

