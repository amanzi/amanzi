/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Transport PK

*/

#include <cmath>
#include <string>

// Amanzi
#include "Point.hh"
#include "Tensor.hh"

#include "MDM_LichtnerKelkarRobinson.hh"

namespace Amanzi {
namespace Transport {

/* ******************************************************************
* Setup fundamental parameters for this model.
****************************************************************** */
MDM_LichtnerKelkarRobinson::MDM_LichtnerKelkarRobinson(Teuchos::ParameterList& plist)
{
  alphaLV_ = plist.get<double>("alpha_lv", 0.0);
  alphaLH_ = plist.get<double>("alpha_lh", 0.0);

  alphaTV_ = plist.get<double>("alpha_tv", 0.0);
  alphaTH_ = plist.get<double>("alpha_th", 0.0);
}


/* ******************************************************************
* Anisotropic tensor defined by four parameters.
****************************************************************** */
WhetStone::Tensor
MDM_LichtnerKelkarRobinson::mech_dispersion(double t,
                                            const AmanziGeometry::Point& xc,
                                            const AmanziGeometry::Point& u,
                                            int axi_symmetry,
                                            double wc,
                                            double phi) const
{
  AMANZI_ASSERT(axi_symmetry >= 0 && axi_symmetry < 3);
  WhetStone::Tensor D(dim_, 2);
  D.PutScalar(0.0);

  // calculate pore velocity
  AmanziGeometry::Point v(u), omega(dim_);
  v /= phi;
  double vnorm = norm(v);

  if (vnorm != 0.0) {
    double a1, a2, a3;
    double theta = v[axi_symmetry] / vnorm; // cosine of angle theta
    double theta2 = theta * theta;

    // define direction orthogonal to symmetry axis
    omega = v * (-theta / vnorm);
    omega[axi_symmetry] += 1.0;

    // we use formula (46) of Lichtner, Water Res. Research, 38 (2002)
    double alphaL = alphaLH_ + theta2 * (alphaLV_ - alphaLH_);
    a1 = alphaTH_ * vnorm;
    a2 = (alphaL - alphaTH_) / vnorm;
    a3 = (alphaTV_ - alphaTH_) * vnorm;

    for (int i = 0; i < dim_; i++) {
      D(i, i) = a1;
      for (int j = i; j < dim_; j++) {
        D(i, j) += a2 * v[i] * v[j] + a3 * omega[i] * omega[j];
        D(j, i) = D(i, j);
      }
    }

    D *= wc;
  }
  return D;
}

} // namespace Transport
} // namespace Amanzi
