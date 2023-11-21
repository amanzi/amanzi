/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_DRIVER_VERIFY_DRIVER_HH
#define AMANZI_DRIVER_VERIFY_DRIVER_HH

#include "CompositeVector.hh"
#include "State.hh"

void VerifyDriverAperture(const Amanzi::State& S, double tol)
{
  double t = 1000.0;
  const auto& a_c = *S.Get<Amanzi::CompositeVector>("fracture-aperture").ViewComponent("cell");

  auto mesh = S.GetMesh("fracture");
  for (int c = 0; c < a_c.MyLength(); ++c) {
    const auto& xc = mesh->getCellCentroid(c);
    double w = (xc[0] + 12.5) * std::sqrt(2000.0 / (t + 1e-10)) - 12.5;
    double tmp = (w > 7.4) ? 1.0e-5 :
                             1.257136447580295e-05 - 5.307442110525289e-07 * w +
                               1.809065105747629e-08 * w * w + 1.184392778715966e-09 * w * w * w -
                               3.620526561612838e-11 * w * w * w * w;
    if (c < 10 && mesh->getComm()->MyPID() == 0)
      std::cout << "aperture: " << a_c[0][c] << ",  exact: " << tmp << std::endl;
    CHECK(std::fabs(a_c[0][c] - tmp) < tol * tmp);
  }
}

#endif

