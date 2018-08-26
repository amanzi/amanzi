/*
  Transport PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Implicit time discretization is unconditionally stable.
*/

#include <string>

#include "DenseMatrix.hh"
#include "DenseVector.hh"

#include "MultiscaleTransportPorosity_DPM.hh"

namespace Amanzi {
namespace Transport {

/* ******************************************************************
* This model is minor extension of the WRM.
****************************************************************** */
MultiscaleTransportPorosity_DPM::MultiscaleTransportPorosity_DPM(const Teuchos::ParameterList& plist)
{
  omega_ = plist.sublist("dual porosity parameters")
                .get<double>("solute transfer coefficient");
}


/* ******************************************************************
* Calculates flux from matrix to fracture and updates matrix node.
****************************************************************** */
double MultiscaleTransportPorosity_DPM::ComputeSoluteFlux(
    double flux_liquid, double& tcc_f, WhetStone::DenseVector& tcc_m, int icomp,
    double dt, double wcf0, double wcf1, double wcm0, double wcm1, double phi)
{
  double bf, bm;
  if (flux_liquid > 0.0) {
    bm = omega_ * dt;
    bf = bm + flux_liquid;
  } else {
    bf = omega_ * dt;
    bm = bf - flux_liquid;
  }

  // solve 2x2 system 
  WhetStone::DenseMatrix A(2, 2);
  WhetStone::DenseVector b(2), x(2);

  A(0, 0) = wcf1 + bf;
  A(0, 1) = -bm;
  A(1, 0) = -bf;
  A(1, 1) = wcm1 + bm;

  b(0) = wcf0 * tcc_f;
  b(1) = wcm0 * tcc_m(0);
 
  A.Inverse();
  A.Multiply(b, x, false);
 
  tcc_f = x(0);
  tcc_m(0) = x(1);

  double tmp = (flux_liquid > 0.0) ? tcc_f : tcc_m(0); 
  return flux_liquid * tmp + omega_ * (tcc_f - tcc_m(0));
}

}  // namespace Transport
}  // namespace Amanzi
  
  
