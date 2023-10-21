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
#include "exceptions.hh"
#include "errors.hh"
#include "Function.hh"
#include "FunctionConstant.hh"
#include "FunctionFactory.hh"
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
  dispersivity_ = false;

  if (plist.isParameter("dispersivity [m]")) {
    ParseAlpha_(plist, "dispersivity [m]");
    dispersivity_ = true;
    count++;
  }
  if (plist.isParameter("dispersion coefficient [m^2 s^-1]")) {
    ParseAlpha_(plist, "dispersion coefficient [m^2 s^-1]");
    dispersivity_ = false;
    count++;
  }

  if (plist.isParameter("dispersion coefficient")) {
    ParseAlpha_(plist, "dispersion coefficient");
    dispersivity_ = false;
    count++;
  }

  if (count == 0 || plist.isParameter("alpha")) {
    ParseAlpha_(plist, "alpha");
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
void
MDM_Isotropic::ParseAlpha_(Teuchos::ParameterList& plist, const std::string& keyword)
{
  if (plist.isSublist(keyword)) {
    FunctionFactory factory;
    alpha_func_ = factory.Create(plist.sublist(keyword));

    std::vector<double> args(4, 0.0);
    alpha_ = (*alpha_func_)(args);
  } else {
    alpha_ = plist.get<double>(keyword, 0.0);
    alpha_func_ = std::make_unique<FunctionConstant>(alpha_);
  }
}


/* ******************************************************************
* Isotropic tensor of rank 1.
****************************************************************** */
WhetStone::Tensor
MDM_Isotropic::mech_dispersion(double t,
                               const AmanziGeometry::Point& xc,
                               const AmanziGeometry::Point& u,
                               int axi_symmetric,
                               double s,
                               double phi) const
{
  std::vector<double> args(1 + dim_);
  args[0] = t;
  for (int i = 0; i < dim_; ++i) args[i + 1] = xc[i];
  alpha_ = (*alpha_func_)(args);

  WhetStone::Tensor D(dim_, 1);
  D(0, 0) = dispersivity_ ? alpha_ * s * phi * norm(u) : alpha_ * s * phi;
  return D;
}

} // namespace Transport
} // namespace Amanzi
