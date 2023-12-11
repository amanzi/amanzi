/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.
*/

/*
  This is the multiphase component of the Amanzi code.

*/

#include <cmath>
#include <string>

#include "MultiphaseDefs.hh"
#include "WRMmp_Custom.hh"

namespace Amanzi {
namespace Multiphase {

/* ******************************************************************
* Setup fundamental parameters for this model.
****************************************************************** */
WRMmp_Custom::WRMmp_Custom(Teuchos::ParameterList& plist)
{
  double coef = plist.get<double>("coefficient", MULTIPHASE_WRM_EXCEPTION);

  Init_(coef);
}

void
WRMmp_Custom::Init_(double coef)
{
  coef_ = 1.0;
}


/* ******************************************************************
* Relative permeability formula.
****************************************************************** */
double
WRMmp_Custom::k_relative(double Sw, int phase)
{
  return 1.0;
}


/* ******************************************************************
* Derivative of relative permeability wrt liquid saturation.
****************************************************************** */
double
WRMmp_Custom::dKdS(double Sw, int phase)
{
  return 0.0;
}

/* ******************************************************************
* Capillary pressure formula.
****************************************************************** */
double
WRMmp_Custom::capillaryPressure(double Sw)
{
  //return Sw/(2.0 * 3.141592653589793 * 3.141592653589793);
  return 0.0;
}


/* ******************************************************************
* Derivative of capillary pressure.
****************************************************************** */
double
WRMmp_Custom::dPc_dS(double Sw)
{
 return 0.0;
 //return 1.0/(2.0 * 3.141592653589793 * 3.141592653589793);
}

} // namespace Multiphase
} // namespace Amanzi
