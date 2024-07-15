/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  EOS

  Tabulated equations of state.
*/

#ifndef AMANZI_LOOKUP_TABLE_HH_
#define AMANZI_LOOKUP_TABLE_HH_

#include "Teuchos_ParameterList.hpp"

namespace Amanzi {
namespace AmanziEOS {

const int EOS_TABLE_UNKNOWN = -1;
const int EOS_TABLE_LIQUID = -2;
const int EOS_TABLE_GAS = -3;
const int EOS_TABLE_SATURATION_CURVE = -4;
const double EOS_TABLE_TOL = 2e-10;

// Equation of State model
class LookupTable {
 public:
  LookupTable(Teuchos::ParameterList& plist){};
  virtual ~LookupTable(){};

  // virtual methods
  virtual double Function(double T, double p, int* ierr);
  virtual double DFunctionDT(double T, double p, int* ierr);
  virtual double DFunctionDp(double T, double p, int* ierr);
  virtual int Location(double T, double p, int* ierr) { return EOS_TABLE_UNKNOWN; }

  // error parsing
  std::string ErrorMessage(double T, double p);

 protected:
  int FindBox_(double T, double p, int* ip, int* jp);

  double DerivativeP_(int i, int j);
  double DerivativeT_(int i, int j);

 protected:
  double shiftP_, shiftT_;
  double scaleP_, scaleT_;
  double scaleF_, shiftF_;

  std::vector<double> axisT_, axisP_;
  std::vector<std::vector<double>> F_;
};

} // namespace AmanziEOS
} // namespace Amanzi

#endif
