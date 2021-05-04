/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  Tabulated equations of state.
*/

#ifndef AMANZI_LOOKUP_TABLE_HH_
#define AMANZI_LOOKUP_TABLE_HH_

#include "Teuchos_ParameterList.hpp"

namespace Amanzi {
namespace AmanziEOS {

// Equation of State model
class LookupTable {
 public:
  LookupTable(Teuchos::ParameterList& plist);
  ~LookupTable() {};

  // Virtual methods that form the Viscosity
  double Function(double T, double p);
  double DFunctionDT(double T, double p);
  double DFunctionDp(double T, double p);

 private:
  void ReadMetaData_(std::ifstream& ifs, const std::string& label,
                     int* n, double* scale, double* shift);

  void ReadBlock_(std::ifstream& ifs, const std::string& field,
                  int nP, int nT, double scale, double shift);

  void FindBox_(double T, double p, int* ip, int* jp);

  double DerivativeP_(int i, int j);
  double DerivativeT_(int i, int j);

 private:
  double shiftP_, shiftT_;
  double scaleP_, scaleT_;
  double scaleF_, shiftF_;

  std::vector<double> axisT_, axisP_;
  std::vector<std::vector<double> > F_; 
};

}  // namespace AmanziEOS
}  // namespace Amanzi

#endif
