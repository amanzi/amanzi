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

#ifndef AMANZI_LOOKUP_TABLE_FEHM_HH_
#define AMANZI_LOOKUP_TABLE_FEHM_HH_

#include "Teuchos_ParameterList.hpp"

#include "LookupTable.hh"

namespace Amanzi {
namespace AmanziEOS {

// Equation of State model
class LookupTable_FEHM : public LookupTable {
 public:
  LookupTable_FEHM(Teuchos::ParameterList& plist);
  ~LookupTable_FEHM(){};

 private:
  std::string ReadBlock_(std::ifstream& ifs, int nP, int nT);

 private:
  double M_;
};

} // namespace AmanziEOS
} // namespace Amanzi

#endif
