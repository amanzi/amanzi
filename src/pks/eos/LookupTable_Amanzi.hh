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

#ifndef AMANZI_LOOKUP_TABLE_AMANZI_HH_
#define AMANZI_LOOKUP_TABLE_AMANZI_HH_

#include "Teuchos_ParameterList.hpp"

#include "LookupTable.hh"

namespace Amanzi {
namespace AmanziEOS {

// Equation of State model
class LookupTable_Amanzi : public LookupTable {
 public:
  LookupTable_Amanzi(Teuchos::ParameterList& plist);
  ~LookupTable_Amanzi(){};

 private:
  void
  ReadMetaData_(std::ifstream& ifs, const std::string& label, int* n, double* scale, double* shift);

  void ReadBlock_(std::ifstream& ifs,
                  const std::string& field,
                  int nP,
                  int nT,
                  double scale,
                  double shift);
};

} // namespace AmanziEOS
} // namespace Amanzi

#endif
