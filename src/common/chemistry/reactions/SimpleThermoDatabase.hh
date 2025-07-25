/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*!

The thermodynamic database has a few sections, some of them are optional.

*/

#ifndef AMANZI_CHEMISTRY_SIMPLETHERMODATABASE_HH_
#define AMANZI_CHEMISTRY_SIMPLETHERMODATABASE_HH_

#include <vector>
#include <string>

#include "Teuchos_ParameterList.hpp"

#include "Beaker.hh"
#include "Species.hh"

namespace Amanzi {
namespace AmanziChemistry {

class SimpleThermoDatabase : public Beaker {
 public:
  SimpleThermoDatabase(Teuchos::RCP<Teuchos::ParameterList> plist, Teuchos::RCP<VerboseObject> vo);
  virtual ~SimpleThermoDatabase() {};

  virtual void Initialize(BeakerState& state, const BeakerParameters& parameters);

 private:
  Teuchos::ParameterList RebuildAqueousComplexes_(const Teuchos::ParameterList& aqlist);

 private:
  Teuchos::RCP<Teuchos::ParameterList> plist_;
};

} // namespace AmanziChemistry
} // namespace Amanzi

#endif
