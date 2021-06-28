/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ben Andre
*/

#ifndef AMANZI_CHEMISTRY_ACTIVITY_MODEL_FACTORY_HH_
#define AMANZI_CHEMISTRY_ACTIVITY_MODEL_FACTORY_HH_

#include <memory>
#include <string>
#include <vector>

#include "ActivityModel.hh"

namespace Amanzi {
namespace AmanziChemistry {

class ActivityModelFactory {
 public:
  ActivityModelFactory() {};
  ~ActivityModelFactory() {};

  std::shared_ptr<ActivityModel> Create(
      const std::string& model,
      const ActivityModel::ActivityModelParameters& parameters,
      const std::vector<Species>& primary_species,
      const std::vector<AqueousEquilibriumComplex>& secondary_species,
      const Teuchos::Ptr<VerboseObject> vo);
  static const std::string debye_huckel;
  static const std::string pitzer_hwm;
  static const std::string unit;
};

}  // namespace AmanziChemistry
}  // namespace Amanzi
#endif
