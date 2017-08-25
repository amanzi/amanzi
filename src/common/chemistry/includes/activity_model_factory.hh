/* -*-  mode: c++; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_ACTIVITY_MODEL_FACTORY_HH_
#define AMANZI_CHEMISTRY_ACTIVITY_MODEL_FACTORY_HH_

#include <string>

#include "activity_model.hh"

namespace Amanzi {
namespace AmanziChemistry {

class ActivityModelFactory {
 public:
  ActivityModelFactory();
  ~ActivityModelFactory();

  ActivityModel* Create(const std::string& model,
                        const ActivityModel::ActivityModelParameters& parameters,
                        const std::vector<Species>& primary_species,
                        const std::vector<AqueousEquilibriumComplex>& secondary_species);
  static const std::string debye_huckel;
  static const std::string pitzer_hwm;
  static const std::string unit;


 protected:

 private:
};

}  // namespace AmanziChemistry
}  // namespace Amanzi
#endif  // AMANZI_CHEMISTRY_ACTIVITY_MODEL_FACTORY_HH_
