/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_ACTIVITY_MODEL_FACTORY_HH_
#define AMANZI_CHEMISTRY_ACTIVITY_MODEL_FACTORY_HH_

#include <string>

#include "activity_model.hh"

namespace amanzi {
namespace chemistry {

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

}  // namespace chemistry
}  // namespace amanzi
#endif  // AMANZI_CHEMISTRY_ACTIVITY_MODEL_FACTORY_HH_
