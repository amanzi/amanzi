/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_ACTIVITY_MODEL_FACTORY_HH_
#define AMANZI_CHEMISTRY_ACTIVITY_MODEL_FACTORY_HH_

#include <string>

#include "activity_model.hh"

class ActivityModelFactory {
 public:
  ActivityModelFactory();
  ~ActivityModelFactory();

  ActivityModel* Create(const std::string& model);

  static const std::string debye_huckel;
  static const std::string unit;

 protected:

 private:
};

#endif  // AMANZI_CHEMISTRY_ACTIVITY_MODEL_FACTORY_HH_
