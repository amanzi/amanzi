/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_ACTIVITY_MODEL_UNIT_HH_
#define AMANZI_CHEMISTRY_ACTIVITY_MODEL_UNIT_HH_

#include "activity_model.hh"

class Species;

class ActivityModelUnit : public ActivityModel {
 public:
  ActivityModelUnit();
  ~ActivityModelUnit();

  double Evaluate(const Species& species);

  void Display(void) const;

 protected:

 private:
};

#endif  // AMANZI_CHEMISTRY_ACTIVITY_MODEL_UNIT_HH_
