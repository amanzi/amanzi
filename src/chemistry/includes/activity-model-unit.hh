/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef __ACTIVITY_MODEL_UNIT_HH__
#define __ACTIVITY_MODEL_UNIT_HH__

#include "activity-model.hh"

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

#endif  // __ACTIVITY_MODEL_UNIT_HPP__

