/* -*-  mode: c++; c-default-style: "google-c-style"; indent-tabs-mode: nil -*- */
#ifndef __ACTIVITY_MODEL_UNIT_HPP__
#define __ACTIVITY_MODEL_UNIT_HPP__

#include "ActivityModel.hpp"

class Species;

class ActivityModelUnit : public ActivityModel {
 public:
  ActivityModelUnit();
  ~ActivityModelUnit();

  double evaluate(const Species& species);

  void display(void) const;

 protected:

 private:
};

#endif  // __ACTIVITY_MODEL_UNIT_HPP__

