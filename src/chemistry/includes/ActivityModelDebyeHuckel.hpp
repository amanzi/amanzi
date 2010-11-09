/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef __ACTIVITY_MODEL_DEBYE_HUCKEL_HPP__
#define __ACTIVITY_MODEL_DEBYE_HUCKEL_HPP__

// Base class for activity calculations

#include <string>
#include <vector>

#include "ActivityModel.hpp"

class Species;

class ActivityModelDebyeHuckel : public ActivityModel {
 public:
  ActivityModelDebyeHuckel();
  ~ActivityModelDebyeHuckel();

  double Evaluate(const Species& species);

  void Display(void) const;

 protected:

 private:
  static const double debyeA;
  static const double debyeB;
  static const double debyeBdot;
};

#endif  // __ACTIVITY_MODEL_DEBYE_HUCKEL_HPP__

