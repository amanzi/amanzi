/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_ACTIVITY_MODEL_DEBYE_HUCKEL_HH_
#define AMANZI_CHEMISTRY_ACTIVITY_MODEL_DEBYE_HUCKEL_HH_

// Base class for activity calculations

#include "activity_model.hh"

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

#endif  // AMANZI_CHEMISTRY_ACTIVITY_MODEL_DEBYE_HUCKEL_HH_

