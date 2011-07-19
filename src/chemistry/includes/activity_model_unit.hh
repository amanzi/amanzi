/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_ACTIVITY_MODEL_UNIT_HH_
#define AMANZI_CHEMISTRY_ACTIVITY_MODEL_UNIT_HH_

#include "activity_model.hh"

namespace amanzi {
namespace chemistry {

class Species;

class ActivityModelUnit : public ActivityModel {
 public:
  ActivityModelUnit();
  ~ActivityModelUnit();

  double Evaluate(const Species& species);

  void EvaluateVector (std::vector<double>& gamma, const std::vector<Species>* prim, const std::vector<AqueousEquilibriumComplex>* sec);

  void Display(void) const;

 protected:

 private:
};

}  // namespace chemistry
}  // namespace amanzi

#endif  // AMANZI_CHEMISTRY_ACTIVITY_MODEL_UNIT_HH_
