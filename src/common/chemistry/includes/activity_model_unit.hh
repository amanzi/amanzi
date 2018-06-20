/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_ACTIVITY_MODEL_UNIT_HH_
#define AMANZI_CHEMISTRY_ACTIVITY_MODEL_UNIT_HH_

#include "activity_model.hh"

namespace Amanzi {
namespace AmanziChemistry {

class Species;

class ActivityModelUnit : public ActivityModel {
 public:
  ActivityModelUnit();
  ~ActivityModelUnit();

  double Evaluate(const Species& species);

  void EvaluateVector(const std::vector<Species>& prim, 
                      const std::vector<AqueousEquilibriumComplex>& sec,
                      std::vector<double>* gamma,
                      double* actw);

  void Display(void) const;

 protected:

 private:
};

}  // namespace AmanziChemistry
}  // namespace Amanzi

#endif  // AMANZI_CHEMISTRY_ACTIVITY_MODEL_UNIT_HH_
