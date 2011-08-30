/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_ACTIVITY_MODEL_DEBYE_HUCKEL_HH_
#define AMANZI_CHEMISTRY_ACTIVITY_MODEL_DEBYE_HUCKEL_HH_

/* 

   Class for activity calculations based on the Debye-Huckel B-dot equation.
   
   TODO(bandre): need to fix the name of this class to be
   DebyeHuckelBdot or something to distinguish it from a pure
   Debye-Huckel. Is it worth worrying about code reuse between
   debye-huckel and debye-huckel b-dot?

 */

#include "activity_model.hh"

namespace amanzi {
namespace chemistry {

class Species;

class ActivityModelDebyeHuckel : public ActivityModel {
 public:
  ActivityModelDebyeHuckel();
  ~ActivityModelDebyeHuckel();

  double Evaluate(const Species& species);

  void EvaluateVector (std::vector<double>& gamma, double& actw, const std::vector<Species>& prim, const std::vector<AqueousEquilibriumComplex>& sec);

  void Display(void) const;

 protected:

 private:
  static const double debyeA;
  static const double debyeB;
  static const double debyeBdot;
};
}  // namespace chemistry
}  // namespace amanzi
#endif  // AMANZI_CHEMISTRY_ACTIVITY_MODEL_DEBYE_HUCKEL_HH_
