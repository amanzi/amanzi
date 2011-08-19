/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_EVALUATOR_HH_
#define AMANZI_CHEMISTRY_EVALUATOR_HH_

//
// Base class for chemistry evaluators (batch/mixing cell, GIA). These
// classes own the function evaluation vector and jacobian. They pass
// them to the reaction network to be populated, then pass them to the
// solver object to do the actual solve. Solution info is then passed
// back to the reaction network for updating species concentrations,
// Q/K, activity coefficients, etc.
//

#include <vector>
#include <ostream>


namespace amanzi {
namespace chemistry {

class Evaluator {
 public:
  Evaluator();
  virtual ~Evaluator();

  virtual void Display(std::ostream& output) const;

 protected:
 private:

};

}  // namespace chemistry
}  // namespace amanzi
#endif  // AMANZI_CHEMISTRY_EVALUATOR_HH_
