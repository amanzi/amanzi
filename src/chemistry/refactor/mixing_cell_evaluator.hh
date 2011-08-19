/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_MIXING_CELL_EVALUATOR_HH_
#define AMANZI_CHEMISTRY_MIXING_CELL_EVALUATOR_HH_

//
// Driver class for batch chemistry calculations
//
// should we just call this BatchEvaluator?
//

#include <vector>
#include <ostream>


namespace amanzi {
namespace chemistry {

class MixingCellEvaluator : public Evaluator {
 public:
  MixingCellEvaluator();
  virtual ~MixingCellEvaluator();

  virtual void Display(std::ostream& output) const;

 protected:
 private:

};

}  // namespace chemistry
}  // namespace amanzi
#endif  // AMANZI_CHEMISTRY_MIXING_CELL_EVALUATOR_HH_
