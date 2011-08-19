/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_GLOBAL_IMPLICIT_EVALUATOR_HH_
#define AMANZI_CHEMISTRY_GLOBAL_IMPLICIT_EVALUATOR_HH_

//
// Class responsible for handling global implicit setup and solve.
//
// How does coupling to the transport solver work?
//

#include <vector>
#include <ostream>


namespace amanzi {
namespace chemistry {

class GlobalImplicitEvaluator : public Evaluator {
 public:
  GlobalImplicitEvaluator();
  virtual ~GlobalImplicitEvaluator();

  virtual void Display(std::ostream& output) const;

 protected:
 private:

};

}  // namespace chemistry
}  // namespace amanzi
#endif  // AMANZI_CHEMISTRY_GLOBAL_IMPLICIT_EVALUATOR_HH_
