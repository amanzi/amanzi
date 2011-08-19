/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_REACTION_RATE_HH_
#define AMANZI_CHEMISTRY_REACTION_RATE_HH_

//
// Base class for reaction rate expressions.
//

#include <vector>
#include <ostream>


namespace amanzi {
namespace chemistry {

class ReactionRate {
 public:
  ReactionRate();
  virtual ~ReactionRate();

  virtual void Evaluate(void) = 0;
  virtual void Display(std::ostream& output) const;

 protected:
 private:

};

}  // namespace chemistry
}  // namespace amanzi
#endif  // AMANZI_CHEMISTRY_REACTION_RATE_HH_
