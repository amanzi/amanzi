/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_TST_REACTION_RATE_HH_
#define AMANZI_CHEMISTRY_TST_REACTION_RATE_HH_

//
// Implementation of TST reaction rate expressions.
//

#include <vector>
#include <ostream>


namespace amanzi {
namespace chemistry {

class TSTReactionRate : public ReactionRate {
 public:
  TSTReactionRate();
  virtual ~TSTReactionRate();

  virtual void Display(std::ostream& output) const;

 protected:
 private:

};

}  // namespace chemistry
}  // namespace amanzi
#endif  // AMANZI_CHEMISTRY_TST_REACTION_RATE_HH_
