/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_MONOD_REACTION_RATE_HH_
#define AMANZI_CHEMISTRY_MONOD_REACTION_RATE_HH_

//
// Implementation of Monod reaction Rate expressions.
//

#include <vector>
#include <ostream>


namespace amanzi {
namespace chemistry {

class MonodReactionRate : public ReactionRate {
 public:
  MonodReactionRate();
  virtual ~MonodReactionRate();

  virtual void Display(std::ostream& output) const;

 protected:
 private:

};

}  // namespace chemistry
}  // namespace amanzi
#endif  // AMANZI_CHEMISTRY_MONOD_REACTIONRATE_HH_
