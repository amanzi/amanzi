/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_REACTION_NETWORK_HH_
#define AMANZI_CHEMISTRY_REACTION_NETWORK_HH_

//
// Class receives the function evaluation vector and jacobian matrix,
// and populates them based on the reactions specified. The solve is
// actually handled by a solver object living inside the
// evaluator. This class will interact with a database class and be
// responsible for basis swapping?
//

#include <vector>
#include <ostream>


namespace amanzi {
namespace chemistry {

class ReactionNetwork {
 public:
  ReactionNetwork();
  virtual ~ReactionNetwork();

  virtual void Display(std::ostream& output) const;

 protected:
 private:

};

}  // namespace chemistry
}  // namespace amanzi
#endif  // AMANZI_CHEMISTRY_REACTION_NETWORK_HH_
