/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_SURFACE_COMPLEXATION_REACTION_HH_
#define AMANZI_CHEMISTRY_SURFACE_COMPLEXATION_REACTION_HH_

//
// Class for representing surface complexation reactions
//

#include <vector>
#include <ostream>


namespace amanzi {
namespace chemistry {

class SurfaceComplexationReaction {
 public:
  SurfaceComplexationReaction();
  virtual ~SurfaceComplexationReaction();

  virtual void Display(std::ostream& output) const;

 protected:
 private:

};

}  // namespace chemistry
}  // namespace amanzi
#endif  // AMANZI_CHEMISTRY_SURFACE_COMPLEXATION_REACTION_HH_
