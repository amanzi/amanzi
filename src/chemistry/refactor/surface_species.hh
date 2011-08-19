/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_SURFACE_SPECIES_HH_
#define AMANZI_CHEMISTRY_SURFACE_SPECIES_HH_

//
// Class for representing surface species
//
// Are we considering a SurfaceSpecies to be a "surface site", "surface complex" or both?  
//

#include <vector>
#include <ostream>


namespace amanzi {
namespace chemistry {

class SurfaceSpecies {
 public:
  SurfaceSpecies();
  virtual ~SurfaceSpecies();

  virtual void Display(std::ostream& output) const;

 protected:
 private:

};

}  // namespace chemistry
}  // namespace amanzi
#endif  // AMANZI_CHEMISTRY_SURFACE_SPECIES_HH_
