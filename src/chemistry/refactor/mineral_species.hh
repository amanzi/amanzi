/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_MINERAL_SPECIES_HH_
#define AMANZI_CHEMISTRY_MINERAL_SPECIES_HH_

//
// Class for describing minerals
//
// Notes: How are we going to deal with sorption onto a mineral? a
// list of surface species associated with this mineral, or a mineral
// associated with a surface species?
//
//

#include <string>
#include <vector>
#include <ostream>

namespace amanzi {
namespace chemistry {

class MineralSpecies {
 public:
  typedef std::vector<MineralSpecies> MineralSpeciesArray;
  typedef std::string MineralSpeciesName;
  typedef int MineralSpeciesId;  // unsigned int?

  MineralSpecies();
  MineralSpecies(const MineralSpeciesId& id, const MineralSpeciesName& name, 
                 const double molecular_weight, const double molar_volume,
                 const double specific_surface_area);
  virtual ~MineralSpecies();

  //
  // public interface
  //

  // does this go here or in MineralReaction?
  void UpdateSurfaceAreaFromVolumeFraction(const double total_volume);
  void Display(std::ostream& output) const;

  //
  // accessor methods
  //

  MineralSpeciesId identifier(void) const {
    return this->identifier_;
  }
  double gram_molecular_weight(void) const {
    return this->gram_molecular_weight_;
  }
  MineralSpeciesName name(void) const {
    return this->name_;
  }


 protected:

 private:


  MineralSpeciesId identifier_;
  double gram_molecular_weight_;  // [grams/mole]
  MineralSpeciesName name_;

  double saturation_index_;
  double molar_volume_;  // [cm^3 / moles]
  double specific_surface_area_;  // [m^2/g]
  double surface_area_;  // [m^2 mineral surface/m^3 mineral]
  double volume_fraction_;  // [m^3 mineral / m^3 bulk]
  double mole_fraction_;  // ???[-]


};

}  // namespace chemistry
}  // namespace amanzi
#endif  // AMANZI_CHEMISTRY_MINERAL_SPECIES_HH_
