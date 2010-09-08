#ifndef __Secondary_Species_hpp__
#define __Secondary_Species_hpp__

#include "Teuchos_RCP.hpp"
#include "State.hpp"
#include "Chemistry_State.hpp"

#include "Species.hpp"

// Base class for species

class SecondarySpecies : public Species {

public:
  SecondarySpecies(const double s_charge, const double s_GMW, 
		   const ActivityCoefficient* s_activityCoefficient, 
		   Species::SpeciesName s_name);

  ~SecondarySpecies();

private:
  std::vector<Species::SpeciesId> components_;
  std::vector<double> stoich_coefficients_;


};

#endif
