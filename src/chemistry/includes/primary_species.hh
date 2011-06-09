/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_PRIMARY_SPECIES_HH_
#define AMANZI_CHEMISTRY_PRIMARY_SPECIES_HH_

#include "Teuchos_RCP.hpp"
#include "State.hpp"
#include "Chemistry_State.hpp"


// Base class for species

class PrimarySpecies : public Species {

  PrimarySpecies(const double s_charge, const double s_GMW, 
		   const ActivityCoefficient* s_activityCoefficient, 
		   Species::SpeciesName s_name);

  ~PrimarySpecies();

private:
  std::vector<Species::SpeciesId> primary_list_;
  std::vector<double> primary_coefficients_;


};

#endif
