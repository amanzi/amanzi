/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Class for mineral reaction, should be written with the mineral
  as the reactant:

    Calcite = 1.0 Ca++ + 1.0 HCO3- -1.0 H+
*/
 
#include <sstream>
#include <iostream>
#include <iomanip>

#include "secondary_species.hh"
#include "matrix_block.hh"

#include "mineral.hh"

namespace Amanzi {
namespace AmanziChemistry {

Mineral::Mineral()
  : SecondarySpecies(),
    molar_volume_(0.0),
    specific_surface_area_(0.0),
    volume_fraction_(0.0)
{}


Mineral::Mineral(int id, const std::string& name,
                 const Teuchos::ParameterList& plist,
                 const std::vector<Species>& primary_species)
  : SecondarySpecies(id, name, plist, primary_species),
    volume_fraction_(0.0)
{
  molar_volume_ = plist.get<double>("molar volume");
  // convert: [cm^3/mole] --> [m^3/mole]
  molar_volume_ /= 1000000.0;

  specific_surface_area_ = plist.get<double>("specific surface area");
  // convert: [cm^2 mineral / cm^3 bulk] --> [m^2 mineral / m^3 bulk]
  specific_surface_area_ *= 100.0;
}


/* *******************************************************************
* NOTE: the rate is a dissolution rate so either need to use -rate
* or vol_frac -= .... inorder to get the correct
* dissolution/precipitation behavior.
*
* TODO(bandre): Right now we are just setting volume fraction to
* zero if they go negative, introducing mass balance errors! Need
* to adjust time step or reaction rate in the N-R solve!
******************************************************************* */
void Mineral::UpdateVolumeFraction(double rate, double delta_time) {
  // delta_vf = [m^3/mole] * [moles/m^3/sec] * [sec]
  volume_fraction_ -= molar_volume_ * rate * delta_time;
  if (volume_fraction_ < 0.0) {
    volume_fraction_ = 0.0;
  }
}


void Mineral::Update(const std::vector<Species>& primary_species, const Species& water_species) {
  double lnQK = -lnK_;
  for (int i = 0; i < ncomp(); i++) {
    lnQK += stoichiometry_.at(i) * primary_species.at(species_ids_.at(i)).ln_activity();
  }
  // Add the contribution of the water activity
  lnQK += SecondarySpecies::h2o_stoich_ * std::log(water_species.act_coef());
  lnQK_ = lnQK;
}


void Mineral::AddContributionToTotal(std::vector<double> *total)
{}


void Mineral::AddContributionToDTotal(const std::vector<Species>& primary_species,
                                      MatrixBlock* dtotal)
{}


void Mineral::Display(const Teuchos::Ptr<VerboseObject> vo) const {
  std::stringstream message;
  message << "    " << name() << " = ";
  for (unsigned int i = 0; i < species_names_.size(); i++) {
    message << std::setprecision(2) << stoichiometry_[i] << " " << species_names_[i];
    if (i < species_names_.size() - 1) {
      message << " + ";
    }
  }

  if (SecondarySpecies::h2o_stoich_!=0.0) {
    message << " + ";
    message << std::setprecision(2) << h2o_stoich_ << " " << "H2O";
  }
  message << std::endl;
  message << std::setw(40) << " "
          << std::setw(10) << std::setprecision(5) << std::fixed << logK_
          << std::setw(13) << std::scientific << molar_volume_
          << std::setw(13) << std::fixed << gram_molecular_weight()
          << std::setw(13) << specific_surface_area()
          << std::setw(13) << std::fixed << volume_fraction()
          << std::endl;
  vo->Write(Teuchos::VERB_HIGH, message.str());
}


void Mineral::DisplayResults(const Teuchos::Ptr<VerboseObject> vo) const {
  std::stringstream message;
  message << std::setw(15) << name()
          << std::scientific << std::setprecision(5)
          << std::setw(15) << Q_over_K()
          << std::fixed << std::setprecision(3)
          << std::setw(15) << saturation_index()
          << std::endl;
  vo->Write(Teuchos::VERB_HIGH, message.str());
}

}  // namespace AmanziChemistry
}  // namespace Amanzi
