/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Class for aqueous equilibrium complexation reaction
*/

#include <sstream>
#include <iostream>
#include <iomanip>

#include "chemistry_utilities.hh"
#include "matrix_block.hh"
#include "surface_complex.hh"

namespace Amanzi {
namespace AmanziChemistry {

namespace acu = Amanzi::AmanziChemistry::utilities;

SurfaceComplex::SurfaceComplex() {
  species_names_.clear();
  species_ids_.clear();
  stoichiometry_.clear();
  logK_array_.clear();
}


SurfaceComplex::SurfaceComplex(const std::string& name,
                               const int id,
                               const std::vector<std::string>& species,
                               const std::vector<double>& stoichiometries,
                               const std::vector<int>& species_ids,
                               const double h2o_stoich,
                               const double free_site_stoich,
                               const double charge,
                               const double logK)
    : name_(name),
      identifier_(id),
      charge_(charge),
      species_names_(species),
      stoichiometry_(stoichiometries),
      species_ids_(species_ids),
      surface_concentration_(0.),
      free_site_name_("Unknown"),
      free_site_stoichiometry_(free_site_stoich),
      free_site_id_(-1),
      h2o_stoichiometry_(h2o_stoich),
      lnK_(acu::log_to_ln(logK)),
      lnQK_(0.),
      logK_(logK) {

  logK_array_.clear();

  set_ncomp(static_cast<int>(stoichiometries.size()));
}


SurfaceComplex::SurfaceComplex(const std::string& name,
                               const int id,
                               const std::vector<std::string>& species,
                               const std::vector<double>& stoichiometries,
                               const std::vector<int>& species_ids,
                               const double h2o_stoich,
                               const std::string& free_site_name,
                               const double free_site_stoich,
                               const int free_site_id,
                               const double charge,
                               const double logK)
    : name_(name),
      identifier_(id),
      charge_(charge),
      species_names_(species),
      stoichiometry_(stoichiometries),
      species_ids_(species_ids),
      surface_concentration_(0.),
      free_site_name_(free_site_name),
      free_site_stoichiometry_(free_site_stoich),
      free_site_id_(free_site_id),
      h2o_stoichiometry_(h2o_stoich),
      lnK_(acu::log_to_ln(logK)),
      lnQK_(0.),
      logK_(logK) {

  logK_array_.clear();

  set_ncomp(static_cast<int>(stoichiometries.size()));
}


void SurfaceComplex::Update(const std::vector<Species>& primarySpecies,
                            const SurfaceSite& surface_site) {
  double lnQK_temp = -lnK_;

  // Need to consider activity of water

  lnQK_temp += free_site_stoichiometry() *
      surface_site.ln_free_site_concentration();

  for (int i = 0; i < ncomp_; i++) {
    lnQK_temp += stoichiometry_[i] *
        primarySpecies[species_ids_[i]].ln_activity();
  }
  set_lnQK(lnQK_temp);
  set_surface_concentration(std::exp(lnQK()));
}


void SurfaceComplex::AddContributionToTotal(std::vector<double> *total) {
  for (int i = 0; i < ncomp_; i++) {
    (*total)[species_ids_[i]] += stoichiometry_[i] * surface_concentration();
  }
}


void SurfaceComplex::AddContributionToDTotal(
    const std::vector<Species>& primarySpecies,
    MatrixBlock* dtotal) {
}


void SurfaceComplex::Display(const Teuchos::Ptr<VerboseObject> vo) const {
  std::stringstream message;
  message << "    " << name() << " = ";
  message << free_site_stoichiometry_ << " " << free_site_name_ << " + ";
  if (h2o_stoichiometry_ > 0) {
    message << h2o_stoichiometry_ << " " << "H2O" << " + ";
  }
  for (unsigned int i = 0; i < species_names_.size(); i++) {
    message << stoichiometry_[i] << " " << species_names_[i];
    if (i < species_names_.size() - 1) {
      message << " + ";
    }
  }
  message << std::endl;
  message << std::setw(40) << " "
            << std::setw(10) << logK_
            << std::setw(10) << charge()
            << std::endl;
  vo->Write(Teuchos::VERB_HIGH, message.str());
}


void SurfaceComplex::display(const Teuchos::Ptr<VerboseObject> vo) const {
  std::stringstream message;
  message << "    " << name() << " = ";
  message << free_site_stoichiometry_ << " " << free_site_name_ << " + ";
  if (h2o_stoichiometry_ > 0) {
    message << h2o_stoichiometry_ << " " << "H2O" << " + ";
  }
  for (unsigned int i = 0; i < species_names_.size(); i++) {
    message << stoichiometry_[i] << " " << species_names_[i];
    if (i < species_names_.size() - 1) {
      message << " + ";
    }
  }
  message << std::endl;
  message << "     log K: " << logK_
            << "\n     charge: " << charge() << std::endl;
  vo->Write(Teuchos::VERB_HIGH, message.str());
}


void SurfaceComplex::DisplayResultsHeader(const Teuchos::Ptr<VerboseObject> vo) const {
  std::stringstream message;
  message << std::setw(15) << "Complex Name"
            << std::setw(15) << "Concentration"
            << std::endl;
  message << std::setw(15) << " "
            << std::setw(15) << "[mol/m^3]"
            << std::endl;
  vo->Write(Teuchos::VERB_HIGH, message.str());
}


void SurfaceComplex::DisplayResults(const Teuchos::Ptr<VerboseObject> vo) const {
  std::stringstream message;
  message << std::setw(15) << name()
          << std::scientific << std::setprecision(5)
          << std::setw(15) << surface_concentration()
          << std::endl;
  vo->Write(Teuchos::VERB_HIGH, message.str());
}

}  // namespace AmanziChemistry
}  // namespace Amanzi
