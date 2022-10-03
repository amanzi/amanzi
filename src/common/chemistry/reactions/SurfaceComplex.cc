/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Class for surface equilibrium complexation reaction
*/

#include <iostream>
#include <iomanip>
#include <sstream>

#include "ChemistryUtilities.hh"
#include "MatrixBlock.hh"
#include "ReactionString.hh"
#include "SurfaceComplex.hh"

namespace Amanzi {
namespace AmanziChemistry {

namespace acu = Amanzi::AmanziChemistry::utilities;

SurfaceComplex::SurfaceComplex(const std::string& name, int id,
                               const std::vector<Species>& primary_species,
                               const std::vector<SurfaceSite>& surface_sites,
                               const Teuchos::ParameterList& plist)
  : name_(name),
    id_(id),
    surface_concentration_(0.0),
    lnQK_(0.0)
{
  std::string reaction = plist.get<std::string>("reaction");
  charge_ = plist.get<int>("charge");

  if (plist.isSublist("equilibrium constant")) {
    auto x = plist.sublist("equilibrium constant").get<Teuchos::Array<double> >("T").toVector();
    auto y = plist.sublist("equilibrium constant").get<Teuchos::Array<double> >("Keq").toVector();
    func_ = Teuchos::rcp(new FunctionTabular(x, y, 0));

    double T = plist.get<double>("temperature");
    logK_ = (*func_)({T});
  } else {
    logK_ = plist.get<double>("equilibrium constant");
  }
  lnK_ = acu::log_to_ln(logK_);

  h2o_stoichiometry_ = 0.0;

  ParseReaction(reaction,
                primary_species, surface_sites,
                &species_names_, &stoichiometry_, &species_ids_,
                &free_site_name_, &free_site_stoichiometry_, &free_site_id_,
                &h2o_stoichiometry_);

  ncomp_ = stoichiometry_.size();
}


/* *******************************************************************
* Recalculates equilibrium constant
******************************************************************* */
void SurfaceComplex::UpdateTemperatureDependentCoefs(double T)
{
  if (func_.get() != nullptr) {
    logK_ = (*func_)({T});
  }
}


/* *******************************************************************
* Recalculate internal data
******************************************************************* */
void SurfaceComplex::Update(const std::vector<Species>& primary_species,
                            const SurfaceSite& surface_site)
{
  double lnQK_tmp = -lnK_;

  // Need to consider activity of water

  lnQK_tmp += free_site_stoichiometry() * surface_site.ln_free_site_concentration();

  for (int i = 0; i < ncomp_; i++) {
    lnQK_tmp += stoichiometry_[i] * primary_species[species_ids_[i]].ln_activity();
  }
  lnQK_ = lnQK_tmp;
  surface_concentration_ = std::exp(lnQK_);
}


void SurfaceComplex::AddContributionToTotal(std::vector<double> *total)
{
  for (int i = 0; i < ncomp_; i++) {
    (*total)[species_ids_[i]] += stoichiometry_[i] * surface_concentration();
  }
}


void SurfaceComplex::AddContributionToDTotal(
    const std::vector<Species>& primary_species,
    MatrixBlock* dtotal) {
}


void SurfaceComplex::Display(const Teuchos::Ptr<VerboseObject> vo) const
{
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
            << std::setw(10) << charge_
            << std::endl;
  vo->Write(Teuchos::VERB_HIGH, message.str());
}


void SurfaceComplex::DisplayResultsHeader(const Teuchos::Ptr<VerboseObject> vo) const
{
  std::stringstream message;
  message << std::setw(15) << "Complex Name"
            << std::setw(15) << "Concentration"
            << std::endl;
  message << std::setw(15) << " "
            << std::setw(15) << "[mol/m^3]"
            << std::endl;
  vo->Write(Teuchos::VERB_HIGH, message.str());
}


void SurfaceComplex::DisplayResults(const Teuchos::Ptr<VerboseObject> vo) const
{
  std::stringstream message;
  message << std::setw(15) << name()
          << std::scientific << std::setprecision(5)
          << std::setw(15) << surface_concentration()
          << std::endl;
  vo->Write(Teuchos::VERB_HIGH, message.str());
}

}  // namespace AmanziChemistry
}  // namespace Amanzi
