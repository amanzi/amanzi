/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ben Andre
           Sergio A Bea
*/

#include <sstream>
#include <string>

#include "ActivityModel.hh"
#include "ActivityModelDebyeHuckel.hh"
#include "ActivityModelPitzerHWM.hh"
#include "ActivityModelUnit.hh"
#include "chemistry_exception.hh"
#include "exceptions.hh"
#include "species.hh"

#include "ActivityModelFactory.hh"

namespace Amanzi {
namespace AmanziChemistry {

const std::string ActivityModelFactory::debye_huckel = "debye-huckel";
const std::string ActivityModelFactory::pitzer_hwm = "pitzer-hwm";
const std::string ActivityModelFactory::unit = "unit";

ActivityModel* ActivityModelFactory::Create(
    const std::string& model,
    const ActivityModel::ActivityModelParameters& parameters,
    const std::vector<Species>& primary_species,
    const std::vector<AqueousEquilibriumComplex>& secondary_species,
    const Teuchos::Ptr<VerboseObject> vo)
{
  ActivityModel* activity_model = nullptr;

  if (model == debye_huckel) {
    activity_model = new ActivityModelDebyeHuckel();
  } else if (model == pitzer_hwm) {
    activity_model = new ActivityModelPitzerHWM();
  } else if (model == unit) {
    activity_model = new ActivityModelUnit();
  } else {
    std::ostringstream oss;
    oss << "Unknown activity model name: " << model << "\n"
        << "  valid names: " << unit << "\n"
        << "               " << debye_huckel << "\n"
        << "               " << pitzer_hwm << "\n" ;
    Exceptions::amanzi_throw(ChemistryInvalidInput(oss.str()));
  }

  activity_model->set_verbosity(vo);
  activity_model->name(model);
  activity_model->Setup(parameters, primary_species, secondary_species);

  return activity_model;
}

}  // namespace AmanziChemistry
}  // namespace Amanzi
