/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include "activity-model-factory.hh"

#include <sstream>
#include <string>

#include "activity-model.hh"
#include "activity-model-debye-huckel.hh"
#include "activity-model-unit.hh"
#include "chemistry-exception.hh"
#include "exceptions.hh"

const std::string ActivityModelFactory::debye_huckel = "debye-huckel";
const std::string ActivityModelFactory::unit = "unit";

ActivityModelFactory::ActivityModelFactory() {
}  // end ActivityModelFactory constructor

ActivityModelFactory::~ActivityModelFactory() {
}  // end ActivityModelFactory destructor

ActivityModel* ActivityModelFactory::Create(const std::string& model) {
  ActivityModel* activity_model = NULL;

  if (model == debye_huckel) {
    activity_model = new ActivityModelDebyeHuckel();
  } else if (model == unit) {
    activity_model = new ActivityModelUnit();
  } else {
    // default type, error...?
    std::ostringstream error_stream;
    error_stream << "ActivityModelFactory::Create(): \n"
                 << "Unknown activity model name: " << model << "\n"
                 << "       valid names: " << unit << "\n"
                 << "                    " << debye_huckel << "\n";
    Exceptions::amanzi_throw(ChemistryInvalidInput(error_stream.str()));
  }

  if (activity_model == NULL) {
    // something went wrong, should throw an exception and exit gracefully....
    std::ostringstream error_stream;
    error_stream << "ActivityModelFactory::Create(): \n"
                 << "Activity model was not created for some reason....\n";
    Exceptions::amanzi_throw(ChemistryException(error_stream.str()));
  } else {
    // finish any additional setup

    // TODO(bandre): set the name in the object constructor so we can
    // verify that the correct object was created.
    activity_model->name(model);
  }
  return activity_model;
}  // end Create()
