/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include <sstream>
#include <string>

#include "ActivityModel.hpp"
#include "ActivityModelDebyeHuckel.hpp"
#include "ActivityModelUnit.hpp"
#include "ActivityModelFactory.hpp"
#include "ChemistryException.hpp"

const std::string ActivityModelFactory::debye_huckel = "debye-huckel";
const std::string ActivityModelFactory::unit = "unit";

ActivityModelFactory::ActivityModelFactory()
{
}  // end ActivityModelFactory constructor

ActivityModelFactory::~ActivityModelFactory()
{
}  // end ActivityModelFactory destructor

ActivityModel* ActivityModelFactory::Create(const std::string& model)
{
  ActivityModel* activity_model = NULL;

  if (model == debye_huckel) {
    activity_model = new ActivityModelDebyeHuckel();
  } else if (model == unit) {
    activity_model = new ActivityModelUnit();
  } else {
    // default type, error...?
    std::ostringstream error_stream;
    error_stream << "ERROR: ActivityModelFactory::Create(): \n";
    error_stream << "ERROR: Unknown activity model name: " << model << "\n"
                 << "       valid names: " << unit << "\n"
                 << "                    " << debye_huckel << "\n";
    throw ChemistryException(error_stream.str(), 
                             ChemistryException::kUnrecoverableError);    
  }

  if (activity_model == NULL) {
    // something went wrong, should throw an exception and exit gracefully....
    std::ostringstream error_stream;
    error_stream << "ERROR: ActivityModelFactory::Create(): \n";
    error_stream << "ERROR: Activity model was not created for some reason....\n";
    throw ChemistryException(error_stream.str(),
                             ChemistryException::kUnrecoverableError);    
  } else {
    // finish any additional setup
    // TODO: set the name in the object constructor so we can verify that the correct object was created.
    activity_model->name(model);
  }
  return activity_model;
}  // end Create()
