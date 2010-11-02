/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

#include <string>

#include "ActivityModel.hpp"
#include "ActivityModelDebyeHuckel.hpp"
#include "ActivityModelUnit.hpp"
#include "ActivityModelFactory.hpp"

const string ActivityModelFactory::debye_huckel = "debye-huckel";
const string ActivityModelFactory::unit = "unit";

ActivityModelFactory::ActivityModelFactory()
{
}  // end ActivityModelFactory constructor

ActivityModelFactory::~ActivityModelFactory()
{
}  // end ActivityModelFactory destructor

ActivityModel* ActivityModelFactory::Create(std::string model)
{
  ActivityModel* activity_model = NULL;

  if (model == debye_huckel) {
    activity_model = new ActivityModelDebyeHuckel();
  } else if (model == unit) {
    activity_model = new ActivityModelUnit();
  } else {
    // default type, error...?
  }

  if (activity_model == NULL) {
    // something went wrong, should throw an exception and exit gracefully....
  } else {
    // finish any additional setup
    activity_model->name(model);
  }
  return activity_model;
}  // end Create()
