/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include <iostream>
#include "stdlib.h"
#include "math.h"
#include "UnitTest++.h"
#include <vector>

#include "SimpleThermoDatabase.hpp"
#include "Beaker.hpp"
#include "ActivityModelFactory.hpp"
#include "Verbosity.hpp"
#include "ChemistryException.hpp"

TEST(CHECK_BAD_COMPONENT_SIZES)
{
  SimpleThermoDatabase chem;

  Beaker::BeakerComponents components;
  components.free_ion.clear();
  components.minerals.clear();
  components.ion_exchange_sites.clear();
  components.total.clear();
  components.total_sorbed.clear();

  components.total.push_back(1.0e-3);  // H+
  components.total.push_back(1.0e-3);  // HCO3-
  components.total.push_back(0.0);
  components.free_ion.push_back(0.0);
  components.minerals.push_back(0.0);
  components.ion_exchange_sites.push_back(0.0);

  Beaker::BeakerParameters parameters = chem.GetDefaultParameters();

  parameters.thermo_database_file = "test_drivers/input/carbonate.bgd";
  parameters.activity_model_name = ActivityModelFactory::unit;

  bool have_exception = false;

  try {
  // should throw an error 
  chem.Setup(components, parameters);
  }
  catch (ChemistryException& e) {
    have_exception = true;
  }

  CHECK_EQUAL(true, have_exception);
}  // end TEST(CHECK_BAD_COMPONENT_SIZES)



TEST(CHECK_BAD_DATABASE_FILE)
{
  SimpleThermoDatabase chem;

  Beaker::BeakerComponents components;
  components.free_ion.clear();
  components.minerals.clear();
  components.ion_exchange_sites.clear();
  components.total.clear();
  components.total_sorbed.clear();

  components.total.push_back(1.0e-3);  // H+
  components.total.push_back(1.0e-3);  // HCO3-

  Beaker::BeakerParameters parameters = chem.GetDefaultParameters();

  parameters.thermo_database_file = "test_drivers/input/bad_database_file.bgd";
  parameters.activity_model_name = ActivityModelFactory::unit;

  bool have_exception = false;

  try {
  // should throw an error 
  chem.Setup(components, parameters);
  }
  catch (ChemistryException& e) {
    have_exception = true;
  }
  CHECK_EQUAL(true, have_exception);


}  // end TEST(CHECK_BAD_DATABASE_FILE)

TEST(CHECK_BAD_ACTIVITY_MODEL)
{
  SimpleThermoDatabase chem;

  Beaker::BeakerComponents components;
  components.free_ion.clear();
  components.minerals.clear();
  components.ion_exchange_sites.clear();
  components.total.clear();
  components.total_sorbed.clear();

  components.total.push_back(1.0e-3);  // H+
  components.total.push_back(1.0e-3);  // HCO3-

  Beaker::BeakerParameters parameters = chem.GetDefaultParameters();

  parameters.thermo_database_file = "test_drivers/input/carbonate.bgd";
  parameters.activity_model_name = "bad activity model name";

  bool have_exception = false;

  try {
  // should throw an error 
  chem.Setup(components, parameters);
  }
  catch (ChemistryException& e) {
    have_exception = true;
  }
  CHECK_EQUAL(true, have_exception);


}  // end TEST(CHECK_BAD_ACTIVITY_MODEL)





// this test will obviously fail
// TEST(EXAMPLE3) {

//   using namespace std;

//   vector<int> x(10);
//   vector<int> y(10);

//   x[1] = 1.0;
//   y[1] = 2.0;

//   CHECK_EQUAL(x[1],y[1]);
//   CHECK_CLOSE(x[1],y[1],0.0001);
// }


