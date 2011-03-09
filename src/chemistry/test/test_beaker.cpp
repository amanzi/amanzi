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

SUITE(BeakerTests)
{

  TEST(CheckBadComponentSizes)
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

    parameters.thermo_database_file = "test/input/carbonate.bgd";
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
  }  // end TEST(CheckBadComponentSizes)



  TEST(CheckBadDatabaseFile)
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


  }  // end TEST(CheckBadDatabaseFile)

  TEST(CheckBadActivityModel)
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


  }  // end TEST(CheckBadActivityModel)


}  // end SUITE(BeakerTests)


