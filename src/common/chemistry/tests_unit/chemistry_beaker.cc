/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include <cstdlib>
#include <cmath>

#include <iostream>
#include <vector>

#include "UnitTest++.h"
#include "Teuchos_RCP.hpp"

#include "VerboseObject.hh"

#include "simple_thermo_database.hh"
#include "beaker.hh"
#include "activity_model_factory.hh"
#include "chemistry_exception.hh"

SUITE(BeakerTests) {
  using Amanzi::AmanziChemistry::Beaker;
  using Amanzi::AmanziChemistry::SimpleThermoDatabase;
  using Amanzi::AmanziChemistry::ActivityModelFactory;
  using Amanzi::AmanziChemistry::ChemistryException;
  using Amanzi::AmanziChemistry::ChemistryUnrecoverableError;
  using Amanzi::AmanziChemistry::ChemistryMemorySizeError;
  using Amanzi::AmanziChemistry::ChemistryInvalidInput;

  TEST(CheckBadComponentSizes) {
    Teuchos::ParameterList plist;
    auto vo = Teuchos::rcp(new Amanzi::VerboseObject("Chemistry PK", plist));

    SimpleThermoDatabase chem(vo);

    Beaker::BeakerComponents components;
    components.free_ion.clear();
    components.mineral_volume_fraction.clear();
    components.ion_exchange_sites.clear();
    components.total.clear();
    components.total_sorbed.clear();

    components.total.push_back(1.0e-3);  // H+
    components.total.push_back(1.0e-3);  // HCO3-
    components.total.push_back(0.0);
    components.free_ion.push_back(0.0);
    components.mineral_volume_fraction.push_back(0.0);
    components.ion_exchange_sites.push_back(0.0);

    Beaker::BeakerParameters parameters = chem.GetDefaultParameters();

    parameters.thermo_database_file = "chemistry_beaker_carbonate.bgd";
    parameters.activity_model_name = ActivityModelFactory::unit;

    bool correct_exception = false;

    try {
      // should throw an error
      chem.Setup(components, parameters);
    } catch (ChemistryMemorySizeError& e) {
      correct_exception = true;
    } catch (ChemistryException& e) {
    } catch (std::exception& e) {
    }

    CHECK(correct_exception);
  }  // end TEST(CheckBadComponentSizes)



  TEST(CheckBadDatabaseFile) {
    Teuchos::ParameterList plist;
    auto vo = Teuchos::rcp(new Amanzi::VerboseObject("Chemistry PK", plist));

    SimpleThermoDatabase chem(vo);

    Beaker::BeakerComponents components;
    components.free_ion.clear();
    components.mineral_volume_fraction.clear();
    components.ion_exchange_sites.clear();
    components.total.clear();
    components.total_sorbed.clear();

    components.total.push_back(1.0e-3);  // H+
    components.total.push_back(1.0e-3);  // HCO3-

    Beaker::BeakerParameters parameters = chem.GetDefaultParameters();

    parameters.thermo_database_file = "test_drivers/input/bad_database_file.bgd";
    parameters.activity_model_name = ActivityModelFactory::unit;

    bool correct_exception = false;

    try {
      // should throw an error
      chem.Setup(components, parameters);
    } catch (ChemistryUnrecoverableError& e) {
    } catch (ChemistryInvalidInput& e) {
      correct_exception = true;
    } catch (ChemistryException& e) {
    } catch (std::exception& e) {
    }

    CHECK(correct_exception);
  }  // end TEST(CheckBadDatabaseFile)

  TEST(CheckBadActivityModel) {
    Teuchos::ParameterList plist;
    auto vo = Teuchos::rcp(new Amanzi::VerboseObject("Chemistry PK", plist));

    SimpleThermoDatabase chem(vo);

    Beaker::BeakerComponents components;
    components.free_ion.clear();
    components.mineral_volume_fraction.clear();
    components.ion_exchange_sites.clear();
    components.total.clear();
    components.total_sorbed.clear();

    components.total.push_back(1.0e-3);  // H+
    components.total.push_back(1.0e-3);  // HCO3-

    Beaker::BeakerParameters parameters = chem.GetDefaultParameters();

    parameters.thermo_database_file = "chemistry_beaker_carbonate.bgd";
    parameters.activity_model_name = "bad activity model name";

    bool correct_exception = false;

    try {
      // should throw an error
      chem.Setup(components, parameters);
    } catch (ChemistryUnrecoverableError& e) {
    } catch (ChemistryInvalidInput& e) {
      correct_exception = true;
    } catch (ChemistryException& e) {
    } catch (std::exception& e) {
    }

    CHECK(correct_exception);
  }  // end TEST(CheckBadActivityModel)
}  // end SUITE(BeakerTests)
