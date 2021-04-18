/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include <cstdlib>
#include <cmath>

#include <iostream>
#include <vector>

#include "UnitTest++.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "exceptions.hh"
#include "VerboseObject.hh"

#include "ActivityModelFactory.hh"
#include "Beaker.hh"
#include "SimpleThermoDatabase.hh"

SUITE(BeakerTests) {
  using Amanzi::AmanziChemistry::Beaker;
  using Amanzi::AmanziChemistry::SimpleThermoDatabase;
  using Amanzi::AmanziChemistry::ActivityModelFactory;

  TEST(CheckBadActivityModel) {
    Teuchos::ParameterList plist;
    auto vo = Teuchos::rcp(new Amanzi::VerboseObject("Chemistry PK", plist));

    auto bgd_list = Teuchos::getParametersFromXmlFile("test/chemistry_beaker_carbonate.xml");
    SimpleThermoDatabase chem(bgd_list, vo);

    Beaker::BeakerState state;
    state.free_ion.clear();
    state.mineral_volume_fraction.clear();
    state.ion_exchange_sites.clear();
    state.total.clear();
    state.total_sorbed.clear();

    state.total.push_back(1.0e-3);  // H+
    state.total.push_back(1.0e-3);  // HCO3-

    Beaker::BeakerParameters parameters;

    parameters.activity_model_name = "bad activity model name";
    parameters.tolerance = 1.0e-12;

    bool correct_exception = false;

    try {
      // should throw an error
      chem.Initialize(parameters);
    } catch (Exceptions::Amanzi_exception& e) {
      correct_exception = true;
    } catch (std::exception& e) {
    }

    CHECK(correct_exception);
  }
}

