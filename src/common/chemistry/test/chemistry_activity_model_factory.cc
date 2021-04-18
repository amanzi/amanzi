/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>
#include <typeinfo>

#include "UnitTest++.h"
#include "Teuchos_RCP.hpp"

#include "exceptions.hh"
#include "VerboseObject.hh"

#include "AqueousEquilibriumComplex.hh"
#include "ActivityModel.hh"
#include "ActivityModelFactory.hh"
#include "ActivityModelDebyeHuckel.hh"
#include "ActivityModelUnit.hh"
#include "Species.hh"

/*!
  @namespace Amanzi::AmanziChemistry::unit_tests::ActivityModelFactory

  @details Test that the activity model factory class returns the
  correct type of object. Use C++ RTTI to determine if the correct
  type of object was returned from the factory, e.g. see "typeid" at
  http://en.wikibooks.org/wiki/C++_Programming/RTTI

  @test Unit tests for: ActivityModelFactory
*/

SUITE(amanzi_chemistry_unit_tests_ActivityModelFactory) {
  /*****************************************************************************
  **  Test for ActivityModelFactory.cpp
  *****************************************************************************/
  namespace ac = Amanzi::AmanziChemistry;

  class ActivityModelFactoryTest {
   protected:
    ActivityModelFactoryTest();
    ~ActivityModelFactoryTest();

    void RunTest(const std::string name);

    ac::ActivityModelFactory amf_;
    ac::ActivityModel* activity_model_;

   private:
    Teuchos::RCP<Amanzi::VerboseObject> vo_;
  };

  ActivityModelFactoryTest::ActivityModelFactoryTest()
    : amf_(),
      activity_model_(NULL) {
    Teuchos::ParameterList plist;
    vo_ = Teuchos::rcp(new Amanzi::VerboseObject("Chemistry", plist));
  }

  ActivityModelFactoryTest::~ActivityModelFactoryTest() {
    delete activity_model_;
  }

  void ActivityModelFactoryTest::RunTest(const std::string name) {
    ac::ActivityModel::ActivityModelParameters parameters;
    parameters.database_filename = "";
    parameters.pitzer_jfunction = "";
    std::vector<ac::Species> primaries;
    std::vector<ac::AqueousEquilibriumComplex> secondaries;

    activity_model_ = amf_.Create(name, parameters, primaries, secondaries, vo_.ptr());
  }

  /*!
    @class Amanzi::AmanziChemistry::unit_tests::ActivityModelFactory::ActivityModelFactory_unit

    @brief ActivityModelFactory_unit

    @details Test that an ActivityModelUnit object is returned from the
    activity model factory when a "unit" model is specified.

    @test ActivityModelFactory::Create()

  */
  TEST_FIXTURE(ActivityModelFactoryTest, ActivityModelFactory_unit) {
    std::string name("unit");
    RunTest(name);
    CHECK_EQUAL(typeid(ac::ActivityModelUnit).name(), typeid(*activity_model_).name());
  }  // end TEST_FIXTURE()

  /*!
    @class Amanzi::AmanziChemistry::unit_tests::ActivityModelFactory::ActivityModelFactory_debyehuckel

    @brief ActivityModelFactory_debyehuckel

    @details Test that an ActivityModelDebyeHuckel object is returned from the
    activity model factory when a "debye-huckel" model is specified.

    @test ActivityModelFactory::Create()

  */
  TEST_FIXTURE(ActivityModelFactoryTest, ActivityModelFactory_debyehuckel) {
    std::string name("debye-huckel");
    RunTest(name);
    CHECK_EQUAL(typeid(ac::ActivityModelDebyeHuckel).name(), typeid(*activity_model_).name());
  }  // end TEST_FIXTURE()

  /*!
    @class Amanzi::AmanziChemistry::unit_tests::ActivityModelFactory::ActivityModelFactory_invalid

    @brief ActivityModelFactory_invalid

    @details Test that a chemistry exception is thrown when an invalid activity
    model is requested.

    @test ActivityModelFactory::Create()

  */
  TEST_FIXTURE(ActivityModelFactoryTest, ActivityModelFactory_invalid) {
    std::string name("invalid-name");
    CHECK_THROW(RunTest(name), Exceptions::Amanzi_exception);
    CHECK(!activity_model_);
  }  // end TEST_FIXTURE()



}  // end SUITE(amanzi_chemistry_unit_tests_ActivityModelFactory)
