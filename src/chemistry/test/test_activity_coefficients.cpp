/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>
#include <typeinfo>

#include <UnitTest++.h>

#include "Species.hpp"
#include "ActivityModelFactory.hpp"
#include "ActivityModelUnit.hpp"
#include "ActivityModelDebyeHuckel.hpp"
#include "ActivityModel.hpp"
#include "ChemistryException.hpp"

SUITE(GeochemistryTestsActivityModels)
{
  /*****************************************************************************
   **
   **  Test for ActivityModelFactory.cpp
   **
   *****************************************************************************/
  
  /*
    if you pass the factory a valid name, it should return the correct
    activity model class. If you pass it an invalid name, it should
    throw an error.
   */

  class ActivityModelFactoryTest
  {
   public:

   protected:
    ActivityModelFactoryTest();
    ~ActivityModelFactoryTest();

    void RunTest(const std::string name);

    ActivityModelFactory amf_;
    ActivityModel* activity_model_;
    
   private:

  };

  ActivityModelFactoryTest::ActivityModelFactoryTest()
      : amf_(),
      activity_model_(NULL)
  {
  }

  ActivityModelFactoryTest::~ActivityModelFactoryTest()
  {
    delete activity_model_;
  }

  void ActivityModelFactoryTest::RunTest(const std::string name)
  {
    activity_model_ = amf_.Create(name);
  }

  // use C++ RTTI to determine if the correct type of object was
  // returned from the factory, e.g. see "typeid" at
  // http://en.wikibooks.org/wiki/C++_Programming/RTTI

  TEST_FIXTURE(ActivityModelFactoryTest, ActivityModelFactory_unit)
  {
    std::string name("unit");
    RunTest(name);
    CHECK_EQUAL(typeid(ActivityModelUnit).name(), typeid(*activity_model_).name());
  } // end TEST_FIXTURE()

  TEST_FIXTURE(ActivityModelFactoryTest, ActivityModelFactory_debyehuckel)
  {
    std::string name("debye-huckel");
    RunTest(name);
    CHECK_EQUAL(typeid(ActivityModelDebyeHuckel).name(), typeid(*activity_model_).name());
  } // end TEST_FIXTURE()

  TEST_FIXTURE(ActivityModelFactoryTest, ActivityModelFactory_invalid)
  {
    std::string name("invalid-name");
    CHECK_THROW(RunTest(name), ChemistryException);
    CHECK( !activity_model_);
  } // end TEST_FIXTURE()


  /*
    Test ionic strength and activity coefficient calculations for
    different activity models. Arbitrary choice of 6 ions with charge
    of +-1, +-2, +-3. Adjust concentrations to yield an ionic strength
    of 0.025.

    - Unit activity coefficients: always 1.0 regardless of species

    - Debye-Huckel: check first two digits of activity coefficients
    at 25C with Langmuir, 1997, Aqueous Environmental Geochemistry,
    Table 4.1 and 4.2, pg 130-131. Note, code uses slightly different
    debyeA and debyeB parameters.

    - source for higher number of sig figs?

    - For temperature dependance, run 5-10 temperature values, then
    store results in a vector and use CHECK_ARRAY_CLOSE(). Source
    for temperature dependance?

  */

  /*****************************************************************************
   **
   **  Common testing code
   **
   *****************************************************************************/
  class ActivityModelTest
  {
   public:
    ActivityModelTest();
    ~ActivityModelTest();

    void RunTest(const std::string name, double* gamma);

    void set_activity_model_name(const std::string name) { activity_model_name_ = name; };
    std::string activity_model_name(void) const { return activity_model_name_; };
    double ionic_strength(void) { return activity_model_->ionic_strength(); };
    double tolerance(void) { return tolerance_; };

   protected:
    ActivityModelFactory amf_;

   private:
    double tolerance_;
    ActivityModel* activity_model_;
    std::string activity_model_name_;
    SpeciesArray species_;
    Species H_p;
    Species OH_m;
    Species Ca_pp;
    Species SO4_mm;
    Species Al_ppp;
    Species PO4_mmm;
    std::vector<AqueousEquilibriumComplex> aqueous_complexes_;
  }; // end class SpeciationTest

  ActivityModelTest::ActivityModelTest()
      : amf_(),
      tolerance_(1.0e-5),
      activity_model_name_(""),
      H_p(0, "H+", 1.0, 1.0079, 9.0),
      OH_m(1, "OH-", -1.0, 17.0073, 3.5),
      Ca_pp(2, "Ca++", 2.0, 40.0780, 6.0),
      SO4_mm(3, "SO4--", -2.0, 96.0636, 4.0),
      Al_ppp(4, "Al+++", 3.0, 26.9815, 9.0),
      PO4_mmm(5, "PO4---", -3.0, 94.9714, 4.0) {
    // set concentrations to get ionic strength of 0.025
    H_p.update(0.0005);
    OH_m.update(0.0015);
    Ca_pp.update(0.001);
    SO4_mm.update(0.002);
    Al_ppp.update(0.003);
    PO4_mmm.update(0.001);
    species_.push_back(H_p);
    species_.push_back(OH_m);
    species_.push_back(Ca_pp);
    species_.push_back(SO4_mm);
    species_.push_back(Al_ppp);
    species_.push_back(PO4_mmm);
    // TODO: should add some aqueous complexes to test ionic strength....
    aqueous_complexes_.clear();
  }

  ActivityModelTest::~ActivityModelTest()
  {
    delete activity_model_;
  }

  void ActivityModelTest::RunTest(const std::string name, double* gamma)
  {
    int index = -1;
    for (std::vector<Species>::iterator primary = species_.begin();
         primary != species_.end(); primary++) {
      if (primary->name() == name) {
        index = primary->identifier();
      }
    }
    *gamma = -1.0; // final value should always be > 0
    activity_model_ = amf_.Create(activity_model_name());
    activity_model_->CalculateIonicStrength(species_, aqueous_complexes_);
    *gamma = activity_model_->Evaluate(species_.at(index));
  }  // end ActivityModelTest::RunTest()

  /*****************************************************************************
   **
   **  Test for ActivityModel.cpp
   **
   *****************************************************************************/
  TEST_FIXTURE(ActivityModelTest, ActivityModelIonicStrength)
  {
    set_activity_model_name(ActivityModelFactory::unit);
    double gamma;
    RunTest("H+", &gamma);
    //std::cout << "ionic strength: " << ionic_strength() << std::endl;
    CHECK_CLOSE(0.025, ionic_strength(), tolerance());
  }  // end TEST_FIXTURE()

  /*****************************************************************************
   **
   **  Test for ActivityModelUnit.cpp
   **
   *****************************************************************************/
  TEST_FIXTURE(ActivityModelTest, ActivityModelUnit_H)
  {
    set_activity_model_name(ActivityModelFactory::unit);
    double gamma;
    RunTest("H+", &gamma);
    CHECK_EQUAL(1.0, gamma);
  }  // end TEST_FIXTURE()

  TEST_FIXTURE(ActivityModelTest, ActivityModelUnit_OH)
  {
    set_activity_model_name(ActivityModelFactory::unit);
    double gamma;
    RunTest("OH-", &gamma);
    CHECK_EQUAL(1.0, gamma);
  }  // end TEST_FIXTURE()

  TEST_FIXTURE(ActivityModelTest, ActivityModelUnit_Ca)
  {
    set_activity_model_name(ActivityModelFactory::unit);
    double gamma;
    RunTest("Ca++", &gamma);
    CHECK_EQUAL(1.0, gamma);
  }  // end TEST_FIXTURE()

  TEST_FIXTURE(ActivityModelTest, ActivityModelUnit_SO4)
  {
    set_activity_model_name(ActivityModelFactory::unit);
    double gamma;
    RunTest("SO4--", &gamma);
    CHECK_EQUAL(1.0, gamma);
  }  // end TEST_FIXTURE()

  TEST_FIXTURE(ActivityModelTest, ActivityModelUnit_Al)
  {
    set_activity_model_name(ActivityModelFactory::unit);
    double gamma;
    RunTest("Al+++", &gamma);
    CHECK_EQUAL(1.0, gamma);
  }  // end TEST_FIXTURE()

  TEST_FIXTURE(ActivityModelTest, ActivityModelUnit_PO4)
  {
    set_activity_model_name(ActivityModelFactory::unit);
    double gamma;
    RunTest("PO4---", &gamma);
    CHECK_EQUAL(1.0, gamma);
  }  // end TEST_FIXTURE()

  /*****************************************************************************
   **
   **  Test for ActivityModelDebyeHuckel.cpp
   **
   *****************************************************************************/
  TEST_FIXTURE(ActivityModelTest, ActivityModelDebyeHuckel_H)
  {
    set_activity_model_name(ActivityModelFactory::debye_huckel);
    double gamma;
    RunTest("H+", &gamma);
    CHECK_CLOSE(0.88, gamma, 1.0e-2);
  }  // end TEST_FIXTURE()

  TEST_FIXTURE(ActivityModelTest, ActivityModelDebyeHuckel_OH)
  {
    set_activity_model_name(ActivityModelFactory::debye_huckel);
    double gamma;
    RunTest("OH-", &gamma);
    CHECK_CLOSE(0.855, gamma, 1.0e-2);
  }  // end TEST_FIXTURE()

  TEST_FIXTURE(ActivityModelTest, ActivityModelDebyeHuckel_Ca)
  {
    set_activity_model_name(ActivityModelFactory::debye_huckel);
    double gamma;
    RunTest("Ca++", &gamma);
    CHECK_CLOSE(0.57, gamma, 1.0e-2);
  }  // end TEST_FIXTURE()

  TEST_FIXTURE(ActivityModelTest, ActivityModelDebyeHuckel_SO4)
  {
    set_activity_model_name(ActivityModelFactory::debye_huckel);
    double gamma;
    RunTest("SO4--", &gamma);
    CHECK_CLOSE(0.545, gamma, 1.0e-2);
  }  // end TEST_FIXTURE()

  TEST_FIXTURE(ActivityModelTest, ActivityModelDebyeHuckel_Al)
  {
    set_activity_model_name(ActivityModelFactory::debye_huckel);
    double gamma;
    RunTest("Al+++", &gamma);
    CHECK_CLOSE(0.325, gamma, 1.0e-2);
  }  // end TEST_FIXTURE()

  TEST_FIXTURE(ActivityModelTest, ActivityModelDebyeHuckel_PO4)
  {
    set_activity_model_name(ActivityModelFactory::debye_huckel);
    double gamma;
    RunTest("PO4---", &gamma);
    CHECK_CLOSE(0.25, gamma, 1.0e-2);
  }  // end TEST_FIXTURE()

}  // end SUITE(GeochemistryTestActivityModels)
