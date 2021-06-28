/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>

#include <UnitTest++.h>

#include "IonExchangeSite.hh"

SUITE(GeochemistryTestsIonExchangeSite) {
  namespace ac = Amanzi::AmanziChemistry;
  /*
    Unit tests for the IonExchangeSite object public interface

    Need better tests for incorrect use of this object
  */
  /*****************************************************************************
   **
   **  Common testing code
   **
   *****************************************************************************/
  class IonExchangeSiteTest {
   public:
    IonExchangeSiteTest();
    ~IonExchangeSiteTest();

   protected:
    ac::SpeciesId exchanger_id_;
    double exchanger_charge_;
    std::string exchanger_location_;
    double molecular_weight_;
    double size_;
    ac::SpeciesName exchanger_name_;
    ac::IonExchangeSite ies_;

   private:
  };  // end class IonExchangeSiteTest

  IonExchangeSiteTest::IonExchangeSiteTest()
      : exchanger_id_(7),
      exchanger_charge_(-2),
      exchanger_location_("Fe(OH)2"),
      molecular_weight_(45.6789),
      size_(9.8),
      exchanger_name_("X--"),
      ies_(exchanger_name_, exchanger_id_, exchanger_charge_,
           exchanger_location_, molecular_weight_, size_) {
  }

  IonExchangeSiteTest::~IonExchangeSiteTest() {
  }

  /*****************************************************************************
   **
   **  Setup test problems
   **
   *****************************************************************************/
  //
  // create a IonExchangeSite object by specifing parameters to the
  // constructor, the primary use of the object
  //

  //
  // check that the parameter data is set correctly
  //
  TEST_FIXTURE(IonExchangeSiteTest, IonExchangeSite_constructor_init_id) {
    CHECK_EQUAL(exchanger_id_, ies_.identifier());
  }

  TEST_FIXTURE(IonExchangeSiteTest, IonExchangeSite_constructor_init_charge) {
    CHECK_EQUAL(exchanger_charge_, ies_.charge());
  }

  TEST_FIXTURE(IonExchangeSiteTest, IonExchangeSite_constructor_init_gmw) {
    CHECK_EQUAL(molecular_weight_, ies_.gram_molecular_weight());
  }

  TEST_FIXTURE(IonExchangeSiteTest, IonExchangeSite_constructor_init_isp) {
    CHECK_EQUAL(size_, ies_.ion_size_parameter());
  }

  TEST_FIXTURE(IonExchangeSiteTest, IonExchangeSite_constructor_init_name) {
    CHECK_EQUAL(exchanger_name_, ies_.name());
  }

  TEST_FIXTURE(IonExchangeSiteTest, IonExchangeSite_constructor_init_location) {
    CHECK_EQUAL(exchanger_location_, ies_.location());
  }

  //
  // check that updating the cation exchange capacity works correctly
  //
  TEST_FIXTURE(IonExchangeSiteTest, IonExchangeSite_default_cation_exchange_capacity) {
    CHECK_EQUAL(0.0, ies_.cation_exchange_capacity());
  }

  TEST_FIXTURE(IonExchangeSiteTest, IonExchangeSite_set_cation_exchange_capacity) {
    double cec = 0.001;
    ies_.set_cation_exchange_capacity(cec);
    CHECK_EQUAL(cec, ies_.cation_exchange_capacity());
  }

  //
  // check that updating the concentrations works correctly, i.e. that
  // nothing happens, act/conc/act_coef always equal to one?
  //
  TEST_FIXTURE(IonExchangeSiteTest, IonExchangeSite_update_molality) {
    double molality = 0.001;
    ies_.update(molality);
    CHECK_EQUAL(1.0, ies_.molality());
  }

  TEST_FIXTURE(IonExchangeSiteTest, IonExchangeSite_update_ln_molality) {
    double molality = 0.001;
    double ln_molality = std::log(molality);
    ies_.update(molality);
    CHECK_EQUAL(0.0, ies_.ln_molality());
  }

  TEST_FIXTURE(IonExchangeSiteTest, IonExchangeSite_update_activity_coefficient) {
    double act_coef = 0.9;
    ies_.act_coef(act_coef);
    CHECK_EQUAL(0.9, ies_.act_coef());
    ies_.update();
    CHECK_EQUAL(1.0, ies_.act_coef());
  }

  TEST_FIXTURE(IonExchangeSiteTest, IonExchangeSite_update_ln_activity_coefficient) {
    double act_coef = 0.9;
    double ln_act_coef = std::log(act_coef);
    ies_.act_coef(act_coef);  // does not update ln_act_coef!
    ies_.update();  // forces the update of ln_act_coef
    CHECK_EQUAL(0.0, ies_.ln_act_coef());
  }

  TEST_FIXTURE(IonExchangeSiteTest, IonExchangeSite_update_activity) {
    double molality = 0.001;
    double act_coef = 0.9;
    double activity = molality * act_coef;
    ies_.act_coef(act_coef);
    ies_.update(molality);
    CHECK_EQUAL(1.0, ies_.activity());
  }


  TEST_FIXTURE(IonExchangeSiteTest, IonExchangeSite_update_ln_activity) {
    double molality = 0.001;
    double act_coef = 0.9;
    double ln_activity = std::log(molality * act_coef);
    ies_.act_coef(act_coef);
    ies_.update(molality);
    CHECK_EQUAL(0.0, ies_.ln_activity());
  }
}  // end SUITE(GeochemistryTestIonExchangeSite)
