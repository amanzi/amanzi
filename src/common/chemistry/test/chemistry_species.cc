#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include <UnitTest++.h>

#include "ChemistryException.hh"
#include "Species.hh"

SUITE(GeochemistryTestsSpecies) {
  namespace ac = Amanzi::AmanziChemistry;
  /*
    Unit tests for the Species object public interface

    Need better tests for incorrect use of this object
  */
  /*****************************************************************************
   **
   **  Common testing code
   **
   *****************************************************************************/
  class SpeciesTest {
   public:
    SpeciesTest();
    ~SpeciesTest();

   protected:
    int id_;
    double charge_;
    double gram_molecular_weight_;
    double ion_size_parameter_;
    std::string name_;

    Teuchos::ParameterList plist_;
  };  // end class SpeciesTest

  SpeciesTest::SpeciesTest()
    : id_(3),
      charge_(2.0),
      gram_molecular_weight_(40.0780),
      ion_size_parameter_(6.0),
      name_("Ca++") {
    plist_.set<int>("charge", charge_)
          .set<double>("gram molecular weight", gram_molecular_weight_)
          .set<double>("ion size parameter", ion_size_parameter_);
  }

  SpeciesTest::~SpeciesTest() {
  }

  /*****************************************************************************
   **
   **  Setup test problems
   **
   *****************************************************************************/
  //
  // create a species object by specifing parameters to the
  // constructor, the primary use of the object
  //

  //
  // check that the parameter data is set correctly
  //
  TEST_FIXTURE(SpeciesTest, Species_constructor_init_id) {
    ac::Species species(id_, name_, plist_);
    CHECK_EQUAL(id_, species.identifier());
  }

  TEST_FIXTURE(SpeciesTest, Species_constructor_init_charge) {
    ac::Species species(id_, name_, plist_);
    CHECK_EQUAL(charge_, species.charge());
  }

  TEST_FIXTURE(SpeciesTest, Species_constructor_init_gmw) {
    ac::Species species(id_, name_, plist_);
    CHECK_EQUAL(gram_molecular_weight_, species.gram_molecular_weight());
  }

  TEST_FIXTURE(SpeciesTest, Species_constructor_init_isp) {
    ac::Species species(id_, name_, plist_);
    CHECK_EQUAL(ion_size_parameter_, species.ion_size_parameter());
  }

  TEST_FIXTURE(SpeciesTest, Species_constructor_init_name) {
    ac::Species species(id_, name_, plist_);
    CHECK_EQUAL(name_, species.name());
  }

  //
  // check that exceptions are thrown for invalid data
  //
  TEST_FIXTURE(SpeciesTest, Species_constructor_invalid_id) {
    CHECK_THROW(ac::Species species(-1, name_, plist_),
                ac::ChemistryException);
  }

  TEST_FIXTURE(SpeciesTest, Species_constructor_invalid_molecular_weight) {
    auto plist = plist_;
    plist.set<double>("gram molecular weight", -45.678);
    CHECK_THROW(ac::Species species(id_, name_, plist),
                ac::ChemistryException);
  }

  TEST_FIXTURE(SpeciesTest, Species_constructor_invalid_ion_size) {
    auto plist = plist_;
    plist.set<double>("ion size parameter", -9.0);
    CHECK_THROW(ac::Species species(id_, name_, plist),
                ac::ChemistryException);
  }

  //
  // check that updating the concentrations works correctly
  //
  TEST_FIXTURE(SpeciesTest, Species_update_molality) {
    ac::Species species(id_, name_, plist_);
    double molality = 0.001;

    species.update(molality);
    CHECK_CLOSE(molality, species.molality(), 1.0e-9);
  }

  TEST_FIXTURE(SpeciesTest, Species_update_ln_molality) {
    ac::Species species(id_, name_, plist_);
    double molality = 0.001;
    double ln_molality = std::log(molality);

    species.update(molality);
    CHECK_CLOSE(ln_molality, species.ln_molality(), 1.0e-9);
  }

  TEST_FIXTURE(SpeciesTest, Species_update_activity_coefficient) {
    ac::Species species(id_, name_, plist_);
    double act_coef = 0.9;

    species.set_act_coef(act_coef);
    CHECK_CLOSE(act_coef, species.act_coef(), 1.0e-9);
  }

  TEST_FIXTURE(SpeciesTest, Species_update_ln_activity_coefficient) {
    ac::Species species(id_, name_, plist_);
    double act_coef = 0.9;
    double ln_act_coef = std::log(act_coef);

    species.set_act_coef(act_coef);  // does not update ln_act_coef!
    species.update();  // forces the update of ln_act_coef
    CHECK_CLOSE(ln_act_coef, species.ln_act_coef(), 1.0e-9);
  }

  TEST_FIXTURE(SpeciesTest, Species_update_activity) {
    ac::Species species(id_, name_, plist_);
    double molality = 0.001;
    double act_coef = 0.9;
    double activity = molality * act_coef;

    species.set_act_coef(act_coef);
    species.update(molality);
    CHECK_CLOSE(activity, species.activity(), 1.0e-9);
  }


  TEST_FIXTURE(SpeciesTest, Species_update_ln_activity) {
    ac::Species species(id_, name_, plist_);
    double molality = 0.001;
    double act_coef = 0.9;
    double ln_activity = std::log(molality * act_coef);

    species.set_act_coef(act_coef);
    species.update(molality);
    CHECK_CLOSE(ln_activity, species.ln_activity(), 1.0e-9);
  }

  //
  // create with the empty constructor, no initialization
  //

  /* We don't really want this interface used, but it appears to be
     necessary for the stl containers, so it should do the "correct"
     thing.... What default behavior do we really want?  The object
     doesn't have another initialization routine or any way to
     verify that it has been initialized, so I guess the default
     behavior needs to be do the correct thing with zero values.
     Species Ca;
     SpeciesId id = 0;
     double charge = 0.0;
     double gram_molecular_weight = 0.0;
     double ion_size_parameter = 0.0;
     SpeciesName name("");
     double molality = 1.0e-9;
     double activity = 1.0;
     double act_coef = 1.0;
     double ln_molality = 0.0;
     double ln_activity = 0.0;
     double ln_act_coef = 0.0;
  */
}  // end SUITE(GeochemistryTestSpecies)
