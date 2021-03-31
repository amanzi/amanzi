#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <typeinfo>

#include <UnitTest++.h>

#include "secondary_species.hh"
#include "chemistry_exception.hh"

SUITE(GeochemistryTestsSecondarySpecies) {
  namespace ac = Amanzi::AmanziChemistry;
  /*****************************************************************************
   **
   **  Test for SecondarySpecies.cpp
   **
   *****************************************************************************/

  class SecondarySpeciesTest {
   public:

   protected:
    SecondarySpeciesTest();
    ~SecondarySpeciesTest();

    std::string name_;
    int secondary_id_;
    double h2o_stoich_;
    double charge_;
    double gram_molecular_weight_;
    double ion_size_parameter_;
    double logK_;

    std::vector<std::string> species_names_;
    std::vector<double> stoichiometry_;
    std::vector<int> species_ids_;
  };

  SecondarySpeciesTest::SecondarySpeciesTest()
    : name_("CO3--"),
      secondary_id_(3),
      h2o_stoich_(0.0),
      charge_(-2),
      gram_molecular_weight_(60.0092),
      ion_size_parameter_(4.5),
      logK_(10.3288) {
    species_names_.clear();
    stoichiometry_.clear();
    species_ids_.clear();

    species_names_.push_back("H+");
    stoichiometry_.push_back(-1.0);
    species_ids_.push_back(0);

    species_names_.push_back("HCO3-");
    stoichiometry_.push_back(1.0);
    species_ids_.push_back(1);
  }


  SecondarySpeciesTest::~SecondarySpeciesTest() {
  }

  //
  // SecondarySpecies has a pure virtual functions, need a dummy
  // object to create and test.
  //
  class MockSecondarySpecies : public ac::SecondarySpecies {
   public:
    MockSecondarySpecies() : SecondarySpecies() {}
    MockSecondarySpecies(const std::string name, const int id,
                         const std::vector<std::string>& species_names,
                         const std::vector<double>& stoichiometry,
                         const std::vector<int>& species_ids,
                         const double h2o_stoich,
                         const double charge,
                         const double mol_wt,
                         const double size,
                         const double logK)
        : SecondarySpecies(name, id, species_names, stoichiometry, species_ids,
                           h2o_stoich, charge, mol_wt, size, logK) {}
    void AddContributionToTotal(std::vector<double> *total) {
      static_cast<void>(total);
    }  // end addContributionToTotal()

    void AddContributionToDTotal(const std::vector<ac::Species>& primary_species,
                                 ac::MatrixBlock* dtotal) {
      static_cast<void>(primary_species);
      static_cast<void>(dtotal);
    }  // end addContributionToDTotal()

   protected:
   private:
  };

  //
  // most of the basic functionality comes from the parent species class.
  // don't bother to repeat those tests here...?
  //

  // make sure we can create an object with the constructor
  TEST_FIXTURE(SecondarySpeciesTest, SecondarySpecies_constructor) {
    MockSecondarySpecies secondary(name_, secondary_id_,
                                   species_names_, stoichiometry_, species_ids_,
                                   h2o_stoich_, charge_, gram_molecular_weight_,
                                   ion_size_parameter_, logK_);
    CHECK_EQUAL(secondary_id_, secondary.identifier());
  }

  // was logK set?
  TEST_FIXTURE(SecondarySpeciesTest, SecondarySpecies_constructor_logK) {
    MockSecondarySpecies secondary(name_, secondary_id_,
                                   species_names_, stoichiometry_, species_ids_,
                                   h2o_stoich_, charge_, gram_molecular_weight_,
                                   ion_size_parameter_, logK_);
    CHECK_EQUAL(logK_, secondary.logK());
  }

  // when logK is set, lnK should also be set
  TEST_FIXTURE(SecondarySpeciesTest, SecondarySpecies_constructor_lnK) {
    MockSecondarySpecies secondary(name_, secondary_id_,
                                   species_names_, stoichiometry_, species_ids_,
                                   h2o_stoich_, charge_, gram_molecular_weight_,
                                   ion_size_parameter_, logK_);
    CHECK_CLOSE(logK_ * 2.30258509299405, secondary.lnK(), 1.0e-12);
  }

  // was species_names set?
  TEST_FIXTURE(SecondarySpeciesTest, SecondarySpecies_constructor_species_names) {
    MockSecondarySpecies secondary(name_, secondary_id_,
                                   species_names_, stoichiometry_, species_ids_,
                                   h2o_stoich_, charge_, gram_molecular_weight_,
                                   ion_size_parameter_, logK_);
    for (unsigned int s = 0; s != species_names_.size(); s++) {
      CHECK_EQUAL(species_names_.at(s), secondary.species_names().at(s));
    }
  }

  // was stoichiometry set?
  TEST_FIXTURE(SecondarySpeciesTest, SecondarySpecies_constructor_stoichiometry) {
    MockSecondarySpecies secondary(name_, secondary_id_,
                                   species_names_, stoichiometry_, species_ids_,
                                   h2o_stoich_, charge_, gram_molecular_weight_,
                                   ion_size_parameter_, logK_);
    for (unsigned int s = 0; s != stoichiometry_.size(); s++) {
      CHECK_EQUAL(stoichiometry_.at(s), secondary.stoichiometry().at(s));
    }
  }

  // was species_ids set?
  TEST_FIXTURE(SecondarySpeciesTest, SecondarySpecies_constructor_species_ids) {
    MockSecondarySpecies secondary(name_, secondary_id_,
                                   species_names_, stoichiometry_, species_ids_,
                                   h2o_stoich_, charge_, gram_molecular_weight_,
                                   ion_size_parameter_, logK_);
    for (unsigned int s = 0; s != species_ids_.size(); s++) {
      CHECK_EQUAL(species_ids_.at(s), secondary.species_ids().at(s));
    }
  }

  //
  // test for invalid input
  //

  // ncomp is taken from the names array
  TEST_FIXTURE(SecondarySpeciesTest, SecondarySpecies_invalid_ncomp) {
    species_names_.clear();
    CHECK_THROW(
        MockSecondarySpecies secondary(name_, secondary_id_,
                                       species_names_, stoichiometry_, species_ids_,
                                       h2o_stoich_, charge_, gram_molecular_weight_,
                                       ion_size_parameter_, logK_),
        ac::ChemistryException);
  }

  TEST_FIXTURE(SecondarySpeciesTest, SecondarySpecies_size_names_stoichiometry) {
    species_names_.push_back("Foo");
    CHECK_THROW(
        MockSecondarySpecies secondary(name_, secondary_id_,
                                       species_names_, stoichiometry_, species_ids_,
                                       h2o_stoich_, charge_, gram_molecular_weight_,
                                       ion_size_parameter_, logK_),
        ac::ChemistryException);
  }

  TEST_FIXTURE(SecondarySpeciesTest, SecondarySpecies_size_names_stoichiometry2) {
    stoichiometry_.push_back(3.4);
    CHECK_THROW(
        MockSecondarySpecies secondary(name_, secondary_id_,
                                       species_names_, stoichiometry_, species_ids_,
                                       h2o_stoich_, charge_, gram_molecular_weight_,
                                       ion_size_parameter_, logK_),
        ac::ChemistryException);
  }

  TEST_FIXTURE(SecondarySpeciesTest, SecondarySpecies_size_names_ids) {
    species_ids_.push_back(5);
    CHECK_THROW(
        MockSecondarySpecies secondary(name_, secondary_id_,
                                       species_names_, stoichiometry_, species_ids_,
                                       h2o_stoich_, charge_, gram_molecular_weight_,
                                       ion_size_parameter_, logK_),
        ac::ChemistryException);
  }
}  // end SUITE(GeochemistryTestSecondarySpecies)
