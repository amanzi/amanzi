#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <typeinfo>

#include <UnitTest++.h>

#include "exceptions.hh"

#include "SecondarySpecies.hh"

SUITE(GeochemistryTestsSecondarySpecies) {
  namespace ac = Amanzi::AmanziChemistry;
  /*****************************************************************************
  **  Test for Secondary Species
  *****************************************************************************/
  class SecondarySpeciesTest {
   protected:
    SecondarySpeciesTest();

    std::string name_;
    int secondary_id_;
    double h2o_stoich_;
    int charge_;
    double gram_molecular_weight_;
    double ion_size_parameter_;
    double logK_;

    std::vector<std::string> species_names_;
    std::vector<double> stoichiometry_;
    std::vector<int> species_ids_;

    Teuchos::ParameterList plist_;
    std::vector<ac::Species> primary_species_;
  };

  SecondarySpeciesTest::SecondarySpeciesTest()
    : name_("CO3--"),
      secondary_id_(3),
      h2o_stoich_(0.0),
      charge_(-2),
      gram_molecular_weight_(60.0092),
      ion_size_parameter_(4.5),
      logK_(10.3288)
  {
    species_names_.clear();
    stoichiometry_.clear();
    species_ids_.clear();

    Teuchos::ParameterList plist;
    plist.set<int>("charge", 1)
         .set<double>("gram molecular weight", 1.0079);
    primary_species_.push_back(ac::Species(0, "H+", plist));

    species_names_.push_back("H+");
    stoichiometry_.push_back(-1.0);
    species_ids_.push_back(0);

    plist.set<int>("charge", -1);
    primary_species_.push_back(ac::Species(1, "HCO3-", plist));

    species_names_.push_back("HCO3-");
    stoichiometry_.push_back(1.0);
    species_ids_.push_back(1);

    plist_.set<int>("charge", charge_)
          .set<double>("ion size parameter", ion_size_parameter_)
          .set<double>("gram molecular weight", gram_molecular_weight_)
          .set<std::string>("reaction", "-1.0 H+  1.0 HCO3-")
          .set<double>("equilibrium constant", logK_);
  }

  //
  // SecondarySpecies has a pure virtual functions, need a dummy
  // object to create and test.
  //
  class MockSecondarySpecies : public ac::SecondarySpecies {
   public:
    MockSecondarySpecies() : SecondarySpecies() {}
    MockSecondarySpecies(const std::string name, const int id,
                         const Teuchos::ParameterList& plist,
                         const std::vector<Species>& primary_species)
        : SecondarySpecies(id, name, plist, primary_species) {}

    virtual void Update(const std::vector<Species>& primary_species, const Species& water_species) {};
    virtual void AddContributionToTotal(std::vector<double> *total) {};

    void AddContributionToDTotal(const std::vector<ac::Species>& primary_species,
                                 ac::MatrixBlock* dtotal) {};
  };

  //
  // most of the basic functionality comes from the parent species class.
  // don't bother to repeat those tests here...?
  //

  // make sure we can create an object with the constructor
  TEST_FIXTURE(SecondarySpeciesTest, SecondarySpecies_constructor) {
    MockSecondarySpecies secondary(name_, secondary_id_, plist_, primary_species_);
    CHECK_EQUAL(secondary_id_, secondary.identifier());
  }

  // was logK set?
  TEST_FIXTURE(SecondarySpeciesTest, SecondarySpecies_constructor_logK) {
    MockSecondarySpecies secondary(name_, secondary_id_, plist_, primary_species_);
    CHECK_EQUAL(logK_, secondary.logK());
  }

  // when logK is set, lnK should also be set
  TEST_FIXTURE(SecondarySpeciesTest, SecondarySpecies_constructor_lnK) {
    MockSecondarySpecies secondary(name_, secondary_id_, plist_, primary_species_);
    CHECK_CLOSE(logK_ * 2.30258509299405, secondary.lnK(), 1.0e-12);
  }

  // was species_names set?
  TEST_FIXTURE(SecondarySpeciesTest, SecondarySpecies_constructor_species_names) {
    MockSecondarySpecies secondary(name_, secondary_id_, plist_, primary_species_);
    for (unsigned int s = 0; s != species_names_.size(); s++) {
      CHECK_EQUAL(species_names_.at(s), secondary.species_names().at(s));
    }
  }

  // was stoichiometry set?
  TEST_FIXTURE(SecondarySpeciesTest, SecondarySpecies_constructor_stoichiometry) {
    MockSecondarySpecies secondary(name_, secondary_id_, plist_, primary_species_);
    for (unsigned int s = 0; s != stoichiometry_.size(); s++) {
      CHECK_EQUAL(stoichiometry_.at(s), secondary.stoichiometry().at(s));
    }
  }

  // was species_ids set?
  TEST_FIXTURE(SecondarySpeciesTest, SecondarySpecies_constructor_species_ids) {
    MockSecondarySpecies secondary(name_, secondary_id_, plist_, primary_species_);
    for (unsigned int s = 0; s != species_ids_.size(); s++) {
      CHECK_EQUAL(species_ids_.at(s), secondary.species_ids().at(s));
    }
  }

  //
  // test for invalid input
  //

  // ncomp is taken from the names array
  TEST_FIXTURE(SecondarySpeciesTest, SecondarySpecies_invalid_ncomp) {
    primary_species_.pop_back();
    CHECK_THROW(
        MockSecondarySpecies secondary(name_, secondary_id_, plist_, primary_species_),
        Exceptions::Amanzi_exception);
  }

  TEST_FIXTURE(SecondarySpeciesTest, SecondarySpecies_size_names_stoichiometry) {
    primary_species_[1] = primary_species_[0];
    CHECK_THROW(
        MockSecondarySpecies secondary(name_, secondary_id_, plist_, primary_species_),
        Exceptions::Amanzi_exception);
  }
}
