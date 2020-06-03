/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>
#include <typeinfo>

#include <UnitTest++.h>

#include "kinetic_rate.hh"
#include "chemistry_exception.hh"

SUITE(GeochemistryTestsKineticRate) {
  namespace ac = Amanzi::AmanziChemistry;
  /*****************************************************************************
   **
   **  Test for KineticRate.cpp
   **
   *****************************************************************************/

  class KineticRateTest {
   public:

   protected:
    KineticRateTest();
    ~KineticRateTest();

    ac::SpeciesArray species_;
   private:
  };  // end class KineticRateTest

  KineticRateTest::KineticRateTest() {
    // set primary species
    ac::Species H_p(0, "H+", 1.0, 1.0079, 9.0);
    H_p.update(0.0005);
    ac::Species OH_m(1, "OH-", -1.0, 17.0073, 3.5);
    OH_m.update(0.0015);
    ac::Species Ca_pp(2, "Ca++", 2.0, 40.0780, 6.0);
    Ca_pp.update(0.001);
    ac::Species CO3_mm(3, "CO3--", -2.0, 96.0636, 4.0);
    CO3_mm.update(0.002);
    ac::Species Al_ppp(4, "Al+++", 3.0, 26.9815, 9.0);
    Al_ppp.update(0.003);
    ac::Species PO4_mmm(5, "PO4---", -3.0, 94.9714, 4.0);
    PO4_mmm.update(0.001);

    species_.push_back(H_p);
    species_.push_back(OH_m);
    species_.push_back(Ca_pp);
    species_.push_back(CO3_mm);
    species_.push_back(Al_ppp);
    species_.push_back(PO4_mmm);
  }


  KineticRateTest::~KineticRateTest() {
  }

  //
  // KineticRate has a pure virtual functions, need a dummy object to
  // create and test. Only test the functions defined in this class,
  // leave the pure virtual functions for the inheriting classes.
  //
  class MockKineticRate : public ac::KineticRate {
   public:
    MockKineticRate() : ac::KineticRate() {
      set_name("abc123");
      set_identifier(456);
    };
    virtual ~MockKineticRate() {}

    void Setup(const ac::SecondarySpecies& reaction,
               const ac::StringTokenizer& reaction_data,
               const ac::SpeciesArray& primary_species) {
      static_cast<void>(reaction);
      static_cast<void>(reaction_data);
      static_cast<void>(primary_species);
    };  // end Setup()

    void Update(const ac::SpeciesArray& primary_species,
                const std::vector<ac::Mineral>& minerals) {
      static_cast<void>(primary_species);
      static_cast<void>(minerals);
    }  // end Update()

    void AddContributionToResidual(const std::vector<ac::Mineral>& minerals,
                                   const double por_den_sat_vol,
                                   std::vector<double> *residual) {
      static_cast<void>(minerals);
      static_cast<void>(por_den_sat_vol);
      static_cast<void>(residual);
    };  // end addContributionToResidual()

    void AddContributionToJacobian(const ac::SpeciesArray& primary_species,
                                   const std::vector<ac::Mineral>& minerals,
                                   const double por_den_sat_vol,
                                   ac::MatrixBlock* J) {
      static_cast<void>(primary_species);
      static_cast<void>(minerals);
      static_cast<void>(por_den_sat_vol);
      static_cast<void>(J);
    };  // end addContributionToJacobian()

    void Display(const Teuchos::Ptr<Amanzi::VerboseObject> vo) const {
      std::cout << this->name() << std::endl;
    };  // end Display()

    void ParseParameters(const ac::StringTokenizer& rate_parameters) {
      static_cast<void>(rate_parameters);
    };  // end ParseParameters()

   protected:
   private:
  };  // end MockKineticRate

  // make sure we can create an object with the constructor
  // can we set the identifier?
  TEST_FIXTURE(KineticRateTest, MockKineticRate_constructor) {
    MockKineticRate rate;
    CHECK_EQUAL(rate.identifier(), 456);
  }

  // can we set the name?
  TEST_FIXTURE(KineticRateTest, MockKineticRate_set_name) {
    MockKineticRate rate;
    CHECK_EQUAL(rate.name(), "abc123");
  }

  // can we set the debug flag?
  TEST_FIXTURE(KineticRateTest, MockKineticRate_set_debug) {
    MockKineticRate rate;
    rate.set_debug(true);
    CHECK_EQUAL(rate.debug(), true);
  }

  // does SetSpeciesIds function work?
  // TODO(bandre): what about testing for dummy/unknown species foo...?
  // TODO(bandre): is SetSpecesIds with a null out_stoichimontery used? needed?
  TEST_FIXTURE(KineticRateTest, MockKineticRate_SetSpeciesIds_test_id) {
    MockKineticRate rate;

    std::string species_type("primary");
    std::vector<ac::SpeciesName> in_names;
    in_names.push_back("Ca++");
    in_names.push_back("OH-");
    in_names.push_back("foo");
    std::vector<double> in_stoichiometry;
    in_stoichiometry.push_back(3.45);
    in_stoichiometry.push_back(0.12);
    in_stoichiometry.push_back(6.78);

    std::vector<ac::SpeciesId> out_ids;
    std::vector<double>* out_stoichiometry = NULL;

    rate.SetSpeciesIds(species_, species_type,
                       in_names, in_stoichiometry,
                       &out_ids, out_stoichiometry);
    // check that the output id's agree with the input
    std::vector<ac::SpeciesId> expeced_ids;
    expeced_ids.push_back(2);
    expeced_ids.push_back(1);
    CHECK_ARRAY_EQUAL(expeced_ids, out_ids, 2);
  }  // end TEST_FIXTURE(MockKineticRate_SetSpeciesIds_test_id)

  TEST_FIXTURE(KineticRateTest, MockKineticRate_SetSpeciesIds_test_stoich) {
    MockKineticRate rate;

    std::string species_type("primary");
    std::vector<ac::SpeciesName> in_names;
    in_names.push_back("Ca++");
    in_names.push_back("OH-");
    in_names.push_back("foo");
    std::vector<double> in_stoichiometry;
    in_stoichiometry.push_back(3.45);
    in_stoichiometry.push_back(0.12);
    in_stoichiometry.push_back(6.78);

    std::vector<ac::SpeciesId> out_ids;
    std::vector<double> out_stoichiometry;

    rate.SetSpeciesIds(species_, species_type,
                       in_names, in_stoichiometry,
                       &out_ids, &out_stoichiometry);

    // check that the output stoichimometry agrees
    std::vector<double> expected_stoich(species_.size());
    expected_stoich[1] = 0.12;
    expected_stoich[2] = 3.45;
    CHECK_ARRAY_EQUAL(expected_stoich, out_stoichiometry, species_.size());
  }  // end TEST_FIXTURE(MockKineticRate_SetSpeciesIds_test_stoich)
}  // end SUITE(GeochemistryTestKineticRate)
