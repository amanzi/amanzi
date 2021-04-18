/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <typeinfo>

#include <UnitTest++.h>

#include "KineticRate.hh"

SUITE(GeochemistryTestsKineticRate) {
  namespace ac = Amanzi::AmanziChemistry;

  /*****************************************************************************
  **  Test for KineticRate.cpp
  *****************************************************************************/

  class KineticRateTest {
   protected:
    KineticRateTest();
    ~KineticRateTest();

    ac::SpeciesArray species_;
  };

  KineticRateTest::KineticRateTest() {
    // set primary species
    Teuchos::ParameterList plist;
    plist.set<int>("charge", 1)
         .set<double>("gram molecular weight", 1.0079)
         .set<double>("ion size parameter", 9.0);
    ac::Species H_p(0, "H+", plist);
    H_p.update(0.0005);

    plist.set<int>("charge", -1)
         .set<double>("gram molecular weight", 17.0073)
         .set<double>("ion size parameter", 3.5);
    ac::Species OH_m(1, "OH-", plist);
    OH_m.update(0.0015);

    plist.set<int>("charge", 2)
         .set<double>("gram molecular weight", 40.0780)
         .set<double>("ion size parameter", 6.0);
    ac::Species Ca_pp(2, "Ca++", plist);
    Ca_pp.update(0.001);

    plist.set<int>("charge", -2)
         .set<double>("gram molecular weight", 96.0636)
         .set<double>("ion size parameter", 4.0);
    ac::Species CO3_mm(3, "CO3--", plist);
    CO3_mm.update(0.002);

    plist.set<int>("charge", 3)
         .set<double>("gram molecular weight", 26.9815)
         .set<double>("ion size parameter", 9.0);
    ac::Species Al_ppp(4, "Al+++", plist);
    Al_ppp.update(0.003);

    plist.set<int>("charge", -3)
         .set<double>("gram molecular weight", 94.9714)
         .set<double>("ion size parameter", 4.0);
    ac::Species PO4_mmm(5, "PO4---", plist);
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
    virtual ~MockKineticRate() {};

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
    std::vector<std::string> in_names;
    in_names.push_back("Ca++");
    in_names.push_back("OH-");
    in_names.push_back("foo");
    std::vector<double> in_stoichiometry;
    in_stoichiometry.push_back(3.45);
    in_stoichiometry.push_back(0.12);
    in_stoichiometry.push_back(6.78);

    std::vector<int> out_ids;
    std::vector<double>* out_stoichiometry = NULL;

    rate.SetSpeciesIds(species_, species_type,
                       in_names, in_stoichiometry,
                       &out_ids, out_stoichiometry);
    // check that the output id's agree with the input
    std::vector<int> expeced_ids;
    expeced_ids.push_back(2);
    expeced_ids.push_back(1);
    CHECK_ARRAY_EQUAL(expeced_ids, out_ids, 2);
  }  // end TEST_FIXTURE(MockKineticRate_SetSpeciesIds_test_id)

  TEST_FIXTURE(KineticRateTest, MockKineticRate_SetSpeciesIds_test_stoich) {
    MockKineticRate rate;

    std::string species_type("primary");
    std::vector<std::string> in_names;
    in_names.push_back("Ca++");
    in_names.push_back("OH-");
    in_names.push_back("foo");
    std::vector<double> in_stoichiometry;
    in_stoichiometry.push_back(3.45);
    in_stoichiometry.push_back(0.12);
    in_stoichiometry.push_back(6.78);

    std::vector<int> out_ids;
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
