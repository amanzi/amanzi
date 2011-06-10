/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>
#include <typeinfo>

#include <UnitTest++.h>

#include "species.hh"
#include "mineral_kinetics_factory.hh"
#include "kinetic_rate.hh"
#include "kinetic_rate_tst.hh"
#include "chemistry_exception.hh"

SUITE(GeochemistryTestsMineralKineticsFactory) {
  /*****************************************************************************
   **
   **  Test for MineralKineticsFactory.cpp
   **
   *****************************************************************************/

  /*
    if you pass the factory a valid name, it should return the correct
    kinetic rate class. If you pass it an invalid name, it should
    throw an error.
  */

  class MineralKineticsFactoryTest {
   public:

   protected:
    MineralKineticsFactoryTest();
    ~MineralKineticsFactoryTest();

    void RunTest(const std::string name);

    MineralKineticsFactory mkf_;
    KineticRate* kinetic_rate_;
    std::vector<Mineral> minerals_;

   private:
    StringTokenizer rate_data_;
    SpeciesArray species_;
    Species H_p;
    Species OH_m;
    Species Ca_pp;
    Species CO3_mm;
    Species Al_ppp;
    Species PO4_mmm;
  };
  MineralKineticsFactoryTest::MineralKineticsFactoryTest()
      : mkf_(),
      kinetic_rate_(NULL),
      rate_data_("log10_rate_constant -9.0 moles_m2_sec", ";"),
      H_p(0, "H+", 1.0, 1.0079, 9.0),
      OH_m(1, "OH-", -1.0, 17.0073, 3.5),
      Ca_pp(2, "Ca++", 2.0, 40.0780, 6.0),
      CO3_mm(3, "CO3--", -2.0, 96.0636, 4.0),
      Al_ppp(4, "Al+++", 3.0, 26.9815, 9.0),
      PO4_mmm(5, "PO4---", -3.0, 94.9714, 4.0) {
    // set primary species
    H_p.update(0.0005);
    OH_m.update(0.0015);
    Ca_pp.update(0.001);
    CO3_mm.update(0.002);
    Al_ppp.update(0.003);
    PO4_mmm.update(0.001);
    species_.push_back(H_p);
    species_.push_back(OH_m);
    species_.push_back(Ca_pp);
    species_.push_back(CO3_mm);
    species_.push_back(Al_ppp);
    species_.push_back(PO4_mmm);

    // setup a list of minerals
    SpeciesName name("Calcite");
    SpeciesId id(6);
    double h2o_stoich(0.0);
    double charge(0.0);
    double gram_molecular_weight(100.0872);
    double ion_size_parameter(0.0);
    double logK(1.8487);
    double molar_volume(36.9340);
    double specific_surface_area(0.987654);
    std::vector<SpeciesName> species_names;
    std::vector<double> stoichiometry;
    std::vector<SpeciesId> species_ids;
    species_names.clear();
    stoichiometry.clear();
    species_ids.clear();

    species_names.push_back("H+");
    stoichiometry.push_back(-1.0);
    species_ids.push_back(0);

    species_names.push_back("HCO3-");
    stoichiometry.push_back(1.0);
    species_ids.push_back(1);

    species_names.push_back("Ca++");
    stoichiometry.push_back(1.0);
    species_ids.push_back(2);
    minerals_.clear();
    // dummy mineral
    minerals_.push_back(Mineral("Foo", 0, species_names, stoichiometry, species_ids,
                                h2o_stoich, gram_molecular_weight, logK, molar_volume, specific_surface_area));
    // "real" mineral
    minerals_.push_back(Mineral(name, 1, species_names, stoichiometry, species_ids,
                                h2o_stoich, gram_molecular_weight, logK, molar_volume, specific_surface_area));
    // dummy mineral
    minerals_.push_back(Mineral("Bar", 2, species_names, stoichiometry, species_ids,
                                h2o_stoich, gram_molecular_weight, logK, molar_volume, specific_surface_area));
  }

  MineralKineticsFactoryTest::~MineralKineticsFactoryTest() {
    delete kinetic_rate_;
  }

  void MineralKineticsFactoryTest::RunTest(const std::string name) {
    kinetic_rate_ = mkf_.Create(name, rate_data_, minerals_.at(1), species_);
  }

  // use C++ RTTI to determine if the correct type of object was
  // returned from the factory, e.g. see "typeid" at
  // http:// en.wikibooks.org/wiki/C++_Programming/RTTI

  TEST_FIXTURE(MineralKineticsFactoryTest, MineralKineticsFactory_create) {
    std::string name("TST");
    RunTest(name);
    CHECK_EQUAL(typeid(KineticRateTST).name(), typeid(*kinetic_rate_).name());
  }  // end TEST_FIXTURE()

  TEST_FIXTURE(MineralKineticsFactoryTest, MineralKineticsFactory_invalid_rate) {
    std::string name("invalid-name");
    CHECK_THROW(RunTest(name), ChemistryException);
    CHECK(!kinetic_rate_);
  }  // end TEST_FIXTURE()

  TEST_FIXTURE(MineralKineticsFactoryTest, MineralKineticsFactory_verify_mineral_valid) {
    SpeciesId mineral_id = mkf_.VerifyMineralName("Calcite", minerals_);
    CHECK_EQUAL(mineral_id, 1);
  }  // end TEST_FIXTURE()

  TEST_FIXTURE(MineralKineticsFactoryTest, MineralKineticsFactory_verify_mineral_invalid) {
    CHECK_THROW(mkf_.VerifyMineralName("Pyrite", minerals_), ChemistryException);
  }  // end TEST_FIXTURE()
}  // end SUITE(GeochemistryTestMineralKineticsFactory)
