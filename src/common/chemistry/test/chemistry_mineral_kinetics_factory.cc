/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <typeinfo>

#include <UnitTest++.h>

#include "exceptions.hh"

#include "KineticRate.hh"
#include "KineticRateTST.hh"
#include "KineticRateFactory.hh"
#include "Species.hh"

SUITE(GeochemistryTestsKineticRateFactory) {
  namespace ac = Amanzi::AmanziChemistry;

  /*****************************************************************************
  **  Test for KineticRateFactory.cpp
  *****************************************************************************/

  /*
    if you pass the factory a valid name, it should return the correct
    kinetic rate class. If you pass it an invalid name, it should
    throw an error.
  */

  class KineticRateFactoryTest {
   protected:
    KineticRateFactoryTest();
    ~KineticRateFactoryTest();

    void RunTest(const std::string& name);

    ac::KineticRateFactory mkf_;
    ac::KineticRate* kinetic_rate_;
    std::vector<ac::Mineral> minerals_;

   private:
    ac::SpeciesArray species_;
    ac::Species H_p;
    ac::Species OH_m;
    ac::Species Ca_pp;
    ac::Species CO3_mm;
    ac::Species Al_ppp;
    ac::Species PO4_mmm;
  };

  KineticRateFactoryTest::KineticRateFactoryTest()
      : mkf_(),
      kinetic_rate_(NULL)
  {
    Teuchos::ParameterList plist;
    plist.set<int>("charge", 1)
         .set<double>("gram molecular weight", 1.0079)
         .set<double>("ion size parameter", 9.0);
    H_p = ac::Species(0, "H+", plist);

    plist.set<int>("charge", -1)
         .set<double>("gram molecular weight", 17.0073)
         .set<double>("ion size parameter", 3.5);
    OH_m = ac::Species(1, "OH-", plist);

    plist.set<int>("charge", 2)
         .set<double>("gram molecular weight", 40.0780)
         .set<double>("ion size parameter", 6.0);
    Ca_pp = ac::Species(2, "Ca++", plist);

    plist.set<int>("charge", -2)
         .set<double>("gram molecular weight", 96.0636)
         .set<double>("ion size parameter", 4.0);
    CO3_mm = ac::Species(3, "CO3--", plist);

    plist.set<int>("charge", 3)
         .set<double>("gram molecular weight", 26.9815)
         .set<double>("ion size parameter", 9.0);
    Al_ppp = ac::Species(4, "Al+++", plist);

    plist.set<int>("charge", -3)
         .set<double>("gram molecular weight", 94.9714)
         .set<double>("ion size parameter", 4.0);
    PO4_mmm = ac::Species(5, "PO4---", plist);

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
    std::string name("Calcite");
    double h2o_stoich(0.0);
    double gram_molecular_weight(100.0872);
    double logK(1.8487);
    double molar_volume(36.9340);
    double specific_surface_area(0.987654);
    std::vector<std::string> species_names;
    std::vector<double> stoichiometry;
    std::vector<int> species_ids;
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

    std::vector<ac::Species> primary_species;

    plist.set<int>("charge", 0)
         .set<double>("gram molecular weight", 100.0)
         .set<double>("equilibrium constant", 10.0)
         .set<std::string>("reaction", "2.0 H2O  -3.0 H+  1.0 HCO3-")
         .set<double>("specific surface area", 1.0)
         .set<double>("molar volume", 1.0);

    primary_species.push_back(ac::Species(0, "H+", plist));
    primary_species.push_back(ac::Species(1, "HCO3-", plist));

    // dummy mineral
    minerals_.push_back(ac::Mineral(0, "Foo", plist, primary_species));
    // "real" mineral
    minerals_.push_back(ac::Mineral(1, name, plist, primary_species));
    // dummy mineral
    minerals_.push_back(ac::Mineral(2, "Bar", plist, primary_species));
  }

  KineticRateFactoryTest::~KineticRateFactoryTest() {
    delete kinetic_rate_;
  }

  void KineticRateFactoryTest::RunTest(const std::string& name) {
    Teuchos::ParameterList plist;
    plist.set<std::string>("rate model", name)
         .set<double>("rate constant", -9.0)
         .set<std::string>("modifiers", "");
    kinetic_rate_ = mkf_.Create(plist, minerals_.at(1), species_);
  }

  TEST_FIXTURE(KineticRateFactoryTest, KineticRateFactory_create) {
    std::string name("TST");
    RunTest(name);
    CHECK_EQUAL(typeid(ac::KineticRateTST).name(), typeid(*kinetic_rate_).name());
  }

  TEST_FIXTURE(KineticRateFactoryTest, KineticRateFactory_invalid_rate) {
    std::string name("invalid-name");
    CHECK_THROW(RunTest(name), Exceptions::Amanzi_exception);
    CHECK(!kinetic_rate_);
  }
}

