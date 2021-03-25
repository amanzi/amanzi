/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>
#include <typeinfo>
#include "species.hh"
#include <UnitTest++.h>

#include "aqueous_equilibrium_complex.hh"
#include "chemistry_exception.hh"
#include "matrix_block.hh"

SUITE(GeochemistryTestsAqueousEquilibriumComplex) {
  namespace ac = Amanzi::AmanziChemistry;
  /*****************************************************************************
  **  Test for Aqueous Equilibrium Complex
  *****************************************************************************/

  class AqueousEquilibriumComplexTest {
   public:

   protected:
    AqueousEquilibriumComplexTest();
    ~AqueousEquilibriumComplexTest();

    ac::SpeciesName name_;
    ac::SpeciesId id_;
    double h2o_stoich_;
    double charge_;
    double gram_molecular_weight_;
    double ion_size_parameter_;
    double logK_;

    std::vector<ac::SpeciesName> species_names_;
    std::vector<double> stoichiometry_;
    std::vector<ac::SpeciesId> species_ids_;

    ac::SpeciesArray primarySpecies_;
    ac::Species water_;

   private:
  };

  AqueousEquilibriumComplexTest::AqueousEquilibriumComplexTest()
    : name_("CO3--"),
      id_(3),
      h2o_stoich_(0.0),
      charge_(-2),
      gram_molecular_weight_(60.0092),
      ion_size_parameter_(4.5),
      logK_(10.3288),
      water_(0,"H2O", 0.0, 18.00, 0.0) {
    species_names_.clear();
    stoichiometry_.clear();
    species_ids_.clear();

    species_names_.push_back("H+");
    stoichiometry_.push_back(-1.0);
    species_ids_.push_back(0);

    species_names_.push_back("HCO3-");
    stoichiometry_.push_back(1.0);
    species_ids_.push_back(1);
    ac::SpeciesId id = 0;
    ac::SpeciesName name = "H+";
    ac::Species H_p = ac::Species(id, name, 1.0, 1.0079, 9.0);
    H_p.update(9.0e-4);
    id = 1;
    name = "HCO3-";
    ac::Species HCO3_m = ac::Species(id, name, -1.0, 61.0171, 4.0);
    HCO3_m.update(1.2e-3);

    primarySpecies_.clear();
    primarySpecies_.push_back(H_p);
    primarySpecies_.push_back(HCO3_m);
  }

  AqueousEquilibriumComplexTest::~AqueousEquilibriumComplexTest() {};

  //
  // most of the basic functionality comes from the parent SecondarySpecies class.
  // only test the virtual methods from the public interface...?
  //

  // make sure we can create an object with the constructor
  TEST_FIXTURE(AqueousEquilibriumComplexTest, AqueousEquilibriumComplex_constructor) {
    ac::AqueousEquilibriumComplex aec(name_, id_,
                                      species_names_, stoichiometry_, species_ids_,
                                      h2o_stoich_, charge_, gram_molecular_weight_,
                                      ion_size_parameter_, logK_);
    CHECK_EQUAL(id_, aec.identifier());
  }

  TEST_FIXTURE(AqueousEquilibriumComplexTest, AqueousEquilibriumComplex_lnK) {
    ac::AqueousEquilibriumComplex aec(name_, id_,
                                      species_names_, stoichiometry_, species_ids_,
                                      h2o_stoich_, charge_, gram_molecular_weight_,
                                      ion_size_parameter_, logK_);
    CHECK_CLOSE(std::log(std::pow(10.0, logK_)), aec.lnK(), 1.0e-10);
  }

  // public methods
  TEST_FIXTURE(AqueousEquilibriumComplexTest, AqueousEquilibriumComplex_Update) {
    ac::AqueousEquilibriumComplex aec(name_, id_,
                                      species_names_, stoichiometry_, species_ids_,
                                      h2o_stoich_, charge_, gram_molecular_weight_,
                                      ion_size_parameter_, logK_);
    aec.Update(primarySpecies_,water_);
    CHECK_CLOSE(aec.lnQK(), -23.4952588360233, 1.0e-10);
  }

  TEST_FIXTURE(AqueousEquilibriumComplexTest, AqueousEquilibriumComplex_AddContributionToTotal) {
    ac::AqueousEquilibriumComplex aec(name_, id_,
                                      species_names_, stoichiometry_, species_ids_,
                                      h2o_stoich_, charge_, gram_molecular_weight_,
                                      ion_size_parameter_, logK_);

    std::vector<double> expected(aec.ncomp(), 0.0);
    expected[0] = -6.25372437900718e-11;
    expected[1] = 6.25372437900718e-11;

    std::vector<double> total(aec.ncomp(), 0.0);
    aec.Update(primarySpecies_,water_);
    aec.AddContributionToTotal(&total);
    CHECK_ARRAY_CLOSE(total, expected, total.size(), 1.0e-15);
  }

  TEST_FIXTURE(AqueousEquilibriumComplexTest, AqueousEquilibriumComplex_AddContributionToDTotal) {
    ac::AqueousEquilibriumComplex aec(name_, id_,
                                      species_names_, stoichiometry_, species_ids_,
                                      h2o_stoich_, charge_, gram_molecular_weight_,
                                      ion_size_parameter_, logK_);

    ac::MatrixBlock expected(aec.ncomp());
    expected.Zero();
    expected(0, 0) = 6.9485826433413e-08;
    expected(0, 1) =-5.21143698250598e-08;
    expected(1, 0) =-6.9485826433413e-08;
    expected(1, 1) = 5.21143698250598e-08;
    double* e = expected.GetValues();

    ac::MatrixBlock dtotal(aec.ncomp());
    dtotal.Zero();
    aec.Update(primarySpecies_,water_);
    aec.AddContributionToDTotal(primarySpecies_, &dtotal);
    double* dt = dtotal.GetValues();
    CHECK_ARRAY_CLOSE(dt, e, aec.ncomp(), 1.0e-15);
  }
}  // end SUITE(GeochemistryTestAqueousEquilibriumComplex)
