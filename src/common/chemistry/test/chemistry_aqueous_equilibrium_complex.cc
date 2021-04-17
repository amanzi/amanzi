#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <typeinfo>
#include "Species.hh"

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

    std::string name_;
    int id_;
    double h2o_stoich_;
    double charge_;
    double gram_molecular_weight_;
    double ion_size_parameter_;
    double logK_;

    std::vector<std::string> species_names_;
    std::vector<double> stoichiometry_;
    std::vector<int> species_ids_;

    ac::SpeciesArray primary_species_;
    ac::Species water_;

    Teuchos::ParameterList plist_;
  };

  AqueousEquilibriumComplexTest::AqueousEquilibriumComplexTest()
    : name_("CO3--"),
      id_(3),
      h2o_stoich_(0.0),
      charge_(-2),
      gram_molecular_weight_(60.0092),
      ion_size_parameter_(4.5),
      logK_(10.3288)
  {
    Teuchos::ParameterList plist;
    plist.set<int>("charge", 0)
         .set<double>("gram molecular weight", 18.0);
    water_ = ac::Species(0,"H2O", plist);

    species_names_.clear();
    stoichiometry_.clear();
    species_ids_.clear();

    species_names_.push_back("H+");
    stoichiometry_.push_back(-1.0);
    species_ids_.push_back(0);

    plist.set<int>("charge", 1)
         .set<double>("gram molecular weight", 1.0079)
         .set<double>("ion size parameter", 9.0);
    int id = 0;
    std::string name = "H+";
    ac::Species H_p = ac::Species(id, name, plist);
    H_p.update(9.0e-4);

    species_names_.push_back("HCO3-");
    stoichiometry_.push_back(1.0);
    species_ids_.push_back(1);

    plist.set<int>("charge", -1)
         .set<double>("gram molecular weight", 61.0171)
         .set<double>("ion size parameter", 4.0);
    id = 1;
    name = "HCO3-";
    ac::Species HCO3_m = ac::Species(id, name, plist);
    HCO3_m.update(1.2e-3);

    primary_species_.clear();
    primary_species_.push_back(H_p);
    primary_species_.push_back(HCO3_m);

    plist_.set<double>("ion size parameter", ion_size_parameter_)
          .set<int>("charge", charge_)
          .set<double>("gram molecular weight", gram_molecular_weight_)
          .set<std::string>("reaction", "-1.0 H+  1.0 HCO3-")
          .set<double>("equilibrium constant", logK_);
  }

  AqueousEquilibriumComplexTest::~AqueousEquilibriumComplexTest() {};

  //
  // most of the basic functionality comes from the parent SecondarySpecies class.
  // only test the virtual methods from the public interface...?
  //

  // make sure we can create an object with the constructor
  TEST_FIXTURE(AqueousEquilibriumComplexTest, AqueousEquilibriumComplex_constructor) {
    ac::AqueousEquilibriumComplex aec(id_, name_, plist_, primary_species_);
    CHECK_EQUAL(id_, aec.identifier());
  }

  TEST_FIXTURE(AqueousEquilibriumComplexTest, AqueousEquilibriumComplex_lnK) {
    ac::AqueousEquilibriumComplex aec(id_, name_, plist_, primary_species_);
    CHECK_CLOSE(std::log(std::pow(10.0, logK_)), aec.lnK(), 1.0e-10);
  }

  // public methods
  TEST_FIXTURE(AqueousEquilibriumComplexTest, AqueousEquilibriumComplex_Update) {
    ac::AqueousEquilibriumComplex aec(id_, name_, plist_, primary_species_);
    aec.Update(primary_species_, water_);
    CHECK_CLOSE(aec.lnQK(), -23.4952588360233, 1.0e-10);
  }

  TEST_FIXTURE(AqueousEquilibriumComplexTest, AqueousEquilibriumComplex_AddContributionToTotal) {
    ac::AqueousEquilibriumComplex aec(id_, name_, plist_, primary_species_);

    std::vector<double> expected(aec.ncomp(), 0.0);
    expected[0] = -6.25372437900718e-11;
    expected[1] = 6.25372437900718e-11;

    std::vector<double> total(aec.ncomp(), 0.0);
    aec.Update(primary_species_, water_);
    aec.AddContributionToTotal(&total);
    CHECK_ARRAY_CLOSE(total, expected, total.size(), 1.0e-15);
  }

  TEST_FIXTURE(AqueousEquilibriumComplexTest, AqueousEquilibriumComplex_AddContributionToDTotal) {
    ac::AqueousEquilibriumComplex aec(id_, name_, plist_, primary_species_);

    ac::MatrixBlock expected(aec.ncomp());
    expected.Zero();
    expected(0, 0) = 6.9485826433413e-08;
    expected(0, 1) =-5.21143698250598e-08;
    expected(1, 0) =-6.9485826433413e-08;
    expected(1, 1) = 5.21143698250598e-08;
    double* e = expected.GetValues();

    ac::MatrixBlock dtotal(aec.ncomp());
    dtotal.Zero();
    aec.Update(primary_species_, water_);
    aec.AddContributionToDTotal(primary_species_, &dtotal);
    double* dt = dtotal.GetValues();
    CHECK_ARRAY_CLOSE(dt, e, aec.ncomp(), 1.0e-15);
  }
}  // end SUITE(GeochemistryTestAqueousEquilibriumComplex)
