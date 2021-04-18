#include <cstdlib>
#include <cmath>
#include <iostream>
#include <typeinfo>
#include <vector>

#include <UnitTest++.h>

#include "exceptions.hh"

#include "ActivityModel.hh"
#include "ActivityModelFactory.hh"
#include "ActivityModelPitzerHWM.hh"
#include "AqueousEquilibriumComplex.hh"
#include "Species.hh"

SUITE(TestPitzer) {
  namespace ac = Amanzi::AmanziChemistry;

  class PitzerTest {
   public:
    PitzerTest();
    ~PitzerTest() {};

    void StorePrimaries();

    ac::ActivityModel::ActivityModelParameters parameters;
    ac::ActivityModelFactory amfac_;
    ac::ActivityModel* am_;
    std::vector<ac::Species> sp_;
    std::vector<ac::AqueousEquilibriumComplex> aqx_;

    ac::Species H;
    ac::Species OH;
    ac::Species Cl;
    ac::Species Na;
    ac::Species K;
    ac::Species Ca;
    ac::Species Mg;
    ac::Species CO3;
    ac::Species CO2;
    ac::Species HCO3;
    ac::Species MgOH;
    ac::Species MgCO3;
    ac::Species CaCO3;
    ac::Species H2O;
    ac::Species HSO4;
    ac::Species SO4;
    ac::Species Br;

    Teuchos::RCP<Amanzi::VerboseObject> vo_;
  };

  PitzerTest::PitzerTest() 
  {
    Teuchos::ParameterList plist;
    plist.set<int>("charge", 1)
         .set<double>("ion size parameters", 9.0)
         .set<double>("gram molecular weight", 1.0079);
    H = ac::Species(0, "H+", plist);

    plist.set<int>("charge", -1);
    OH = ac::Species(1, "OH-", plist);
    Cl = ac::Species(2, "Cl-", plist);

    plist.set<int>("charge", 1);
    Na = ac::Species(3, "Na+", plist);
    K = ac::Species(4, "K+", plist);

    plist.set<int>("charge", 2);
    Ca = ac::Species(5, "Ca+2", plist);
    Mg = ac::Species(6, "Mg+2", plist);

    plist.set<int>("charge", -2);
    CO3 = ac::Species(7, "CO3-2", plist);

    plist.set<int>("charge", 0);
    CO2 = ac::Species(8, "CO2", plist);

    plist.set<int>("charge", -1);
    HCO3 = ac::Species(9, "HCO3-", plist);

    plist.set<int>("charge", 1);
    MgOH = ac::Species(10, "MgOH+", plist);

    plist.set<int>("charge", 0);
    MgCO3 = ac::Species(11, "MgCO3", plist),
    CaCO3 = ac::Species(12, "CaCO3", plist),
    H2O = ac::Species(13, "H2O", plist);

    plist.set<int>("charge", -1);
    HSO4 = ac::Species(13, "HSO4-", plist);

    plist.set<int>("charge", -2);
    SO4 = ac::Species(14, "SO4-2", plist);

    plist.set<int>("charge", -1);
    Br = ac::Species(16, "Br-", plist);

    // parameters.database_filename = "phreeqc_pitzer";
    parameters.database_filename = "test/chemistry_pitzer";
    parameters.pitzer_jfunction = "pitzer1975";

    vo_ = Teuchos::rcp(new Amanzi::VerboseObject("Chemistry", plist));
  }

  void PitzerTest::StorePrimaries() {
    sp_.clear();
    sp_.push_back(Cl);
    sp_.push_back(Na);
    sp_.push_back(H);
    sp_.push_back(K);
    sp_.push_back(H2O);
    sp_.push_back(Ca);
    sp_.push_back(Mg);
    sp_.push_back(CO3);
    sp_.push_back(HCO3);
    sp_.push_back(CO2);
    sp_.push_back(MgCO3);
    sp_.push_back(CaCO3);
    sp_.push_back(MgOH);
    sp_.push_back(OH);
  }

  /*!
    @class Amanzi::AmanziChemistry::unit_tests::ActivityModelPitzer::TestComputeActivityCoeff_System_I

    @brief TestComputeActivityCoeff_System_I

    @details Test the calculation of activity coefficients.
    Results are compared with PHREEQC.

    @test ActivityModelPitzer::EvaluateVector()
  */
  TEST_FIXTURE(PitzerTest, TestComputeActivityCoeff_System_I) {
    H.update(4.776e-08);
    OH.update(2.570e-07);
    Cl.update(3.0e0);
    Na.update(3.0e0);
    Ca.update(9.516e-02);
    Mg.update(9.711e-02);
    CO3.update(3.339e-03);
    HCO3.update(2.760e-01);
    CO2.update(1.299e-02);
    MgOH.update(1.131e-6);
    MgCO3.update(2.885e-03);
    CaCO3.update(4.838e-3);
    K.update(0.1e0);
    H2O.update(1.0);
    StorePrimaries();
    aqx_.clear();
    std::vector<double> gamma(sp_.size(), 1.0);

    am_ = amfac_.Create("pitzer-hwm", parameters, sp_, aqx_, vo_.ptr());

    if (vo_->getVerbLevel() >= Teuchos::VERB_EXTREME) {
      am_->Display();
    }
    am_->CalculateActivityCoefficients(&sp_, &aqx_, &H2O);
    double actw = log10(H2O.act_coef());
    for (int i = 0; i < sp_.size(); i++) {
      gamma[i] = log10(sp_[i].act_coef());
    }
    // std::cout << "Testing coeff. 1" << std::endl;
    // Results are compared with PHREEQC
    CHECK_CLOSE(-0.243, gamma[0], 1.0e-2);  // Cl-
    CHECK_CLOSE(-0.027, gamma[1], 1.0e-2);  // Na+
    CHECK_CLOSE(0.321, gamma[2], 1.0e-2);  // H+
    CHECK_CLOSE(-0.177, gamma[3], 1.0e-2);  // K+
    CHECK_CLOSE(-0.055, gamma[4], 1.0e-2);  // H2O
    CHECK_CLOSE(-0.109, gamma[5], 1.0e-2);  // Ca
    CHECK_CLOSE(-0.119, gamma[6], 1.0e-2);  // Mg
    CHECK_CLOSE(-1.859, gamma[7], 1.0e-2);  // CO3
    CHECK_CLOSE(-0.05, gamma[12], 1.0e-2);  // MgOH
    CHECK_CLOSE(0.283, gamma[9], 1.0e-2);  // CO2
    CHECK_CLOSE(-0.463, gamma[13], 1.0e-2);  // OH
    CHECK_CLOSE(-0.437, gamma[8], 1.0e-2);  // HCO3
    CHECK_CLOSE(-0.055, actw, 1.0e-2);  // H2O
  }

  /*!
    @class Amanzi::AmanziChemistry::unit_tests::ActivityModelPitzer::TestComputeActivityCoeff_System_II

    @brief TestComputeActivityCoeff_System_II

    @details Test the calculation of activity coefficients.
    Results are compared with PHREEQC.

    @test ActivityModelPitzer::EvaluateVector()
  */
  TEST_FIXTURE(PitzerTest, TestComputeActivityCoeff_System_II) {
    H.update(4.98e-08);
    OH.update(2.574e-07);
    Cl.update(3.0e0);
    Na.update(3.0e0);
    Ca.update(9.545e-02);
    Mg.update(9.717e-02);
    CO3.update(3.382e-03);
    HCO3.update(2.764e-01);
    CO2.update(1.283e-02);
    MgOH.update(1.224e-6);
    MgCO3.update(2.829e-03);
    CaCO3.update(4.55e-3);
    HSO4.update(3.0e-8);
    SO4.update(0.1e0);
    K.update(0.1e0);
    H2O.update(1.0);
    StorePrimaries();
    sp_.push_back(HSO4);
    sp_.push_back(SO4);
    aqx_.clear();
    std::vector<double> gamma(sp_.size(), 1.0);

    am_ = amfac_.Create("pitzer-hwm", parameters, sp_, aqx_, vo_.ptr());
    am_->Display();

    am_->CalculateActivityCoefficients(&sp_, &aqx_, &H2O);
    double actw = log10(H2O.act_coef());
    for (int i = 0; i < sp_.size(); i++) {
      gamma[i] = log10(sp_[i].act_coef());
      // std::cout << sp_[i].name() << "  " <<gamma[i] << std::endl;
    }
    // Results are compared with PHREEQC
    CHECK_CLOSE(-0.241, gamma[0], 1.0e-2);  // Cl-
    CHECK_CLOSE(-0.039, gamma[1], 1.0e-2);  // Na+
    CHECK_CLOSE(0.303, gamma[2], 1.0e-2);  // H+
    CHECK_CLOSE(-0.192, gamma[3], 1.0e-2);  // K+
    CHECK_CLOSE(-0.055, gamma[4], 1.0e-2);  // H2O
    CHECK_CLOSE(-0.14, gamma[5], 1.0e-2);  // Ca
    CHECK_CLOSE(-0.131, gamma[6], 1.0e-2);  // Mg
    CHECK_CLOSE(-1.862, gamma[7], 1.0e-2);  // CO3
    CHECK_CLOSE(-0.096, gamma[12], 1.0e-2);  // MgOH
    CHECK_CLOSE(0.291, gamma[9], 1.0e-2);  // CO2
    CHECK_CLOSE(-0.464, gamma[13], 1.0e-2);  // OH
    CHECK_CLOSE(-0.435, gamma[8], 1.0e-2);  // HCO3
    CHECK_CLOSE(-0.321, gamma[14], 1.0e-2);  // HSO4
    CHECK_CLOSE(-1.823, gamma[15], 1.0e-2);  // SO4
    CHECK_CLOSE(-0.055, actw, 1.0e-2);  // H2O
  }

  /*!
    @class Amanzi::AmanziChemistry::unit_tests::ActivityModelPitzer::TestComputeActivityCoeff_System_III

    @brief TestComputeActivityCoeff_System_III

    @details Test the calculation of activity coefficients.
    This a solution in equilibrium with halite and gypsum.
    Results are compared with PHREEQC.

    @test ActivityModelPitzer::EvaluateVector()
  */
  TEST_FIXTURE(PitzerTest, TestComputeActivityCoeff_System_III) {
    H.update(9.326e-010);
    OH.update(4.644e-006);
    Cl.update(5.948e+000);
    Na.update(5.178e0);
    Ca.update(1.108e-02);
    Mg.update(6.282e-001);
    CO3.update(1.100e-006);
    HCO3.update(3.402e-006);
    CO2.update(3.906e-009);
    MgOH.update(9.303e-004);
    MgCO3.update(1.242e-005);
    CaCO3.update(1.964e-007);
    HSO4.update(2.052e-009);
    SO4.update(2.250e-001);
    K.update(1.208e-001);
    Br.update(1.0e-9);
    H2O.update(1.0);
    aqx_.clear();
    StorePrimaries();
    sp_.push_back(HSO4);
    sp_.push_back(SO4);
    sp_.push_back(Br);
    std::vector<double> gamma(sp_.size(), 1.0);

    am_ = amfac_.Create("pitzer-hwm", parameters, sp_, aqx_, vo_.ptr());
    am_->Display();

    am_->CalculateActivityCoefficients(&sp_, &aqx_, &H2O);
    double actw = log10(H2O.act_coef());
    for (int i = 0; i < sp_.size(); i++) {
      gamma[i] = log10(sp_[i].act_coef());
      // std::cout << sp_[i].name() << "  " << gamma[i] << std::endl;
    }
    CHECK_CLOSE(-0.192, gamma[0], 1.0e-2);  // Cl-
    CHECK_CLOSE(0.274, gamma[1], 1.0e-2);  // Na+
    CHECK_CLOSE(0.931, gamma[2], 1.0e-2);  // H+
    CHECK_CLOSE(-0.065, gamma[3], 1.0e-2);  // K+
    CHECK_CLOSE(-0.13, gamma[4], 1.0e-2);  // H2O
    CHECK_CLOSE(0.563, gamma[5], 1.0e-2);  // Ca
    CHECK_CLOSE(0.834, gamma[6], 1.0e-2);  // Mg
    CHECK_CLOSE(-2.507, gamma[7], 1.0e-2);  // CO3
    CHECK_CLOSE(-0.176, gamma[12], 1.0e-2);  // MgOH
    CHECK_CLOSE(0.55, gamma[9], 1.0e-2);  // CO2
    CHECK_CLOSE(-0.695, gamma[13], 1.0e-2);  // OH
    CHECK_CLOSE(-0.758, gamma[8], 1.0e-2);  // HCO3
    CHECK_CLOSE(-0.361, gamma[14], 1.0e-2);  // HSO4
    CHECK_CLOSE(-2.28, gamma[15], 1.0e-2);  // SO4
    CHECK_CLOSE(-0.038, gamma[16], 1.0e-2);  // Br
    CHECK_CLOSE(-0.13, actw, 1.0e-2);  // H2O
  }

  /*!
    @class Amanzi::AmanziChemistry::unit_tests::ActivityModelPitzer::TestInvalidActivityModel

    @brief TestInvalidActivityModel

    @details Test that a chemistry exception is thrown when an invalid activity model
    is provided.

    @test ActivityModelPitzer::Create()
  */
  TEST_FIXTURE(PitzerTest, TestInvalidActivityModel) {
    Teuchos::ParameterList plist;
    plist.set<int>("charge", -1)
         .set<double>("gram molecular weight", 0.0);
    ac::Species Cl(0, "Cl-", plist);

    plist.set<int>("charge", 1);
    ac::Species Na(1, "Na+", plist);

    Cl.update(1.0);
    Na.update(1.0);
    aqx_.clear();
    sp_.clear();
    sp_.push_back(Cl);
    sp_.push_back(Na);
    CHECK_THROW(am_ = amfac_.Create("invalid activity model", parameters, sp_, aqx_, vo_.ptr()), Exceptions::Amanzi_exception);
  }

  /*!
    @class Amanzi::AmanziChemistry::unit_tests::ActivityModelPitzer::TestZeroConcentrations

    @brief TestZeroConcentrations

    @details Test that a chemistry exception is thrown when zero concentrations are provided.

    @test ActivityModelPitzer::EvaluateVector()
  */
  TEST_FIXTURE(PitzerTest, TestZeroConcentrations) {
    Teuchos::ParameterList plist;
    plist.set<int>("charge", -1)
         .set<double>("gram molecular weight", 0.0);
    ac::Species Cl(0, "Cl-", plist);

    plist.set<int>("charge", 1);
    ac::Species Na(1, "Na+", plist);

    plist.set<int>("charge", 0);
    ac::Species H2O(2, "H2O", plist);

    Cl.update(0.0);
    Na.update(0.0);
    H2O.update(0.0);
    aqx_.clear();
    sp_.clear();
    sp_.push_back(Cl);
    sp_.push_back(Na);
    std::vector<double> gamma(sp_.size(), 1.0);

    am_ = amfac_.Create("pitzer-hwm", parameters, sp_, aqx_, vo_.ptr());
    am_->Display();
    CHECK_THROW(am_->CalculateActivityCoefficients(&sp_, &aqx_, &H2O), Exceptions::Amanzi_exception);
  }

  /*!
    @class Amanzi::AmanziChemistry::unit_tests::ActivityModelPitzer::TestZeroConcentrations

    @brief TestNumberSpecies

    @details Test the number of species.

    @test ActivityModelPitzer::EvaluateVector()
  */
  TEST_FIXTURE(PitzerTest, TestNumberSpecies) {
    Teuchos::ParameterList plist;
    plist.set<int>("charge", -1)
         .set<double>("gram molecular weight", 0.0);
    ac::Species Cl(0, "Cl-", plist);

    plist.set<int>("charge", 1);
    ac::Species Na(1, "Na+", plist);

    plist.set<int>("charge", 0);
    ac::Species H2O(2, "H2O", plist);

    plist.set<int>("charge", 2);
    ac::Species Ca(3, "Ca+2", plist);

    Cl.update(3.0);
    Na.update(1.0);
    Ca.update(1.0);
    H2O.update(1.0);
    aqx_.clear();
    sp_.clear();
    sp_.push_back(Cl);
    sp_.push_back(Na);
    sp_.push_back(Ca);
    std::vector<double> gamma(sp_.size(), 1.0);
    am_ = amfac_.Create("pitzer-hwm", parameters, sp_, aqx_, vo_.ptr());
    am_->Display();

    sp_.pop_back();
    CHECK_THROW(am_->CalculateActivityCoefficients(&sp_, &aqx_, &H2O), Exceptions::Amanzi_exception);
  }
}  // end SUITE(TestPitzer)

