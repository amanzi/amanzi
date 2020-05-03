/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>
#include <typeinfo>
#include "species.hh"
#include <UnitTest++.h>

#include "mineral.hh"
#include "chemistry_exception.hh"

SUITE(GeochemistryTestsMineral) {
  namespace ac = Amanzi::AmanziChemistry;
  /*****************************************************************************
   **
   **  Test for Mineral.cpp
   **
   *****************************************************************************/

  class MineralTest {
   public:

   protected:
    MineralTest();
    ~MineralTest();

    ac::SpeciesName name_;
    ac::SpeciesId id_;
    double h2o_stoich_;
    double charge_;
    double gram_molecular_weight_;
    double ion_size_parameter_;
    double logK_;
    double molar_volume_;
    double specific_surface_area_;

    std::vector<ac::SpeciesName> species_names_;
    std::vector<double> stoichiometry_;
    std::vector<ac::SpeciesId> species_ids_;

    ac::SpeciesArray primarySpecies_;
    ac::Species water_;
    ac::Mineral* mineral_;

   private:
  };  // end class MineralTest

  MineralTest::MineralTest()
      : name_("Calcite"),
      id_(3),
      h2o_stoich_(0.0),
      charge_(0.0),
      gram_molecular_weight_(100.0872),
      ion_size_parameter_(0.0),
      logK_(1.8487),
      molar_volume_(36.9340),
      specific_surface_area_(0.987654),
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

    species_names_.push_back("Ca++");
    stoichiometry_.push_back(1.0);
    species_ids_.push_back(2);

    mineral_ = new ac::Mineral(name_, id_,
                               species_names_, stoichiometry_, species_ids_,
                               h2o_stoich_, gram_molecular_weight_, logK_,
                               molar_volume_, specific_surface_area_);

    ac::SpeciesId id = 0;
    ac::SpeciesName name = "H+";
    ac::Species H_p = ac::Species(id, name, 1.0, 1.0079, 9.0);
    H_p.act_coef(1.0);
    H_p.update(2.74965e-9);

    id = 1;
    name = "HCO3-";
    ac::Species HCO3_m = ac::Species(id, name, -1.0, 61.0171, 4.0);
    HCO3_m.act_coef(1.0);
    HCO3_m.update(9.71848e-4);

    id = 2;
    name = "Ca++";
    ac::Species Ca_pp = ac::Species(id, name, 2.0, 40.0780, 6.0);
    Ca_pp.act_coef(1.0);
    Ca_pp.update(9.61037e-5);

    water_.act_coef(1.0);
    water_.update(1.0);

    primarySpecies_.clear();
    primarySpecies_.push_back(H_p);
    primarySpecies_.push_back(HCO3_m);
    primarySpecies_.push_back(Ca_pp);
  }


  MineralTest::~MineralTest() {
    delete mineral_;
  }

  class MineralTestHydrated {
     public:

     protected:
      MineralTestHydrated();
      ~MineralTestHydrated();

      ac::SpeciesName name_;
      ac::SpeciesId id_;
      double h2o_stoich_;
      double charge_;
      double gram_molecular_weight_;
      double ion_size_parameter_;
      double logK_;
      double molar_volume_;
      double specific_surface_area_;

      std::vector<ac::SpeciesName> species_names_;
      std::vector<double> stoichiometry_;
      std::vector<ac::SpeciesId> species_ids_;

      ac::SpeciesArray primarySpecies_;
      ac::Species water_;
      ac::Mineral* mineral_;

     private:
    };  // end class MineralTest

    MineralTestHydrated::MineralTestHydrated()
        : name_("Alunogen"),
        id_(2),
        h2o_stoich_(17.0),
        charge_(0.0),
        gram_molecular_weight_(0.0),
        ion_size_parameter_(0.0),
        logK_(-6.204e0),
        molar_volume_(0.0),
        specific_surface_area_(0.0),
        water_(0,"H2O", 0.0, 18.00, 0.0) {
      species_names_.clear();
      stoichiometry_.clear();
      species_ids_.clear();

      species_names_.push_back("Al+3");
      stoichiometry_.push_back(2.0);
      species_ids_.push_back(0);

      species_names_.push_back("SO4-2");
      stoichiometry_.push_back(3.0);
      species_ids_.push_back(1);

      mineral_ = new ac::Mineral(name_, id_,
                                 species_names_, stoichiometry_, species_ids_,
                                 h2o_stoich_, gram_molecular_weight_, logK_,
                                 molar_volume_, specific_surface_area_);

      ac::SpeciesId id = 0;
      ac::SpeciesName name = "Al+3";
      ac::Species Al = ac::Species(id, name, 3.0, 0.0, 0.0);
      Al.act_coef(1.0);
      Al.update(0.5);

      id = 1;
      name = "SO4-2";
      ac::Species SO4 = ac::Species(id, name, -2.0, 0.0, 0.0);
      SO4.act_coef(1.0);
      SO4.update(0.5);

      water_.act_coef(0.5);
      water_.update(1.0);

      primarySpecies_.clear();
      primarySpecies_.push_back(Al);
      primarySpecies_.push_back(SO4);
    }


    MineralTestHydrated::~MineralTestHydrated() {
      delete mineral_;
    }

  //
  // most of the basic functionality comes from the parent SecondarySpecies class.
  // only test the virtual methods from the public interface...?
  //

  // make sure we can create an object with the constructor
  TEST_FIXTURE(MineralTest, Mineral_constructor) {
    CHECK_EQUAL(id_, mineral_->identifier());
  }

  //
  // virtual public methods from parent class
  //
  TEST_FIXTURE(MineralTest, Mineral_Update) {
    mineral_->Update(primarySpecies_,water_);
    CHECK_CLOSE(mineral_->lnQK(), -0.731390740816932, 1.0e-10);
  }

  // not currently used...
  // TEST_FIXTURE(MineralTest, Mineral_AddContributionToTotal)
  // TEST_FIXTURE(MineralTest, Mineral_AddContributionToDTotal)

  //
  // local public methods
  //
  TEST_FIXTURE(MineralTest, Mineral_QoverK) {
    mineral_->Update(primarySpecies_,water_);
    CHECK_CLOSE(mineral_->Q_over_K(), 0.481239245416214, 1.0e-10);
  }

  TEST_FIXTURE(MineralTest, Mineral_saturation_index) {
    mineral_->Update(primarySpecies_,water_);
    CHECK_CLOSE(mineral_->saturation_index(), -0.317638962851925, 1.0e-10);
  }

  TEST_FIXTURE(MineralTest, Mineral_molar_volume) {
    CHECK_CLOSE(mineral_->molar_volume(), 36.9340, 1.0e-4);
  }

  TEST_FIXTURE(MineralTest, Mineral_specific_surface_area) {
    CHECK_CLOSE(mineral_->specific_surface_area(), 0.987654, 1.0e-6);
  }
  TEST_FIXTURE(MineralTest, Mineral_volume_fraction) {
    mineral_->set_volume_fraction(0.23456);
    CHECK_CLOSE(mineral_->volume_fraction(), 0.23456, 1.0e-10);
  }

/*
  TEST_FIXTURE(MineralTest, Mineral_UpdateVolumeFraction) {
    double rate = 1.0;
    double delta_t = 1.0;
    mineral_->UpdateVolumeFraction(rate, delta_time);
    CHECK_CLOSE(mineral_->volume_fraction(), 100.0, 1.0e-10);
  }
*/
/*
  TEST_FIXTURE(MineralTest, Mineral_UpdateSpecificSurfaceArea) {
    // code currently has a hard coded value of 100
    mineral_->set_volume_fraction(0.2);
    mineral_->UpdateSpecificSurfaceArea(void);
    CHECK_CLOSE(mineral_->specific_surface_area(), 100.0, 1.0e-10);
  }
*/

  TEST_FIXTURE(MineralTest, Mineral_verbosity) {
    mineral_->set_verbosity(ac::kDebugMineralKinetics);
    CHECK_EQUAL(mineral_->verbosity(), ac::kDebugMineralKinetics);
  }
  //---------------------------------------------------------------
  // Test the log10(IAP/K) for Alunogen
  //---------------------------------------------------------------
  TEST_FIXTURE(MineralTestHydrated, Mineral_saturation_index) {
      mineral_->Update(primarySpecies_,water_);
      CHECK_CLOSE(mineral_->saturation_index(),-0.418659904607586, 1.0e-10);
  }
}  // end SUITE(GeochemistryTestMineral)
