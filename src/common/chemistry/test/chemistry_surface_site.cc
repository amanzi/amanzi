/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include <cstdlib>
#include <cmath>

#include <iostream>
#include <vector>

#include <UnitTest++.h>

#include "Species.hh"
#include "SurfaceSite.hh"

SUITE(GeochemistryTestsSurfaceSite)
{
  namespace ac = Amanzi::AmanziChemistry;
  /*
    Unit tests for the SurfaceSite object public interface

    Need better tests for incorrect use of this object
  */
  /*****************************************************************************
   **
   **  Common testing code
   **
   *****************************************************************************/
  class SurfaceSiteTest {
   public:
    SurfaceSiteTest();
    ~SurfaceSiteTest();

   protected:
    std::string name_;
    int id_;
    double molar_density_;
    ac::SurfaceSite site_;

   private:
  }; // end class SurfaceSiteTest

  SurfaceSiteTest::SurfaceSiteTest() : name_(">FeOH"), id_(5), molar_density_(1.23456)
  {
    Teuchos::ParameterList plist;
    plist.set<double>("density", molar_density_);
    site_ = ac::SurfaceSite(name_, id_, plist);
  }

  SurfaceSiteTest::~SurfaceSiteTest() {}

  /*****************************************************************************
  *  Setup test problems
  *****************************************************************************/
  //
  // create a SurfaceSite object by specifing parameters to the
  // constructor, the primary use of the object. check that the
  // parameter data is set correctly by the constructor
  //
  TEST_FIXTURE(SurfaceSiteTest, SurfaceSite_constructor_init_id)
  {
    CHECK_EQUAL(id_, site_.identifier());
  }

  TEST_FIXTURE(SurfaceSiteTest, SurfaceSite_constructor_init_molar_density)
  {
    CHECK_EQUAL(molar_density_, site_.molar_density());
  }

  TEST_FIXTURE(SurfaceSiteTest, SurfaceSite_constructor_init_name)
  {
    CHECK_EQUAL(name_, site_.name());
  }

  TEST_FIXTURE(SurfaceSiteTest, SurfaceSite_constructor_init_change)
  {
    CHECK_EQUAL(0.0, site_.charge());
  }

  TEST_FIXTURE(SurfaceSiteTest, SurfaceSite_constructor_init_molar_surface_density)
  {
    CHECK_EQUAL(0.0, site_.molar_surface_density());
  }

  TEST_FIXTURE(SurfaceSiteTest, SurfaceSite_constructor_init_free_site_conc)
  {
    CHECK_EQUAL(0.1 * site_.molar_density(), site_.free_site_concentration());
  }

  TEST_FIXTURE(SurfaceSiteTest, SurfaceSite_constructor_init_ln_free_site_conc)
  {
    CHECK_EQUAL(std::log(0.1 * site_.molar_density()), site_.ln_free_site_concentration());
  }

  //
  // check that updating the individual parameters works correctly
  //
  TEST_FIXTURE(SurfaceSiteTest, SurfaceSite_set_molar_density)
  {
    double molar_density(0.12345);
    site_.set_molar_density(molar_density);
    CHECK_EQUAL(molar_density, site_.molar_density());
  }

  TEST_FIXTURE(SurfaceSiteTest, SurfaceSite_set_molar_surface_density)
  {
    double molar_surface_density(-0.9876);
    site_.set_molar_surface_density(molar_surface_density);
    CHECK_EQUAL(molar_surface_density, site_.molar_surface_density());
  }

  TEST_FIXTURE(SurfaceSiteTest, SurfaceSite_set_free_site_conc)
  {
    // set free_site_concentration and update ln_free_site_concentration!
    double fsc(0.0076543);
    site_.set_free_site_concentration(fsc);
    CHECK_EQUAL(fsc, site_.free_site_concentration());
    CHECK_CLOSE(std::log(fsc), site_.ln_free_site_concentration(), 1.0e-9);
  }

  TEST_FIXTURE(SurfaceSiteTest, SurfaceSite_set_ln_free_site_conc)
  {
    // set ln_free_site_concentration and update free_site_concentration
    double ln_fsc(-1.234);
    site_.set_ln_free_site_concentration(ln_fsc);
    CHECK_EQUAL(ln_fsc, site_.ln_free_site_concentration());
    CHECK_CLOSE(std::exp(ln_fsc), site_.free_site_concentration(), 1.0e-9);
  }

  //
  // check that public interface works correctly
  //
  TEST_FIXTURE(SurfaceSiteTest, SurfaceSite_SiteDensity)
  {
    // current behavior, should be different when minerals added....
    CHECK_EQUAL(molar_density_, site_.SiteDensity());
  }
} // end SUITE(GeochemistryTestSurfaceSite)
