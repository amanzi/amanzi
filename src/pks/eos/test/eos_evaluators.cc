/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*
  EOS

*/

#include <iostream>
#include <cstdio>
#include <cmath>

#include "UnitTest++.h"

#include "errors.hh"

#include "EOS_Density.hh"
#include "EOS_SaturatedVaporPressure.hh"
#include "EOSFactory.hh"
#include "H2O_Density.hh"
#include "H2O_DensityFEHM.hh"
#include "H2O_ViscosityFEHM.hh"
#include "LookupTable_Amanzi.hh"
#include "LookupTable_FEHM.hh"


TEST(DensityEOS)
{
  using namespace Amanzi::AmanziEOS;

  Teuchos::ParameterList plist;
  H2O_Density eos(plist);

  double p(101325.0), T0(273.15), d0, d1;
  for (double T = 0; T < 30; T += 0.1) {
    d1 = eos.Density(T + T0, p);
    CHECK_CLOSE(d1, 998.0, 2.0);

    if (T > 10.0) {
      CHECK(d1 < d0);
      double der = eos.DDensityDT(T + T0, p);
      CHECK(der < 0.0);
    }
    d0 = d1;
  }

  // next EOS
  H2O_DensityFEHM eos_fehm(plist);

  p = 101325.0;
  for (int loop = 0; loop < 5; loop++) {
    p *= 2.0;
    d0 = 1200.0;
    for (double T = 15; T < 300; T += 1.5) {
      d1 = eos_fehm.Density(T + T0, p);
      CHECK(d1 < 1001.0 && d1 > 550.0 && d1 < d0);
      d0 = d1;

      double der = eos_fehm.DDensityDT(T + T0, p);
      CHECK(der < 0.0);
    }
  }
}


TEST(ViscosityEOS)
{
  using namespace Amanzi::AmanziEOS;

  Teuchos::ParameterList plist;
  H2O_ViscosityFEHM eos(plist);

  double p(101325.0), T0(273.15), nu0, nu1;
  for (double T = 1; T < 300; T += 0.5) {
    nu1 = eos.Viscosity(T + T0, p);
    if (T < 100.0) CHECK_CLOSE(nu1, 0.0012, 0.001);

    if (T > 10.0) {
      CHECK(nu1 < nu0);
      double der = eos.DViscosityDT(T + T0, p);
      CHECK(der < 0.0);
    }
    nu0 = nu1;
  }
}


TEST(TabularEOS_Amanzi)
{
  using namespace Amanzi::AmanziEOS;

  int ierr;
  Teuchos::ParameterList plist;
  plist.set<std::string>("table name", "test/h2o.eos").set<std::string>("field name", "density");

  LookupTable_Amanzi eos(plist);

  for (double T = 283.15; T < 320; T += 10.0) {
    for (double p = 1e+5; p < 1.4e+5; p += 1.0e+4) {
      double d = eos.Function(T, p, &ierr);
      CHECK_CLOSE(d, 995.0, 5.0);
    }
  }

  // verify derivatives
  double T(293.15), p(101325), dT(10.0), dp(100.0);
  double dFdT = eos.DFunctionDT(T + dT / 2, p, &ierr);
  CHECK_CLOSE(
    (eos.Function(T + dT, p, &ierr) - eos.Function(T, p, &ierr)) / dT, dFdT, fabs(dFdT) * 5e-2);

  double dFdP = eos.DFunctionDp(T, p + dp / 2, &ierr);
  CHECK_CLOSE(
    (eos.Function(T, p + dp, &ierr) - eos.Function(T, p, &ierr)) / dp, dFdP, fabs(dFdP) * 5e-2);

  std::cout << "water density at 20C  and 1atm  = " << eos.Function(293.15, 101325.0, &ierr)
            << ", derivatives: " << eos.DFunctionDT(293.15, 101325.0, &ierr) << "  "
            << eos.DFunctionDp(293.15, 101325.0, &ierr) << std::endl;
}

TEST(TabularEOS_FEHM)
{
  using namespace Amanzi::AmanziEOS;

  std::vector<std::string> fields({ "density", "viscosity", "internal_energy" });
  for (int loop = 0; loop < 3; ++loop) {
    int ierr;
    Teuchos::ParameterList plist;
    plist.set<std::string>("table name", "test/air.eos")
         .set<std::string>("field name", fields[loop])
         .set<double>("molar weight", 0.02896);

    LookupTable_FEHM eos(plist);

    for (double T = 293.15; T < 320; T += 10.0) {
      for (double p = 1e+5; p < 1.4e+5; p += 1.0e+4) {
         double d = eos.Function(T, p, &ierr);
      }
    }

    // verify derivatives
    double T(273.15), p(101325), dT(20.0), dp(100.0);
    double dFdT = eos.DFunctionDT(T + dT / 2, p, &ierr);
    CHECK_CLOSE(
      (eos.Function(T + dT, p, &ierr) - eos.Function(T, p, &ierr)) / dT, dFdT, fabs(dFdT) * 5e-2);

    double dFdP = eos.DFunctionDp(T, p + dp / 2, &ierr);
    CHECK_CLOSE(
      (eos.Function(T, p + dp, &ierr) - eos.Function(T, p, &ierr)) / dp, dFdP, fabs(dFdP) * 5e-2);

    std::cout << "air " << fields[loop] << " at 0C and 1atm  = " << eos.Function(T, p, &ierr)
              << ", P/T derivatives: " << eos.DFunctionDT(T, p, &ierr) << "  "
              << eos.DFunctionDp(T, p, &ierr) << std::endl;
  }
  std::cout << std::endl;
}


TEST(FactoryEOS)
{
  using namespace Amanzi::AmanziEOS;

  std::vector<std::string> names = { "liquid water 0-30C", "liquid water FEHM", "lookup table" };

  for (auto& name : names) {
    Teuchos::ParameterList plist;
    plist.set<std::string>("table name", "test/h2o.eos")
      .set<std::string>("field name", "density")
      .set<std::string>("eos type", name);

    EOSFactory<EOS_Density> factory;
    auto eos = factory.Create(plist);

    std::cout << name << ": water density at 20C and 1atm = " << eos->Density(293.15, 101325.0)
              << std::endl;
  }
  std::cout << std::endl;

  // viscosity
  for (auto& name : names) {
    Teuchos::ParameterList plist;
    plist.set<std::string>("table name", "test/h2o.eos")
      .set<std::string>("field name", "viscosity")
      .set<std::string>("eos type", name);

    EOSFactory<EOS_Viscosity> factory;
    auto eos = factory.Create(plist);

    std::cout << name << ": water viscosity at 20C and 1atm = " << eos->Viscosity(293.15, 101325.0)
              << std::endl;
  }
}


TEST(Exceptions)
{
  using namespace Amanzi::AmanziEOS;

  std::vector<std::string> names = { "liquid water 0-30C", "liquid water FEHM", "lookup table" };

  std::cout << std::endl;
  for (auto& name : names) {
    Teuchos::ParameterList plist;
    plist.set<std::string>("table name", "test/h2o.eos")
      .set<std::string>("field name", "density")
      .set<std::string>("eos type", name);

    // density
    {
      EOSFactory<EOS_Density> factory;
      auto eos = factory.Create(plist);

      try {
        eos->Density(193.15, -10.0);
        std::cout << name << ": density (err=" << eos->error_code() << ", \"" << eos->error_msg()
                  << "\")\n";
      } catch (const Errors::CutTimeStep& e) {
        std::cout << name << ": density exception: " << e.what() << std::endl;
      } catch (...) {
        AMANZI_ASSERT(false);
      }
    }

    // viscosity
    {
      EOSFactory<EOS_Viscosity> factory;
      auto eos = factory.Create(plist);

      try {
        eos->Viscosity(193.15, 1.0e+10);
        std::cout << name << ": viscosity (err=" << eos->error_code() << ", \"" << eos->error_msg()
                  << "\")\n";
      } catch (const Errors::CutTimeStep& e) {
        std::cout << name << ": viscosity exception: " << e.what() << std::endl;
      } catch (...) {
        AMANZI_ASSERT(false);
      }
    }
  }

  // saturated vapor pressure
  Teuchos::ParameterList plist;
  plist.set<std::string>("eos type", "water vapor over water/ice");
  EOSFactory<EOS_SaturatedVaporPressure> factory;
  auto eos = factory.Create(plist);

  try {
    eos->Pressure(93.15);
    std::cout << "saturated vapor pressure (err=" << eos->error_code() << ", \"" << eos->error_msg()
              << "\")\n";
  } catch (const Errors::CutTimeStep& e) {
    std::cout << "saturated vapor pressure exception: " << e.what() << std::endl;
  } catch (...) {
    AMANZI_ASSERT(false);
  }
}
