/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Energy PK
*/

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// TPLs
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "UnitTest++.h"

// Amanzi
#include "CommonDefs.hh"
#include "CompositeVector.hh"
#include "EnergyPressureTemperature_PK.hh"
#include "evaluators_reg.hh"
#include "IAPWS95.hh"
#include "MeshFactory.hh"
#include "PK_Physical.hh"
#include "State.hh"
#include "IAPWS95_StateEvaluators.hh"
#include "VerboseObject.hh"

TEST(GIBBS_ENERGY_IAPWS95_SINGLE_PHASE)
{
  using namespace Amanzi;
  std::cout << "Test: Gibbs energy, single phase" << std::endl;

  Teuchos::ParameterList plist;
  AmanziEOS::IAPWS95 eos(plist);

  int n = 100;
  double scale(50.0 / n), eps(2e-6);
  double dp(0.05 * scale), dT(1.0 * scale), p, pl, pr, p1, T, Tl, Tr, T1, der, der1, der2, der3;

  for (int i = -n; i < n; ++i) {
    for (int j = -n; j < n; ++j) {
      p = eos.PC + (i + 0.5) * dp; // Mpa
      T = eos.TC + (j + 0.5) * dT;

      pl = p * (1.0 - eps);
      pr = p * (1.0 + eps);

      Tl = T * (1.0 - eps);
      Tr = T * (1.0 + eps);

      auto [prop, liquid, vapor] = eos.ThermodynamicsPT(p, T);
      CHECK_CLOSE(prop.helmholtz, prop.u - T * prop.s, eps * std::fabs(prop.helmholtz));
      CHECK(prop.cv > 0.0);
      // std::cout << p << " " << T << " " << prop.helmholtz << std::endl;

      auto [prop1, liquid1, vapor1] = eos.ThermodynamicsPT(pr, T);
      auto [prop2, liquid2, vapor2] = eos.ThermodynamicsPT(pl, T);

      auto [prop3, liquid3, vapor3] = eos.ThermodynamicsPT(p, Tr);
      auto [prop4, liquid4, vapor4] = eos.ThermodynamicsPT(p, Tl);

      der = (prop1.helmholtz - prop.helmholtz) / (eps * p);
      CHECK_CLOSE(der, 1000.0 * p * prop.kt * prop.v, 3000 * eps * std::fabs(der));

      // Gibbs energy is concave in T and p 
      der1 = (prop4.gibbs - 2 * prop.gibbs + prop3.gibbs) / (eps * eps * T * T);
      CHECK(der1 < 0.0);
      der2 = (prop1.gibbs - 2 * prop.gibbs + prop2.gibbs) / (eps * eps * p * p);
      CHECK(der2 < 0.0);
      der3 = (prop1.gibbs + prop3.gibbs - prop2.gibbs - prop4.gibbs) / (eps * eps * p * T) / 4;

      // determinant of Hessian is positive
      double der3(0.0);
      for (int k = -1; k < 2; k+=2) {
        p1 = p + k * dp;
        for (int l = -1; l < 2; l+=2) {
          T1 = T + l * dT;
          auto [prop1, liquid1, vapor1] = eos.ThermodynamicsPT(p1, T1);
          der3 += prop1.gibbs * k * l; 
        }      
      }      
      der3 /= 4 * p * T;
      CHECK(der1 * der2 - der3 * der3 > 0.0);
    }
  }
}


TEST(GIBBS_ENERGY_IAPWS95_TWO_PHASES)
{
  using namespace Amanzi;
  std::cout << "Test: Helmholtz energy, two phases" << std::endl;

  Teuchos::ParameterList plist;
  AmanziEOS::IAPWS95 eos(plist);

  int n(200), count(0);
  double eps(1e-6);
  double dT(1.0), T, rhol0, rhov0, rhol, rhov, rho, p; 

  for (int i = 0; i < n; ++i) {
    T = eos.TC - (i + 0.5) * dT;
    rhol0 = eos.DensityLiquid(T);
    rhov0 = eos.DensityVapor(T);
    std::tie(rhol, rhov, p) = eos.SaturationLine(T, rhol0, rhov0);

    rho = (rhol + rhov) / 2;
    auto [prop, liquid, vapor] = eos.ThermodynamicsRhoT(rho, T);
    CHECK_CLOSE(liquid.gibbs, vapor.gibbs, eps * std::fabs(liquid.gibbs));
    if (liquid.gibbs == 0.0) count++;
  }

  CHECK(count == 0);
}


TEST(HELMHOLTZ_ENERGY_IAPWS95_SINGLE_PHASE)
{
  using namespace Amanzi;
  std::cout << "Test: Helmholtz energy, single phase" << std::endl;

  Teuchos::ParameterList plist;
  AmanziEOS::IAPWS95 eos(plist);

  int n = 100;
  double scale(50.0 / n), eps(2e-6);
  double drho(6.0 * scale), dT(1.0 * scale), rho, T, Tl, Tr, der;

  for (int i = -n; i < n; ++i) {
    for (int j = -n; j < n; ++j) {
      rho = 2 * eos.RHOC + (i + 0.5) * drho;
      T = eos.TC + (j + 0.5) * dT;

      Tl = T * (1.0 - eps);
      Tr = T * (1.0 + eps);

      auto [prop, liquid, vapor] = eos.ThermodynamicsRhoT(rho, T);

      auto [prop3, liquid3, vapor3] = eos.ThermodynamicsRhoT(rho, Tr);
      auto [prop4, liquid4, vapor4] = eos.ThermodynamicsRhoT(rho, Tl);

      // Helmholtz energy is concave in T
      der = (prop4.helmholtz - 2 * prop.helmholtz + prop3.helmholtz) / (eps * eps * T * T);
      CHECK(der < 0.0);
    }
  }
}


// TEST(HELMHOLTZ_IAPWS95_SPLINE)
// {
// }

