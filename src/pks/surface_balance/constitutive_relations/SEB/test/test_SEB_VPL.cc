#include "UnitTest++.h"
#include "TestReporterStdout.h"

#include <iomanip>
#include <limits>
#include <cmath>
#include <fstream>
#include <vector>

#include "SnowEnergyBalance_VPL.hh"

#define DO_ALL 1

using namespace SurfaceEnergyBalance_VPL;

struct TestSEB {
  double poro;
  double temp;
  double pressure;
  double snow_ht;
  double QswIn;
  LocalData seb;

  TestSEB(double poro_, double temp_, double pressure_, double snow_ht_, double QswIn_) :
      poro(poro_), temp(temp_), pressure(pressure_), snow_ht(snow_ht_), QswIn(QswIn_) {

    double uf = temp > 273.15 ? 1. : temp < 272.15 ? 0. : temp - 272.15;
    double pd = pressure > 101325 ? (pressure - 101325.) / 10000. : 0.;

    double dt = 86400.;
    seb.st_energy.dt = dt; // one day
    seb.st_energy.AlbedoTrans = 0.02;
    seb.vp_ground.relative_humidity = 1.;

    seb.st_energy.Zo=0.005;
    double air_temp = 274.15;
    if (air_temp > 270){// Little ditty I wrote for the roughness lenght ~ AA 1/10/14
      double Zsmooth = 0.005;
      double Zrough = 0.04;
      double Zfraction = -0.1*air_temp + 28;
      if (air_temp>=280){
        Zfraction = 0;
      }
      seb.st_energy.Zo=(Zsmooth*Zfraction) + (Zrough*(1-Zfraction));
    }

    seb.st_energy.water_depth = pd;
    seb.st_energy.surface_pressure = pressure;
    seb.st_energy.water_fraction = uf;
    seb.st_energy.temp_ground = temp;
    seb.vp_ground.temp = temp;
    seb.vp_ground.actual_vaporpressure = 0.;
    seb.st_energy.surface_porosity = poro;
    seb.st_energy.temp_air = 274.15;
    seb.st_energy.QswIn = QswIn;
    seb.st_energy.Us = 5.;
    seb.st_energy.Pr = (snow_ht > 0. ? 0. : 1.e-8) * dt;
    seb.st_energy.Ps = (snow_ht > 0. ? 1.e-8 : 0.) * dt;
    seb.vp_air.temp = air_temp;
    seb.vp_air.relative_humidity = 0.8;
    seb.st_energy.ht_snow = snow_ht;
    seb.st_energy.density_snow = 100.;
    seb.st_energy.age_snow = 1.;
    seb.st_energy.stored_surface_pressure = pressure;
  }


  void static WriteFile(std::string prefix,
                 std::vector<double>& indep,
                 std::vector<double>& e_src,
                 std::vector<double>& m_src_surf,
                 std::vector<double>& m_src_subsurf,
                 std::vector<double>& temp_m_src) {
    std::ofstream fid;
    fid.open(("SEBVPL_"+prefix+std::string(".dat")).c_str());
    for (unsigned int i=0; i!=indep.size(); ++i) {
      fid << std::setprecision(15)
          << indep[i] << " "
          << e_src[i] << " "
          << m_src_surf[i] << " "
          << m_src_subsurf[i] << " "
          << temp_m_src[i] << std::endl;
    }
  };

};


SUITE(SEB_VPL) {

#if DO_ALL
  TEST(TEST_TEMP_UNSATURATED_SNOW) {
    std::vector<double> temps;
    std::vector<double> e_src;
    std::vector<double> m_src_surf;
    std::vector<double> m_src_subsurf;
    std::vector<double> temp_m_src;

    double T0 = 268.15;
    double dT = 0.1;
    for (int i=0; i!=100; ++i) {
      double T = T0 + i * dT;
      temps.push_back(T);

      // initialize
      TestSEB seb(0.5, T, 30000., 0.3, 65);

      SnowEnergyBalance(seb.seb);

      CHECK(!std::isnan(seb.seb.st_energy.fQc));
      e_src.push_back(seb.seb.st_energy.fQc);

      CHECK(!std::isnan(seb.seb.st_energy.Mr));
      m_src_surf.push_back(seb.seb.st_energy.Mr);

      CHECK(!std::isnan(seb.seb.st_energy.SurfaceVaporFlux));
      m_src_subsurf.push_back(seb.seb.st_energy.SurfaceVaporFlux);

      CHECK(!std::isnan(seb.seb.st_energy.Trw));
      temp_m_src.push_back(seb.seb.st_energy.Trw);
    }

    TestSEB::WriteFile("temp_unsat_snow", temps, e_src, m_src_surf, m_src_subsurf, temp_m_src);
  }

#endif

#if DO_ALL
  TEST(TEST_TEMP_UNSATURATED_NOSNOW) {
    std::vector<double> temps;
    std::vector<double> e_src;
    std::vector<double> m_src_surf;
    std::vector<double> m_src_subsurf;
    std::vector<double> temp_m_src;

    double T0 = 268.15;
    double dT = 0.1;
    for (int i=0; i!=100; ++i) {
      double T = T0 + i * dT;
      temps.push_back(T);

      // initialize
      TestSEB seb(0.5, T, 30000., 0., 65.);

      SnowEnergyBalance(seb.seb);

      CHECK(!std::isnan(seb.seb.st_energy.fQc));
      e_src.push_back(seb.seb.st_energy.fQc);

      CHECK(!std::isnan(seb.seb.st_energy.Mr));
      m_src_surf.push_back(seb.seb.st_energy.Mr);

      CHECK(!std::isnan(seb.seb.st_energy.SurfaceVaporFlux));
      m_src_subsurf.push_back(seb.seb.st_energy.SurfaceVaporFlux);

      CHECK(!std::isnan(seb.seb.st_energy.Trw));
      temp_m_src.push_back(seb.seb.st_energy.Trw);
    }

    TestSEB::WriteFile("temp_unsat_nosnow", temps, e_src, m_src_surf, m_src_subsurf, temp_m_src);
  }


#endif

#if DO_ALL
  TEST(TEST_TEMP_SATURATED_SNOW) {
    std::vector<double> temps;
    std::vector<double> e_src;
    std::vector<double> m_src_surf;
    std::vector<double> m_src_subsurf;
    std::vector<double> temp_m_src;

    double T0 = 268.15;
    double dT = 0.1;
    for (int i=0; i!=100; ++i) {
      double T = T0 + i * dT;
      temps.push_back(T);

      // initialize
      TestSEB seb(0.5, T, 102325., 0.3, 65.);

      SnowEnergyBalance(seb.seb);

      CHECK(!std::isnan(seb.seb.st_energy.fQc));
      e_src.push_back(seb.seb.st_energy.fQc);

      CHECK(!std::isnan(seb.seb.st_energy.Mr));
      m_src_surf.push_back(seb.seb.st_energy.Mr);

      CHECK(!std::isnan(seb.seb.st_energy.SurfaceVaporFlux));
      m_src_subsurf.push_back(seb.seb.st_energy.SurfaceVaporFlux);

      CHECK(!std::isnan(seb.seb.st_energy.Trw));
      temp_m_src.push_back(seb.seb.st_energy.Trw);
    }

    TestSEB::WriteFile("temp_sat_snow", temps, e_src, m_src_surf, m_src_subsurf, temp_m_src);
  }


#endif

#if DO_ALL
  TEST(TEST_TEMP_SATURATED_NOSNOW) {
    std::vector<double> temps;
    std::vector<double> e_src;
    std::vector<double> m_src_surf;
    std::vector<double> m_src_subsurf;
    std::vector<double> temp_m_src;

    double T0 = 268.15;
    double dT = 0.1;
    for (int i=0; i!=100; ++i) {
      double T = T0 + i * dT;
      temps.push_back(T);

      // initialize
      TestSEB seb(0.5, T, 102325., 0., 65.);

      SnowEnergyBalance(seb.seb);

      CHECK(!std::isnan(seb.seb.st_energy.fQc));
      e_src.push_back(seb.seb.st_energy.fQc);

      CHECK(!std::isnan(seb.seb.st_energy.Mr));
      m_src_surf.push_back(seb.seb.st_energy.Mr);

      CHECK(!std::isnan(seb.seb.st_energy.SurfaceVaporFlux));
      m_src_subsurf.push_back(seb.seb.st_energy.SurfaceVaporFlux);

      CHECK(!std::isnan(seb.seb.st_energy.Trw));
      temp_m_src.push_back(seb.seb.st_energy.Trw);
    }

    TestSEB::WriteFile("temp_sat_nosnow", temps, e_src, m_src_surf, m_src_subsurf, temp_m_src);
  }


#endif

#if DO_ALL
  TEST(TEST_PRES_UNFROZEN_SNOW) {
    std::vector<double> pres;
    std::vector<double> e_src;
    std::vector<double> m_src_surf;
    std::vector<double> m_src_subsurf;
    std::vector<double> temp_m_src;

    double p0 = 0.;
    double dp = 1000.;
    for (int i=0; i!=110; ++i) {
      double p = p0 + i * dp;
      pres.push_back(p);

      // initialize
      TestSEB seb(0.5, 275., p, 0.3, 65.);

      SnowEnergyBalance(seb.seb);

      CHECK(!std::isnan(seb.seb.st_energy.fQc));
      e_src.push_back(seb.seb.st_energy.fQc);

      CHECK(!std::isnan(seb.seb.st_energy.Mr));
      m_src_surf.push_back(seb.seb.st_energy.Mr);

      CHECK(!std::isnan(seb.seb.st_energy.SurfaceVaporFlux));
      m_src_subsurf.push_back(seb.seb.st_energy.SurfaceVaporFlux);

      CHECK(!std::isnan(seb.seb.st_energy.Trw));
      temp_m_src.push_back(seb.seb.st_energy.Trw);
    }

    TestSEB::WriteFile("pres_unfrozen_snow", pres, e_src, m_src_surf, m_src_subsurf, temp_m_src);
  }

#endif

#if DO_ALL
  TEST(TEST_PRES_UNFROZEN_NOSNOW) {
    std::vector<double> pres;
    std::vector<double> e_src;
    std::vector<double> m_src_surf;
    std::vector<double> m_src_subsurf;
    std::vector<double> temp_m_src;

    double p0 = 90000.;
    double dp = 10.;
    for (int i=0; i!=2100; ++i) {
      double p = p0 + i * dp;
      pres.push_back(p);

      // initialize
      TestSEB seb(0.5, 275., p, 0., 65.);

      SnowEnergyBalance(seb.seb);

      CHECK(!std::isnan(seb.seb.st_energy.fQc));
      e_src.push_back(seb.seb.st_energy.fQc);

      CHECK(!std::isnan(seb.seb.st_energy.Mr));
      m_src_surf.push_back(seb.seb.st_energy.Mr);

      CHECK(!std::isnan(seb.seb.st_energy.SurfaceVaporFlux));
      m_src_subsurf.push_back(seb.seb.st_energy.SurfaceVaporFlux);

      CHECK(!std::isnan(seb.seb.st_energy.Trw));
      temp_m_src.push_back(seb.seb.st_energy.Trw);
    }

    TestSEB::WriteFile("pres_unfrozen_nosnow", pres, e_src, m_src_surf, m_src_subsurf, temp_m_src);
  }


#endif

#if DO_ALL
  TEST(TEST_QSWIN_NOSNOW) {
    std::vector<double> QswIns;
    std::vector<double> e_src;
    std::vector<double> m_src_surf;
    std::vector<double> m_src_subsurf;
    std::vector<double> temp_m_src;

    double Q0 = 10.;
    double dQ = 1.;
    for (int i=0; i!=100; ++i) {
      double Q = Q0 + i * dQ;
      QswIns.push_back(Q);

      // initialize
      TestSEB seb(0.5, 275., 101300., 0., Q);

      SnowEnergyBalance(seb.seb);

      CHECK(!std::isnan(seb.seb.st_energy.fQc));
      e_src.push_back(seb.seb.st_energy.fQc);

      CHECK(!std::isnan(seb.seb.st_energy.Mr));
      m_src_surf.push_back(seb.seb.st_energy.Mr);

      CHECK(!std::isnan(seb.seb.st_energy.SurfaceVaporFlux));
      m_src_subsurf.push_back(seb.seb.st_energy.SurfaceVaporFlux);

      CHECK(!std::isnan(seb.seb.st_energy.Trw));
      temp_m_src.push_back(seb.seb.st_energy.Trw);
    }

    TestSEB::WriteFile("QswIn_nosnow", QswIns, e_src, m_src_surf, m_src_subsurf, temp_m_src);
  }

#endif

#if DO_ALL
  TEST(TEST_QSWIN_SNOW) {
    std::vector<double> QswIns;
    std::vector<double> e_src;
    std::vector<double> m_src_surf;
    std::vector<double> m_src_subsurf;
    std::vector<double> temp_m_src;

    double Q0 = 10.;
    double dQ = 1.;
    for (int i=0; i!=100; ++i) {
      double Q = Q0 + i * dQ;
      QswIns.push_back(Q);

      // initialize
      TestSEB seb(0.5, 275., 101300., 0.3, Q);

      SnowEnergyBalance(seb.seb);

      CHECK(!std::isnan(seb.seb.st_energy.fQc));
      e_src.push_back(seb.seb.st_energy.fQc);

      CHECK(!std::isnan(seb.seb.st_energy.Mr));
      m_src_surf.push_back(seb.seb.st_energy.Mr);

      CHECK(!std::isnan(seb.seb.st_energy.SurfaceVaporFlux));
      m_src_subsurf.push_back(seb.seb.st_energy.SurfaceVaporFlux);

      CHECK(!std::isnan(seb.seb.st_energy.Trw));
      temp_m_src.push_back(seb.seb.st_energy.Trw);
    }

    TestSEB::WriteFile("QswIn_snow", QswIns, e_src, m_src_surf, m_src_subsurf, temp_m_src);
  }
#endif
}


