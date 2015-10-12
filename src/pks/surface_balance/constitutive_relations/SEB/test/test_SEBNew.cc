#include "UnitTest++.h"
#include "TestReporterStdout.h"

#include <iomanip>
#include <fstream>
#include <vector>

#include "seb_physics_defs.hh"
#include "seb_physics_funcs.hh"

#define DO_ALL 1

using namespace Amanzi::SurfaceBalance::SEBPhysics;

struct TestSEB {
  double poro;
  double temp;
  double pressure;
  double snow_ht;
  double QswIn;
  SEB seb;

  TestSEB(double poro_, double temp_, double pressure_, double snow_ht_, double QswIn_) :
      poro(poro_), temp(temp_), pressure(pressure_), snow_ht(snow_ht_), QswIn(QswIn_) {

    double uf = temp > 273.15 ? 1. : temp < 272.15 ? 0. : temp - 272.15;
    double pd = pressure > 101325 ? (pressure - 101325.) / 10000. : 0.;

    seb.in.dt = 86400.; // one day

    // ground properties
    seb.in.vp_ground.temp = temp;
    seb.in.vp_ground.pressure = pressure;

    // snow properties
    seb.in.snow_old.ht = snow_ht;
    seb.in.snow_old.density = seb.params.density_freshsnow;
    seb.in.snow_old.age = 1.;
    seb.out.snow_new = seb.in.snow_old;
    seb.in.vp_snow.temp = 270.15;

    // met data
    seb.in.met.Us = 5.;
    seb.in.met.QswIn = QswIn;
    seb.in.met.Ps = snow_ht > 0. ? 1.e-8 : 0.;
    seb.in.met.Pr = snow_ht > 0. ? 0. : 1.e-8;
    seb.in.met.vp_air.temp = 274.15;
    seb.in.met.vp_air.relative_humidity = 0.8;

    // smoothed/interpolated surface properties
    SurfaceParams surf_pars;

    Partition al_part = Partitioner().CalcPartition(snow_ht,pd,uf);
    seb.in.surf.albedo = al_part.Interpolate(CalcAlbedoSnow(seb.in.snow_old.density),
            surf_pars.a_water, surf_pars.a_ice, surf_pars.a_tundra);

    Partition other_part = Partitioner(0.02,0.02).CalcPartition(snow_ht,pd,uf);
    seb.in.surf.emissivity = other_part.Interpolate(surf_pars.e_snow,
            surf_pars.e_water, surf_pars.e_ice, surf_pars.e_tundra);
    seb.in.vp_ground.porosity = other_part.Interpolate(1., 1., 1., poro);

    // roughness factor
    seb.in.surf.Zo = CalcRoughnessFactor(seb.in.met.vp_air.temp);
  }

  void static WriteFile(std::string prefix,
                 std::vector<double>& indep,
                 std::vector<double>& e_src,
                 std::vector<double>& m_src_surf,
                 std::vector<double>& m_src_subsurf,
                 std::vector<double>& temp_m_src) {
    std::ofstream fid;
    fid.open(("SEBNew_"+prefix+std::string(".dat")).c_str());
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


SUITE(NEW_SEB) {

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

      CalculateSurfaceBalance(seb.seb, true);

      CHECK(!std::isnan(seb.seb.out.eb.fQc));
      e_src.push_back(seb.seb.out.eb.fQc);

      CHECK(!std::isnan(seb.seb.out.mb.MWg));
      m_src_surf.push_back(seb.seb.out.mb.MWg);

      CHECK(!std::isnan(seb.seb.out.mb.MWg_subsurf));
      m_src_subsurf.push_back(seb.seb.out.mb.MWg_subsurf);

      CHECK(!std::isnan(seb.seb.out.mb.MWg_temp));
      temp_m_src.push_back(seb.seb.out.mb.MWg_temp);
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

      CalculateSurfaceBalance(seb.seb, true);

      CHECK(!std::isnan(seb.seb.out.eb.fQc));
      e_src.push_back(seb.seb.out.eb.fQc);

      CHECK(!std::isnan(seb.seb.out.mb.MWg));
      m_src_surf.push_back(seb.seb.out.mb.MWg);

      CHECK(!std::isnan(seb.seb.out.mb.MWg_subsurf));
      m_src_subsurf.push_back(seb.seb.out.mb.MWg_subsurf);

      CHECK(!std::isnan(seb.seb.out.mb.MWg_temp));
      temp_m_src.push_back(seb.seb.out.mb.MWg_temp);
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

      CalculateSurfaceBalance(seb.seb, true);

      CHECK(!std::isnan(seb.seb.out.eb.fQc));
      e_src.push_back(seb.seb.out.eb.fQc);

      CHECK(!std::isnan(seb.seb.out.mb.MWg));
      m_src_surf.push_back(seb.seb.out.mb.MWg);

      CHECK(!std::isnan(seb.seb.out.mb.MWg_subsurf));
      m_src_subsurf.push_back(seb.seb.out.mb.MWg_subsurf);

      CHECK(!std::isnan(seb.seb.out.mb.MWg_temp));
      temp_m_src.push_back(seb.seb.out.mb.MWg_temp);
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

      CalculateSurfaceBalance(seb.seb, true);

      CHECK(!std::isnan(seb.seb.out.eb.fQc));
      e_src.push_back(seb.seb.out.eb.fQc);

      CHECK(!std::isnan(seb.seb.out.mb.MWg));
      m_src_surf.push_back(seb.seb.out.mb.MWg);

      CHECK(!std::isnan(seb.seb.out.mb.MWg_subsurf));
      m_src_subsurf.push_back(seb.seb.out.mb.MWg_subsurf);

      CHECK(!std::isnan(seb.seb.out.mb.MWg_temp));
      temp_m_src.push_back(seb.seb.out.mb.MWg_temp);
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

      CalculateSurfaceBalance(seb.seb, true);

      CHECK(!std::isnan(seb.seb.out.eb.fQc));
      e_src.push_back(seb.seb.out.eb.fQc);

      CHECK(!std::isnan(seb.seb.out.mb.MWg));
      m_src_surf.push_back(seb.seb.out.mb.MWg);

      CHECK(!std::isnan(seb.seb.out.mb.MWg_subsurf));
      m_src_subsurf.push_back(seb.seb.out.mb.MWg_subsurf);

      CHECK(!std::isnan(seb.seb.out.mb.MWg_temp));
      temp_m_src.push_back(seb.seb.out.mb.MWg_temp);
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

      CalculateSurfaceBalance(seb.seb, true);

      CHECK(!std::isnan(seb.seb.out.eb.fQc));
      e_src.push_back(seb.seb.out.eb.fQc);

      CHECK(!std::isnan(seb.seb.out.mb.MWg));
      m_src_surf.push_back(seb.seb.out.mb.MWg);

      CHECK(!std::isnan(seb.seb.out.mb.MWg_subsurf));
      m_src_subsurf.push_back(seb.seb.out.mb.MWg_subsurf);

      CHECK(!std::isnan(seb.seb.out.mb.MWg_temp));
      temp_m_src.push_back(seb.seb.out.mb.MWg_temp);
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

      CalculateSurfaceBalance(seb.seb, true);

      CHECK(!std::isnan(seb.seb.out.eb.fQc));
      e_src.push_back(seb.seb.out.eb.fQc);

      CHECK(!std::isnan(seb.seb.out.mb.MWg));
      m_src_surf.push_back(seb.seb.out.mb.MWg);

      CHECK(!std::isnan(seb.seb.out.mb.MWg_subsurf));
      m_src_subsurf.push_back(seb.seb.out.mb.MWg_subsurf);

      CHECK(!std::isnan(seb.seb.out.mb.MWg_temp));
      temp_m_src.push_back(seb.seb.out.mb.MWg_temp);
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

      CalculateSurfaceBalance(seb.seb, true);

      CHECK(!std::isnan(seb.seb.out.eb.fQc));
      e_src.push_back(seb.seb.out.eb.fQc);

      CHECK(!std::isnan(seb.seb.out.mb.MWg));
      m_src_surf.push_back(seb.seb.out.mb.MWg);

      CHECK(!std::isnan(seb.seb.out.mb.MWg_subsurf));
      m_src_subsurf.push_back(seb.seb.out.mb.MWg_subsurf);

      CHECK(!std::isnan(seb.seb.out.mb.MWg_temp));
      temp_m_src.push_back(seb.seb.out.mb.MWg_temp);
    }

    TestSEB::WriteFile("QswIn_snow", QswIns, e_src, m_src_surf, m_src_subsurf, temp_m_src);
  }
#endif

}


