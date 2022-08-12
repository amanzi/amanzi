#include <unistd.h>

//#define ABORT_ON_FLOATING_POINT_EXCEPTIONS
#ifdef __APPLE__
  #include <xmmintrin.h>
#endif

#include "boost/detail/fenv.hpp"

#include <cstdlib>
#include <cctype>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

// TPLs
#include "Teuchos_XMLParameterListHelpers.hpp"
#include <UnitTest++.h>

// Chemistry
#include "Beaker.hh"
#include "ActivityModelFactory.hh"
#include "ChemistryUtilities.hh"
#include "SimpleThermoDatabase.hh"

namespace ac = Amanzi::AmanziChemistry;

int CompareFiles(const std::string& file1, const std::string& file2)
{
  std::ifstream ifs1(file1.c_str(), std::ios::in);
  std::ifstream ifs2(file2.c_str(), std::ios::in);
  if(!ifs1.good() || !ifs2.good()) return 1;

  do {
    std::string word1, word2;
    ifs1 >> word1;
    ifs2 >> word2;
    if (ifs1.fail() || ifs2.fail()) break;
    if (ifs1.eof() || ifs2.eof()) break;

    // first check that the words match
    if (std::memcmp(word1.c_str(), word2.c_str(), word1.size()) != 0) {
      double val1 = std::atof(word1.c_str());
      double val2 = std::atof(word2.c_str());
      if (std::fabs(val1 - val2) > 1e-12 * std::max(1.0, std::fabs(val1))) return 3;
    }
  } while(true);

  return 0;
}


void RunBatchNative(const std::string& filexml,
                    const std::string& filetest,
                    const std::string& activity_model,
                    const std::vector<double>& ic_total,
                    const std::vector<double>& ic_mineral,
                    const std::vector<double>& ic_ion_exchange,
                    const std::vector<double>& ic_free_ion,
                    double porosity, double saturation, double volume,
                    double dt = 0.0, int max_dt_steps = 0, int frequency = 1)
{
#ifdef ABORT_ON_FLOATING_POINT_EXCEPTIONS
#ifdef __APPLE__
  // Make floating point exceptions abort the program. runtime error
  // message isn't helpful, but running in gdb will stop at the
  // correct line. This may code may not be apple specific....
  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);
  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_DENORM);
  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_DIV_ZERO);
  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_OVERFLOW);
  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_UNDERFLOW);
#endif
#endif

  // default i/o level
  double dt0 = dt;
  auto plist = Teuchos::getParametersFromXmlFile(filexml);
  auto vo = Teuchos::rcp(new Amanzi::VerboseObject("Beaker", *plist));

  ac::BeakerState state;
  ac::BeakerParameters parameters;

  ac::Beaker* chem = new ac::SimpleThermoDatabase(plist, vo);

  parameters.tolerance = 1e-12;
  parameters.max_iterations = 250;
  parameters.update_activity_newton = false;
  parameters.activity_model_name = activity_model;
  if (activity_model == "pitzer-hwm") {
    parameters.pitzer_database = "test/chemistry_pitzer";
    parameters.pitzer_jfunction = "pitzer1975";
  }

  state.porosity = porosity;
  state.saturation = saturation;
  state.water_density = 997.16;
  state.volume = volume;

  chem->Initialize(state, parameters);

  state.mineral_volume_fraction = ic_mineral;
  state.ion_exchange_sites = ic_ion_exchange;

  chem->CopyStateToBeaker(state);

  // we do not have external state in this test, we need to initialize 
  // chemistry state from beaker's data
  chem->CopyBeakerToState(&state);

  int ncomp = chem->primary_species().size();
  int nmineral = chem->minerals().size();
  int nsorbed = chem->total_sorbed().size();
  int nion_site = chem->ion_exchange_rxns().size();
  int nisotherm = chem->sorption_isotherm_rxns().size();

  state.total = ic_total;
  if (ic_free_ion.size() > 0) {
    state.free_ion = ic_free_ion;
  } else { 
    // state.free_ion.resize(ncomp, 1.0e-9);
    for (int i = 0; i < ncomp; ++i)
      state.free_ion[i] = std::max(0.1 * state.total[i], 1e-9);
  }

  if (nmineral > 0) {
    state.mineral_volume_fraction.resize(nmineral, 0.0);
    state.mineral_specific_surface_area.resize(nmineral, 0.0);
  }
  if (nsorbed > 0)
    state.total_sorbed.resize(nsorbed, 0.0);

  if (nion_site > 0)
    state.ion_exchange_sites.resize(nion_site, 0.0);

  if (nisotherm > 0) {
    // state.surface_site_density.resize(nisotherm, 0.0);
    state.isotherm_kd.resize(ncomp, 0.0);
    state.isotherm_freundlich_n.resize(ncomp, 0.0);
    state.isotherm_langmuir_b.resize(ncomp, 0.0);
  }

  // io
  chem->Display();
  chem->DisplayComponents(state);

  // solve for free-ion concentrations
  int itrs = chem->Speciate(&state);
  std::cout << "Speciation: " << filetest << " " << itrs << " itrs" << std::endl;

  chem->CopyBeakerToState(&state);
  chem->DisplayResults();

  // kinetics
  if (dt > 0.0) {
    double time(0.0);

    for (int n = 0; n < max_dt_steps; ++n) {
      try {
        chem->ReactionStep(&state, parameters, dt);
        // chem->CopyBeakerToState(&state);
        time += dt;
        if ((n + 1) % frequency == 0) chem->DisplayTotalColumns(time, state, false);
        dt = std::min(dt0, dt * 1.1);

        // non-conservative clipping of small values
        // for (int i = 0; i < ncomp; ++i)
        //   state.total[i] = std::max(1.0e-200, state.total[i]);
      } catch (...) {
        chem->CopyStateToBeaker(state);
        dt /= 2;
      }
    }
    chem->Speciate(&state);
    chem->DisplayResults();
  }

  vo = Teuchos::null;  // closing the stream
  std::string tmp = plist->sublist("verbose object").get<std::string>("output filename");
  int ok = CompareFiles(tmp, filetest);
  if (ok != 0) std::cout << "Error = " << ok << std::endl;
  CHECK(ok == 0);

  // cleanup memory
  delete chem;
}

TEST(NATIVE_CA_DEBYE_HUCKEL) {
  std::vector<double> ict = {3.0e-3, 1.0e-3, 1.0e-3};
  std::vector<double> icm, icie, icfi, icsd, icssa;
  RunBatchNative("test/native/ca-carbonate.xml",
                 "test/native/ca-carbonate-debye-huckel.test",
                 "debye-huckel",
                 ict, icm, icie, icfi, // initial conditions
                 0.5, 1.0, 1.0);  // porosity, saturation, cell volume
}

TEST(NATIVE_CA_UNIT) {
  std::vector<double> ict = {3.0e-3, 1.0e-3, 1.0e-3};
  std::vector<double> icm, icie, icfi, icsd, icssa;
  RunBatchNative("test/native/ca-carbonate.xml",
                 "test/native/ca-carbonate-unit.test",
                 "unit",
                 ict, icm, icie, icfi,  // initial conditions
                 0.5, 1.0, 1.0);  // porosity, saturation, cell volume
}


TEST(NATIVE_CALCITE_KINETICS) {
  std::vector<double> ict = {-1.0e-5, 1.0e-5, 1.0e-5};
  std::vector<double> icm = {0.2};
  std::vector<double> icie, icfi;
  RunBatchNative("test/native/calcite.xml",
                 "test/native/calcite-kinetics.test",
                 "debye-huckel",
                 ict, icm, icie, icfi,  // initial conditions
                 0.5, 1.0, 1.0);  // porosity, saturation, cell volume
}

TEST(NATIVE_CALCITE_KINETICS_VOLUME_FRACTIONS) {
  std::vector<double> ict = {1.0e-2, 1.0e-2, 1.0e-19};
  std::vector<double> icm = {0.2};
  std::vector<double> icie, icfi;
  RunBatchNative("test/native/calcite.xml",
                 "test/native/calcite-kinetics-volume-fractions.test",
                 "debye-huckel",
                 ict, icm, icie, icfi,  // initial conditions
                 0.5, 1.0, 1.0,  // porosity, saturation, cell volume, dt
                 2592000.0, 60);  // dt, max time steps
}


TEST(NATIVE_CARBONATE_DEBYE_HUCKEL) {
  std::vector<double> ict = {1.0e-3, 1.0e-3};
  std::vector<double> icm, icie, icfi;
  RunBatchNative("test/native/carbonate.xml",
                 "test/native/carbonate-debye-huckel.test",
                 "debye-huckel",
                 ict, icm, icie, icfi,  // initial conditions
                 0.5, 1.0, 1.0);  // porosity, saturation, cell volume
}

TEST(NATIVE_CARBONATE_UNIT) {
  std::vector<double> ict = {1.0e-3, 1.0e-3};
  std::vector<double> icm, icie, icfi;
  RunBatchNative("test/native/carbonate.xml",
                 "test/native/carbonate-unit.test",
                 "unit",
                 ict, icm, icie, icfi,  // initial conditions
                 0.5, 1.0, 1.0);  // porosity, saturation, cell volume
}


TEST(NATIVE_FAREA17_UNIT) {
  std::vector<double> ict = {3.4363E-02, 1.2475E-05, 3.0440E-05, 1.7136E-05,
                             2.8909E-05, 3.6351E-03, 1.3305E-03, 3.4572E-02,
                             2.1830E-03, 3.3848E-05, 6.2463E-04, 7.1028E-05,
                             7.8954E-05, 2.5280E-04, 3.5414E-05, 2.6038E-04, 3.5414E-05};
  std::vector<double> icm = { 0.0, 0.21, 0.15,
                              0.0, 0.1,  0.15,
                              0.0, 0.0,  0.0,
                              0.0, 0.0};
  std::vector<double> icie, icfi;
  // [total_sorbed] - all zeros 
  // std::vector<double> icfi = { 9.9969E-06, 9.9746E-06, 2.2405E-18, 1.8874E-04,
  //                              5.2970E-16, 3.2759E-08, 1.0000E-05, 1.0000E-05,
  //                              1.9282E-04, 9.9999E-06, 9.9860E-07, 9.9886E-07,
  //                              1.0000E-06, 1.8703E-04, 1.7609E-20, 2.5277E-04, 1.0000E-15 };
  RunBatchNative("test/native/fbasin-17.xml",
                 "test/native/fbasin-17-unit.test",
                 "unit",
                 ict, icm, icie, icfi,  // initial conditions
                 0.5, 1.0, 1.0,
                 2592000.0, 12);  // porosity, saturation, cell volume
}

TEST(NATIVE_FAREA17_DEBYE_HUCKEL) {
  std::vector<double> ict = { 1.3132E-04, 1.0000E-05, 1.0000E-12, 1.0000E-06,
                              1.0000E-12, 1.0716E-05, 1.0000E-05, 1.0000E-05,
                              6.0081E-05, 1.0000E-05, 1.0000E-06, 1.0000E-06,
                              7.8954E-05, 1.0000E-05, 1.0000E-15, 2.5277E-04, 1.0000E-15 };
  std::vector<double> icm = { 0.0, 0.21, 0.15,
                              0.0, 0.1,  0.15,
                              0.0, 0.0,  0.0,
                              0.0, 0.0};
  std::vector<double> icie, icfi;
  RunBatchNative("test/native/fbasin-17.xml",
                 "test/native/fbasin-17-debye-huckel.test",
                 "debye-huckel",
                 ict, icm, icie, icfi,  // initial conditions
                 0.5, 1.0, 1.0,
                 2592000.0, 12);  // porosity, saturation, cell volume
}

TEST(NATIVE_FAREA17_DEBYE_HUCKEL_10PERCENT_FREE_ION) {
  std::vector<double> ict = { 1.3132E-04, 1.0000E-05, 1.0000E-12, 1.0000E-06,
                              1.0000E-12, 1.0716E-05, 1.0000E-05, 1.0000E-05,
                              6.0081E-05, 1.0000E-05, 1.0000E-06, 1.0000E-06,
                              7.8954E-05, 1.0000E-05, 1.0000E-15, 2.5277E-04, 1.0000E-15 };
  std::vector<double> icm = { 0.0, 0.21, 0.15,
                              0.0, 0.1,  0.15,
                              0.0, 0.0,  0.0,
                              0.0, 0.0};
  std::vector<double> icie;
  // [total_sorbed] - all zeros 
  std::vector<double> icfi = { 9.9969E-06, 9.9746E-06, 2.2405E-18, 1.8874E-04,
                               5.2970E-16, 3.2759E-08, 1.0000E-05, 1.0000E-05,
                               1.9282E-04, 9.9999E-06, 9.9860E-07, 9.9886E-07,
                               1.0000E-06, 1.8703E-04, 1.7609E-20, 2.5277E-04, 1.0000E-15 };  // free ion
  RunBatchNative("test/native/fbasin-17.xml",
                 "test/native/fbasin-17-debye-huckel.test",
                 "debye-huckel",
                 ict, icm, icie, icfi,  // initial conditions
                 0.5, 1.0, 1.0,
                 2592000.0, 12);  // porosity, saturation, cell volume
}

TEST(NATIVE_FAREA17_DEBYE_HUCKEL_CONSTANT_FREE_ION) {
  std::vector<double> ict = { 1.3132E-04, 1.0000E-05, 1.0000E-12, 1.0000E-06,
                              1.0000E-12, 1.0716E-05, 1.0000E-05, 1.0000E-05,
                              6.0081E-05, 1.0000E-05, 1.0000E-06, 1.0000E-06,
                              7.8954E-05, 1.0000E-05, 1.0000E-15, 2.5277E-04, 1.0000E-15 };
  std::vector<double> icm = { 0.0, 0.21, 0.15,
                              0.0, 0.1,  0.15,
                              0.0, 0.0,  0.0,
                              0.0, 0.0};
  std::vector<double> icie;
  std::vector<double> icfi = { 1.0e-9, 1.0e-9, 1.0e-9, 1.0e-9, 1.0e-9, 1.0e-9,
                               1.0e-9, 1.0e-9, 1.0e-9, 1.0e-9, 1.0e-9, 1.0e-9,
                               1.0e-9, 1.0e-9, 1.0e-9, 1.0e-9, 1.0e-9 };
  RunBatchNative("test/native/fbasin-17.xml",
                 "test/native/fbasin-17-debye-huckel.test",
                 "debye-huckel",
                 ict, icm, icie, icfi,  // initial conditions
                 0.5, 1.0, 1.0,
                 2592000.0, 12);  // porosity, saturation, cell volume
}

TEST(NATIVE_FAREA17_PITZEL) {
  std::vector<double> ict = { 1.3132E-04, 1.0000E-05, 1.0000E-12, 1.0000E-06,
                              1.0000E-12, 1.0716E-05, 1.0000E-05, 1.0000E-05,
                              6.0081E-05, 1.0000E-05, 1.0000E-06, 1.0000E-06,
                              7.8954E-05, 1.0000E-05, 1.0000E-15, 2.5277E-04, 1.0000E-15 };
  std::vector<double> icm = { 0.0, 0.21, 0.15,
                              0.0, 0.1,  0.15,
                              0.0, 0.0,  0.0,
                              0.0, 0.0};
  std::vector<double> icie, icfi;
  RunBatchNative("test/native/fbasin-17.xml",
                 "test/native/fbasin-17-pitzel-hwm.test",
                 "pitzer-hwm",
                 ict, icm, icie, icfi,  // initial conditions
                 0.5, 1.0, 1.0);  // porosity, saturation, cell volume
}

TEST(NATIVE_GENERAL_KINETICS) {
  std::vector<double> ict = { 1.0e-4, 2.0e-5};
  std::vector<double> icm, icie, icfi;
  RunBatchNative("test/native/general-reaction.xml",
                 "test/native/general-reaction.test",
                 "unit",
                 ict, icm, icie, icfi,  // initial conditions
                 0.25, 1.0, 1.0,  // porosity, saturation, cell volume
                 8640.0, 500, 5);
}

/*
TEST(NATIVE_GENERAL_KINETICSi_QUADRATIC) {
  std::vector<double> ict = { 1.0e-4, 2.0e-5, 1e-20};
  std::vector<double> icm, icie, icfi;
  RunBatchNative("test/native/general-reaction-quadratic.xml",
                 "test/native/general-reaction-quadrtic.test",
                 "unit",
                 ict, icm, icie, icfi,  // initial conditions
                 0.25, 1.0, 1.0,  // porosity, saturation, cell volume
                 8640.0, 500, 5);
}
*/


TEST(NATIVE_VALOCCHI_INITIAL) {
  std::vector<double> ict = { 8.65e-02, 1.82e-02, 1.11e-02, 0.1451 };
  std::vector<double> icm = { 1.0e-5 };
  std::vector<double> icts = { 1.607889E+02, 1.415668E+02, 1.530388E+02, 0.000000E+00 };  // total sorbed
  std::vector<double> icie = { 750.0 };  // ion exchange
  std::vector<double> icfi, icsd, icssa;  // free ion, side density, specific surface area
  RunBatchNative("test/native/ion-exchange-valocchi.xml",
                 "test/native/ion-exchange-valocchi-initial.test",
                 "unit",
                 ict, icm, icie, icfi,  // initial conditions
                 0.25, 1.0, 1.0,  // porosity, saturation, cell volume
                 86400.0, 10);
}


TEST(NATIVE_SURFACE_COMPLEXATION_1) {
  std::vector<double> ict = { 1.1973159693031387E-05, 9.9987826840306965E-02, 0.1, 1.0e-7 };
  std::vector<double> icm, icie, icfi, icsd, icssa;
  std::vector<double> icts = { 7.585367E+04, 0.0, 0.0, 9.695984E-03 }; // total sorbed
  RunBatchNative("test/native/surface-complexation-1.xml",
                 "test/native/surface-complexation-1.test",
                 "debye-huckel",
                 ict, icm, icie, icfi,  // initial conditions
                 0.9, 1.0, 1.0,  // porosity, saturation, cell volume
                 3600.0, 10);
}

TEST(NATIVE_SURFACE_COMPLEXATION_2) {
  std::vector<double> ict = { 1.1973159693031387E-05, 9.9987826840306965E-02, 0.1, 1.0e-7 };
  std::vector<double> icm, icie, icfi, icsd, icssa;
  std::vector<double> icts = { 7.774868E+04, 0.0, 0.0, 2.410517E-01 }; // total sorbed
  RunBatchNative("test/native/surface-complexation-2.xml",
                 "test/native/surface-complexation-2.test",
                 "debye-huckel",
                 ict, icm, icie, icfi,  // initial conditions
                 0.9, 1.0, 1.0,  // porosity, saturation, cell volume
                 1.0, 10);
}


TEST(NATIVE_RADIOACTIVE_DECAY_AQUEOUS) {
  std::vector<double> ict = { 1.0e-01, 1.0e-20, 1.0e-20, 1.0e-01 };
  std::vector<double> icm, icie, icfi, icts;
  RunBatchNative("test/native/radioactive-decay-aqueous.xml",
                 "test/native/radioactive-decay-aqueous.test",
                 "unit",
                 ict, icm, icie, icfi,  // initial conditions
                 0.25, 1.0, 1.0,  // porosity, saturation, cell volume
                 864.0, 500, 10);
}

TEST(NATIVE_RADIOACTIVE_DECAY_SORBED_TINY) {
  std::vector<double> ict = { 1.0e-20 };
  std::vector<double> icm, icie, icfi;
  std::vector<double> icts = { 1.0e-20 };  // total sorbed
  RunBatchNative("test/native/radioactive-decay-sorbed-tiny.xml",
                 "test/native/radioactive-decay-sorbed-tiny.test",
                 "unit",
                 ict, icm, icie, icfi,  // initial conditions
                 1.0, 1.0, 1.0,  // porosity, saturation, cell volume
                 86400.0, 2000, 50);
}

TEST(NATIVE_RADIOACTIVE_DECAY_SORBED) {
  std::vector<double> ict = { 1.0e-01, 1.0e-20, 1.0e-20 };
  std::vector<double> icm, icie, icfi;
  std::vector<double> icts = { 1.0e-20, 1.0e-20, 1.0e-20 };  // total sorbed
  RunBatchNative("test/native/radioactive-decay-sorbed.xml",
                 "test/native/radioactive-decay-sorbed.test",
                 "unit",
                 ict, icm, icie, icfi,  // initial conditions
                 0.25, 1.0, 1.0,  // porosity, saturation, cell volume
                 864000.0, 2000, 50);
}

TEST(NATIVE_RADIOACTIVE_DECAY_BRANCHES) {
  std::vector<double> ict = { 8.7086e-05, 3.2169e-03, 3.4698e-03, 5.1930e-04,
                              1.7495e-04, 2.5079e-06, 5.5912e-10, 1.9627e-10 };

  std::vector<double> icm, icie, icfi, icts;
  RunBatchNative("test/native/radioactive-decay-branches.xml",
                 "test/native/radioactive-decay-branches.test",
                 "unit",
                 ict, icm, icie, icfi,  // initial conditions
                 0.25, 1.0, 1.0,  // porosity, saturation, cell volume
                 1.0e-1, 1000, 10);
}

TEST(NATIVE_SORPTION_ISOTHERMS) {
  std::vector<double> ict = { 1.0e-4, 1.0e-4, 1.0e-4 };
  std::vector<double> icm, icie, icfi;
  std::vector<double> icts = { 1.0e-4, 1.0e-4, 1.0e-4 };  // total sorbed
  RunBatchNative("test/native/sorption-isotherms.xml",
                 "test/native/sorption-isotherms.test",
                 "unit",
                 ict, icm, icie, icfi,  // initial conditions
                 0.25, 1.0, 1.0,  // porosity, saturation, cell volume
                 86400.0, 500, 10);
}


TEST(NATIVE_FAREA5_INITIAL) {
  std::vector<double> ict = { 6.1426E-09, 2.2923e-08, 1.0000E-06, 1.8703E-04, 1.0000E-15 };
  std::vector<double> icm = { 0.15, 0.21, 0.0 };
  std::vector<double> icie;
  // [total_sorbed] - all zeros 
  std::vector<double> icfi = { 5.2970E-16, 3.2759E-08, 9.9886E-07, 1.8703E-04, 1.7609E-20 };
  RunBatchNative("test/native/uo2-5-component.xml",
                 "test/native/uo2-5-component-initial.test",
                 "debye-huckel",
                 ict, icm, icie, icfi,  // initial conditions
                 0.5, 1.0, 1.0,  // porosity, saturation, cell volume
                 2592000.0, 12);
}

TEST(NATIVE_FAREA5_OUTLET) {
  std::vector<double> ict = { 1.0000E-12, -1.1407E-09, 1.0000E-06, 1.0000E-05, 1.0000E-15 };
  std::vector<double> icm = { 0.15, 0.21, 0.0 };
  std::vector<double> icie;
  // [total_sorbed] - all zeros 
  std::vector<double> icfi = { 5.2970E-16, 3.2759E-08, 9.9886E-07, 1.8703E-04, 1.7609E-20 };
  RunBatchNative("test/native/uo2-5-component.xml",
                 "test/native/uo2-5-component-outlet.test",
                 "debye-huckel",
                 ict, icm, icie, icfi,  // initial conditions
                 0.5, 1.0, 1.0,  // porosity, saturation, cell volume
                 2592000.0, 12);
}

TEST(NATIVE_FAREA5_SOURCE) {
  std::vector<double> ict = { 2.8909E-05, 1.2786E-03, 7.1028E-05, 2.5280E-04, 3.5414E-05 };
  std::vector<double> icm = { 0.15, 0.21, 0.0 };
  std::vector<double> icie;
  // [total_sorbed] - all zeros 
  std::vector<double> icfi = { 5.2970E-16, 3.2759E-08, 9.9886E-07, 1.8703E-04, 1.7609E-20 };
  RunBatchNative("test/native/uo2-5-component-source.xml",
                 "test/native/uo2-5-component-source.test",
                 "debye-huckel",
                 ict, icm, icie, icfi,  // initial conditions
                 0.5, 1.0, 1.0,  // porosity, saturation, cell volume
                 2592000.0, 12);
}


