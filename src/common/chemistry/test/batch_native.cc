#include <unistd.h>

//#define ABORT_ON_FLOATING_POINT_EXCEPTIONS
#ifdef __APPLE__
  #include <xmmintrin.h>
#endif

#include <cstdlib>
#include <cctype>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>

// TPLs
#include "boost/algorithm/string.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include <UnitTest++.h>

// Chemistry
#include "simple_thermo_database.hh"
#include "beaker.hh"
#include "activity_model_factory.hh"
#include "chemistry_utilities.hh"
#include "chemistry_exception.hh"

namespace ac = Amanzi::AmanziChemistry;
int BUFFER_SIZE = 100000;

int CompareFiles(const std::string& file1, const std::string& file2)
{
  std::ifstream ifs1(file1.c_str(), std::ios::in | std::ios::binary);
  std::ifstream ifs2(file2.c_str(), std::ios::in | std::ios::binary);
  if(!ifs1.good() || !ifs2.good()) return 1;

  char *buffer1 = new char[BUFFER_SIZE];
  char *buffer2 = new char[BUFFER_SIZE];

  do {
    ifs1.read(buffer1, BUFFER_SIZE);
    ifs2.read(buffer2, BUFFER_SIZE);
    std::streamsize count1 = ifs1.gcount();
    std::streamsize count2 = ifs2.gcount();

    if (count1 != count2 || std::memcmp(buffer1, buffer2, count1) != 0) return 2;
  } while (ifs1.good() || ifs2.good());

  return 0;
}


void RunBatchNative(const std::string& filexml,
                    const std::string& filetest,
                    const std::string& activity_model,
                    const std::vector<double>& ict,
                    const std::vector<double>& icm,
                    double porosity, double saturation, double volume,
                    double dt = 0.0, int max_dt_steps = 0)
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
  auto plist = Teuchos::getParametersFromXmlFile(filexml);
  auto vo = Teuchos::rcp(new Amanzi::VerboseObject("Beaker", *plist));

  ac::Beaker::BeakerState state;

  ac::Beaker* chem = new ac::SimpleThermoDatabase(plist, vo);

  ac::Beaker::BeakerParameters parameters;
  parameters.tolerance = 1e-12;
  parameters.max_iterations = 100;
  parameters.activity_model_name = activity_model;

  state.porosity = porosity;
  state.saturation = saturation;
  state.water_density = 997.16;
  state.volume = volume;

  chem->Initialize(parameters);

  state.mineral_volume_fraction = icm;

  chem->CopyStateToBeaker(state);

  // we do not have external state in this test, we need to initialize 
  // chemistry state from beaker's data
  chem->CopyBeakerToState(&state);

  int ncomp = chem->primary_species().size();
  int nmineral = chem->minerals().size();
  int nsorbed = chem->total_sorbed().size();
  int nion_site = chem->ion_exchange_rxns().size();
  int nisotherm = chem->sorption_isotherm_rxns().size();

  state.total = ict;
  state.free_ion.resize(ncomp, 1.0e-9);

  if (nmineral > 0) {
    state.mineral_volume_fraction.resize(nmineral, 0.0);
    state.mineral_specific_surface_area.resize(nmineral, 0.0);
  }
  if (nsorbed > 0)
    state.total_sorbed.resize(nsorbed, 0.0);

  if (nion_site > 0)
    state.ion_exchange_sites.resize(nion_site, 0.0);

  if (nisotherm > 0) {
    state.surface_site_density.resize(nisotherm, 0.0);
    state.isotherm_kd.resize(ncomp, 0.0);
    state.isotherm_freundlich_n.resize(ncomp, 0.0);
    state.isotherm_langmuir_b.resize(ncomp, 0.0);
  }

  // io
  chem->Display();
  chem->DisplayComponents(state);

  // solve for free-ion concentrations
  chem->Speciate(&state);

  chem->CopyBeakerToState(&state);
  chem->DisplayResults();

  // kinetics
  if (dt > 0.0) {
    double time(0.0);

    for (int n = 0; n < max_dt_steps; ++n) {
      chem->ReactionStep(&state, dt);
      // chem->CopyBeakerToState(&state);
      time += dt;
      chem->DisplayTotalColumns(time, state, false);
    }
    chem->Speciate(&state);
    chem->DisplayResults();
  }

  vo = Teuchos::null;  // closing the stream
  std::string tmp = plist->sublist("verbose object").get<std::string>("output filename");
  int ok = CompareFiles(tmp, filetest);
  CHECK(ok == 0);

  // cleanup memory
  delete chem;
}


TEST(NATIVE_CA_DEBYE_HUCKEL) {
  std::vector<double> ict = {3.0e-3, 1.0e-3, 1.0e-3};
  std::vector<double> icm;
  RunBatchNative("test/native/ca-carbonate.xml",
                 "test/native/ca-carbonate-debye-huckel.test",
                 "debye-huckel",
                 ict, icm,  // initial conditions
                 0.5, 1.0, 1.0);  // porosity, saturation, cell volume
}

TEST(NATIVE_CA_UNIT) {
  std::vector<double> ict = {3.0e-3, 1.0e-3, 1.0e-3};
  std::vector<double> icm;
  RunBatchNative("test/native/ca-carbonate.xml",
                 "test/native/ca-carbonate-unit.test",
                 "unit",
                 ict, icm,  // initial conditions
                 0.5, 1.0, 1.0);  // porosity, saturation, cell volume
}


TEST(NATIVE_CALCITE_KINETICS) {
  std::vector<double> ict = {-1.0e-5, 1.0e-5, 1.0e-5};
  std::vector<double> icm = {0.2};
  RunBatchNative("test/native/calcite.xml",
                 "test/native/calcite-kinetics.test",
                 "debye-huckel",
                 ict, icm,  // initial conditions
                 0.5, 1.0, 1.0);  // porosity, saturation, cell volume
}

TEST(NATIVE_CALCITE_KINETICS_VOLUME_FRACTIONS) {
  std::vector<double> ict = {1.0e-2, 1.0e-2, 1.0e-19};
  std::vector<double> icm = {0.2};
  RunBatchNative("test/native/calcite.xml",
                 "test/native/calcite-kinetics-volume-fractions.test",
                 "debye-huckel",
                 ict, icm,  // initial conditions
                 0.5, 1.0, 1.0,  // porosity, saturation, cell volume, dt
                 2592000.0, 60);  // dt, max time steps
}


TEST(NATIVE_CARBONATE_DEBYE_HUCKEL) {
  std::vector<double> ict = {1.0e-3, 1.0e-3};
  std::vector<double> icm;
  RunBatchNative("test/native/carbonate.xml",
                 "test/native/carbonate-debye-huckel.test",
                 "debye-huckel",
                 ict, icm,  // initial conditions
                 0.5, 1.0, 1.0);  // porosity, saturation, cell volume
}

TEST(NATIVE_CARBONATE_UNIT) {
  std::vector<double> ict = {1.0e-3, 1.0e-3};
  std::vector<double> icm;
  RunBatchNative("test/native/carbonate.xml",
                 "test/native/carbonate-debye-huckel.test",
                 "unit",
                 ict, icm,  // initial conditions
                 0.5, 1.0, 1.0);  // porosity, saturation, cell volume
}

