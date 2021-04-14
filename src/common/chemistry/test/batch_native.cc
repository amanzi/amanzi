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


void RunBatchNative(const std::string& filename,
                    const std::string& activity_model,
                    const std::vector<double>& ic,
                    double porosity, double saturation, double volume)
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
  auto plist = Teuchos::getParametersFromXmlFile(filename + ".xml");
  auto vo = Teuchos::rcp(new Amanzi::VerboseObject("Beaker", *plist));

  ac::Beaker::BeakerState state;

  ac::Beaker* chem = new ac::SimpleThermoDatabase(plist, vo);

  ac::Beaker::BeakerParameters parameters = chem->GetDefaultParameters();
  parameters.thermo_database_file = "";
  parameters.activity_model_name = activity_model;
  parameters.max_iterations = 100;
  parameters.tolerance = 1e-12;
  parameters.porosity = porosity;
  parameters.saturation = saturation;
  parameters.water_density = 997.16;
  parameters.volume = volume;

  chem->Setup(state, parameters);

  int ncomp = chem->primary_species().size();
  int nmineral = chem->minerals().size();
  int nsorbed = chem->total_sorbed().size();
  int nion_site = chem->ion_exchange_rxns().size();
  int nisotherm = chem->sorption_isotherm_rxns().size();

  state.total = ic;
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
  chem->Speciate(&state, parameters);

  chem->CopyBeakerToState(&state);
  chem->DisplayResults();

  vo = Teuchos::null;  // closing the stream
  std::string tmp = plist->sublist("verbose object").get<std::string>("output filename");
  int ok = CompareFiles(tmp, filename + ".test");
  CHECK(ok == 0);

  // cleanup memory
  delete chem;
}

TEST(NATIVE_CHEMISTRY) {
  std::vector<double> ic = {3.0e-3, 1.0e-3, 1.0e-3};
  RunBatchNative("test/native/ca-carbonate",
                 "debye-huckel",
                 ic,  // initial conditions
                 0.5, 1.0, 1.0);  // porosity, saturation, cell volume
}

