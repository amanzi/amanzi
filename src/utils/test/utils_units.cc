#include "UnitTest++.h"

#include "Units.hh"

using namespace Amanzi;
using namespace Amanzi::Utils;

using namespace boost::units;

TEST(UNITS_TIME) 
{
  Units units("molar");
  bool flag;

  double t = units.ConvertTime(1.0, "y", "s", flag);
  std::cout << "1 y = " << t << " s, flag=" << flag << std::endl;
  CHECK_CLOSE(t, 31557600.0, 0.1);

  t = units.ConvertTime(1.0, "d", "s", flag);
  std::cout << "1 d = " << t << " s, flag=" << flag << std::endl;
  CHECK_CLOSE(t, 86400.0, 0.1);

  t = units.ConvertTime(1.0, "y", "d", flag);
  std::cout << "1 y = " << t << " d, flag=" << flag << std::endl;
  CHECK_CLOSE(t, 365.25, 1e-3);
}


TEST(UNITS_LENGTH) 
{
  Units units("molar");
  bool flag;

  double len = units.ConvertLength(1.0, "m", "cm", flag);
  std::cout << "1 m = " << len << " cm, flag=" << flag << std::endl;
  CHECK_CLOSE(len, 100.0, 1e-8);

  len = units.ConvertLength(1.0, "ft", "m", flag);
  std::cout << "1 ft = " << len << " m, flag=" << flag << std::endl;
  CHECK_CLOSE(len, 0.3048, 2e-5);

  len = units.ConvertLength(1.0, "ft", "in", flag);
  std::cout << "1 ft = " << len << " in, flag=" << flag << std::endl;
  CHECK_CLOSE(len, 12.0, 1e-5);
}


TEST(UNITS_CONCENTRATION) 
{
  Units units("molar");
  bool flag;

  double conc = units.ConvertConcentration(1.0, "mol/m^3", "molar", 1.0, flag);
  std::cout << "1 mol/m^3 = " << conc << " molar, flag=" << flag << std::endl;
  CHECK_CLOSE(conc, 1e-3, 1e-4);

  conc = units.ConvertConcentration(1.0, "molar", "mol/m^3", 1.0, flag);
  std::cout << "1 molar = " << conc << " mol/m^3, flag=" << flag << std::endl;
  CHECK_CLOSE(conc, 1e+3, 1e-4);

  conc = units.ConvertConcentration(1.0, "ppm", "mol/m^3", 51.9961e-3, flag);
  std::cout << "1 ppm = " << conc << " mol/m^3, flag=" << flag << std::endl;
  // CHECK_CLOSE(conc, 1e+3, 1e-4);

  conc = units.ConvertConcentration(1.0, "molar", "ppb", 51.9961e-3, flag);
  std::cout << "1 molar = " << conc << " ppb, flag=" << flag << std::endl;
  // CHECK_CLOSE(conc, 1e+3, 1e-4);
}


TEST(UNITS_DERIVED) 
{
  Units units("molar");
  bool flag;

  std::string derived("m/s");
  AtomicUnitForm form = units.ComputeAtomicUnitForm_(derived, &flag);
 
  for (AtomicUnitForm::iterator it = form.begin(); it != form.end(); ++it) {
    std::cout << "unit: " << it->first << " degree: " << it->second << std::endl;
  }
}

