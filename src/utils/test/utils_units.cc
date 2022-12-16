/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include "UnitTest++.h"

#include "Units.hh"

using namespace Amanzi;
using namespace Amanzi::Utils;

TEST(UNITS_TIME)
{
  Units units("molar");
  bool flag;

  double t = units.ConvertTime(1.0, "y", "s", flag);
  std::cout << "Time tests:\n  1 y = " << t << " s, flag=" << flag << std::endl;
  CHECK_CLOSE(t, 31557600.0, 0.1);

  t = units.ConvertTime(1.0, "d", "s", flag);
  std::cout << "  1 d = " << t << " s, flag=" << flag << std::endl;
  CHECK_CLOSE(t, 86400.0, 0.1);

  t = units.ConvertTime(1.0, "y", "d", flag);
  std::cout << "  1 y = " << t << " d, flag=" << flag << std::endl;
  CHECK_CLOSE(t, 365.25, 1e-3);

  t = units.ConvertTime(1.0, "noleap", "d", flag);
  std::cout << "  1 noleap = " << t << " d, flag=" << flag << std::endl;
  CHECK_CLOSE(t, 365., 1e-3);

  CHECK(units.IsValidTime("y"));
  CHECK(units.IsValidTime("noleap"));
  CHECK(!units.IsValidTime("yr"));
  CHECK(!units.IsValidTime("m"));
  std::cout << "Valid times are: " << units.ValidTimeStrings() << std::endl;
}


TEST(UNITS_LENGTH)
{
  Units units("molar");
  bool flag;

  double len = units.ConvertLength(1.0, "m", "cm", flag);
  std::cout << "Length tests:\n  1 m = " << len << " cm, flag=" << flag << std::endl;
  CHECK_CLOSE(len, 100.0, 1e-8);

  len = units.ConvertLength(1.0, "ft", "m", flag);
  std::cout << "  1 ft = " << len << " m, flag=" << flag << std::endl;
  CHECK_CLOSE(len, 0.3048, 2e-5);

  len = units.ConvertLength(1.0, "ft", "in", flag);
  std::cout << "  1 ft = " << len << " in, flag=" << flag << std::endl;
  CHECK_CLOSE(len, 12.0, 1e-5);

  CHECK(units.IsValidLength("ft"));
  CHECK(units.IsValidLength("m"));
  CHECK(!units.IsValidLength("meter"));
  std::cout << "Valid lengths are: " << units.ValidLengthStrings() << std::endl;
}


TEST(UNITS_CONCENTRATION)
{
  Units units("molar");
  bool flag;

  double conc = units.ConvertConcentration(1.0, "SI", "molar", 1.0, flag);
  std::cout << "Concentration tests:\n  1 mol/m^3 = " << conc << " molar, flag=" << flag
            << std::endl;
  CHECK_CLOSE(conc, 1e-3, 1e-4);

  conc = units.ConvertConcentration(1.0, "molar", "SI", 1.0, flag);
  std::cout << "  1 molar = " << conc << " mol/m^3, flag=" << flag << std::endl;
  CHECK_CLOSE(conc, 1e+3, 1e-4);

  conc = units.ConvertConcentration(1.0, "ppm", "SI", 51.9961e-3, flag);
  std::cout << "  1 ppm = " << conc << " mol/m^3, flag=" << flag << std::endl;
  // CHECK_CLOSE(conc, 1e+3, 1e-4);

  conc = units.ConvertConcentration(1.0, "molar", "ppb", 51.9961e-3, flag);
  std::cout << "  1 molar = " << conc << " ppb, flag=" << flag << std::endl;
  // CHECK_CLOSE(conc, 1e+3, 1e-4);

  CHECK(units.IsValidConcentration("molar"));
  CHECK(units.IsValidConcentration("SI"));
  CHECK(!units.IsValidConcentration("y"));
  std::cout << "Valid concentrations are: " << units.ValidConcentrationStrings() << std::endl;
}


TEST(UNITS_DERIVED_DOUBLE)
{
  Units units("molar");
  bool flag;

  std::string in_unit("m/d"), out_unit("ft/y");
  double tmp = units.ConvertUnitD(1.0, in_unit, out_unit, 51.9961e-3, flag);
  std::cout << "Derived tests:\n  1 m/d = " << tmp << " ft/y, flag=" << flag << std::endl;
  CHECK_CLOSE(tmp, 1198.326, 1e-3);

  out_unit = "ft/y/m";
  tmp = units.ConvertUnitD(1.0, in_unit, out_unit, 51.9961e-3, flag);
  CHECK(!flag);

  in_unit = "g*m/s^2";
  out_unit = "kg*in/h^2";
  tmp = units.ConvertUnitD(1.0, in_unit, out_unit, 51.9961e-3, flag);
  std::cout << "  1 g*m/s^2 = " << tmp << " kg*in/h^2, flag=" << flag << std::endl;
  CHECK_CLOSE(tmp, 5.10236e+05, 1.0);

  in_unit = "g*s/m^3";
  out_unit = "kg*h/L";
  tmp = units.ConvertUnitD(1.0, in_unit, out_unit, 51.9961e-3, flag);
  std::cout << "  1 g*s/m^3 = " << tmp << " kg*h/L, flag=" << flag << std::endl;
  CHECK_CLOSE(tmp, 2.77778e-10, 1.0e-15);

  in_unit = "Pa/s";
  out_unit = "kg/m/s^3";
  tmp = units.ConvertUnitD(2.0, in_unit, out_unit, 51.9961e-3, flag);
  std::cout << "  2 Pa/s = " << tmp << " kg/m/s^3, flag=" << flag << std::endl;
  CHECK_CLOSE(tmp, 2.0, 1.0e-10);

  in_unit = "J/s";
  out_unit = "kg*m^2/s^3";
  tmp = units.ConvertUnitD(3.0, in_unit, out_unit, 1.0, flag);
  std::cout << "  3 J/s = " << tmp << " kg*m^2/s^3, flag=" << flag << std::endl;
  CHECK_CLOSE(tmp, 3.0, 1.0e-10);
}


TEST(UNITS_DERIVED_STRING)
{
  Units units("molar");
  std::string out_unit;

  {
    UnitsSystem system("h", "kg", "m", "molar", "mol", "K");
    std::string in_unit = "g/s";
    out_unit = units.ConvertUnitS(in_unit, system);
    std::cout << "Derived tests:\n  g/s -> " << out_unit << std::endl;
    CHECK(out_unit == "h^-1*kg");
  }

  {
    UnitsSystem system("s", "kg", "m", "molar", "mol", "K");
    std::string in_unit = "Pa/s";
    out_unit = units.ConvertUnitS(in_unit, system);
    std::cout << "  Pa/s -> " << out_unit << std::endl;
    CHECK(out_unit == "kg*m^-1*s^-3");
  }
}


TEST(UNITS_FANCY_OUTPUT)
{
  Units units("molar");

  for (double val = 1e-3; val < 1e+8; val *= 4) {
    std::cout << val << " res=" << units.OutputTime(val) << std::endl;
  }
}
