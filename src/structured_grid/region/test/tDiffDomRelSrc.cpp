/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include <iostream>
#include <fstream>
#include <iomanip>

// closing DSO objects
#include "VerboseObject_objs.hh"

using std::cout;
using std::endl;
#include <ParmParse.H>
#include <VisMF.H>

#include <RegionManager.H>
#include <DiffDomRelSrc.H>

int
main (int   argc,
      char* argv[])
{
  BoxLib::Initialize(argc,argv);

  RegionManager rm;

  Array<Real> plo(BL_SPACEDIM,0);
  Array<Real> dx(BL_SPACEDIM,1);
  const Array<const Region*> regions = rm.RegionPtrArray();

  Real mixingLength = 2;
  Real Deff = 9;
  Real totalInventory = 10;
  Real timeScale = 1.e3;
  Real startTime = 0;
  Real endTime = 2*timeScale;

  DiffDomRelSrc rd("myLabel",regions,"all",mixingLength,Deff,totalInventory,startTime,endTime,timeScale);

  Real Pi = 2*std::asin(1.);

  Real t1 = timeScale;
  Real t2 = 2 * t1;
  Real t3 = 4 * t1;

  Real Qav01 = 2.*totalInventory/mixingLength*std::sqrt(Deff/(timeScale*Pi));
  Real Q1 = 2 * (totalInventory/mixingLength) * std::sqrt(Deff*t1/Pi);
  Real Q2 = 2 * (totalInventory/mixingLength) * std::sqrt(Deff*t2/Pi);
  Real Qav12 = (Q2 - Q1)/(t2 - t1);

  Real eps = 1.e-8 * (Qav01 + Qav12);

  Real diff01 = std::abs(Qav01 - rd(0,timeScale)[0]);
  Real diff12 = std::abs(Qav12 - rd(t1,t2)[0]);
  Real diff02 = std::abs(0.5*(Qav12 + Qav01) - rd(0,t2)[0]);
  Real diff03 = std::abs((Qav12 + Qav01)/4. - rd(0,t3)[0]);

  bool ret = (diff01 > eps
              || diff12 > eps
              || diff02 > eps
              || diff03 > eps);

  return ret;
}
