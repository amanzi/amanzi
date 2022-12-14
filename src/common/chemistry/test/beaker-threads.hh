/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#ifndef __Beaker_Threads_hpp__
#define __Beaker_Threads_hpp__

#include <string>

#include "Beaker.hh"

int
CommandLineOptions(int argc, char** argv);

void
fbasin_source(Amanzi::AmanziChemistry::Beaker::BeakerState* state);
void
fbasin_aqueous_source(Amanzi::AmanziChemistry::Beaker::BeakerState* state);
void
fbasin_free_ions(Amanzi::AmanziChemistry::Beaker::BeakerState* state);
void
fbasin_minerals(Amanzi::AmanziChemistry::Beaker::BeakerState* state);
void
fbasin_sorbed(Amanzi::AmanziChemistry::Beaker::BeakerState* state);

#endif
