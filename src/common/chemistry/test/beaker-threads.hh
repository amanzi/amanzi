/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
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
