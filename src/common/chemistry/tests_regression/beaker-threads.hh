/* -*-  mode: c++; indent-tabs-mode: nil -*- */
#ifndef __Beaker_Threads_hpp__
#define __Beaker_Threads_hpp__

#include <string>

#include "beaker.hh"
#include "chemistry_verbosity.hh"

int CommandLineOptions(int argc, char** argv,
                       Amanzi::AmanziChemistry::Verbosity* verbosity);

void fbasin_source(Amanzi::AmanziChemistry::Beaker::BeakerComponents* components);
void fbasin_aqueous_source(Amanzi::AmanziChemistry::Beaker::BeakerComponents* components);
void fbasin_free_ions(Amanzi::AmanziChemistry::Beaker::BeakerComponents* components);
void fbasin_minerals(Amanzi::AmanziChemistry::Beaker::BeakerComponents* components);
void fbasin_sorbed(Amanzi::AmanziChemistry::Beaker::BeakerComponents* components);

#endif
