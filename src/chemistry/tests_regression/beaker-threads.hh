/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef __Beaker_Threads_hpp__
#define __Beaker_Threads_hpp__

#include <string>

#include "beaker.hh"
#include "verbosity.hh"

int CommandLineOptions(int argc, char** argv,
                       amanzi::chemistry::Verbosity* verbosity);

void fbasin_source(amanzi::chemistry::Beaker::BeakerComponents* components);
void fbasin_aqueous_source(amanzi::chemistry::Beaker::BeakerComponents* components);
void fbasin_free_ions(amanzi::chemistry::Beaker::BeakerComponents* components);
void fbasin_minerals(amanzi::chemistry::Beaker::BeakerComponents* components);
void fbasin_sorbed(amanzi::chemistry::Beaker::BeakerComponents* components);

#endif
