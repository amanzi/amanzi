/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef __Beaker_Threads_hpp__
#define __Beaker_Threads_hpp__

#include <string>

#include "beaker.hh"
#include "verbosity.hh"

int CommandLineOptions(int argc, char **argv,
                       Verbosity* verbosity);

void fbasin_source(Beaker::BeakerComponents* components);
void fbasin_aqueous_source(Beaker::BeakerComponents* components);
void fbasin_free_ions(Beaker::BeakerComponents* components);
void fbasin_minerals(Beaker::BeakerComponents* components);
void fbasin_sorbed(Beaker::BeakerComponents* components);

#endif 
