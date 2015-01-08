/*
  This is the energy component of the ATS and Amanzi codes. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Self-registering factory for TC implementations.
*/

#include "twophase_thermal_conductivity_factory.hh"

// explicity instantitate the static data of Factory<EOS>
template<> 
Amanzi::Utils::Factory<Amanzi::Energy::ThermalConductivityTwoPhase>::map_type* 
Amanzi::Utils::Factory<Amanzi::Energy::ThermalConductivityTwoPhase>::map_;

