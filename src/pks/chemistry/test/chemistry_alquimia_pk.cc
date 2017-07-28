/*
  Chemistry PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.
*/

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>
#include <typeinfo>

#include <UnitTest++.h>

#include "Teuchos_RCP.hpp"
#include "XMLParameterListWriter.hh"

// Amanzi
#include "dbc.hh"
#include "errors.hh"
#include "exceptions.hh"
#include "State.hh"

// Amanzi::Chemistry
#include "Alquimia_PK.hh"


TEST(INTERFACE_LIBRARY) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziChemistry;

  auto engine = Teuchos::rcp(new AmanziChemistry::ChemistryEngine("PFloTran", "test/chemistry_alquimia_pk.in"));
  CHECK(engine->NumPrimarySpecies() == 14);
  CHECK(engine->NumAqueousComplexes() == 37);
  CHECK(engine->NumSorbedSpecies() == 14);
  CHECK(engine->NumSurfaceSites() == 1);
  CHECK(engine->NumIonExchangeSites() == 1);
  CHECK(engine->NumIsothermSpecies() == 0);
  CHECK(engine->NumAqueousKinetics() == 1);

  std::vector<std::string> species;
  engine->GetPrimarySpeciesNames(species);
  CHECK(species.size() == 14);
  CHECK(species[0] == "H+" && species[13] == "UO2++");

  std::vector<std::string> minerals;
  engine->GetMineralNames(minerals);
  CHECK(minerals.size() == 8);
  CHECK(minerals[0] == "Quartz" && minerals[7] == "Opal");

  std::vector<std::string> sites;
  engine->GetSurfaceSiteNames(sites);
  CHECK(sites.size() == 1);
  CHECK(sites[0] == ">davis_OH");

  std::vector<std::string> ions;
  engine->GetIonExchangeNames(ions);
  CHECK(ions.size() == 1);
  CHECK(ions[0] == "X-");

  engine->GetIsothermSpeciesNames(species);
  CHECK(species.size() == 0);

  std::vector<std::string> aux;
  engine->GetAuxiliaryOutputNames(aux);
  CHECK(aux.size() == 119);
  CHECK(aux[0] == "pH");

  std::vector<std::string> names;
  engine->GetAqueousKineticNames(names);
  CHECK(names.size() == 1);
  std::cout << "aqueous kinetic names: " << names[0] << " size=" << names[0].size() << std::endl;

  // engine->CreateCondition("background"); 
}  
