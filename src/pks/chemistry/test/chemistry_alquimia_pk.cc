/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*
  Chemistry PK

*/

#include <cmath>
#include <vector>

// TPLs
#include "Teuchos_RCP.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include <UnitTest++.h>
#include "XMLParameterListWriter.hh"

// Amanzi
#include "dbc.hh"
#include "errors.hh"
#include "exceptions.hh"
#include "MeshFactory.hh"
#include "State.hh"

// Amanzi::Chemistry
#include "Alquimia_PK.hh"


TEST(INTERFACE_LIBRARY_INIT)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziChemistry;

  auto engine =
    Teuchos::rcp(new AmanziChemistry::ChemistryEngine("PFloTran", "test/chemistry_alquimia_pk.in"));
  CHECK(engine->Name() == "PFloTran");

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
  std::vector<std::vector<std::string>> aux_subfields;
  engine->GetAuxiliaryOutputNames(aux, aux_subfields);
  CHECK(aux.size() == 7);
  int count = 0;
  for (const auto& sf : aux_subfields) count += sf.size();
  CHECK(count == 119);
  CHECK(aux[0] == "pH");

  std::vector<std::string> names;
  engine->GetAqueousKineticNames(names);
  CHECK(names.size() == 1);
}


TEST(INTERFACE_LIBRARY_ADVANCE)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziChemistry;

  auto engine =
    Teuchos::rcp(new AmanziChemistry::ChemistryEngine("PFloTran", "test/chemistry_alquimia_pk.in"));

  AlquimiaState state, state_tmp;
  AlquimiaProperties mat_props;
  AlquimiaAuxiliaryData aux_data;
  AlquimiaAuxiliaryOutputData aux_output;

  engine->InitState(mat_props, state_tmp, aux_data, aux_output);
  engine->InitState(mat_props, state, aux_data, aux_output);

  mat_props.volume = 0.1;
  mat_props.saturation = 1.0;
  // mat_props.aqueous_kinetic_rate_cnst.data[0] = 1.78577e-09;
  {
    std::vector<double> data(
      { 1.3345e+05, -7.94e+04, -1.2967e+05, 2.0e+04, -1.15e+05, -8.0e+04, -8.0e+04, -1.2135e+05 });
    for (int i = 0; i < 8; ++i) mat_props.mineral_rate_cnst.data[i] = data[i];
  }

  state.water_density = 998.0;
  state.porosity = 0.2;
  state.cation_exchange_capacity.data[0] = 2.750;
  state.surface_site_density.data[0] = 1.56199e-01;

  {
    std::vector<double> data({ 2.95786307706898224e-06,
                               2.21343191783434494e-08,
                               1.0e-05,
                               9.98922451603534434e-03,
                               2.52312817040197963e-16,
                               1.21461876452330652e-05,
                               3.32e-05,
                               5.35e-03,
                               2.78e-04,
                               1.77282139025972417e-04,
                               2.25e-05,
                               1.0e-15,
                               1.0e-03,
                               1.25e-10 });
    for (int i = 0; i < 14; ++i) state.total_mobile.data[i] = data[i];
    for (int i = 0; i < 14; ++i) state.total_immobile.data[i] = 0.0;
  }

  {
    std::vector<double> data({ 8.8e-01, 1.6e-02, 1.1e-01, 0.0, 0.0, 0.0, 0.0 });
    for (int i = 0; i < 8; ++i) state.mineral_volume_fraction.data[i] = data[i];
  }

  {
    std::vector<double> data(
      { 3.2623e+05, 1.10763e+06, 5.90939e+06, 10.0, 10.0, 10.0, 10.0, 10.0 });
    for (int i = 0; i < 8; ++i) state.mineral_specific_surface_area.data[i] = data[i];
  }

  CopyAlquimiaState(&state, &state_tmp);
  engine->EnforceCondition("background", 0.0, mat_props, state, aux_data, aux_output);

  // fancy output
  std::vector<std::string> species;
  engine->GetPrimarySpeciesNames(species);

  for (int i = 0; i < 14; ++i) {
    double v0 = state_tmp.total_mobile.data[i];
    double v1 = state.total_mobile.data[i];
    double diff = 200 * fabs(v0 - v1) / (fabs(v0) + fabs(v1) + 1e-30);

    CHECK(diff < 0.1);
    printf("%10s  %14.6g  -> %12.6g   diff:%5.2f\n", species[i].c_str(), v0, v1, diff);
  }
}


TEST(INITIALIZE_CRUNCH)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziChemistry;
  using namespace Amanzi::AmanziMesh;

  // read parameter list
  std::string xmlFileName = "test/chemistry_alquimia_pk.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  // create a mesh framework
  Comm_ptr_type comm = Amanzi::getDefaultComm();
  Teuchos::ParameterList regions_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, regions_list, *comm));

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  Teuchos::RCP<const Mesh> mesh = meshfactory.create(0.0, 0.0, 0.0, 100.0, 1.0, 1.0, 100, 1, 1);

  Teuchos::ParameterList state_list = plist->sublist("state");
  auto S = Teuchos::rcp(new State(state_list));
  S->RegisterDomainMesh(Teuchos::rcp_const_cast<Mesh>(mesh));

  auto pks_list = plist->sublist("PK tree").sublist("chemistry");
  Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());
  auto CPK = Teuchos::rcp(new Alquimia_PK(pks_list, plist, S, soln));
  CPK->Setup();

  S->Setup();
  S->InitializeFields();
  S->InitializeEvaluators();
  S->set_time(0.0);

  CPK->Initialize();

  // Now create second chemistry PK
  auto CPK2 = Teuchos::rcp(new Alquimia_PK(pks_list, plist, S, soln));
  CPK2->Setup();
  CPK2->Initialize();
}
