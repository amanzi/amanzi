/*
  Multi-Process Coordinator

  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy(dasvyat@lanl.gov)
           Ethan Coon (ecoon@lanl.gov)

  A factory for creating observations.
*/

#include <iostream>
#include <sstream>

#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_StrUtils.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "dbc.hh"
#include "errors.hh"

#include "ObservableAmanzi.hh"

#include "ObservableFactory.hh"
#include "ObservableAqueous.hh"
#include "ObservableLineSegment.hh"
#include "ObservableLineSegmentAqueous.hh"
#include "ObservableLineSegmentSolute.hh"
#include "ObservableSolute.hh"

namespace Amanzi {

Teuchos::RCP<Observable>
CreateObservable(Teuchos::ParameterList& coord_plist,
                 Teuchos::ParameterList& observable_plist,
                 Teuchos::ParameterList& units_plist,
                 Teuchos::RCP<const AmanziMesh::Mesh> mesh)
{
  Teuchos::RCP<Observable> observe;
  std::string var = observable_plist.get<std::string>("variable");
  std::string func = observable_plist.get<std::string>("functional");
  std::string region = observable_plist.get<std::string>("region");

  Teuchos::Array<std::string> comp_names;
  int num_liquid = 0;

  if (coord_plist.isParameter("component names")) {
    comp_names = coord_plist.get<Teuchos::Array<std::string>>("component names");
  }
  num_liquid = coord_plist.get<int>("number of liquid components", comp_names.size());

  // check if observation of solute was requested
  bool obs_solute_liquid(false), obs_solute_gas(false), obs_aqueous(true);
  int tcc_index(-1);
  for (tcc_index = 0; tcc_index != comp_names.size(); ++tcc_index) {
    int pos = var.find(comp_names[tcc_index]);
    if (pos == 0) {
      (tcc_index < num_liquid) ? obs_solute_liquid = true : obs_solute_gas = true;
      obs_aqueous = false;
      break;
    }
  }

  bool obs_solute = obs_solute_liquid || obs_solute_gas;

  Teuchos::RCP<const AmanziGeometry::GeometricModel> gm_ptr = mesh->geometric_model();
  Teuchos::RCP<const AmanziGeometry::Region> reg_ptr = gm_ptr->FindRegion(region);

  if (obs_solute) {
    if (reg_ptr->get_type() == AmanziGeometry::RegionType::LINE_SEGMENT) {
      observe = Teuchos::rcp(
        new ObservableLineSegmentSolute(var, region, func, observable_plist, units_plist, mesh));
    } else {
      observe =
        Teuchos::rcp(new ObservableSolute(var, region, func, observable_plist, units_plist, mesh));
    }
    Teuchos::rcp_dynamic_cast<ObservableSolute>(observe)->RegisterComponentNames(
      comp_names.toVector(), num_liquid, tcc_index);

  } else if (obs_aqueous) {
    if (reg_ptr->get_type() == AmanziGeometry::RegionType::LINE_SEGMENT) {
      observe = Teuchos::rcp(
        new ObservableLineSegmentAqueous(var, region, func, observable_plist, units_plist, mesh));
    } else {
      observe =
        Teuchos::rcp(new ObservableAqueous(var, region, func, observable_plist, units_plist, mesh));
    }
  }

  int n = observe->ComputeRegionSize();

  if (n == 0) {
    Errors::Message msg;
    msg << "Observation: variable \"" << var << "\" has no assigned mesh objects.";
    Exceptions::amanzi_throw(msg);
  }

  return observe;
}

} // namespace Amanzi
