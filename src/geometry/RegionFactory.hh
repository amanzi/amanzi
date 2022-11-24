/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/*
  Nonmember function for creating regions.

  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Rao Garimella
           Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_REGION_FACTORY_HH_
#define AMANZI_REGION_FACTORY_HH_

#include "AmanziTypes.hh"
#include "Teuchos_ParameterList.hpp"

#include "GeometryDefs.hh"

namespace Amanzi {
namespace AmanziGeometry {

class Region;

Teuchos::RCP<Region>
createRegion(const std::string reg_name,
             int reg_id,
             Teuchos::ParameterList& reg_spec,
             const Comm_type& comm);

} // namespace AmanziGeometry
} // namespace Amanzi

#endif
