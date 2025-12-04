/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_MODEL_ASSUMPTIONS_HH_
#define AMANZI_MODEL_ASSUMPTIONS_HH_

#include "Teuchos_ParameterList.hpp"

namespace Amanzi {

struct ModelAssumptions {
 public:
  ModelAssumptions()
    : flow_on_manifold(false),
      vapor_diffusion(false),
      use_ppm(false),
      mu_model(""),
      msm_name(""),
      msm_porosity(false),
      undrained_split(false),
      fixed_stress_split(false),
      poroelasticity(false),
      thermoelasticity(false),
      external_aperture(false),
      use_bulk_modulus(false),
      use_overburden_stress(false),
      use_transport_porosity(false),
      use_effective_diffusion(false),
      use_dispersion(true),
      abs_perm(false),
      vol_flowrate(""),
      sat_liquid(""),
      lookup_table(""),
      use_gravity(true),
      L_scheme(false),
      L_scheme_key("") {};

  void Init(Teuchos::ParameterList plist, const AmanziMesh::Mesh& mesh)
  {
    flow_on_manifold = plist.get<bool>("flow and transport in fractures", false);
    flow_on_manifold &= (mesh.getManifoldDimension() != mesh.getSpaceDimension());

    vapor_diffusion = plist.get<bool>("vapor diffusion", false);

    use_ppm = plist.get<bool>("permeability porosity model", false);
    mu_model = plist.get<std::string>("viscosity model", "constant viscosity");

    msm_name = plist.get<std::string>("multiscale model", "single continuum");
    msm_porosity = (msm_name != "single continuum");

    undrained_split = plist.get<bool>("biot scheme: undrained split", false);
    fixed_stress_split = plist.get<bool>("biot scheme: fixed stress split", false);
    poroelasticity = undrained_split || fixed_stress_split;
    thermoelasticity = plist.get<bool>("thermoelasticity", false);

    external_aperture = plist.get<bool>("external aperture", false);
    use_bulk_modulus = plist.get<bool>("use bulk modulus", false);
    use_overburden_stress = plist.get<bool>("use overburden stress", false);

    use_transport_porosity = plist.get<bool>("effective transport porosity", false);
    use_effective_diffusion = plist.get<bool>("effective transport diffusion", false);
    use_dispersion = plist.get<bool>("use dispersion solver", true);
    abs_perm = plist.get<bool>("permeability field is required", false);

    vol_flowrate = plist.get<std::string>("volumetric flow rate key", "volumetric_flow_rate");
    sat_liquid = plist.get<std::string>("saturation key", "saturation_liquid");

    lookup_table = plist.get<std::string>("eos lookup table", "");

    use_gravity = plist.get<bool>("use gravity", true);

    L_scheme = plist.get<bool>("L-scheme stabilization", false);
    L_scheme_key = plist.get<std::string>("L-scheme key", "");
  }

 public:
  bool flow_on_manifold; // true for a DFN model
  bool vapor_diffusion;

  bool use_ppm;
  std::string mu_model;

  std::string msm_name;
  bool msm_porosity;

  bool undrained_split, fixed_stress_split;
  bool poroelasticity, thermoelasticity;

  bool external_aperture;
  bool use_bulk_modulus;
  bool use_overburden_stress;

  bool use_transport_porosity;
  bool use_effective_diffusion;
  bool use_dispersion;
  bool abs_perm;

  std::string vol_flowrate;
  std::string sat_liquid;

  std::string lookup_table;

  bool use_gravity;

  bool L_scheme;
  std::string L_scheme_key;
};

} // namespace Amanzi

#endif
