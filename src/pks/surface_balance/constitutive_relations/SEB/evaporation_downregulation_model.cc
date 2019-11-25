/*
  The evaporation downregulation via soil resistance model is an algebraic model with dependencies.

  Generated via evaluator_generator with:
Downregulates evaporation from a potential.

    
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "Teuchos_ParameterList.hpp"
#include "dbc.hh"
#include "evaporation_downregulation_model.hh"
#include "seb_physics_funcs.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

// Constructor from ParameterList
EvaporationDownregulationModel::EvaporationDownregulationModel(Teuchos::ParameterList& plist)
{
  InitializeFromPlist_(plist);
}


// Initialize parameters
void
EvaporationDownregulationModel::InitializeFromPlist_(Teuchos::ParameterList& plist)
{
  dess_dz_ = plist.get<double>("dessicated zone thickness [m]", 0.1);
  Clapp_Horn_b_ = plist.get<double>("Clapp and Hornberger b of surface soil [-]", 1.0);
}


// main method
double
EvaporationDownregulationModel::Evaporation(double sg, double poro, double pot_evap) const
{
  return pot_evap / (1. + SEBPhysics::EvaporativeResistanceCoef(sg, poro, dess_dz_, Clapp_Horn_b_));
}

double
EvaporationDownregulationModel::DEvaporationDSaturationGas(double sg, double poro, double pot_evap) const
{
  return 0.;
}

double
EvaporationDownregulationModel::DEvaporationDPorosity(double sg, double poro, double pot_evap) const
{
  return 0.;
}

double
EvaporationDownregulationModel::DEvaporationDPotentialEvaporation(double sg, double poro, double pot_evap) const
{
  return 1. / (1. + SEBPhysics::EvaporativeResistanceCoef(sg, poro, dess_dz_, Clapp_Horn_b_));
}

} //namespace
} //namespace
} //namespace
  
