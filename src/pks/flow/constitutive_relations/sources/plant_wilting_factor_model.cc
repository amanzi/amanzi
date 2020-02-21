/*
  The plant wilting factor model is an algebraic model with dependencies.

  Generated via evaluator_generator with:
Wilting factor.

Beta, or the water availability factor, or the plant wilting factor.

Beta =  (p_closed - p) / (p_closed - p_open)

where p is the capillary pressure or soil mafic potential, and closed
and open indicate the values at which stomates are fully open or fully
closed (the wilting point).


    
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "Teuchos_ParameterList.hpp"
#include "dbc.hh"
#include "plant_wilting_factor_model.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

// Constructor from ParameterList
PlantWiltingFactorModel::PlantWiltingFactorModel(Teuchos::ParameterList& plist)
{
  InitializeFromPlist_(plist);
}


// Initialize parameters
void
PlantWiltingFactorModel::InitializeFromPlist_(Teuchos::ParameterList& plist)
{
  pc_o_ = plist.get<double>("capillary pressure at fully open stomates [Pa]");
  pc_c_ = plist.get<double>("capillary pressure at wilting point [Pa]");
}


// main method
double
PlantWiltingFactorModel::PlantWiltingFactor(double pc) const
{
  return pc_c_ < pc ? 0. : (pc < pc_o_ ? 1. : ((-pc + pc_c_)/(pc_c_ - pc_o_)));
}

double
PlantWiltingFactorModel::DPlantWiltingFactorDCapillaryPressureGasLiq(double pc) const
{
  return pc_c_ < pc ? 0. : (pc < pc_o_ ? 0. : (-1/(pc_c_ - pc_o_)));
}

} //namespace
} //namespace
} //namespace
  
