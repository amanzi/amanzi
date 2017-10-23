/*
  The overland subgrid water content model is an algebraic model with dependencies.

  Generated via evaluator_generator with:
Subgrid water content.
    
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "Teuchos_ParameterList.hpp"
#include "dbc.hh"
#include "overland_subgrid_water_content_model.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

// Constructor from ParameterList
OverlandSubgridWaterContentModel::OverlandSubgridWaterContentModel(Teuchos::ParameterList& plist)
{
  InitializeFromPlist_(plist);
}


// Initialize parameters
void
OverlandSubgridWaterContentModel::InitializeFromPlist_(Teuchos::ParameterList& plist)
{
  delta_max_ = plist.get<double>("max ponded depth");
  delta_ex_ = plist.get<double>("excluded volume");
}


// main method
double
OverlandSubgridWaterContentModel::OverlandSubgridWaterContent(double pd) const
{
  return pow(pd, 2)*(-3*delta_ex + 2*delta_max)/pow(delta_max, 2) + pow(pd, 3)*(2*delta_ex - delta_max)/pow(delta_max, 3);
}

double
OverlandSubgridWaterContentModel::DOverlandSubgridWaterContentDHeight(double pd) const
{
  return 2*pd*(-3*delta_ex + 2*delta_max)/pow(delta_max, 2) + 3*pow(pd, 2)*(2*delta_ex - delta_max)/pow(delta_max, 3);
}

} //namespace
} //namespace
} //namespace
  