/*
  The overland subgrid water content model is an algebraic model with dependencies.

  Generated via evaluator_generator with:
Subgrid water content.
    
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOW_OVERLAND_SUBGRID_WATER_CONTENT_MODEL_HH_
#define AMANZI_FLOW_OVERLAND_SUBGRID_WATER_CONTENT_MODEL_HH_

namespace Amanzi {
namespace Flow {
namespace Relations {

class OverlandSubgridWaterContentModel {

 public:
  explicit
  OverlandSubgridWaterContentModel(Teuchos::ParameterList& plist);

  double OverlandSubgridWaterContent(double pd) const;

  double DOverlandSubgridWaterContentDHeight(double pd) const;
  
 protected:
  void InitializeFromPlist_(Teuchos::ParameterList& plist);

 protected:

  double delta_max_;
  double delta_ex_;

};

} //namespace
} //namespace
} //namespace

#endif