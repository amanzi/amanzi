/*
  The overland subgrid water content evaluator is an algebraic evaluator of a given model.
Subgrid water content.  
  Generated via evaluator_generator.
*/

#include "overland_subgrid_water_content_evaluator.hh"
#include "overland_subgrid_water_content_model.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

// Constructor from ParameterList
OverlandSubgridWaterContentEvaluator::OverlandSubgridWaterContentEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist){ 

  M_ = plist_.get<double>("molar mass", 0.0180153);
  bar_ = plist_.get<bool>("water content bar", false);
  rollover_ = plist_.get<double>("water content rollover", 0.);
  
  
  Key domain;

  if(!my_key_.empty())
    domain = getDomain(my_key_);

  if (my_key_.empty()) {
    domain = plist.get<std::string>("domain name", "surface");
    my_key_ = getKey(domain, "water_content");
    if (bar_) my_key_ += std::string("_bar");
    my_key_ = plist_.get<std::string>("water content key", my_key_);
  }


  delta_max_key_ = plist_.get<std::string>("maximum ponded depth key", getKey(domain,"maximum_ponded_depth"));
  dependencies_.insert(delta_max_key_);
  delta_ex_key_ = plist_.get<std::string>("excluded volume key", getKey(domain,"excluded_volume"));
  dependencies_.insert(delta_ex_key_);
  // my dependencies
  pres_key_ = plist_.get<std::string>("pressure key", getKey(domain,"pressure"));
  dependencies_.insert(pres_key_);
  // vpd_key_ = plist_.get<std::string>("volumetric height key", getKey(domain,"volumetric_ponded_depth"));
  //dependencies_.insert(vpd_key_);
  cv_key_ = plist_.get<std::string>("cell volume key", getKey(domain,"cell_volume"));
  dependencies_.insert(cv_key_);

  
}


// Copy constructor
OverlandSubgridWaterContentEvaluator::OverlandSubgridWaterContentEvaluator(const OverlandSubgridWaterContentEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    pres_key_(other.pres_key_),
    model_(other.model_),
    M_(other.M_),
    delta_max_key_(other.delta_max_key_),
    delta_ex_key_(other.delta_ex_key_),
    bar_(other.bar_),
    rollover_(other.rollover_),
    cv_key_(other.cv_key_)
 {}


// Virtual copy constructor
Teuchos::RCP<FieldEvaluator>
OverlandSubgridWaterContentEvaluator::Clone() const
{
  return Teuchos::rcp(new OverlandSubgridWaterContentEvaluator(*this));
}


void
OverlandSubgridWaterContentEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{
  Epetra_MultiVector& res = *result->ViewComponent("cell",false);
  const Epetra_MultiVector& pres = *S->GetFieldData(pres_key_)->ViewComponent("cell",false);
  //  const Epetra_MultiVector& vpd = *S->GetFieldData(vpd_key_)->ViewComponent("cell",false);

  const Epetra_MultiVector& cv = *S->GetFieldData(cv_key_)->ViewComponent("cell",false);

  const double& p_atm = *S->GetScalarData("atmospheric_pressure");
  const Epetra_Vector& gravity = *S->GetConstantVectorData("gravity");
  double gz = -gravity[2];  // check this
  
  Teuchos::RCP<const CompositeVector> max_pd = S->GetFieldData(delta_max_key_);
  Teuchos::RCP<const CompositeVector> ex_vol = S->GetFieldData(delta_ex_key_);
  // cell values
  const Epetra_MultiVector& max_pd_v = *max_pd->ViewComponent("cell", false);
  const Epetra_MultiVector& ex_vol_v = *ex_vol->ViewComponent("cell", false);

 
  //assert(max_pd_v[0][3]);
  //assert(ex_vol_v[0][3]);

  int ncells = res.MyLength();
  if (bar_) {
    for (int c=0; c!=ncells; ++c) {
      double delta_max = max_pd_v[0][c]*1000./M_;
      double delta_ex = ex_vol_v[0][c]*1000. / M_;
    
      double pd = (pres[0][c] - p_atm)/ (gz * M_);
      if (pd <=delta_max){
	res[0][c] = std::pow(pd,2)*(2*delta_max - 3*delta_ex)/std::pow(delta_max,2) + std::pow(pd,3)*(2*delta_ex - delta_max)/std::pow(delta_max,3);
	res[0][c] *= cv[0][c];
      }
      else if( pd > delta_max)
	res[0][c] = cv[0][c] *(pd - delta_ex);  
    }
  } else if (rollover_ > 0.) {
    std::cout<<"Debuggin....ROLLOVER_ in the subgrid water content evaluator: \n"; abort();
    for (int c=0; c!=ncells; ++c) {
      double dp = pres[0][c] - p_atm;
      double dp_eff = dp < 0. ? 0. :
          dp < rollover_ ?
            dp*dp/(2*rollover_) :
            dp - rollover_/2.;
      res[0][c] = cv[0][c] * dp_eff / (gz * M_);
    }
  } else {

    for (int c=0; c!=ncells; ++c) {
      double pd = pres[0][c] < p_atm ? 0.0 : (pres[0][c] - p_atm)/ (gz * M_);
      double delta_max = max_pd_v[0][c]*1000./M_;
      double delta_ex = ex_vol_v[0][c]*1000. / M_;
      if (0 <= pd && pd <=delta_max){
	res[0][c] = std::pow(pd,2)*(2*delta_max - 3*delta_ex)/std::pow(delta_max,2) + std::pow(pd,3)*(2*delta_ex - delta_max)/std::pow(delta_max,3);
	res[0][c] *= cv[0][c];
      }
      else if( pd > delta_max)
	res[0][c] = cv[0][c] *(pd - delta_ex);  
    }
  }
  
}


void
OverlandSubgridWaterContentEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{

  ASSERT(wrt_key == pres_key_);

  Epetra_MultiVector& res = *result->ViewComponent("cell",false);
  const Epetra_MultiVector& pres = *S->GetFieldData(pres_key_)
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& cv = *S->GetFieldData(cv_key_)
      ->ViewComponent("cell",false);

  const double& p_atm = *S->GetScalarData("atmospheric_pressure");
  const Epetra_Vector& gravity = *S->GetConstantVectorData("gravity");
  double gz = -gravity[2];  // check this
 

  Teuchos::RCP<const CompositeVector> max_pd = S->GetFieldData(delta_max_key_);
  Teuchos::RCP<const CompositeVector> ex_vol = S->GetFieldData(delta_ex_key_);
  // cell values
  const Epetra_MultiVector& max_pd_v = *max_pd->ViewComponent("cell", false);
  const Epetra_MultiVector& ex_vol_v = *ex_vol->ViewComponent("cell", false);
  
  

  int ncells = res.MyLength();

  
  if (bar_) {
    for (int c=0; c!=ncells; ++c) {
      double pd = (pres[0][c] - p_atm)/ (gz * M_);
      double delta_max = max_pd_v[0][c]*1000./M_;
      double delta_ex = ex_vol_v[0][c]*1000. / M_;
      if (pd <=delta_max){
	res[0][c] = 2*pd*(2*delta_max - 3*delta_ex )/ std::pow(delta_max,2) + 3*pd*pd*(2*delta_ex - delta_max)/std::pow(delta_max,3);
	res[0][c] *= cv[0][c] / (gz * M_);
      }
      else if( pd > delta_max)
	res[0][c] = cv[0][c]/ (gz * M_);
    }
  } else if (rollover_ > 0.) {
    std::cout<<"Debugging...ROLLOVER-DERIVATIVE in the subgrid water content. \n"; abort();
    for (int c=0; c!=ncells; ++c) {
      double dp = pres[0][c] - p_atm;
      double ddp_eff = dp < 0. ? 0. :
          dp < rollover_ ? dp/rollover_ : 1.;
      res[0][c] = cv[0][c] * ddp_eff / (gz * M_);
    }
  } else {
    
     for (int c=0; c!=ncells; ++c) {
       double pd = pres[0][c] < p_atm ? 0.0 : (pres[0][c] - p_atm)/ (gz * M_);
       double delta_max =  max_pd_v[0][c]*1000./M_;
       double delta_ex = ex_vol_v[0][c]*1000. / M_;
      if (0 <= pd && pd <=delta_max){
	res[0][c] = 2*pd*(2*delta_max - 3*delta_ex )/ std::pow(delta_max,2) + 3*pd*pd*(2*delta_ex - delta_max)/std::pow(delta_max,3);
	res[0][c] *= cv[0][c] / (gz * M_);
      }
      else if( pd > delta_max)
	res[0][c] = cv[0][c]/ (gz * M_);
       
    }
     
  }


}


} //namespace
} //namespace
} //namespace
