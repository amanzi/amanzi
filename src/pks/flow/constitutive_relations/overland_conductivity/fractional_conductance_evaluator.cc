/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  The elevation evaluator gets the surface elevation, slope, and updates pres + elev.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "fractional_conductance_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

FractionalConductanceEvaluator::FractionalConductanceEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {

  Key domain = getDomain(my_key_);

  if (my_key_.empty())
    my_key_ = plist_.get<std::string>("fractional conductance key", getKey(domain,"fractional_conductance"));


  pd_key_ = plist_.get<std::string>("height key", getKey(domain,"ponded_depth"));
  dependencies_.insert(pd_key_);
  
  vpd_key_ = plist_.get<std::string>("volumetric height key", getKey(domain,"volumetric_ponded_depth"));
  dependencies_.insert(vpd_key_); 

  pdd_key_ = plist_.get<std::string>("ponded depression depth key", getKey(domain,"ponded_depression_depth"));
  dependencies_.insert(pdd_key_);
  
  delta_max_key_ = plist_.get<std::string>("maximum ponded depth key", getKey(domain,"maximum_ponded_depth"));
  dependencies_.insert(delta_max_key_);
  delta_ex_key_ = plist_.get<std::string>("excluded volume key", getKey(domain,"excluded_volume"));
  dependencies_.insert(delta_ex_key_);
  depr_depth_key_ = plist_.get<std::string>("depression depth key", getKey(domain,"depression_depth"));
  dependencies_.insert(depr_depth_key_);
  //  depr_depth_ = plist_.get<double>("depression depth");//, 0.043);
  //delta_max_ = plist_.get<double>("maximum ponded depth");//, 0.483);
  //delta_ex_ = plist_.get<double>("excluded volume");//,0.23);

}


FractionalConductanceEvaluator::FractionalConductanceEvaluator(const FractionalConductanceEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    pdd_key_(other.pdd_key_),
    pd_key_(other.pd_key_),
    vpd_key_(other.vpd_key_),
    delta_ex_key_(other.delta_ex_key_),
    delta_max_key_(other.delta_max_key_),
    depr_depth_key_(other.delta_max_key_)
{};

Teuchos::RCP<FieldEvaluator>
FractionalConductanceEvaluator::Clone() const {
  return Teuchos::rcp(new FractionalConductanceEvaluator(*this));
}


void FractionalConductanceEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {

  Epetra_MultiVector& res = *result->ViewComponent("cell",false);
  //  const Epetra_MultiVector& pd_depth = *S->GetFieldData(pdd_key_)->ViewComponent("cell",false);
  const Epetra_MultiVector& depth = *S->GetFieldData(pd_key_)->ViewComponent("cell",false);
  const Epetra_MultiVector& vpd = *S->GetFieldData(vpd_key_)->ViewComponent("cell",false);
  
  const Epetra_MultiVector& delta_max_v = *S->GetFieldData(delta_max_key_)->ViewComponent("cell", false);
  const Epetra_MultiVector& delta_ex_v = *S->GetFieldData(delta_ex_key_)->ViewComponent("cell", false);
  const Epetra_MultiVector& depr_depth_v = *S->GetFieldData(depr_depth_key_)->ViewComponent("cell", false);

 
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  int ncells = res.MyLength();
  for (int c=0; c!=ncells; ++c) {
    double depr_depth = depr_depth_v[0][c];
    double delta_max = delta_max_v[0][c];
    double delta_ex = delta_ex_v[0][c];
    const double fixed_depth = std::pow(depr_depth,2)*(2*delta_max - 3*delta_ex)/std::pow(delta_max,2) + std::pow(depr_depth,3)*(2*delta_ex - delta_max)/std::pow(delta_max,3);

    if (depth[0][c] <= depr_depth)
      res[0][c] = 0;
    else{
      //      double pd = std::pow(depth[0][c],2)*(2*delta_max_ - 3*delta_ex_)/std::pow(delta_max_,2) + std::pow(depth[0][c],3)*(2*delta_ex_ - delta_max_)/std::pow(delta_max_,3);
      // res[0][c] = (pd - depr) / depth[0][c];
      res[0][c] = (vpd[0][c] - fixed_depth) / depth[0][c];
    }
  }
}

// This is hopefully never called?
void FractionalConductanceEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {
  std::cout<<"Fractional Conductance Evaluator Partial Derivative -- NEVER CALLED??: \n"; abort();
  result->PutScalar(1.0/double(maxpd_depth_ - depr_depth_));
}

} //namespace
} //namespace
} //namespace
