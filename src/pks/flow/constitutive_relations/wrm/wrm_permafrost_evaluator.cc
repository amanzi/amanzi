/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  This WRM evaluator evaluates saturation of gas, liquid, and ice from the
  constituents, A and B in the permafrost notes.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "wrm_permafrost_evaluator.hh"


namespace Amanzi {
namespace Flow {
namespace FlowRelations {

#define DEBUG_FLAG 0

WRMPermafrostEvaluator::WRMPermafrostEvaluator(Teuchos::ParameterList& wrm_plist) :
    wrm_plist_(wrm_plist) {

  // my keys are for saturation
  s_l_key_ = wrm_plist_.get<string>("liquid saturation key", "saturation_liquid");
  my_keys_.push_back(s_l_key_);
  my_keys_.push_back(wrm_plist_.get<string>("ice saturation key", "saturation_ice"));
  my_keys_.push_back(wrm_plist_.get<string>("gas saturation key", "saturation_gas"));

  // 1/A is the ice-liquid
  one_on_A_key_ = wrm_plist_.get<string>("1/A key", "wrm_permafrost_one_on_A");
  dependencies_.insert(one_on_A_key_);

  // 1/B is the gas-liquid
  one_on_B_key_ = wrm_plist_.get<string>("1/B key", "wrm_permafrost_one_on_B");
  dependencies_.insert(one_on_B_key_);
}

WRMPermafrostEvaluator::WRMPermafrostEvaluator(const WRMPermafrostEvaluator& other) :
    SecondaryVariablesFieldEvaluator(other),
    wrm_plist_(other.wrm_plist_),
    one_on_A_key_(other.one_on_A_key_),
    one_on_B_key_(other.one_on_B_key_),
    s_l_key_(other.s_l_key_) {}

Teuchos::RCP<FieldEvaluator>
WRMPermafrostEvaluator::Clone() const {
  return Teuchos::rcp(new WRMPermafrostEvaluator(*this));
}


void WRMPermafrostEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const std::vector<Teuchos::Ptr<CompositeVector> >& results) {
  Teuchos::Ptr<CompositeVector> sat = results[0];
  Teuchos::Ptr<CompositeVector> sat_i = results[1];
  Teuchos::Ptr<CompositeVector> sat_g = results[2];

  Teuchos::RCP<const CompositeVector> one_on_A = S->GetFieldData(one_on_A_key_);
  Teuchos::RCP<const CompositeVector> one_on_B = S->GetFieldData(one_on_B_key_);

  // Loop over names in the target and then owned entities in that name,
  // evaluating the evaluator to calculate saturations.
  for (CompositeVector::name_iterator comp=sat->begin();
       comp!=sat->end(); ++comp) {
    for (int id=0; id!=sat->size(*comp); ++id) {
      double s_l = 1.0 / (1.0/(*one_on_A)(*comp, id) + 1.0/(*one_on_B)(*comp, id) - 1.0);
      (*sat)(*comp, id) = s_l;
      (*sat_i)(*comp, id) = s_l * ( 1.0/(*one_on_A)(*comp, id) - 1.0);
      (*sat_g)(*comp, id) = s_l * ( 1.0/(*one_on_B)(*comp, id) - 1.0);

#if DEBUG_FLAG
      if (id==0) {
        std::cout << " sat_l( 0) = " << s_l << ",    sat_i( 0) = " << (*sat_i)(*comp, id) << ",    sat_g( 0) = " << (*sat_g)(*comp, id) << std::endl;
      }
      if (id==99) {
        std::cout << " sat_l(99) = " << s_l << ",    sat_i(99) = " << (*sat_i)(*comp, id) << ",    sat_g(99) = " << (*sat_g)(*comp, id) << std::endl;
      }
#endif
    }
  }
}


void WRMPermafrostEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const std::vector<Teuchos::Ptr<CompositeVector> > & results) {
  Teuchos::Ptr<CompositeVector> dsat = results[0];
  Teuchos::Ptr<CompositeVector> dsat_i = results[1];
  Teuchos::Ptr<CompositeVector> dsat_g = results[2];

  Teuchos::RCP<const CompositeVector> one_on_A = S->GetFieldData(one_on_A_key_);
  Teuchos::RCP<const CompositeVector> one_on_B = S->GetFieldData(one_on_B_key_);
  Teuchos::RCP<const CompositeVector> sat = S->GetFieldData(s_l_key_);

  if (wrt_key == one_on_A_key_) {
    for (CompositeVector::name_iterator comp=sat->begin();
         comp!=sat->end(); ++comp) {
      for (int id=0; id!=sat->size(*comp); ++id) {
        double Ainv = (*one_on_A)(*comp, id);
        double A = 1.0 / Ainv;
        double dA = - A * A;
        double B = 1.0 / (*one_on_B)(*comp, id);
        double sl = (*sat)(*comp, id);

        (*dsat)(*comp, id) = (- sl * sl) * dA;
        (*dsat_i)(*comp, id) = sl*dA + (A - 1.0)*(*dsat)(*comp, id);
        (*dsat_g)(*comp, id) = (B - 1.0)*(*dsat)(*comp, id);
      }
    }
  } else if (wrt_key == one_on_B_key_) {
    for (CompositeVector::name_iterator comp=sat->begin();
         comp!=sat->end(); ++comp) {
      for (int id=0; id!=sat->size(*comp); ++id) {
        double Binv = (*one_on_B)(*comp, id);
        double B = 1.0 / Binv;
        double dB = - B * B;
        double A = 1.0 / (*one_on_A)(*comp, id);
        double sl = (*sat)(*comp, id);

        (*dsat)(*comp, id) = (- sl * sl) * dB;
        (*dsat_i)(*comp, id) = (A - 1.0)*(*dsat)(*comp, id);
        (*dsat_g)(*comp, id) = sl*dB + (B - 1.0)*(*dsat)(*comp, id);
      }
    }
  } else {
    ASSERT(0);
  }
}


} // namespace
} // namespace
} // namespace



