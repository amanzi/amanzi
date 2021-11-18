/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Interface for a heat flux at the surface of lake model.

  License: BSD
  Authors: Svetlana Tokareva (tokareva@lanl.gov)
*/

#include "heat_flux_bc_evaluator.hh"

namespace Amanzi {
namespace LakeThermo {

HeatFluxBCEvaluator::HeatFluxBCEvaluator(
      Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {
  if (my_key_ == std::string("")) {
    my_key_ = plist_.get<std::string>("heat flux bc key",
            "surface-heat_flux_bc");
  }

  Key domain = Keys::getDomain(my_key_);

  // Set up my dependencies.
  std::string domain_name = Keys::getDomain(my_key_);

  // -- temperature
  temperature_key_ = Keys::readKey(plist_, domain_name, "temperature", "temperature");
  dependencies_.insert(temperature_key_);

//  AMANZI_ASSERT(plist_.isSublist("heat flux bc parameters"));
//  Teuchos::ParameterList sublist = plist_.sublist("heat flux bc parameters");

  // later: read these parameters from xml
  SS = 0.;      // solar radiation (read from met data)
  alpha = 0.07; // water albedo
  E_a = 0.;     // atmospheric downward radiation (read from met data)
  E_s = 0.;     // surface radiation
  H = 0.;       // "sensible" heat
  LE = 0.;      // latent heat

}


HeatFluxBCEvaluator::HeatFluxBCEvaluator(
      const HeatFluxBCEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    SS(other.SS),
    alpha(other.alpha),
    E_a(other.E_a),
    E_s(other.E_s),
    H(other.H),
    LE(other.LE),
    temperature_key_(other.temperature_key_){}


Teuchos::RCP<FieldEvaluator>
HeatFluxBCEvaluator::Clone() const {
  return Teuchos::rcp(new HeatFluxBCEvaluator(*this));
}

void HeatFluxBCEvaluator::EvaluateField_(
      const Teuchos::Ptr<State>& S,
      const Teuchos::Ptr<CompositeVector>& result) {

  ice_cover_ = false; // first always assume that there is no ice

  // get temperature
  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(temperature_key_);

  // read parameters from the met data
  Teuchos::ParameterList& param_list = plist_.sublist("parameters");
  Amanzi::FunctionFactory fac;
  Teuchos::RCP<Amanzi::Function> SS_func_ = Teuchos::rcp(fac.Create(param_list.sublist("solar radiation")));
  Teuchos::RCP<Amanzi::Function> E_a_func_ = Teuchos::rcp(fac.Create(param_list.sublist("atmospheric downward radiation")));

  std::vector<double> args(1);
  args[0] = S->time();
  SS = (*SS_func_)(args);
  E_a = (*E_a_func_)(args);

  for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {

      const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp,false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);

      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = SS*(1.-alpha) + E_a - E_s - H - LE;
      } // i

    }

}

void HeatFluxBCEvaluator::EvaluateFieldPartialDerivative_(
      const Teuchos::Ptr<State>& S, Key wrt_key,
      const Teuchos::Ptr<CompositeVector>& result) {
  std::cout<<"HEAT FLUX BC: Derivative not implemented yet!"<<wrt_key<<"\n";
  AMANZI_ASSERT(0); // not implemented, not yet needed
  result->Scale(1.e-6); // convert to MJ
}

} //namespace
} //namespace
