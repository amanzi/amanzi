/*
  Drainage rate.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "drainage_evaluator.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

DrainageEvaluator::DrainageEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {

  std::string domain = plist_.get<std::string>("layer name");
  
  // Set up my dependencies.
  // -- the extent of material, LAI for 
  ai_key_ = plist_.get<std::string>("area index key", Keys::getKey(domain, "area_index"));
  dependencies_.insert(ai_key_);

  // -- water content of the layer drained
  wc_key_ = plist_.get<std::string>("water content key", Keys::getKey(domain, "water_content"));
  dependencies_.insert(wc_key_);

  // -- uptake, i.e. rewetting of layer from surface water
  is_uptake_ = plist_.get<bool>("wet layer from surface water", false);

  pd_key_ = "";

  if (is_uptake_) {
    pd_key_ = plist_.get<std::string>("ponded depth key",
				      "surface-ponded_depth");
    dependencies_.insert(pd_key_);
  }

  // -- source of water into the layer
  source_key_ = plist_.get<std::string>("source key", Keys::getKey(domain, "interception"));
  dependencies_.insert(source_key_);


  // 
  dependencies_.insert("surface-cell_volume");

  // parameters for the drainage model
  tau_ = plist_.get<double>("drainage timescale [s]");
  wc_sat_ = plist_.get<double>("saturated moisture content [m^3 H20 / m^2 biomass area]"); // this is somehow related to a LAI?
  n_liq_ = plist_.get<double>("density of liquid water [mol/m^3]", 1000. / 0.0180153);
};


DrainageEvaluator::DrainageEvaluator(const DrainageEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    wc_key_(other.wc_key_),
    ai_key_(other.ai_key_),
    pd_key_(other.pd_key_),
    source_key_(other.source_key_),
    tau_(other.tau_),
    wc_sat_(other.wc_sat_),
    n_liq_(other.n_liq_),
    is_uptake_(other.is_uptake_)
{}


Teuchos::RCP<FieldEvaluator> DrainageEvaluator::Clone() const {
  return Teuchos::rcp(new DrainageEvaluator(*this));
}


void DrainageEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
                         const Teuchos::Ptr<CompositeVector>& result) {
  // Pull dependencies out of state.
  const Epetra_MultiVector& wc =
      *S->GetFieldData(wc_key_)->ViewComponent("cell",false);
  const Epetra_MultiVector& ai =
      *S->GetFieldData(ai_key_)->ViewComponent("cell",false);
  const Epetra_MultiVector& source =
      *S->GetFieldData(source_key_)->ViewComponent("cell",false);  
  const Epetra_MultiVector& cv =
      *S->GetFieldData("surface-cell_volume")->ViewComponent("cell",false);
  Epetra_MultiVector& res_c = *result->ViewComponent("cell",false);

  Teuchos::RCP<const Epetra_MultiVector> pd;
  if (is_uptake_)
    pd = S->GetFieldData(pd_key_)->ViewComponent("cell",false);

  
  // evaluate the model
  for (int c=0; c!=res_c.MyLength(); ++c) {
    res_c[0][c] = 0.0;

    if (ai[0][c] == 0.0) {
      res_c[0][c] = source[0][c] * cv[0][c];
    } else {
      double wc_sat = n_liq_ * ai[0][c] * cv[0][c] * wc_sat_;
      if (wc[0][c] > wc_sat) {

        //  is oversaturated and draining
        res_c[0][c] = (wc[0][c] - wc_sat) / tau_;
      } else if (is_uptake_) {
        //  is undersaturated and there is surface water to be absorbed
        double wetting = std::min((*pd)[0][c] / ai[0][c], 1.0);
        res_c[0][c] = wetting * (wc[0][c] - wc_sat) / tau_;

      }
    }

    // turn back into mols / (m^2 s)
    res_c[0][c] /= cv[0][c];
  }
}


void DrainageEvaluator::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S, Key wrt_key,
    const Teuchos::Ptr<CompositeVector>& result) {

  // Pull dependencies out of state.
  const Epetra_MultiVector& wc =
      *S->GetFieldData(wc_key_)->ViewComponent("cell",false);
  const Epetra_MultiVector& ai =
      *S->GetFieldData(ai_key_)->ViewComponent("cell",false);
  const Epetra_MultiVector& cv =
      *S->GetFieldData("surface-cell_volume")->ViewComponent("cell",false);
  const Epetra_MultiVector& source =
      *S->GetFieldData(source_key_)->ViewComponent("cell",false);  
  Epetra_MultiVector& res_c = *result->ViewComponent("cell",false);

  Teuchos::RCP<const Epetra_MultiVector> pd;
  if (is_uptake_)
    pd = S->GetFieldData(pd_key_)->ViewComponent("cell",false);
  
  if (wrt_key == wc_key_) {
    for (int c=0; c!=res_c.MyLength(); ++c) {
      if (ai[0][c] == 0.) {

        res_c[0][c] = 0.;
      } else {
        double wc_sat = n_liq_ * ai[0][c] * cv[0][c] * wc_sat_;
        if (wc[0][c] > wc_sat) {
          //  is oversaturated and draining
          res_c[0][c] = 1.0 / tau_;
        } else if (is_uptake_) {
          //  is undersaturated and there is surface water to be absorbed
          double wetting = std::min((*pd)[0][c] / ai[0][c], 1.0);
          res_c[0][c] = wetting / tau_;
        }

      }
      // turn back into mols / (m^2 s)
      res_c[0][c] /= cv[0][c];
    }

  } else if (wrt_key == ai_key_) {
    for (int c=0; c!=res_c.MyLength(); ++c) {
      if (ai[0][c] == 0.) {

        res_c[0][c] = 0.;
      } else {
        double wc_sat = n_liq_ * ai[0][c] * cv[0][c] * wc_sat_;
        if (wc[0][c] > wc_sat) {
          //  is oversaturated and draining
          res_c[0][c] = -n_liq_ * cv[0][c] * wc_sat_ / tau_;
        } else if (is_uptake_) {
          double wetting = std::min((*pd)[0][c] / ai[0][c], 1.0);
          res_c[0][c] = -wetting * n_liq_ * cv[0][c] * wc_sat_ / tau_;
        }

      }
      // turn back into mols / (m^2 s)
      res_c[0][c] /= cv[0][c];
    }


  } else if (wrt_key == source_key_) {
    for (int c=0; c!=res_c.MyLength(); ++c) {
      if (ai[0][c] == 0.) {
        res_c[0][c] = 1.;
      } else {
        res_c[0][c] = 0.;

      }
      // turn back into mols / (m^2 s)
      res_c[0][c] /= cv[0][c];
    }


  } else if (is_uptake_ && wrt_key == pd_key_) {
    for (int c=0; c!=res_c.MyLength(); ++c) {
      if (ai[0][c] == 0.) {
        res_c[0][c] = 0.;
      } else {
        double wc_sat = n_liq_ * ai[0][c] * cv[0][c] * wc_sat_;
        if (wc[0][c] > wc_sat) {
          //  is oversaturated and draining
          res_c[0][c] = 0.;
        } else if (is_uptake_) {
          //  is undersaturated and there is surface water to be absorbed
          double wetting = (*pd)[0][c] > ai[0][c] ? 0. : 1./ai[0][c];
          res_c[0][c] = wetting * (wc[0][c] - wc_sat) / tau_;
        }

      }
      // turn back into mols / (m^2 s)
      res_c[0][c] /= cv[0][c];
    }

  } else {
    AMANZI_ASSERT(0);
  }

}

} // namespace
} // namespace
} // namespace
