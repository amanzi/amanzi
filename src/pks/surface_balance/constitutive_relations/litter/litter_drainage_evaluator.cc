/*
  Litter drainage rate.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "litter_drainage_evaluator.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

LitterDrainageEvaluator::LitterDrainageEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {

  // Set up my dependencies.
  litter_thickness_key_ = plist_.get<std::string>("litter thickness key",
          "litter_thickness");
  dependencies_.insert(litter_thickness_key_);

  litter_wc_key_ = plist_.get<std::string>("litter water content key",
          "litter_water_content");
  dependencies_.insert(litter_wc_key_);

  pd_key_ = plist_.get<std::string>("ponded depth key",
          "ponded_depth");
  dependencies_.insert(pd_key_);  

  dependencies_.insert("surface_cell_volume");

  tau_ = plist_.get<double>("litter drainage timescale [s]");
  wc_sat_ = plist_.get<double>("litter moisture (saturated) [-]"); // this is somehow related to a LAI?
  n_liq_ = plist_.get<double>("density of liquid water [mol/m^3]", 1000. / 0.0180153);
  rewetting_ = plist_.get<bool>("wet litter from surface water", true);

  source_key_ = plist_.get<std::string>("litter source key");
  source_coef_ = plist_.get<double>("litter source coefficient", 1.0);
};


LitterDrainageEvaluator::LitterDrainageEvaluator(const LitterDrainageEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    litter_wc_key_(other.litter_wc_key_),
    litter_thickness_key_(other.litter_thickness_key_),
    pd_key_(other.pd_key_),
    tau_(other.tau_),
    wc_sat_(other.wc_sat_),
    n_liq_(other.n_liq_),
    rewetting_(other.rewetting_),
    source_key_(other.source_key_),
    source_coef_(other.source_coef_)
{}


Teuchos::RCP<FieldEvaluator> LitterDrainageEvaluator::Clone() const {
  return Teuchos::rcp(new LitterDrainageEvaluator(*this));
}


void LitterDrainageEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
                         const Teuchos::Ptr<CompositeVector>& result) {
  // Pull dependencies out of state.
  const Epetra_MultiVector& wc =
      *S->GetFieldData(litter_wc_key_)->ViewComponent("cell",false);
  const Epetra_MultiVector& ld =
      *S->GetFieldData(litter_thickness_key_)->ViewComponent("cell",false);
  const Epetra_MultiVector& pd =
      *S->GetFieldData(pd_key_)->ViewComponent("cell",false);  
  const Epetra_MultiVector& source =
      *S->GetFieldData(source_key_)->ViewComponent("cell",false);  
  const Epetra_MultiVector& cv =
      *S->GetFieldData("surface_cell_volume")->ViewComponent("cell",false);
  Epetra_MultiVector& res_c = *result->ViewComponent("cell",false);

  // evaluate the model
  for (int c=0; c!=res_c.MyLength(); ++c) {
    res_c[0][c] = 0.0;

    if (ld[0][c] == 0.0) {
      res_c[0][c] = source[0][c] * source_coef_;
    } else {
      double wc_sat = n_liq_ * ld[0][c] * cv[0][c] * wc_sat_;
      if (wc[0][c] > wc_sat) {
	// litter is oversaturated and draining
	res_c[0][c] = (wc[0][c] - wc_sat) / tau_;
      } else if (rewetting_) {
	// litter is undersaturated and there is surface water to be absorbed
	double litter_wetting = ld[0][c] > 0 ? std::min(pd[0][c] / ld[0][c], 1.0) : 1.0;
	res_c[0][c] = litter_wetting * (wc[0][c] - wc_sat) / tau_;
      }
    }
  }
}


void LitterDrainageEvaluator::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S, Key wrt_key,
    const Teuchos::Ptr<CompositeVector>& result) {

  // Pull dependencies out of state.
  const Epetra_MultiVector& wc =
      *S->GetFieldData(litter_wc_key_)->ViewComponent("cell",false);
  const Epetra_MultiVector& ld =
      *S->GetFieldData(litter_thickness_key_)->ViewComponent("cell",false);
  const Epetra_MultiVector& pd =
      *S->GetFieldData(pd_key_)->ViewComponent("cell",false);  
  const Epetra_MultiVector& cv =
      *S->GetFieldData("surface_cell_volume")->ViewComponent("cell",false);
  const Epetra_MultiVector& source =
      *S->GetFieldData(source_key_)->ViewComponent("cell",false);  
  Epetra_MultiVector& res_c = *result->ViewComponent("cell",false);

  if (wrt_key == litter_wc_key_) {
    for (int c=0; c!=res_c.MyLength(); ++c) {
      if (ld[0][c] == 0.) {
	res_c[0][c] = 0.;
      } else {
	double wc_sat = n_liq_ * ld[0][c] * cv[0][c] * wc_sat_;
	if (wc[0][c] > wc_sat) {
	  // litter is oversaturated and draining
	  res_c[0][c] = 1.0 / tau_;
	} else if (rewetting_) {
	  // litter is undersaturated and there is surface water to be absorbed
	  double litter_wetting = std::min(pd[0][c] / ld[0][c], 1.0);
	  res_c[0][c] = litter_wetting / tau_;
	}
      }
    }

  } else if (wrt_key == litter_thickness_key_) {
    for (int c=0; c!=res_c.MyLength(); ++c) {
      if (ld[0][c] == 0.) {
	res_c[0][c] = 0.;
      } else {
	double wc_sat = n_liq_ * ld[0][c] * cv[0][c] * wc_sat_;
	if (wc[0][c] > wc_sat) {
	  // litter is oversaturated and draining
	  res_c[0][c] = -n_liq_ * cv[0][c] * wc_sat_ / tau_;
	} else if (rewetting_) {
	  double litter_wetting = std::min(pd[0][c] / ld[0][c], 1.0);
	  res_c[0][c] = -litter_wetting * n_liq_ * cv[0][c] * wc_sat_ / tau_;
	}
      }
    }

  } else if (wrt_key == "surface_cell_volume") {
    for (int c=0; c!=res_c.MyLength(); ++c) {
      if (ld[0][c] == 0.) {
	res_c[0][c] = 0.;
      } else {
	double wc_sat = n_liq_ * ld[0][c] * cv[0][c] * wc_sat_;
	if (wc[0][c] > wc_sat) {
	  // litter is oversaturated and draining
	  res_c[0][c] = -n_liq_ * ld[0][c] * wc_sat_ / tau_;
	} else if (rewetting_) {
	  // litter is undersaturated and there is surface water to be absorbed
	  double litter_wetting = std::min(pd[0][c] / ld[0][c], 1.0);
	  res_c[0][c] = -litter_wetting * n_liq_ * ld[0][c] * wc_sat_ / tau_;
	}
      }
    }

  } else if (wrt_key == pd_key_) {
    for (int c=0; c!=res_c.MyLength(); ++c) {
      if (ld[0][c] == 0.) {
	res_c[0][c] = 0.;
      } else {
	double wc_sat = n_liq_ * ld[0][c] * cv[0][c] * wc_sat_;
	if (wc[0][c] > wc_sat) {
	  // litter is oversaturated and draining
	  res_c[0][c] = 0.;
	} else if (rewetting_) {
	  // litter is undersaturated and there is surface water to be absorbed
	  double dlitter_wetting = pd[0][c] > ld[0][c] ? 0. : 1./ld[0][c];
	  res_c[0][c] = dlitter_wetting * (wc[0][c] - wc_sat) / tau_;
	}
      }
    }

  } else if (wrt_key == source_key_) {
    for (int c=0; c!=res_c.MyLength(); ++c) {
      if (ld[0][c] == 0.) {
	res_c[0][c] = source_coef_;
      } else {
	res_c[0][c] = 0.;
      }
    }
    
  } else {
    ASSERT(0);
  }  
}

} // namespace
} // namespace
} // namespace
