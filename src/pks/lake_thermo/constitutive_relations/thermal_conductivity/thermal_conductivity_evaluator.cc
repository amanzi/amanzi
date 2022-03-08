/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Interface for a thermal conductivity of lake model.

  License: BSD
  Authors: Svetlana Tokareva (tokareva@lanl.gov)
 */

#include "thermal_conductivity_evaluator.hh"

namespace Amanzi {
namespace LakeThermo {

ThermalConductivityEvaluator::ThermalConductivityEvaluator(
    Teuchos::ParameterList& plist) :
        SecondaryVariableFieldEvaluator(plist) {
  if (my_key_ == std::string("")) {
    my_key_ = plist_.get<std::string>("thermal conductivity key",
        "surface-thermal_conductivity");
  }

  Key domain = Keys::getDomain(my_key_);

  // Set up my dependencies.
  std::string domain_name = Keys::getDomain(my_key_);

  // -- temperature
  temperature_key_ = Keys::readKey(plist_, domain_name, "temperature", "temperature");
  dependencies_.insert(temperature_key_);

  AMANZI_ASSERT(plist_.isSublist("thermal conductivity parameters"));
  Teuchos::ParameterList sublist = plist_.sublist("thermal conductivity parameters");
  //  K_liq_ = sublist.get<double>("thermal conductivity of water [W/(m-K)]", 0.58);
  //  K_ice_ = sublist.get<double>("thermal conductivity of ice [W/(m-K)]", 2.18);
  //  min_K_ = sublist.get<double>("minimum thermal conductivity", 1.e-14);

  // later: read these parameters from xml
  K_max_    = 150; // [W/(m * K)]
  K_0_      = 50.;
  V_wind_   = 10.; // [m/s]
  V_wind_0_ = 20.;
}


ThermalConductivityEvaluator::ThermalConductivityEvaluator(
    const ThermalConductivityEvaluator& other) :
        SecondaryVariableFieldEvaluator(other),
        K_max_(other.K_max_),
        K_0_(other.K_0_),
        V_wind_(other.V_wind_),
        V_wind_0_(other.V_wind_0_),
        temperature_key_(other.temperature_key_){}
//    uf_key_(other.uf_key_),
//    height_key_(other.height_key_),
//    K_liq_(other.K_liq_),
//    K_ice_(other.K_ice_),
//    min_K_(other.min_K_) {}


Teuchos::RCP<FieldEvaluator>
ThermalConductivityEvaluator::Clone() const {
  return Teuchos::rcp(new ThermalConductivityEvaluator(*this));
}

void ThermalConductivityEvaluator::EvaluateField_(
    const Teuchos::Ptr<State>& S,
    const Teuchos::Ptr<CompositeVector>& result) {

  ice_cover_ = false; // first always assume that there is no ice

  double z_ice = 1.0;
  double z_w   = 1.0;

  double lambda_ice = 2.2;
  double lambda_w   = 1.5; //3.*0.561; //1.5
  //  lambda_ice = lambda_w;

  // get temperature
  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(temperature_key_);

  // get mesh
  Teuchos::RCP<const AmanziMesh::Mesh> mesh = result->Mesh();

  for (CompositeVector::name_iterator comp=result->begin();
      comp!=result->end(); ++comp) {
    // much more efficient to pull out vectors first
    //      const Epetra_MultiVector& eta_v = *eta->ViewComponent(*comp,false);
    //      const Epetra_MultiVector& height_v = *height->ViewComponent(*comp,false);
    const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp,false);
    Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

    int ncomp = result->size(*comp, false);

    //      std::cout << "ncomp in lambda = " << ncomp << std::endl;

    int i_ice_max = 0;
    int i_ice_min = ncomp-1;


    for (int i=0; i!=ncomp; ++i) {
      if (temp_v[0][i] < 273.15) { // check if there is ice cover
        //          std::cout << "temp_v[0][" << i << "] = " << temp_v[0][i] << std::endl;
        ice_cover_ = true;
        i_ice_max = i;
      }
    } // i

    for (int i=ncomp-2; i!=0; --i) {
      if (temp_v[0][i] < 273.15) { // check if there is ice cover
        //          std::cout << "temp_v[0][" << i << "] = " << temp_v[0][i] << std::endl;
        ice_cover_ = true;
        i_ice_min = i;
      }
    } // i

    //      i_ice_max = i_ice_min;
    //      std::cout << "i_ice_max/min = " << i_ice_max << " " << i_ice_min << std::endl;

    if (ice_cover_) {
      const AmanziGeometry::Point& zci = mesh->cell_centroid(i_ice_max+1);
      z_ice = zci[2];

      const AmanziGeometry::Point& zcw = mesh->cell_centroid(i_ice_max-1);
      z_w = zcw[2];
    }

    //      std::cout << "z_ice = " << z_ice << ", z_w = " << z_w << std::endl;

    for (int i=0; i!=ncomp; ++i) {
      if (temp_v[0][i] < 273.15) { // this cell is in ice layer
        result_v[0][i] = lambda_ice;
      }
      else {
        if (ice_cover_) {
          result_v[0][i] = lambda_w;
        } else {
          result_v[0][i] = lambda_w; //10.*K_0_ + V_wind_/V_wind_0_*(K_max_ - K_0_);
        }
      }
    } // i

    int d_thawed = ncomp-1-i_ice_max;
    int d_ice = i_ice_max-i_ice_min+1;

    std::vector<double> lambda_new(ncomp);

    for (int i=0; i!=ncomp; ++i) {
      lambda_new[i] = result_v[0][i];
    }

    if (ice_cover_ && d_thawed > 0) {
      // if thawing occured at the top, swap cells
      for (int i=0; i < std::min(d_thawed,d_ice); ++i) {
        lambda_new[ncomp-1-i] = result_v[0][i_ice_max-i];
        lambda_new[ncomp-d_ice-1-i] = result_v[0][ncomp-1-i];
      } // i

      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = lambda_new[i];
      }

    }


    //  	// if melting occured at the top, swap cells
    //  	for (int i=ncomp-1; i!=1; --i) {
    //  	  if (ice_cover_ && i > i_ice_max && temp_v[0][i] >= 273.15 ) {
    //  		result_v[0][i] = result_v[0][i-1];
    //  		result_v[0][i_ice_min] = result_v[0][i_ice_min-1];
    //  	  }
    //  	} // i

    /*
      // continuous conductivity between ice and water, no ice movement or changes depending on temperature
      for (int i=0; i!=ncomp; ++i) {

        const AmanziGeometry::Point& zc = mesh->cell_centroid(i);

        if (zc[2] >= z_ice) {
          result_v[0][i] = lambda_ice;
        } else {
            if (zc[2] < z_w) {
              result_v[0][i] = lambda_w;
            } else {
                result_v[0][i] = lambda_w + (lambda_ice - lambda_w)/(z_ice - z_w)*(zc[2] - z_w);
            }
        }
        std::cout << "z = " << zc[2] << ", lambda = " << result_v[0][i] << std::endl;
      } // i
     */

    // simple interpolation between neighboring cells
    std::vector<double> lambda(ncomp);

    lambda[0] = result_v[0][0];

    for (int i=1; i!=ncomp; ++i) {

      const AmanziGeometry::Point& zc = mesh->cell_centroid(i);

      lambda[i] = 0.5*(result_v[0][i-1]+result_v[0][i]);

    } // i

    for (int i=0; i!=ncomp; ++i)  result_v[0][i] = lambda[i];

    //      // if melting occured at the top, swap cells
    //      while (ice_cover_ && temp_v[0][0] >= 273.15 ) {
    //        if (ice_cover_ && temp_v[0][0] >= 273.15) { // check temperature at top cell
    //            double tmp = temp_v[0][0];
    //            std::cout << "temp_v[0][0] = " << temp_v[0][0] << std::endl;
    //            for (int i=0; i!=ncomp; ++i) std::cout << "i = " << i << " result_v[0][i] = " << result_v[0][i] << std::endl;
    //            for (int i=0; i!=ncomp-1; ++i) { // check temperature in other cells
    //                if (temp_v[0][i] <= 273.15) {
    //                    result_v[0][i] = result_v[0][i+1];
    //                    temp_v[0][i]   = temp_v[0][i+1];
    //                }
    //                if (temp_v[0][i] > 273.15) {
    //                    result_v[0][i] = result_v[0][i+1];
    //                    temp_v[0][i]   = temp_v[0][i+1];
    //                }
    //            } //i
    //            for (int i=0; i!=ncomp; ++i) std::cout << "i = " << i << " result_v[0][i] = " << result_v[0][i] << std::endl;
    //        }
    //      }

    //      for (int i=0; i!=ncomp; ++i) {
    //        if (ice_cover_ && i < i_ice_max && temp_v[0][i] >= 273.15 ) {
    //          result_v[0][i] = lambda_ice;
    //          result_v[0][i+i_ice_max] = lambda_w;
    //        }
    //      } // i

    //      // no ice
    //      for (int i=0; i!=ncomp; ++i) {
    ////          result_v[0][i] = 1.5;
    //          result_v[0][i] = lambda_w;
    //      }

  }

}


void ThermalConductivityEvaluator::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S, Key wrt_key,
    const Teuchos::Ptr<CompositeVector>& result) {
  std::cout<<"THERMAL CONDUCITIVITY: Derivative not implemented yet!"<<wrt_key<<"\n";
  AMANZI_ASSERT(0); // not implemented, not yet needed
  result->Scale(1.e-6); // convert to MJ
}

} //namespace
} //namespace
