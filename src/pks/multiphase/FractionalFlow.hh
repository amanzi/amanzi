/*
  This is the multiphase component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Quan Bui (mquanbui@math.umd.edu)
*/

#ifndef AMANZI_MULTIPHASE_FRACTIONAL_FLOW_HH_
#define AMANZI_MULTIPHASE_FRACTIONAL_FLOW_HH_

#include <string>
#include <vector>

#include "errors.hh"
#include "Epetra_IntVector.h"
#include "Epetra_Vector.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Mesh.hh"
#include "CompositeVector.hh"
#include "ParallelCommunication.hh"
//#include "Tensor.hh"
#include "VerboseObject.hh"

#include "WaterRetentionModel.hh"
#include "RelativePermeability.hh"
#include "MultiphaseTypeDefs.hh"

namespace Amanzi {
namespace Multiphase {

class FractionalFlow {
 public:
  FractionalFlow(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                 double mu_w, const Teuchos::RCP<const RelativePermeability>& rel_perm_w,
                 double mu_nw, const Teuchos::RCP<const RelativePermeability>& rel_perm_nw) : 
      mesh_(mesh), mu_w_(mu_w), mu_nw_(mu_nw), 
      rel_perm_w_(rel_perm_w), rel_perm_nw_(rel_perm_nw) {};

  ~FractionalFlow() { delete vo_; }

  // main methods
  void Init(Teuchos::ParameterList& plist);
  void Compute(const CompositeVector& Sw);

  double Value(int c, double Sw) const { 
      double lambda_wc, lambda_nwc;
      lambda_wc = WRM_[(*map_c2mb_)[c]]->k_relative(Sw, "wetting")/mu_w_;
      lambda_nwc = WRM_[(*map_c2mb_)[c]]->k_relative(Sw, "non wetting")/mu_nw_;
      return lambda_wc/(lambda_wc + lambda_nwc); 
  } 
  double Derivative(int c, double Sw) const {
    double lambda_wc, lambda_nwc, dlamdaw_ds, dlamdanw_ds;
    lambda_wc = WRM_[(*map_c2mb_)[c]]->k_relative(Sw, "wetting")/mu_w_;
    lambda_nwc = WRM_[(*map_c2mb_)[c]]->k_relative(Sw, "non wetting")/mu_nw_;
    dlamdaw_ds = WRM_[(*map_c2mb_)[c]]->dKdS(Sw, "wetting")/mu_w_;
    dlamdanw_ds = WRM_[(*map_c2mb_)[c]]->dKdS(Sw, "non wetting")/mu_nw_;
    double df_dS = (dlamdaw_ds*(lambda_wc + lambda_nwc) - lambda_wc*(dlamdaw_ds + dlamdanw_ds)) 
                    / (pow(lambda_wc + lambda_nwc, 2.0));
    return df_dS;
  } 
  double Value(int c, double Sw, const std::string name) const { 
    if (name == "fractional_flow") {
      double lambda_wc, lambda_nwc;
      lambda_wc = WRM_[(*map_c2mb_)[c]]->k_relative(Sw, "wetting")/mu_w_;
      lambda_nwc = WRM_[(*map_c2mb_)[c]]->k_relative(Sw, "non wetting")/mu_nw_;
      return lambda_wc/(lambda_wc + lambda_nwc);
    } else if (name == "df_dS") {
      double lambda_wc, lambda_nwc, dlamdaw_ds, dlamdanw_ds;
      lambda_wc = WRM_[(*map_c2mb_)[c]]->k_relative(Sw, "wetting")/mu_w_;
      lambda_nwc = WRM_[(*map_c2mb_)[c]]->k_relative(Sw, "non wetting")/mu_nw_;
      dlamdaw_ds = WRM_[(*map_c2mb_)[c]]->dKdS(Sw, "wetting")/mu_w_;
      dlamdanw_ds = WRM_[(*map_c2mb_)[c]]->dKdS(Sw, "non wetting")/mu_nw_;
      double df_dS = (dlamdaw_ds*(lambda_wc + lambda_nwc) - lambda_wc*(dlamdaw_ds + dlamdanw_ds)) 
                      / (pow(lambda_wc + lambda_nwc, 2.0));
      return df_dS;
    } else {
      Errors::Message msg;
      msg << "Multiphase Saturation PK: Upwind for fractional flow only supports <fractional_flow> and <df_dS>\n";
      Exceptions::amanzi_throw(msg); 
    }
    return 0.0;
  }

  // access methods
  int Method() { return method_; }

  std::vector<Teuchos::RCP<WaterRetentionModel> >& WRM() { return WRM_; }

  Teuchos::RCP<CompositeVector> Frac_Flow() { return fractional_flow_w_; }
  Teuchos::RCP<CompositeVector> dF_dS() { return dfw_dS_; }

  //int method() { return method_; }
  const Epetra_IntVector& map_c2mb() { return *map_c2mb_; }

 private:
  void ProcessParameterList_(Teuchos::ParameterList& plist);
  void PopulateMapC2MB_();

 private:
  int method_;

  double mu_w_, mu_nw_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::RCP<ParallelCommunication> pp_;

  std::vector<Teuchos::RCP<WaterRetentionModel> > WRM_;
  std::string phase_;

  //int method_;  // method for calculating relative permeability
  const Teuchos::RCP<const RelativePermeability> rel_perm_w_;  // realitive permeability wetting phase
  const Teuchos::RCP<const RelativePermeability> rel_perm_nw_;  // realitive permeability non wetting phase
  
  Teuchos::RCP<CompositeVector> fractional_flow_w_;
  Teuchos::RCP<CompositeVector> dfw_dS_;
  Teuchos::RCP<CompositeVector> total_mobility_;

  Teuchos::RCP<Epetra_IntVector> map_c2mb_;  // cell->model map

 protected:
  VerboseObject* vo_;
};

typedef double(FractionalFlow::*FractionalFlowUpwindFn)(int c, double Sw) const;
}  // namespace Flow
}  // namespace Amanzi

#endif

