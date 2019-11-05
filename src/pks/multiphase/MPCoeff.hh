/*
This is the multiphase component of the Amanzi code. 

Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided in the top-level COPYRIGHT file.

Authors: Quan Bui (mquanbui@math.umd.edu)
*/

#ifndef AMANZI_MULTIPHASE_MP_COEF_HH_
#define AMANZI_MULTIPHASE_MP_COEF_HH_

#include <string>
#include <vector>

#include "Epetra_IntVector.h"
#include "Epetra_Vector.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Mesh.hh"
#include "CompositeVector.hh"
#include "ParallelCommunication.hh"
#include "Tensor.hh"
#include "VerboseObject.hh"

#include "WaterRetentionModel.hh"
#include "MultiphaseTypeDefs.hh"

namespace Amanzi {
namespace Multiphase {

class MPCoeff {
 public:
  MPCoeff(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) : mesh_(mesh) {};
  ~MPCoeff() { delete vo_; }

  // main methods
  void Init(std::string phase_name, Teuchos::ParameterList& plist);
  void Init(std::string phase_name, Teuchos::ParameterList& plist, double Cg);
  void Compute(const CompositeVector& p, const CompositeVector& s,
               const std::vector<int>& bc_model, const std::vector<double>& bc_value);
  void Compute(const CompositeVector& s,
               const std::vector<int>& bc_model, const std::vector<double>& bc_value);

  double ValueKrel(int c, double Sw) const { return WRM_[(*map_c2mb_)[c]]->k_relative(Sw, phase_); } 
  double ValuePrimaryVar(int c, double PrimaryVar) const { return PrimaryVar; }
  double ValueRhoKrel(int c, double PrimaryVar, double Sw) const { 
    if (phase_ == "non wetting")
      return Cg_ * (PrimaryVar + WRM_[(*map_c2mb_)[c]]->capillaryPressure(Sw)) * WRM_[(*map_c2mb_)[c]]->k_relative(Sw, phase_);
    else
      return PrimaryVar * WRM_[(*map_c2mb_)[c]]->k_relative(Sw, phase_);
  }
  double ValueRhoDerivKrel(int c, double PrimaryVar, double Sw) const { 
    if (phase_ == "non wetting")
      return Cg_ * (PrimaryVar + WRM_[(*map_c2mb_)[c]]->capillaryPressure(Sw)) * WRM_[(*map_c2mb_)[c]]->dKdS(Sw, phase_);
    else
      return PrimaryVar * WRM_[(*map_c2mb_)[c]]->dKdS(Sw, phase_);
  }

  // wrt Sw
  double DerivativeKrel(int c, double Sw) const { return WRM_[(*map_c2mb_)[c]]->dKdS(Sw, phase_); } 

  // wrt Sw
  double DerivativePc(int c, double Sw) const { return WRM_[(*map_c2mb_)[c]]->dPc_dS(Sw); } 

  // access methods
  std::vector<Teuchos::RCP<WaterRetentionModel> >& WRM() { return WRM_; }
  
  Teuchos::RCP<const AmanziMesh::Mesh> Mesh() { return mesh_; }
  Teuchos::RCP<CompositeVector> dKdS() { return dKdS_; }
  Teuchos::RCP<CompositeVector> Krel() { return Krel_; }
  Teuchos::RCP<CompositeVector> Coeff() { return mpCoeff_; }
  Teuchos::RCP<CompositeVector> RhoDerivKrel() { return rhoDerivKrel_; }
  Teuchos::RCP<CompositeVector> DPc() { return dPc_; }

  //int method() { return method_; }
  const Epetra_IntVector& map_c2mb() { return *map_c2mb_; }

 private:
  void ProcessParameterList_(Teuchos::ParameterList& plist);
  void PopulateMapC2MB_();

 private:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::RCP<ParallelCommunication> pp_;

  std::vector<Teuchos::RCP<WaterRetentionModel> > WRM_;
  std::string phase_;
  //double atm_pressure;

  //int method_;  // method for calculating relative permeability
  Teuchos::RCP<CompositeVector> Krel_;  // realitive permeability 
  Teuchos::RCP<CompositeVector> rho_;
  Teuchos::RCP<CompositeVector> mpCoeff_;
  Teuchos::RCP<CompositeVector> rhoDerivKrel_;
  Teuchos::RCP<CompositeVector> dKdS_;  // derivative of realitive permeability wrt saturation
  Teuchos::RCP<CompositeVector> dPc_;

  double Cg_;

  Teuchos::RCP<Epetra_IntVector> map_c2mb_;  // cell->model map

 protected:
  VerboseObject* vo_;
};

typedef double(MPCoeff::*CoefUpwindFn1)(int c, double s) const; 
typedef double(MPCoeff::*CoefUpwindFn2)(int c, double x, double s) const; 

}  // namespace Flow
}  // namespace Amanzi

#endif

