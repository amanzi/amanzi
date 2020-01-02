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

  double ValuePrimaryVar(int c, double PrimaryVar) const { return PrimaryVar; }

  // access methods
  Teuchos::RCP<WRMmpPartition>& wrm() { return wrm_; }
  
  Teuchos::RCP<const AmanziMesh::Mesh> Mesh() { return mesh_; }
  Teuchos::RCP<CompositeVector> dKdS() { return dKdS_; }
  Teuchos::RCP<CompositeVector> Krel() { return Krel_; }
  Teuchos::RCP<CompositeVector> Coeff() { return mpCoeff_; }
  Teuchos::RCP<CompositeVector> RhoDerivKrel() { return rhoDerivKrel_; }
  Teuchos::RCP<CompositeVector> DPc() { return dPc_; }

 private:
  void ProcessParameterList_(Teuchos::ParameterList& plist);

 private:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;

  Teuchos::RCP<WRMmpPartition> wrm_;
  std::string phase_;

  //int method_;  // method for calculating relative permeability
  Teuchos::RCP<CompositeVector> Krel_;  // realitive permeability 
  Teuchos::RCP<CompositeVector> rho_;
  Teuchos::RCP<CompositeVector> mpCoeff_;
  Teuchos::RCP<CompositeVector> rhoDerivKrel_;
  Teuchos::RCP<CompositeVector> dKdS_;  // derivative of realitive permeability wrt saturation
  Teuchos::RCP<CompositeVector> dPc_;

  double Cg_;

 protected:
  VerboseObject* vo_;
};

typedef double(MPCoeff::*CoefUpwindFn1)(int c, double s) const; 
typedef double(MPCoeff::*CoefUpwindFn2)(int c, double x, double s) const; 

}  // namespace Flow
}  // namespace Amanzi

#endif

