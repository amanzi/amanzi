/*
  Flow PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)

  Reliative permeability as a function of capillary pressure, k=k(pc).
*/

#ifndef AMANZI_FLOW_REL_PERM_HH_
#define AMANZI_FLOW_REL_PERM_HH_

#include "Teuchos_ParameterList.hpp"

#include "Mesh.hh"
#include "secondary_variable_field_evaluator.hh"
#include "Tensor.hh"
#include "WRM.hh"
#include "WRMPartition.hh"

namespace Amanzi {
namespace Flow {

class RelPerm {
 public:
  RelPerm(Teuchos::ParameterList& plist,
          Teuchos::RCP<const AmanziMesh::Mesh> mesh,
          double patm,
          const Teuchos::RCP<WRMPartition>& wrm);

  // main members
  void Compute(Teuchos::RCP<const CompositeVector> p,
               const std::vector<int>& bc_model,
               const std::vector<double>& bc_value,
               const Teuchos::RCP<CompositeVector>& krel);
  void ComputeDerivative(Teuchos::RCP<const CompositeVector> p,
                         const std::vector<int>& bc_model,
                         const std::vector<double>& bc_value,
                         const Teuchos::RCP<CompositeVector>& dKdP);

  // WRM relations
  double Compute(int c, double p) const;
  double ComputeDerivative(int c, double p) const;

  void Compute_dSdP(const Epetra_MultiVector& p, Epetra_MultiVector& dSdP);

  // other
  void PlotWRMcurves();

 protected:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::RCP<WRMPartition> wrm_;

  double patm_;
};

typedef double(RelPerm::*RelPermUpwindFn)(int c, double p) const; 

}  // namespace Flow
}  // namespace Amanzi

#endif
