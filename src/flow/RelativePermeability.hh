/*
This is the flow component of the Amanzi code. 

Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided in the top-level COPYRIGHT file.

Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef __RELATIVE_PERMEABILITY_HH__
#define __RELATIVE_PERMEABILITY_HH__

#include <string>
#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Mesh.hh"
#include "WaterRetentionModel.hh"

namespace Amanzi {
namespace AmanziFlow {

class RelativePermeability {
 public:
  RelativePermeability(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                       Teuchos::ParameterList& list)
      : mesh_(mesh), list_(list) {};
  ~RelativePermeability() {};

  // main methods
  void Init();

  void VerifyWRMparameters(double m, double alpha, double sr, double pc0);
  void VerifyStringMualemBurdine(const std::string name);

  // access methods
  std::vector<Teuchos::RCP<WaterRetentionModel> >& WRM() { return WRM_; }

 private:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::ParameterList list_;
  std::vector<Teuchos::RCP<WaterRetentionModel> > WRM_;
};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif

