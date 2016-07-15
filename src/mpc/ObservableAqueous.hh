/*
  Multi-Process Coordinator

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Markus Berndt
           Konstantin Lipnikov
           Daniil Svyatsky
*/

#ifndef AMANZI_OBSERVABLE_AQUEOUS_HH
#define AMANZI_OBSERVABLE_AQUEOUS_HH

#include "Observable.hh"
#include "Units.hh"

namespace Amanzi{

class ObservableAqueous : public virtual Observable{
 public:
  ObservableAqueous(std::string variable,
                    std::string region,
                    std::string functional,
                    Teuchos::ParameterList& plist,
                    Teuchos::ParameterList& units_plist,
                    Teuchos::RCP<AmanziMesh::Mesh> mesh);

  virtual void ComputeObservation(State& S, double* value, double* volume);
  virtual int ComputeRegionSize();

 protected:
  double CalculateWaterTable_(State& S, AmanziMesh::Entity_ID_List& ids);

  bool obs_boundary_, obs_planar_;
  std::vector<double> vofs_;
  AmanziGeometry::Point reg_normal_;
};

}  // namespace Amanzi

#endif
