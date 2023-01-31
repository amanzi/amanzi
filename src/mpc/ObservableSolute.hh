/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Markus Berndt
           Konstantin Lipnikov
           Daniil Svyatsky
*/

/*
  Multi-Process Coordinator

*/

#ifndef AMANZI_OBSERVABLE_SOLUTE_HH
#define AMANZI_OBSERVABLE_SOLUTE_HH

#include "ObservableAmanzi.hh"
#include "Point.hh"

namespace Amanzi {

class ObservableSolute : public virtual Observable {
 public:
  ObservableSolute(std::string variable,
                   std::string region,
                   std::string functional,
                   Teuchos::ParameterList& plist,
                   Teuchos::ParameterList& units_plist,
                   Teuchos::RCP<const AmanziMesh::Mesh> mesh);

  virtual void
  ComputeObservation(State& S, double* value, double* volume, std::string& unit, double dt);
  virtual int ComputeRegionSize();
  void RegisterComponentNames(std::vector<std::string> comp_names, int num_liquid, int tcc_index)
  {
    comp_names_ = comp_names;
    num_liquid_ = num_liquid;
    tcc_index_ = tcc_index;
  }

 protected:
  bool obs_boundary_, obs_planar_;
  AmanziMesh::Double_View vofs_;
  AmanziGeometry::Point reg_normal_;
  Teuchos::Array<std::string> comp_names_;
  int num_liquid_, tcc_index_;
};

} // namespace Amanzi

#endif
