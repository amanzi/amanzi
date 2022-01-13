/*
  This is the mpc_pk component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Process kernel for coupling of Flow_PK and ShallowWater_PK.
*/


#ifndef AMANZI_SURFACESUBSURFACE_PK_HH_
#define AMANZI_SURFACESUBSURFACE_PK_HH_

#include "Teuchos_RCP.hpp"

#include "PK.hh"
#include "PK_BDF.hh"
#include "Flow_PK.hh"
#include "PK_Factory.hh"
#include "PK_MPCSubcycled.hh"

namespace Amanzi {

class SurfaceSubsurface_PK : public PK_MPCSubcycled {
 public:
  SurfaceSubsurface_PK(Teuchos::ParameterList& pk_tree,
                       const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                       const Teuchos::RCP<State>& S,
                       const Teuchos::RCP<TreeVector>& soln);

  ~SurfaceSubsurface_PK() {};

  // PK methods
  // -- dt is the minimum of the sub pks
  virtual double get_dt();
  virtual void set_dt(double dt);
  
  // Setup and initialization
  virtual void Initialize();

  // -- advance each sub pk from t_old to t_new.
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit = false);

  virtual void CommitStep(double t_old, double t_new, const Tag& tag);

  std::string name() { return "surface subsurface";}
  
 protected:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_domain_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_surface_;

 private:
  // factory registration
  static RegisteredPKFactory<SurfaceSubsurface_PK> reg_;
};

}  // namespace Amanzi
#endif
