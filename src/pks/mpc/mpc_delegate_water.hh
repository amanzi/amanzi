// Delegate for heuristic corrections based upon coupled surface/subsurface water.


#ifndef AMANZI_MPC_DELEGATE_WATER_HH_
#define AMANZI_MPC_DELEGATE_WATER_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "VerboseObject.hh"
#include "TreeVector.hh"
#include "CompositeVector.hh"

namespace Amanzi {

class MPCDelegateWater {

 public:
  MPCDelegateWater(const Teuchos::RCP<Teuchos::ParameterList>& plist,
                   int i_domain, int i_surf);

  bool
  ModifyCorrection(double h, Teuchos::RCP<const TreeVector> res,
          Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> du);

 protected:
  int
  ModifyCorrection_WaterFaceLimiter_(double h, Teuchos::RCP<const TreeVector> res,
          Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> du);

  int
  ModifyCorrection_WaterSpurt_(double h, Teuchos::RCP<const TreeVector> res,
          Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> du);

 protected:
  Teuchos::RCP<Teuchos::ParameterList> plist_;

  // Preconditioned correction alteration
  // -- control
  bool cap_the_spurt_;
  bool damp_the_spurt_;

  // -- parameters
  double cap_size_;
  double face_limiter_;

  // Predictor alteration
  bool modify_predictor_heuristic_;

  // indices into the TreeVector
  int i_surf_;
  int i_domain_;

  // Verbose Object
  Teuchos::RCP<VerboseObject> vo_;
};

} // namespace



#endif 
