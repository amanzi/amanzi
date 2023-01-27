/*
  Pipe Flow PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Pipe flow model inherited from shallow water model. 
*/

#ifndef AMANZI_PIPE_FLOW_PK_HH_
#define AMANZI_PIPE_FLOW_PK_HH_

#include "ShallowWater_PK.hh"

namespace Amanzi {
namespace ShallowWater {

class PipeFlow_PK : public ShallowWater_PK {

 public:
  PipeFlow_PK(Teuchos::ParameterList& pk_tree,
                  const Teuchos::RCP<Teuchos::ParameterList>& glist,
                  const Teuchos::RCP<State>& S,
                  const Teuchos::RCP<TreeVector>& soln);

  ~PipeFlow_PK() {};

  virtual void Initialize() override;

 private:
  static RegisteredPKFactory<PipeFlow_PK> reg_;

  //TODO: why is gravity private in the shallow water PK?
  double g_;

  std::string hydrostatic_pressure_force_type_;

  double pipe_diameter_;

};


}  // namespace ShallowWater
}  // namespace Amanzi

#endif
