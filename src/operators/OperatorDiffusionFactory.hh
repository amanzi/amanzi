/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_OPERATOR_DIFFUSION_FACTORY_HH_
#define AMANZI_OPERATOR_DIFFUSION_FACTORY_HH_

#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "OperatorDiffusion.hh"
#include "OperatorDiffusionWithGravity.hh"

namespace Amanzi {
namespace Operators {

class BCs;

struct OperatorDiffusionFactory {
  // Diffusion operators with optional gravity.
  // Decision is made based on data in the parameter list.
  Teuchos::RCP<OperatorDiffusion>
  Create(Teuchos::ParameterList& oplist,
         const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
         const Teuchos::RCP<BCs>& bc,
         double rho,
         const AmanziGeometry::Point& g);

  Teuchos::RCP<OperatorDiffusion>
  Create(Teuchos::ParameterList& oplist,
         const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
         const Teuchos::RCP<BCs>& bc,
         const Teuchos::RCP<const CompositeVector>& rho,
         const AmanziGeometry::Point& g);

  // Diffusion operators without gravity.
  Teuchos::RCP<OperatorDiffusion>
  Create(Teuchos::ParameterList& oplist,
         const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
         const Teuchos::RCP<BCs>& bc);

  Teuchos::RCP<OperatorDiffusion>
  Create(Teuchos::ParameterList& oplist,
         const Teuchos::RCP<const AmanziMesh::Mesh>& mesh);
  
  Teuchos::RCP<OperatorDiffusion>
  Create(Teuchos::ParameterList& oplist,
         const Teuchos::RCP<Operator>& global_op);

  // Diffusion operators with gravity.
  Teuchos::RCP<OperatorDiffusionWithGravity>
  CreateWithGravity(Teuchos::ParameterList& oplist,
                    const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                    const Teuchos::RCP<BCs>& bc);
                    
  Teuchos::RCP<OperatorDiffusionWithGravity>
  CreateWithGravity(Teuchos::ParameterList& oplist,
                    const Teuchos::RCP<Operator>& global_op,
                    const Teuchos::RCP<BCs>& bc);

  Teuchos::RCP<OperatorDiffusionWithGravity>
  CreateWithGravity(Teuchos::ParameterList& oplist,
                    const Teuchos::RCP<const AmanziMesh::Mesh>& mesh);
                    
  Teuchos::RCP<OperatorDiffusionWithGravity>
  CreateWithGravity(Teuchos::ParameterList& oplist,
                    const Teuchos::RCP<Operator>& global_op);
  
 private:
  inline void SetCellSchema_(Teuchos::ParameterList& oplist);
  inline void SetCellFaceSchema_(Teuchos::ParameterList& oplist);
};

}  // namespace Operators
}  // namespace Amanzi

#endif
