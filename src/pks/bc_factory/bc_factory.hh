#ifndef AMANZI_BC_FACTORY_HH_
#define AMANZI_BC_FACTORY_HH_

/* -------------------------------------------------------------------------
ATS

Author: ...
    Ethan Coon (ATS version) (ecoon@lanl.gov)
 ------------------------------------------------------------------------- */

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Point.hh"
#include "Mesh.hh"
#include "boundary_function.hh"

namespace Amanzi {

class BCFactory {

public:
  BCFactory(const Teuchos::RCP<const AmanziMesh::Mesh> &mesh,
                const Teuchos::ParameterList& plist)
     : mesh_(mesh), plist_(plist) {}

  Teuchos::RCP<Functions::BoundaryFunction>
  CreateWithFunction(std::string list_name, std::string function_name) const;

  Teuchos::RCP<Functions::BoundaryFunction>
  CreateWithoutFunction(std::string list_name) const;

 private:

  void ProcessListWithFunction_(const Teuchos::ParameterList&,
          std::string function_name,
          const Teuchos::RCP<Functions::BoundaryFunction>&) const;

  void ProcessListWithoutFunction_(const Teuchos::ParameterList&,
          const Teuchos::RCP<Functions::BoundaryFunction>&) const;

  void ProcessSpecWithFunction_(const Teuchos::ParameterList&,
          std::string function_name,
          const Teuchos::RCP<Functions::BoundaryFunction>&) const;

  void ProcessSpecWithoutFunction_(const Teuchos::ParameterList&,
          const Teuchos::RCP<Functions::BoundaryFunction>&) const;

 private:
  const Teuchos::RCP<const AmanziMesh::Mesh>& mesh_;
  Teuchos::ParameterList plist_;
};

}  // namespace

#endif // AMANZI_BC_FACTORY_HH_
