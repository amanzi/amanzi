/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Flow PK

  A collection of porosity models along with a mesh partition.
*/

#include "dbc.hh"

#include "PorosityModel_Constant.hh"
#include "PorosityModel_Compressible.hh"
#include "PorosityModelPartition.hh"

namespace Amanzi {
namespace Flow {

/* ******************************************************************
* Non-member factory.
****************************************************************** */
Teuchos::RCP<PorosityModelPartition>
CreatePorosityModelPartition(Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                             Teuchos::RCP<Teuchos::ParameterList> plist)
{
  std::vector<Teuchos::RCP<PorosityModel>> pom_list;
  std::vector<std::vector<std::string>> region_list;

  for (auto lcv = plist->begin(); lcv != plist->end(); ++lcv) {
    std::string name = lcv->first;
    if (plist->isSublist(name)) {
      Teuchos::ParameterList sublist = plist->sublist(name);
      region_list.push_back(sublist.get<Teuchos::Array<std::string>>("regions").toVector());

      std::string model = sublist.get<std::string>("porosity model");
      if (model == "constant") {
        pom_list.push_back(Teuchos::rcp(new PorosityModel_Constant(sublist)));
      } else if (model == "compressible") {
        pom_list.push_back(Teuchos::rcp(new PorosityModel_Compressible(sublist)));
      } else {
        AMANZI_ASSERT(0);
      }
    } else {
      AMANZI_ASSERT(0);
    }
  }

  auto partition = Teuchos::rcp(new Functions::MeshPartition());
  partition->Initialize(mesh, AmanziMesh::Entity_kind::CELL, region_list, -1);
  partition->Verify();

  return Teuchos::rcp(new PorosityModelPartition(partition, pom_list));
}

} // namespace Flow
} // namespace Amanzi
