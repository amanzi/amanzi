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

  A collection of fractrue aperture models along with a mesh partition.
*/

#include "dbc.hh"

#include "ApertureModel_BartonBandis.hh"
#include "ApertureModel_Constant.hh"
#include "ApertureModel_ExponentialLaw.hh"
#include "ApertureModelPartition.hh"

namespace Amanzi {
namespace Flow {

/* ******************************************************************
* Non-member factory.
****************************************************************** */
Teuchos::RCP<ApertureModelPartition>
CreateApertureModelPartition(Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                             Teuchos::RCP<Teuchos::ParameterList> plist)
{
  std::vector<Teuchos::RCP<ApertureModel>> fam_list;
  std::vector<std::vector<std::string>> region_list;

  for (auto lcv = plist->begin(); lcv != plist->end(); ++lcv) {
    std::string name = lcv->first;
    if (plist->isSublist(name)) {
      Teuchos::ParameterList sublist = plist->sublist(name);
      region_list.push_back(sublist.get<Teuchos::Array<std::string>>("regions").toVector());

      std::string model = sublist.get<std::string>("fracture aperture model");
      if (model == "Barton-Bandis") {
        fam_list.push_back(Teuchos::rcp(new ApertureModel_BartonBandis(sublist)));
      } else if (model == "exponential law") {
        fam_list.push_back(Teuchos::rcp(new ApertureModel_ExponentialLaw(sublist)));
      } else if (model == "constant") {
        fam_list.push_back(Teuchos::rcp(new ApertureModel_Constant(sublist)));
      } else {
        AMANZI_ASSERT(0);
      }
    } else {
      AMANZI_ASSERT(0);
    }
  }

  auto partition = Teuchos::rcp(new Functions::MeshPartition());
  partition->Initialize(mesh, AmanziMesh::CELL, region_list, -1);
  partition->Verify();

  return Teuchos::rcp(new ApertureModelPartition(partition, fam_list));
}

} // namespace Flow
} // namespace Amanzi
