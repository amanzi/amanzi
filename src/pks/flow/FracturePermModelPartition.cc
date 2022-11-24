/*
  Flow PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)

  A collection of fracture permeability models along with a mesh partition.
*/

#include "dbc.hh"

#include "FracturePermModel_Constant.hh"
#include "FracturePermModel_CubicLaw.hh"
#include "FracturePermModel_Linear.hh"
#include "FracturePermModelPartition.hh"

namespace Amanzi {
namespace Flow {

/* ******************************************************************
* Non-member factory.
****************************************************************** */
Teuchos::RCP<FracturePermModelPartition>
CreateFracturePermModelPartition(Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                                 Teuchos::RCP<Teuchos::ParameterList> plist)
{
  std::vector<Teuchos::RCP<FracturePermModel>> fpm_list;
  std::vector<std::string> region_list;

  for (auto lcv = plist->begin(); lcv != plist->end(); ++lcv) {
    std::string name = lcv->first;
    if (plist->isSublist(name)) {
      Teuchos::ParameterList sublist = plist->sublist(name);
      region_list.push_back(sublist.get<std::string>("region"));

      std::string model = sublist.get<std::string>("model");
      if (model == "cubic law") {
        fpm_list.push_back(Teuchos::rcp(new FracturePermModel_CubicLaw(sublist)));
      } else if (model == "linear") {
        fpm_list.push_back(Teuchos::rcp(new FracturePermModel_Linear(sublist)));
      } else if (model == "constant") {
        fpm_list.push_back(Teuchos::rcp(new FracturePermModel_Constant(sublist)));
      } else {
        AMANZI_ASSERT(0);
      }
    } else {
      AMANZI_ASSERT(0);
    }
  }

  auto partition = Teuchos::rcp(new Functions::MeshPartition(AmanziMesh::CELL, region_list));
  partition->Initialize(mesh, -1);
  partition->Verify();

  return Teuchos::rcp(new FracturePermModelPartition(partition, fpm_list));
}

} // namespace Flow
} // namespace Amanzi
