
#include "MeshWrapper.hpp"
#include "Epetra_Vector.h"

#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/fem/FieldTraits.hpp>
#include <stk_mesh/fem/FieldDeclarations.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>
#include <stk_mesh/fem/EntityRanks.hpp>

MeshWrapper::MeshWrapper( Teuchos::RCP<STKMesh1D> mesh1D_, 
			  Teuchos::RCP<STK1DDisc> disc1D_, 
			  Teuchos::RCP<DataLayout> data_layout_):
  mesh1D(mesh1D_),
  disc1D(disc1D_),
  data_layout(data_layout_)
{
  // compute the element volumes
  
  element_volumes = Teuchos::rcp(new Epetra_Vector( *data_layout->get_element_map() ));

  // here comes the mesh specific code that computed element volumes

  stk::mesh::MetaData *metaData = disc1D->getMesh()->metaData;
  stk::mesh::BulkData *bulkData = disc1D->getMesh()->bulkData;

  stk::mesh::Selector select_owned_in_part =
    stk::mesh::Selector( metaData->universal_part() ) &
    stk::mesh::Selector( metaData->locally_owned_part() );
  
  std::vector< stk::mesh::Entity * > cells ;
  stk::mesh::get_selected_entities( select_owned_in_part ,
                                    bulkData->buckets( stk::mesh::Element ),
                                    cells );  

  Teuchos::ArrayRCP<double> coords = disc1D->getCoordinates();
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<int> > elNodeID = disc1D->getElNodeID(); 
 
  int indices[2];

  for (int i=0; i<cells.size(); i++) {

    stk::mesh::Entity& e = *cells[i];
    stk::mesh::PairIterRelation rel = e.relations();    

    stk::mesh::Entity& lnode =  * rel[0].entity();
    stk::mesh::Entity& rnode =  * rel[1].entity();
    
    int lnode_gid = lnode.identifier() - 1;
    int rnode_gid = rnode.identifier() - 1;

    indices[0] = lnode_gid;
    indices[1] = rnode_gid;
    
    double h;
    (*element_volumes)[i] = coords[ elNodeID[i][1] ] - coords[ elNodeID[i][0] ];
  }
}
  
