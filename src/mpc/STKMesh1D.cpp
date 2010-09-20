#include <iostream>

#include "STKMesh1D.hpp"

#include <Shards_BasicTopologies.hpp>

#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/Selector.hpp>

#include <stk_mesh/fem/FieldDeclarations.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>
#include <stk_mesh/fem/EntityRanks.hpp>

enum { field_data_chunk_size = 1001 };

STKMesh1D::STKMesh1D(const Teuchos::RCP<const Epetra_Comm>& comm,
		     const Teuchos::RCP<Teuchos::ParameterList>& params) 
{

  // Create global mesh
  const int nelem = params->get<int>("1D Elements");
  const double xleft  = params->get("Left Coordinate",     0.0);
  const double xright = params->get("Right Coordinate",     1.0);
  std::vector<double> x(nelem+1);
  double h = (xright-xleft)/nelem;
  for (unsigned int i=0; i<=nelem; i++) x[i] = xleft+h*i;

  // Distribute the elements equally among processors
  Teuchos::RCP<Epetra_Map> elem_map = Teuchos::rcp(new Epetra_Map(nelem, 0, *comm));
  int numMyElements = elem_map->NumMyElements();

  //Start STK stuff, from UseCase_2 constructor
  metaData = new stk::mesh::MetaData(stk::mesh::fem_entity_rank_names() );
  bulkData = new stk::mesh::BulkData(*metaData , MPI_COMM_WORLD , field_data_chunk_size );
  coordinates_field = & metaData->declare_field< VectorFieldType >( "coordinates" );

  partVec.resize(1);
  partVec[0] = &  metaData->declare_part( "Block_1", stk::mesh::Element );

  nsPartVec.resize(2);
  nsPartVec[0] = & metaData->declare_part( "NodeSet0", stk::mesh::Node );
  nsPartVec[1] = & metaData->declare_part( "NodeSet1", stk::mesh::Node );

  stk::mesh::set_cell_topology< shards::Line<2> >(*partVec[0]);
  stk::mesh::put_field( *coordinates_field , stk::mesh::Node , metaData->universal_part() , 1 );
  metaData->commit();

  // Finished with metaData, now work on bulk data

  // STK
  bulkData->modification_begin(); // Begin modifying the mesh
  std::vector<stk::mesh::Part*> noPartVec;

  int rightNode=0;
  // Create elements and node IDs
  for (unsigned int i=0; i<numMyElements; i++) {
    const unsigned int elem_GID = elem_map->GID(i);
    
    const unsigned int left_node  = elem_GID;
    unsigned int right_node = left_node+1;
    if (rightNode < right_node) rightNode = right_node;

//Start STK stuff, from UseCase_2 populate
    stk::mesh::EntityId elem_id = (stk::mesh::EntityId) elem_GID;
// ELEM ID CAN NOT BE 0 ?????
    stk::mesh::Entity& edge  = bulkData->declare_entity(stk::mesh::Element, 1+elem_id, partVec);

    stk::mesh::Entity& lnode = bulkData->declare_entity(stk::mesh::Node, 1+left_node, noPartVec);
    stk::mesh::Entity& rnode = bulkData->declare_entity(stk::mesh::Node, 1+right_node, noPartVec);

    bulkData->declare_relation(edge, lnode, 0);
    bulkData->declare_relation(edge, rnode, 1);

    double* lnode_coord = stk::mesh::field_data(*coordinates_field, lnode);
    lnode_coord[0] = x[elem_GID];
    double* rnode_coord = stk::mesh::field_data(*coordinates_field, rnode);
    rnode_coord[0] = x[elem_GID+1];
    // Set node sets
    std::vector<stk::mesh::Part*> tmpPartVec(1);

    if (elem_GID==0) {
       tmpPartVec[0] = nsPartVec[0];
       bulkData->change_entity_parts(lnode, tmpPartVec);
    }
    if ((elem_GID+1)==elem_map->NumGlobalElements()) {
      tmpPartVec[0] = nsPartVec[1];
      bulkData->change_entity_parts(rnode, tmpPartVec);
    }
  }

  // STK
  bulkData->modification_end();
}




STKMesh1D::~STKMesh1D()
{ 
  delete metaData;
  delete bulkData;
}
