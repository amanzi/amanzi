#include "STKMesh1D.hpp"
#include "STK1DDisc.hpp"
#include "Epetra_Export.h"

// Start of STK stuff
#include <stk_util/parallel/Parallel.hpp>

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

STK1DDisc::STK1DDisc (Teuchos::RCP<STKMesh1D>& stkMesh_,
				const Teuchos::RCP<const Epetra_Comm>& comm) :
  stkMesh(stkMesh_)
{

  // get the STK mesh data from stkMesh
  stk::mesh::MetaData& metaData = *stkMesh->metaData;
  stk::mesh::BulkData& bulkData = *stkMesh->bulkData;

  std::vector<stk::mesh::Part*>& nsPartVec = stkMesh->nsPartVec;

  STKMesh1D::VectorFieldType& coordinates_field
    = *stkMesh->coordinates_field;

  stk::mesh::Part& universalPart = metaData.universal_part();

  nodes_per_element = 2;

  // this is hardwired for 1D, so edges are the top entity rank
  stk::mesh::EntityRankEnum topEntityRank = stk::mesh::Edge;
  
 // Constructs overlap_map
  int row, col;
  //STK version
  // maps for owned nodes elements and unknowns
  stk::mesh::Selector select_owned_in_part =
    stk::mesh::Selector( universalPart ) &
    stk::mesh::Selector( metaData.locally_owned_part() );

  // the node map
  std::vector< stk::mesh::Entity * > nodes ;
  stk::mesh::get_selected_entities( select_owned_in_part ,
				    bulkData.buckets( stk::mesh::Node ) ,
				    nodes );

  std::vector<int> node_indices(nodes.size());
  for (int i=0; i < nodes.size(); i++) 
    node_indices[i] = nodes[i]->identifier() - 1;
  
  node_map = Teuchos::rcp(new Epetra_Map(-1, nodes.size(),
					 &(node_indices[0]), 0, *comm));


  // the element map
  std::vector< stk::mesh::Entity * > elements;
  stk::mesh::get_selected_entities( select_owned_in_part,
				   bulkData.buckets( stk::mesh::Element),
				   elements );

  std::vector<int> element_indices( elements.size() );
  
  for (int i=0; i < elements.size(); i++) 
    element_indices[i] = elements[i]->identifier() - 1;  
  element_map = Teuchos::rcp(new Epetra_Map(-1, elements.size(),
					    &(element_indices[0]), 0, *comm));


  // the unknown map
  node_indices.resize(nodes.size());
  for (int i=0; i < nodes.size(); i++)
    node_indices[i] = getnodeDOF(*nodes[i]);

  map = Teuchos::rcp(new Epetra_Map(-1, node_indices.size(),
                              &(node_indices[0]), 0, *comm));


  // maps for overlap unknowns
  stk::mesh::Selector select_overlap_in_part =
    stk::mesh::Selector( universalPart ) &
    ( stk::mesh::Selector( metaData.locally_owned_part() )
      | stk::mesh::Selector( metaData.globally_shared_part() ) );
  

  //  overlapnodes used for overlap map -- stored for changing coords
  stk::mesh::get_selected_entities( select_overlap_in_part ,
				    bulkData.buckets( stk::mesh::Node ) ,
				    overlapnodes );

  node_indices.resize(overlapnodes.size());
  for (int i=0; i < overlapnodes.size(); i++)
    node_indices[i] = getnodeDOF(*overlapnodes[i]);

  overlap_map = Teuchos::rcp(new Epetra_Map(-1, node_indices.size(),
					    &(node_indices[0]), 0, *comm));
  
  coordinates.resize(overlapnodes.size());
  
  // create overlap graph
  overlap_graph =
    Teuchos::rcp(new Epetra_CrsGraph(Copy, *overlap_map,
                                     nodes_per_element, false));
  

  std::vector< stk::mesh::Entity * > cells ;
  stk::mesh::get_selected_entities( select_owned_in_part ,
				    bulkData.buckets( stk::mesh::Element ) ,
				    cells );

  for (int i=0; i < cells.size(); i++) {
    stk::mesh::Entity& e = *cells[i];
    stk::mesh::PairIterRelation rel = e.relations();

    // loop over local nodes
    for (int j=0; j < rel.size(); j++) {
      stk::mesh::Entity& rowNode = * rel[j].entity();

      row = getnodeDOF(rowNode);
      for (int l=0; l < rel.size(); l++) {
	stk::mesh::Entity& colNode = * rel[l].entity();
	col = getnodeDOF(colNode);
	overlap_graph->InsertGlobalIndices(row, 1, &col);
      }
    }
  }
  overlap_graph->FillComplete();

  // Fill  elNodeID(el_LID, local node) => node_LID
  elNodeID.resize(cells.size());
  for (int i=0; i < cells.size(); i++) {
    stk::mesh::Entity& e = *cells[i];
    int  el_lid = i;

    // Q: Does (i==el_lid) always???
    TEST_FOR_EXCEPTION(el_lid<0 || el_lid >= cells.size(), std::logic_error,
		       "STK1D_Disc: el_lid out of range " << el_lid << " ("<<cells.size()<<")"<<endl);
    stk::mesh::PairIterRelation rel = e.relations();

    elNodeID[el_lid].resize(nodes_per_element);
    // loop over local nodes
    for (int j=0; j < rel.size(); j++) {
      stk::mesh::Entity& rowNode = * rel[j].entity();
      int node_gid = rowNode.identifier() - 1;
      int node_lid = overlap_map->LID(node_gid);
      TEST_FOR_EXCEPTION(node_lid<0, std::logic_error,
			 "STK1D_Disc: node_lid out of range " << node_lid << endl);
      elNodeID[el_lid][j] = node_lid;
    }
  }
  
  int  estNonzeroesPerRow=3;
  
  graph = Teuchos::rcp(new Epetra_CrsGraph(Copy, *map, estNonzeroesPerRow, false));
  
  // Create non-overlapped matrix using two maps and export object
  Epetra_Export exporter(*overlap_map, *map);
  graph->Export(*overlap_graph, exporter, Insert);
  graph->FillComplete();

  // STK: NodeSets

  //MB nodeSets.resize(nsPartVec.size());
  for (int ns=0; ns<nsPartVec.size(); ns++) {
    stk::mesh::Selector select_owned_in_nspart =
      stk::mesh::Selector( *nsPartVec[ns]) &
      stk::mesh::Selector( metaData.locally_owned_part() );

    stk::mesh::get_selected_entities( select_owned_in_nspart ,
                             bulkData.buckets( stk::mesh::Node ) ,
                             nodes );

    //MB nodeSets[ns].resize(nodes.size());
    //MB cout << "STKDisc: nodeset "<<ns<<" has size " << nodes.size() << endl;
    //MB for (int i=0; i < nodes.size(); i++) nodeSets[ns][i] = nodes[i]->identifier() - 1;
    
    
    
  }
}


STK1DDisc::~STK1DDisc () {};


inline int STK1DDisc::getnodeDOF(stk::mesh::Entity& node) const
{ return (node.identifier()-1); }


Teuchos::ArrayRCP<double>&  STK1DDisc::getCoordinates() const
{
  for (int i=0; i < overlapnodes.size(); i++)  {
    int node_gid = overlapnodes[i]->identifier() - 1;
    int node_lid = overlap_map->LID(node_gid);
    double* x = stk::mesh::field_data(*stkMesh->coordinates_field, *overlapnodes[i]);
    coordinates[node_lid ] = x[0];
  }

  return coordinates;
}


