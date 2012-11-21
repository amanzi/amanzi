#ifndef _ENTITY_MAPS_H_
#define _ENTITY_MAPS_H_

#include "dbc.hh"
#include "MeshDefs.hh"

#include <stk_mesh/fem/FEMMetaData.hpp>
#include <stk_mesh/base/GetEntities.hpp>

#include <map>


namespace Amanzi {
namespace AmanziMesh {
namespace STK {

/* Class Entity_map
 *
 * Models the relationship between STK_mesh entity labels and
 * Mesh_data labels and a consecutive indexed set of Mesh_data labels
 * (to facilitate looping over entity kinds).
 */

class Entity_map
{

  const stk::mesh::fem::FEMMetaData *meta_data_;

  unsigned int dimension_;

  stk::mesh::EntityRank ranks_  [3];
  Entity_kind kinds_ [3];
  
  std::map<stk::mesh::EntityRank, Entity_kind> rank_to_kind_;
  std::map<Entity_kind, stk::mesh::EntityRank> kind_to_rank_;
  
  inline void create_maps_ ();

public:
  
  Entity_map (const stk::mesh::fem::FEMMetaData *meta_data) : meta_data_ (meta_data),
                                                              dimension_(meta_data_->spatial_dimension())
  {
    create_maps_ ();
  }
  
  bool valid_dimension_(unsigned int dimension);

  unsigned int dimension () const { return dimension_; }
  
  stk::mesh::EntityRank  kind_to_rank (Entity_kind kind) const;
  Entity_kind rank_to_kind (stk::mesh::EntityRank   rank) const;
  
};

void Entity_map::create_maps_ ()
{

  rank_to_kind_ [meta_data_->node_rank()] = NODE;
  rank_to_kind_ [meta_data_->edge_rank()] = (dimension_ == 2) ? FACE : EDGE;
  rank_to_kind_ [meta_data_->face_rank()] = (dimension_ == 2) ? CELL : FACE;
  rank_to_kind_ [meta_data_->volume_rank()] = CELL;


  kind_to_rank_ [NODE] = meta_data_->node_rank();
  kind_to_rank_ [EDGE] = meta_data_->edge_rank();
  kind_to_rank_ [FACE] = (dimension_ == 2) ? meta_data_->edge_rank() : meta_data_->face_rank();
  kind_to_rank_ [CELL] = (dimension_ == 2) ? meta_data_->face_rank() : meta_data_->volume_rank();

}


} // close namespace STK 
} // close namespace Mesh 
} // close namespace Amanzi 




#endif /* _ENTITY_MAPS_H_ */
