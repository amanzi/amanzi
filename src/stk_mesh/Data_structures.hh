#ifndef _DATA_STRUCTURES_H_
#define _DATA_STRUCTURES_H_

#include <stk_mesh/fem/FieldTraits.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <Teuchos_RCPDecl.hpp>

#include "Mesh.hh"


namespace STK_mesh
{

class Mesh;

typedef stk::mesh::Field<double, stk::mesh::Cartesian> Vector_field_type;
typedef stk::mesh::Field<double>                       Scalar_field_type;
typedef std::vector<stk::mesh::Entity*>                Entity_vector;
typedef std::vector<stk::mesh::EntityId>               Entity_Ids;

typedef std::pair<stk::mesh::EntityRank, unsigned int> Rank_and_id;
typedef std::map<Rank_and_id, stk::mesh::Part*> Id_map;

typedef Teuchos::RCP<STK_mesh::Mesh> Mesh_p;

typedef std::map<unsigned int, unsigned int> Index_map;
}

#endif /* _DATA_STRUCTURES_H_ */
