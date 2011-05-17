#ifndef _DATA_STRUCTURES_H_
#define _DATA_STRUCTURES_H_

#include <stk_mesh/fem/FieldTraits.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <Teuchos_RCPDecl.hpp>

namespace Amanzi {
namespace AmanziMesh {
namespace STK {

class Mesh_STK_Impl;

typedef stk::mesh::Field<double, stk::mesh::Cartesian> Vector_field_type;
typedef stk::mesh::Field<double>                       Scalar_field_type;
typedef stk::mesh::Field<unsigned>                     Id_field_type;
typedef std::vector<stk::mesh::Entity*>                Entity_vector;
typedef std::vector<stk::mesh::EntityId>               Entity_Ids;

typedef std::pair<stk::mesh::EntityRank, unsigned int> Rank_and_id;
typedef std::map<Rank_and_id, stk::mesh::Part*> Id_map;

typedef Teuchos::RCP<Mesh_STK_Impl> Mesh_STK_Impl_p;

typedef std::map<unsigned int, unsigned int> Index_map;

} // close namespace STK 
} // close namespace Mesh 
} // close namespace Amanzi 

#endif /* _DATA_STRUCTURES_H_ */
