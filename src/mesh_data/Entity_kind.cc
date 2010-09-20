#include "Entity_kind.hh"

namespace Mesh_data 
{


bool valid_entity_kind (Entity_kind kind)
{
    return (kind >= NODE && kind <= CELL);
}

}
