#ifndef _ENTITY_KINDS_H_
#define _ENTITY_KINDS_H_


// A stk-mesh independent enumeration of mesh elements.
// In 2D EDGE and FACE are synonymous. In 3D, all of these are unique.

namespace Mesh_data
{


enum Entity_kind
{
    NODE,
    EDGE,
    FACE,
    CELL
};

bool valid_entity_kind (Entity_kind kind);
}

#endif /* _ENTITY_KINDS_H_ */
