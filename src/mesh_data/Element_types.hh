#ifndef _ELEMENT_TYPES_HH_
#define _ELEMENT_TYPES_HH_

#include <string>

namespace Mesh_data
{

enum ELEMENT_TYPE
{
    CIRCLE,
    SPHERE,
    TRUSS,
    BEAM,
    TRIANGLE,
    QUAD,
    SHELL,
    TETRA,
    WEDGE,
    HEX,
    UNKNOWN
};

std::string type_to_name (ELEMENT_TYPE type);

bool ok_type (ELEMENT_TYPE);

}

#endif
