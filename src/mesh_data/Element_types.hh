#ifndef _ELEMENT_TYPES_HH_
#define _ELEMENT_TYPES_HH_

#include <string>

namespace Mesh_data
{

enum ELEMENT_TYPE
{
    UNKNOWN = -1,
    CIRCLE = 0,
    SPHERE,
    TRUSS,
    BEAM,
    TRIANGLE,
    QUAD,
    SHELL,
    TETRA,
    PYRAMID,
    WEDGE,
    HEX
};

std::string type_to_name (ELEMENT_TYPE type);

bool ok_type (ELEMENT_TYPE);

}

#endif
