#include "Element_types.hh"

#include <string>

#include "dbc.hh"

namespace Mesh_data
{

std::string type_to_name (Mesh_data::ELEMENT_TYPE type)
{

    ASSERT (ok_type (type));

    switch (type)
    {
    case Mesh_data::CIRCLE:
        return "circle";
    case Mesh_data::SPHERE:
        return "sphere";
    case Mesh_data::TRUSS:
        return "truss";
    case Mesh_data::BEAM:
        return "beam";
    case Mesh_data::TRIANGLE:
        return "triangle";
    case Mesh_data::QUAD:
        return "quad";
    case Mesh_data::SHELL:
        return "shell";
    case Mesh_data::TETRA:
        return "tetrahedron";
    case Mesh_data::WEDGE:
        return "wedge";
    case Mesh_data::HEX:
        return "hexahedron";
    default:
        return "unknown";
    }
}

bool ok_type (Mesh_data::ELEMENT_TYPE type)
{

    bool result = (type >= Mesh_data::CIRCLE);
    result &= (type <= Mesh_data::UNKNOWN);

    return result;
}

}
