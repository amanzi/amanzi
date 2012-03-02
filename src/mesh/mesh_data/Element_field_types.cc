#include "Element_field_types.hh"

#include "dbc.hh"

namespace Amanzi {
namespace AmanziMesh {
namespace Data {

// bool ok_field_location (Entity_types type)
// {
//     return (type>=NODE) && (type<=ELEMENT);
// }

// std::string location_to_name (FIELD_LOCATION type)
// {
//     ASSERT (ok_field_location (type));

//     switch (type)
//     {
//     case VERTEX:
//         return "vertex";
//     case EDGE:
//         return "edge";
//     case FACE:
//         return "face";
//     case ELEMENT:
//         return "element";
//     }

// }

bool ok_field_type (FIELD_TYPE type)
{
    return (type>=SCALAR) && (type<=VECTOR);
}


std::string type_to_name (FIELD_TYPE type)
{
    ASSERT (ok_field_type (type));

    switch (type)
    {
    case SCALAR:
        return "scalar";
    case VECTOR:
        return "vector";
    default:
      return "unknown";
    }
}

} // namespace Data
} // namespace Mesh
} // namespace Amanzi
