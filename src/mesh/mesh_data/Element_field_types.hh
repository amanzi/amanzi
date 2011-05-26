#ifndef _ELEMENT_FIELD_TYPES_HH_
#define _ELEMENT_FIELD_TYPES_HH_

#include <string>

namespace Amanzi {
namespace AmanziMesh {
namespace Data {

// enum FIELD_LOCATION
// {
//     VERTEX,
//     EDGE,
//     FACE,
//     ELEMENT
// };

// bool ok_field_location (FIELD_LOCATION type);
// std::string location_to_name (FIELD_LOCATION type);


enum FIELD_TYPE
{
    SCALAR,
    VECTOR
};

bool ok_field_type (FIELD_TYPE type);
std::string type_to_name (FIELD_TYPE type);


} // namespace Data
} // namespace Mesh
} // namespace Amanzi
#endif
