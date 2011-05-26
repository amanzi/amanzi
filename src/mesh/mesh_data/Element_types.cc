#include "Element_types.hh"

#include <string>

#include "dbc.hh"

namespace Amanzi {
namespace AmanziMesh {
namespace Data {

std::string type_to_name (Cell_type type)
{

  ASSERT (cell_valid_type(type));

  switch (type)
  {
    case TRI:
      return "triangle";
    case QUAD:
      return "quad";
    case POLYGON:
      return "polygon";
    case TET:
      return "tetrahedron";
    case PYRAMID:
      return "pyramid";
    case PRISM:
      return "prism";
    case HEX:
      return "hexahedron";
    case POLYHED:
      return "polyhedron";
    default:
      return "unknown";
  }
}

} // namespace Data
} // namespace Mesh
} // namespace Amanzi
