#ifndef MESH_WRITER_HH_
#define MESH_WRITER_HH_

#include "Mesh.hh"


namespace Amanzi {
namespace AmanziGeometry {

void writeExodus(const Mesh3D& m, const std::string& filename);

}
}


#endif
