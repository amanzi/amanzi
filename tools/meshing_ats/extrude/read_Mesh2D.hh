#ifndef MESH_READER_HH_
#define MESH_READER_HH_

#include "Mesh.hh"

namespace Amanzi {
namespace AmanziGeometry {


struct PointFactory {
  PointFactory() {}

  bool addPoint(const Point& p, int& id) {
    auto loc = std::find(points.begin(), points.end(), p);
    if (loc == points.end()) {
      id = points.size();
      points.push_back(p);
      return true;
    } else {
      id = loc - points.begin();
      return false;
    }
  }

  std::vector<Point> points;
};

Mesh2D readFile(const std::string& filename,
                std::vector<int>& soil_type,
                std::vector<int>& bedrock_type,
                std::vector<double>& depth_to_bedrock);

}
}


#endif
