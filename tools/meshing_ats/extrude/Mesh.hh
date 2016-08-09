#ifndef MESH_2D_HH_
#define MESH_2D_HH_

#include <vector>
#include <utility>
#include <algorithm>

#include "Point.hh"

namespace Amanzi {
namespace AmanziGeometry {

struct Mesh2D {
  Mesh2D(std::vector<Point> coords_,
         std::vector<std::vector<int> > cell2node_) :
      coords(std::move(coords_)),
      cell2node(std::move(cell2node_))
  {
    for (auto& c : cell2node) {
      Point v1(2), v2(2);
      v1.set(coords[c[1]][0] - coords[c[0]][0], coords[c[1]][1] - coords[c[0]][1]);
      v2.set(coords[c[2]][0] - coords[c[0]][0], coords[c[2]][1] - coords[c[0]][1]);
      Point cross = v1^v2;
      if (cross[0] < 0) {
        std::reverse(c.begin(), c.end());
      }
    }
  }
  
  std::vector<Point> coords;
  std::vector<std::vector<int> > cell2node;
};


struct Mesh3D {
  Mesh3D(const Mesh2D& m_, int n_layers);

  void extrude(double dz, const std::vector<int>& cell_set) {
    std::vector<double> dzs(m.coords.size(), dz);
    extrude(dzs, cell_set);
  }

  void extrude(const std::vector<double>& dzs, int my_cell_set) {
    std::vector<int> cell_set(m.coords.size(), my_cell_set);
    extrude(dzs, cell_set);
  }


  void extrude(double dz, int my_cell_set) {
    std::vector<double> dzs(m.coords.size(), dz);
    std::vector<int> cell_set(m.coords.size(), my_cell_set);
    extrude(dzs, cell_set);
  }

  int node_structure(int n_2d, int layer) {
    return n_2d + layer*m.coords.size();
  }
  int cell_structure(int c_2d, int layer) {
    return c_2d + layer*m.cell2node.size();
  }

  // assumes sorted!
  bool equal(const std::vector<int>& n1,
             const std::vector<int>& n2) {
    if (n1 == n2) return true;
    return false;
  }
  
  int face_constructor(const std::vector<int>& nodes,
                       int face_in_cell,
                       int cell,
                       bool guaranteed_new=false);  
  void extrude(const std::vector<double>& dz,
               const std::vector<int>& cell_set);
  void finish_sets(); 

  const Mesh2D& m;

  std::vector<Point> coords;
  std::vector<std::vector<int> > cell2face;
  std::vector<int> cell_sets;
  std::vector<std::vector<int> > face2node;
  std::vector<std::vector<int> > face2node_sorted;
  std::vector<int> side_face_counts;
  std::vector<int> face_in_cell_when_created;
  std::vector<int> face_cell_when_created;

  std::vector<std::pair<std::vector<int>, std::vector<int> > > face_sets;


  int current_layer;
  int total_layers;

};

}
}

#endif
