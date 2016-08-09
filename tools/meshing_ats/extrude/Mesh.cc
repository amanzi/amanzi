#include "Mesh.hh"

namespace Amanzi {
namespace AmanziGeometry {

Mesh3D::Mesh3D(const Mesh2D& m_, int n_layers) :
    m(m_),
    current_layer(0),
    total_layers(n_layers) {
  int n_coords_2d = m.coords.size();
  
  Point d(3);
  coords.resize(n_coords_2d*(n_layers+1), d);

  cell2face.resize(n_layers * m.cell2node.size());
  for (int i=0; i!=cell2face.size(); ++i) {
    cell2face[i].resize(5,-1);
  }
}


int
Mesh3D::face_constructor(const std::vector<int>& nodes,
                         int index_in_cell,
                         int cell,
                         bool guaranteed_new) {

  std::vector<int> nodes_sorted(nodes);
  std::sort(nodes_sorted.begin(), nodes_sorted.end());

  if (!guaranteed_new) {
    for (int f=0; f!=face2node.size(); ++f) {
      if (equal(nodes_sorted, face2node_sorted[f])) {
        side_face_counts[f]++;
        return f;
      }
    }
  }

  // not already created
  face2node.emplace_back(nodes);
  face2node_sorted.emplace_back(std::move(nodes_sorted));
  face_in_cell_when_created.push_back(index_in_cell);
  face_cell_when_created.push_back(cell);
  if (guaranteed_new) {
    side_face_counts.push_back(0);
  } else {
    side_face_counts.push_back(1);
  }
  return face2node.size()-1;
}
  
void
Mesh3D::extrude(const std::vector<double>& dz,
                const std::vector<int>& cell_set) {
  std::cout << "PRE-Extruding: currently " << cell2face.size() << " cells and " << face2node.size() << " faces." << std::endl;

  if (current_layer == 0) {
    // surface set
    std::vector<int> surface_c(m.cell2node.size(), -1);

    // copy the top layer of coords
    std::copy(m.coords.begin(), m.coords.end(), coords.begin());
    
    // create the top layer of faces
    for (int c=0; c!=m.cell2node.size(); ++c) {
      int f = face_constructor(m.cell2node[c], 0, c, true);
      cell2face[c][0] = f;
      surface_c[c] = c;
    }

    std::vector<int> surface_f(surface_c.size(), 0);
    face_sets.emplace_back(std::make_pair(surface_c, surface_f));
  }

  // shift the coordinates
  for (int n=0; n!=m.coords.size(); ++n) {
    coords[node_structure(n, current_layer+1)] =
        coords[node_structure(n, current_layer)];
    coords[node_structure(n,current_layer+1)][2] -= dz[n];
  }

  // potentially the list of bottom faces
  std::vector<int> bottom(m.cell2node.size(), -1);
    
  for (int c=0; c!=m.cell2node.size(); ++c) {
    int c3 = cell_structure(c, current_layer);
    
    std::vector<int> top_nodes(3,-1), bottom_nodes(3,-1);
    for (int n=0; n!=3; ++n) {
      top_nodes[n] = node_structure(m.cell2node[c][n], current_layer);
      bottom_nodes[n] = node_structure(m.cell2node[c][n], current_layer+1);
    }
    
    // add the bottom face
    int f = face_constructor(bottom_nodes, 1, c3, true);
    cell2face[c3][1] = f;
    if (current_layer != total_layers-1) {
      cell2face[cell_structure(c, current_layer+1)][0] = f;
    } else {
      bottom[c] = c3;
    }
    
    // add the side faces
    for (int n=0; n!=3; ++n) {
      std::vector<int> nodes = {top_nodes[n], top_nodes[(n+1)%3],
                                bottom_nodes[(n+1)%3], bottom_nodes[n]};
      int f = face_constructor(nodes, n+2, c3);
      cell2face[c3][n+2] = f;
    }
  }

  // create the bottom set
  if (current_layer == total_layers - 1) {
    std::vector<int> bottom_f(bottom.size(), 1);
    face_sets.emplace_back(std::make_pair(bottom,bottom_f));
  }

  // copy over the cell sets
  cell_sets.insert(cell_sets.end(), cell_set.begin(), cell_set.end());

  // increment
  current_layer++;

  std::cout << "POST-Extruding: currently " << cell2face.size() << " cells and " << face2node.size() << " faces." << std::endl;

}

void
Mesh3D::finish_sets() {
  std::vector<int> sides_c;
  std::vector<int> sides_f;
  
  for (int f=0; f!=side_face_counts.size(); ++f) {
    if (side_face_counts[f] == 1) {
      sides_c.push_back(face_cell_when_created[f]);
      sides_f.push_back(face_in_cell_when_created[f]);
    }
  }
  face_sets.emplace_back(std::make_pair(sides_c, sides_f));
}


}
}
