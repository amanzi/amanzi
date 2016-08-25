#include <set>
#include "Mesh.hh"

namespace Amanzi {
namespace AmanziGeometry {

Mesh2D::Mesh2D(std::vector<Point>&& coords_,
               std::vector<std::vector<int> >&& cell2node_,
               std::vector<std::vector<int> >&& cell_sets_) :
    coords(std::move(coords_)),
    cell2node(std::move(cell2node_)),
    cell_sets(std::move(cell_sets_)),
    datum(2)
{
  double area_eps = 1.e-6;
  
  for (auto& set : cell_sets) {
    ASSERT(set.size() == cell2node.size());
  }

  nnodes = coords.size();
  ncells = cell2node.size();
  
  for (auto& c : cell2node) {
    Point v1(2), v2(2);
    v1.set(coords[c[1]][0] - coords[c[0]][0], coords[c[1]][1] - coords[c[0]][1]);
    v2.set(coords[c[2]][0] - coords[c[0]][0], coords[c[2]][1] - coords[c[0]][1]);
    Point cross = v1^v2;
    if (cross[0] < 0) {
      std::reverse(c.begin(), c.end());
    }
    if (std::abs(cross[0]) < area_eps) {
      std::cout << "Zero area triangle:" << std::endl
                << " " << coords[c[0]] << std::endl
                << " " << coords[c[1]] << std::endl
                << " " << coords[c[2]] << std::endl;
      ASSERT(false);
    }
  }

  for (int c=0; c!=cell2node.size(); ++c) {
    std::vector<int> c2f(3);
    
    std::vector<int> f1 = {cell2node[c][0], cell2node[c][1]};
    c2f[0] = face_constructor(f1, 0, c);
    std::vector<int> f2 = {cell2node[c][1], cell2node[c][2]};
    c2f[1] = face_constructor(f2, 1, c);
    std::vector<int> f3 = {cell2node[c][2], cell2node[c][0]};
    c2f[2] = face_constructor(f3, 2, c);      
    cell2face.emplace_back(c2f);
  }

  // set sizes
  nfaces = face2node.size();
  ASSERT(ncells == cell2face.size());

  // set boundaries
  std::vector<int> boundary_c, boundary_f;
  for (int lcv=0; lcv!=side_face_counts.size(); ++lcv) {
    if (side_face_counts[lcv] == 1) {
      boundary_c.push_back(face_cell_when_created[lcv]);
      boundary_f.push_back(face_in_cell_when_created[lcv]);
    }
  }
  boundary_faces = std::make_pair(std::move(boundary_c),
          std::move(boundary_f));

  // normalize the nodes
  int64_t x=0., y=0.;
  for (auto& p : coords) {
    x += p[0];
    y += p[1];
  }
  datum[0] = x / coords.size();
  datum[1] = y / coords.size();
  for (auto& p : coords) {
    p[0] -= datum[0];
    p[1] -= datum[1];
  }
}

int
Mesh2D::face_constructor(const std::vector<int>& nodes,
                         int index_in_cell,
                         int cell) {
  
  auto h = nodes[0] > nodes[1] ? hash(nodes[1], nodes[0]) :
      hash(nodes[0],nodes[1]);

  auto match = faces_sorted.find(h);
  if (match != faces_sorted.end()) {
    int f = match->second;
    side_face_counts[f]++;
    return f;
  }

  // not already created
  int f = face2node.size();
  faces_sorted[h] = f;
  face2node.emplace_back(nodes);
  face_in_cell_when_created.push_back(index_in_cell);
  face_cell_when_created.push_back(cell);
  side_face_counts.push_back(1);
  return f;
}


Mesh3D::Mesh3D(const Mesh2D& m_, int n_layers) :
    m(m_),
    current_layer(0),
    total_layers(n_layers),
    datum(m_.datum)
{

  // reserve/allocate space
  Point d(3);
  int n_nodes = m.nnodes*(n_layers+1);
  coords.resize(n_nodes, d);

  int n_cells = n_layers * m.ncells;
  cell2face.reserve(n_cells);

  int n_faces = n_layers * m.nfaces
      + (n_layers+1)*m.ncells;
  face2node.reserve(n_faces);

  // copy the top surface coords
  std::copy(m.coords.begin(), m.coords.end(), coords.begin());

  // create the top layer of faces
  face2node.insert(face2node.end(), m.cell2node.begin(),
                   m.cell2node.end());

  std::vector<int> top_c(m.ncells, -1);
  std::vector<int> top_f(m.ncells, 0);
  for (int i=0; i!=m.ncells; ++i)
    top_c[i] = i;
  face_sets.emplace_back(std::make_pair(std::move(top_c),
          std::move(top_f)));
  face_sets_id.push_back(1);

  // move the 2d cell sets to face sets on the surface
  std::set<int> set_ids;
  for (auto& part : m.cell_sets) {
    set_ids.insert(part.begin(), part.end());
  }
  for (int sid : set_ids) {
    std::vector<int> set_cells;
    for (auto& part : m.cell_sets) {
      for (int c=0; c!=part.size(); ++c) {
        if (part[c] == sid) {
          set_cells.push_back(c);
        }
      }
    }
    std::vector<int> set_faces(set_cells.size(), 0);
    face_sets.emplace_back(std::make_pair(std::move(set_cells),
            std::move(set_faces)));
    face_sets_id.push_back(sid);
  }
}


void
Mesh3D::extrude(const std::vector<double>& dz,
                const std::vector<int>& block_ids_) {
  ASSERT(dz.size() == coords.size());
  ASSERT(block_ids_.size() == coords.size());

  // shift the coordinates
  for (int n=0; n!=m.coords.size(); ++n) {
    coords[node_structure(n, current_layer+1)] =
        coords[node_structure(n, current_layer)];
    coords[node_structure(n,current_layer+1)][2] -= dz[n];
  }

  // add side faces
  int nc = cell2face.size();
  int nf = face2node.size();

  // add the side faces
  for (int f=0; f!=m.nfaces; ++f) {
    std::vector<int> nodes =
        { node_structure(m.face2node[f][1],current_layer),
          node_structure(m.face2node[f][0],current_layer),
          node_structure(m.face2node[f][0],current_layer+1),
          node_structure(m.face2node[f][1],current_layer+1)
        };
    face2node.emplace_back(std::move(nodes));
  }

  // add the bottom faces
  for (int c=0; c!=m.ncells; ++c) {
    std::vector<int> nodes =
        { node_structure(m.cell2node[c][0],current_layer+1),
          node_structure(m.cell2node[c][1],current_layer+1),
          node_structure(m.cell2node[c][2],current_layer+1)
        };
    if (current_layer+1 == total_layers) {
      nodes = { node_structure(m.cell2node[c][2],current_layer+1),
          node_structure(m.cell2node[c][1],current_layer+1),
          node_structure(m.cell2node[c][0],current_layer+1)
        };
    }
    face2node.emplace_back(std::move(nodes));
  }

  // add the cell
  for (int c=0; c!=m.ncells; ++c) {
    // cell2face
    std::vector<int> c2f =
        {
          nf - m.ncells + c,
          nf + m.nfaces + c,
          nf + m.cell2face[c][0],
          nf + m.cell2face[c][1],
          nf + m.cell2face[c][2]
        };
    cell2face.emplace_back(c2f);
  }
  
  // copy over the cell sets
  block_ids.insert(block_ids.end(), block_ids_.begin(), block_ids_.end());

  // increment
  current_layer++;

  std::cout << "POST-Extruding: currently " << cell2face.size() << " cells and " << face2node.size() << " faces." << std::endl;

}

void
Mesh3D::finish_sets() {
  // create the bottom set
  std::vector<int> bottom_c(m.ncells);
  int bottom_c_begin = cell2face.size() - m.ncells;
  for (int i=0; i!=m.ncells; ++i) bottom_c[i] = i + bottom_c_begin;
                                      
  std::vector<int> bottom_f(bottom_c.size(), 1);
  face_sets.emplace_back(std::make_pair(std::move(bottom_c),
          std::move(bottom_f)));
  face_sets_id.push_back(2);
  
  // side sets
  int n_boundary_faces = total_layers * m.boundary_faces.first.size();
  std::vector<int> sides_c, sides_f;
  sides_c.reserve(n_boundary_faces);
  sides_f.reserve(n_boundary_faces);
  auto prism_boundary_f(m.boundary_faces.second);
  for (auto& f: prism_boundary_f) f += 2;
  for (int layer=0; layer!=total_layers; ++layer) {
    int c_start = layer*m.ncells;
    for (auto c : m.boundary_faces.first)
      sides_c.push_back(c_start + c);

    sides_f.insert(sides_f.begin(),
                   prism_boundary_f.begin(),
                   prism_boundary_f.end());
  }
  
  face_sets.emplace_back(std::make_pair(std::move(sides_c),
          std::move(sides_f)));
  face_sets_id.push_back(3);
}


}
}
