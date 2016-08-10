#include "Mesh.hh"
#include "writeMesh3D.hh"
#include "readMesh2D.hh"


int main() {
  using namespace Amanzi::AmanziGeometry;

  std::string mesh_in = "Mesh.txt";
  std::string mesh_out = "Mesh3D.exo";

  std::vector<int> soil_type;
  std::vector<int> bedrock_type;
  std::vector<double> depths;

  Mesh2D m = readMesh2D_text(mesh_in, soil_type, bedrock_type, depths);
  // ASSERT(m.coords.size() == 4);
  // ASSERT(m.cell2node.size() == 2);
  
  Mesh3D m3(m, 3);
  m3.extrude(depths, soil_type);
  m3.extrude(10, bedrock_type);
  m3.extrude(10, bedrock_type);
  m3.finish_sets();

  std::cout << "NNodes on the surf = " << m.coords.size() << std::endl;
  std::cout << "Ncells on the surf = " << m.cell2node.size() << std::endl;
  std::cout << "NNodes on 3D = " << m3.coords.size() << std::endl;
  std::cout << "Ncells on 3D = " << m3.cell2face.size() << std::endl;
  ASSERT(m3.coords.size() == (4*m.coords.size()));
  ASSERT(m3.cell2face.size() == (3*m.cell2node.size()));

  writeMesh3D_exodus(m3, mesh_out);
  return 0;
}
