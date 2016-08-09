#include "Mesh.hh"
#include "write_Mesh.hh"
#include "read_Mesh2d.hh"


int main() {
  using namespace Amanzi::AmanziGeometry;

  std::string mesh_in = "Mesh.txt";
  std::string mesh_out = "Mesh3D.exo";

  std::vector<int> soil_type;
  std::vector<int> bedrock_type;
  std::vector<double> depths;

  Mesh2D m = readFile(mesh_in, soil_type, bedrock_type, depths);
  ASSERT(m.coords.size() == 4);
  ASSERT(m.cell2node.size() == 2);
  
  Mesh3D m3(m, 3);
  m3.extrude(depths, soil_type);
  m3.extrude(10, bedrock_type);
  m3.extrude(10, bedrock_type);
  m3.finish_sets();


  ASSERT(m3.coords.size() == 4*4);
  ASSERT(m3.cell2face.size() == 2*3);
  ASSERT(m3.face2node.size() == 5*3+2*4);
  
  
  writeExodus(m3, mesh_out);
  return 0;
}
