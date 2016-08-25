#include "Mesh.hh"
#include "writeMesh3D.hh"
#include "readMesh2D.hh"


int main() {
  using namespace Amanzi::AmanziGeometry;

  std::string mesh_in = "Mesh.txt";
  std::string mesh_out = "Mesh3D_VariableSoil.exo";
  
  std::vector<double> ref_soil_mlay_dz = {1.0e-2, 2.0e-2, 5.0e-2, 1.2e-1, 2.5e-1, 5.0e-1, 1.0e1};
  std::vector<double> ref_bedrock_mlay_dz = {1.0e1, 1.0e1};
  double eps_dz = 1.0e-3;
  
  int nsoil_lay = ref_soil_mlay_dz.size();
  int nbedrock_lay = ref_bedrock_mlay_dz.size();

  std::cout << "Extruding: " << mesh_in << " and writing to: " << mesh_out << std::endl;

  std::vector<int> soil_type;
  std::vector<int> bedrock_type;
  std::vector<double> depths;

  Mesh2D m = readMesh2D_text(mesh_in, soil_type, bedrock_type, depths);

  int nsnodes = m.coords.size();
  std::vector<double> dzs(nsnodes, 0.0);
  std::vector<double> rem_soil = depths;
  
  Mesh3D m3(&m, nsoil_lay + nbedrock_lay);
  
  for (int ilay = 0; ilay < nsoil_lay; ilay++) {
    for (int inode = 0; inode < nsnodes; inode++) {
      if (rem_soil[inode] < eps_dz) {
        dzs[inode] = 0.0;
        continue;
      }
      if (rem_soil[inode] > 2*ref_soil_mlay_dz[ilay])
        dzs[inode] = ref_soil_mlay_dz[ilay];
      else if (rem_soil[inode] > ref_soil_mlay_dz[ilay] + eps_dz)
        dzs[inode] = 0.5*rem_soil[inode];
      else if (rem_soil[inode] < ref_soil_mlay_dz[ilay] - eps_dz)
        dzs[inode] = rem_soil[inode];
      else
        dzs[inode] = ref_soil_mlay_dz[ilay];
      
      rem_soil[inode] -= dzs[inode];
    }
    m3.extrude(dzs, soil_type);
  }
  
  for (int ilay = 0; ilay < nbedrock_lay; ilay++)
    m3.extrude(ref_bedrock_mlay_dz[ilay], bedrock_type);
  
  m3.finish_sets();

  std::cout << "NNodes on the surf = " << m.coords.size() << std::endl;
  std::cout << "Ncells on the surf = " << m.cell2node.size() << std::endl;
  std::cout << "NNodes on 3D = " << m3.coords.size() << std::endl;
  std::cout << "Ncells on 3D = " << m3.cell2face.size() << std::endl;

  writeMesh3D_exodus(m3, mesh_out);
  return 0;
}
