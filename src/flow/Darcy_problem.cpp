#include "Darcy_problem.hpp"
#include "Epetra_FECrsGraph.h"

void Darcy_problem::ComputeF(const Epetra_Vector & x, Epetra_Vector & f)
{
  Teuchos::RCP<const STK_mesh::Mesh_maps> mesh = FS->get_mesh_maps();
  const int ncell_use = mesh->count_entities(Mesh_data::CELL, STK_mesh::USED);
  const int ncell_own = mesh->count_entities(Mesh_data::CELL, STK_mesh::OWNED);
  const int nface_use = mesh->count_entities(Mesh_data::FACE, STK_mesh::USED);
  const int nface_own = mesh->count_entities(Mesh_data::FACE, STK_mesh::OWNED);
  
  // The cell and face-based DoF are packed together into the X and F Epetra
  // vectors: cell-based DoF in the first part, followed by the face-based DoF.
  // In addition, only the owned DoF belong to the vectors; for computation we
  // need to expand these to include the ghost values.
  
  // Cell segment of the output Epetra F vector.
  double *f_cell_own;
  f.ExtractView(&f_cell_own); // f_cell_own = &f[0]
  double *f_cell_use = f_cell_own; // assume no ghost cells
  
  // Face segment of the output Epetra F vector.
  double *f_face_own = f_cell_own + ncell_own; // = &f[ncell_own]
  double *f_face_use = new double[nface_use];
  
  // Cell segment of the input Epetra X vector -- pressures.
  double *p_cell_own;
  x.ExtractView(&p_cell_own); // p_cell_own = &x[0]
  double *p_cell_use = p_cell_own;  // assume no ghost cells
  
  // Face segment of the input Epetra X vector -- pressures.
  double *p_face_own = p_cell_own + ncell_own; // = &x[ncell_own]
  double *p_face_use = new double[nface_use];
  for (int j = 0; j < nface_own; ++j) p_face_use[j] = p_face_own[j];
  
  // Apply initial BC fixups to P_FACE_USE.
  FBC_initial(p_face_use);
  
  // Gather the ghost P_FACE_USE values.
  // CODE NEEDED HERE -- epetra import/export?
  
  int face[6];
  double aux1[6], aux2[6];
  
  // Compute the diffusion operator contribution to F.
  for (int j = 0; j < nface_use; ++j) f_face_use[j] = 0.0;
  for (int j = 0; j < ncell_use; ++j) {
     // Get the list of process-local face indices for the cell.
     mesh->cell_to_faces(j, face, face+6);
     // Gather the local face pressures.
     for (int i = 0; i < 6; ++i) aux1[i] = p_face_use[face[i]];
     // Compute the local value of the diffusion operator.
     MD[j].diff_op(K[j], p_cell_use[j], aux1, f_cell_use[j], aux2);
     // Scatter local face result into F_FACE_USE.
     for (int i = 0; i < 6; ++i) f_face_use[face[i]] += aux2[i];
  }
  
  // Parallel assembly: sum ghost F_FACE_USE values to the master value.
  // CODE NEEDED HERE -- epetra import/export?
  
  // Apply final BC fixups to F_FACE_USE.
  FBC_final(f_face_use);
  
  // Copy owned part of F_FACE_USE into face segment of F.
  for (int j = 0; j < nface_own; ++j) f_face_own[j] = f_face_use[j];
  
  delete [] p_face_use;
  delete [] f_face_use;
}


// BC fixups for F computation: initial pass.
void FBC_initial(double p_face[])
{
  //std::vector<flow_bc> bc = fbcs.get_BCs();
}


// BC fixups for F computation: final pass.
void FBC_final(double f_face[])
{
}



void Darcy_problem::initialize()
{

  // initialize the crs matrix

  // first get the mesh maps
  
  Epetra_Map face_map = FS->get_mesh_maps()->face_map(false);
  Epetra_Map face_map_ovl = FS->get_mesh_maps()->face_map(true);

  Epetra_Map cell_map = FS->get_mesh_maps()->cell_map(false);

  // now make an Epetra_FECrsGraph
  Epetra_FECrsGraph PrecMat_Graph(Copy, face_map, 11);

  // loop over the cells in the mesh

  std::vector<int> cface(6);

  for (int icell=cell_map.MinLID(); icell < cell_map.MaxLID(); icell++)
    {
      // get the local face indices of faces connected to the 
      // current cell
      FS->get_mesh_maps()->cell_to_faces(icell,cface.begin(),cface.end());
      
      // translate to global indices in a data structure that
      // Epetra_FECrsGraph.InsertGlobalIndices will understand
      int col[6];
      
      for (int i=0; i<6; i++) 
	{
	  col[i] = face_map.GID(cface[i]);
	}
	
      // for each face in this list of attached faces, enter 
      // its connectivity into the graph

      for (int i=0; i<6; i++)
	{
	  PrecMat_Graph.InsertGlobalIndices(1, &col[i], 6, col); 
	}
    }

  PrecMat_Graph.GlobalAssemble();

  PrecMat = Teuchos::rcp(new Epetra_FECrsMatrix(Copy, PrecMat_Graph));
  
}
