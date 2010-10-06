#include "Darcy_problem.hpp"

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
  double *f_cell_own; // = &f[0]; // address of initial ghost cell dof
  double *f_cell_use = f_cell_own; // assume no ghost cells
  
  // Face segment of the output Epetra F vector.
  double *f_face_own; // = &f[ncell_own]; // address of initial ghost face dof
  double *f_face_use = new double[nface_use];
  
  // Cell segment of the input Epetra X vector -- pressures.
  double *p_cell_own; // = &x[0]; address of initial cell dof
  double *p_cell_use = p_cell_own;  // assume no ghost cells
  
  // Face segment of the input Epetra X vector -- pressures.
  double *p_face_own; // = &x[ncell_own]; address of initial face dof
  double *p_face_use = new double[nface_use];
  for (int j = 0; j < nface_own; ++j) p_face_use[j] = p_face_own[j];
  
  // Apply initial BC fixups to P_FACE_USE.
  // CODE NEEDED HERE
  
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
  // CODE NEEDED HERE
  
  // Copy owned part of F_FACE_USE into face segment of F.
  for (int j = 0; j < nface_own; ++j) f_face_own[j] = f_face_use[j];
  
  delete [] p_face_use;
  delete [] f_face_use;
}



void Darcy_problem::initialize() 
{

  // here we intialize things such as the structure of the CRS matrix


  // inititialize the CRS matrix...

  const Epetra_Map face_map = FS->get_mesh_maps()->face_map(false);
  const Epetra_Map face_map_ovl = FS->get_mesh_maps()->face_map(true);
  const Epetra_Map cell_map = FS->get_mesh_maps()->cell_map(false);

  // first we construct an Epetra_FECrsGraph
  Epetra_FECrsGraph PrecFECrsGraph(View, face_map, 11);
  
  std::vector<unsigned int> local_faces (6);
  for ( int icellLID = cell_map.MinLID(); icellLID < cell_map.MaxLID(); icellLID++)  
    {
      // get an array of local face indices for the current cell
      FS->get_mesh_maps()->cell_to_faces(icellLID, local_faces.begin(), local_faces.end());
      
      // figure out the columns (ie. the global face indices) for all the 
      // faces in the current cell
      
      int cols[6];
      
      for (int icol=0; icol<6; icol++) 
	{ 
	  cols[icol] = face_map_ovl.GID( local_faces[icol] );
	}
      
      // each face in the cell must be connected to the other faces
      // in the cell
      for (int icol = 0; icol<6; icol++) 
	{ 
	  PrecFECrsGraph.InsertGlobalIndices(1, &cols[icol], 6, cols );
	}
    }
	  
	  
  PrecMat = Teuchos::rcp(new Epetra_CrsMatrix(Copy, PrecFECrsGraph) );
}


Darcy_problem::Darcy_problem(Teuchos::RCP<Flow_State> FS_ ): 
  FS(FS_)
{

  initialize();

}
