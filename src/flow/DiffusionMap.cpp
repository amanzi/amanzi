#include "DiffusionMap.hpp"
#include "dbc.hh"

DiffusionMap::DiffusionMap(const Epetra_Map &cell_map, const Epetra_Map &face_map)
    : cell_map_(cell_map), face_map_(face_map)
{
  int ncell_tot = cell_map.NumGlobalElements();
  int ndof_tot = ncell_tot + face_map.NumGlobalElements();
  int ncell = cell_map.NumMyElements();
  int ndof = ncell + face_map.NumMyElements();
  int *gids = new int[ndof];
  cell_map.MyGlobalElements(&(gids[0]));
  face_map.MyGlobalElements(&(gids[ncell]));
  for (int i = ncell; i < ndof; ++i) gids[i] += ncell_tot;
  dof_map_ = new Epetra_Map(ndof_tot, ndof, gids, 0, cell_map_.Comm());
  delete [] gids;
}

DiffusionMap::~DiffusionMap()
{
  delete dof_map_;
}


Epetra_Vector* DiffusionMap::CreateCellView(const Epetra_Vector &X) const
{
  ASSERT( X.Map().SameAs(Map()) );
  double *data;
  X.ExtractView(&data);
  return new Epetra_Vector(View, CellMap(), data);
}


Epetra_Vector* DiffusionMap::CreateFaceView(const Epetra_Vector &X) const
{
  ASSERT( X.Map().SameAs(Map()) );
  double *data;
  X.ExtractView(&data);
  int ncell = CellMap().NumMyElements();
  return new Epetra_Vector(View, FaceMap(), data+ncell);
}


Epetra_MultiVector* DiffusionMap::CreateCellView(const Epetra_MultiVector &X) const
{
  ASSERT( X.Map().SameAs(Map()) );
  double **data_ptrs = X.Pointers();
  return new Epetra_MultiVector(View, CellMap(), data_ptrs, X.NumVectors());
}


Epetra_MultiVector* DiffusionMap::CreateFaceView(const Epetra_MultiVector &X) const
{
  ASSERT( X.Map().SameAs(Map()) );
  double **data_ptrs = X.Pointers();
  double **fdat_ptrs = new double*[X.NumVectors()];
  int ncell = CellMap().NumMyElements();
  for (int i = 0; i < X.NumVectors(); ++i) fdat_ptrs[i] = data_ptrs[i] + ncell;
  Epetra_MultiVector *fvec = new Epetra_MultiVector(View, FaceMap(), fdat_ptrs, X.NumVectors());
  delete [] fdat_ptrs;
  return fvec;
}
