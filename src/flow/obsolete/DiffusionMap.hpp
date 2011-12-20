#ifndef __DIFFUSIONMAP_H__
#define __DIFFUSIONMAP_H__

#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"

class DiffusionMap
{
public:
      
  DiffusionMap(const Epetra_Map &cell_map, const Epetra_Map &face_map);
  ~DiffusionMap();
  
  const Epetra_Map& Map() const { return *dof_map_; }
  const Epetra_Map& CellMap() const { return cell_map_; }
  const Epetra_Map& FaceMap() const { return face_map_; }
      
  Epetra_Vector* CreateCellView(const Epetra_Vector&) const;
  Epetra_Vector* CreateFaceView(const Epetra_Vector&) const;
  
  Epetra_MultiVector* CreateCellView(const Epetra_MultiVector&) const;
  Epetra_MultiVector* CreateFaceView(const Epetra_MultiVector&) const;
  
private:
    
  Epetra_Map cell_map_;
  Epetra_Map face_map_;
  Epetra_Map *dof_map_;
  
};

#endif
