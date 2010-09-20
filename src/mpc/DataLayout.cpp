#include "Epetra_Map.h"
#include "Teuchos_RCP.hpp"
#include "DataLayout.hpp"



DataLayout::DataLayout (Teuchos::RCP<const Epetra_Map> vertex_map_,
			Teuchos::RCP<const Epetra_Map> vertex_overlap_map_,
			Teuchos::RCP<const Epetra_Map> face_map_,
			Teuchos::RCP<const Epetra_Map> face_overlap_map_,
			Teuchos::RCP<const Epetra_Map> element_map_,
			Teuchos::RCP<const Epetra_Map> element_overlap_map_):
  vertex_map(vertex_map_),
  vertex_overlap_map(vertex_overlap_map_),
  face_map(face_map_),
  face_overlap_map(face_overlap_map_),
  element_map(element_map_),
  element_overlap_map(element_overlap_map_)
{
  // create exporter objects
  
  if ( (!is_null(vertex_map)) && (!is_null(vertex_overlap_map))) 
    vertex_exporter = Teuchos::rcp(new Epetra_Export(*vertex_overlap_map, *vertex_map));
  if ( (!is_null(face_map)) && (!is_null(face_overlap_map)))
    face_exporter = Teuchos::rcp(new Epetra_Export(*face_overlap_map, *face_map));
  if ( (!is_null(element_map)) && (!is_null(element_overlap_map)))  
    element_exporter = Teuchos::rcp(new Epetra_Export(*element_overlap_map, *element_map));
  
}
