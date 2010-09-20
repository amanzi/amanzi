#ifndef __DataLayout_hpp__
#define __DataLayout_hpp__

#include "Teuchos_RCP.hpp"
#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_Export.h"

class DataLayout {

public:
  DataLayout(Teuchos::RCP<const Epetra_Map> vertex_map,
	     Teuchos::RCP<const Epetra_Map> vertex_overlap_map,
	     Teuchos::RCP<const Epetra_Map> face_map,
	     Teuchos::RCP<const Epetra_Map> face_overlap_map,
	     Teuchos::RCP<const Epetra_Map> element_map,
	     Teuchos::RCP<const Epetra_Map> element_overlap_map);

  ~DataLayout() {};

  // access methods for Epetra_Map objects
  Teuchos::RCP<const Epetra_Map> get_vertex_map () const { return vertex_map; };
  Teuchos::RCP<const Epetra_Map> get_vertex_overlap_map() const { return vertex_overlap_map; } 
  
  Teuchos::RCP<const Epetra_Map> get_face_map () const { return face_map; };
  Teuchos::RCP<const Epetra_Map> get_face_overlap_map () const { return face_overlap_map; };   

  Teuchos::RCP<const Epetra_Map> get_element_map () const { return element_map; }; 
  Teuchos::RCP<const Epetra_Map> get_element_overlap_map ()  const{ return element_overlap_map; };

  // we will also need access methods for the Epetra_CrsGraph objects

  
  // access methods for Epetra_Export objects
  Teuchos::RCP<const Epetra_Export> get_vertex_exporter () const { return vertex_exporter; };   
  Teuchos::RCP<const Epetra_Export> get_face_exporter () const { return face_exporter; };   
  Teuchos::RCP<const Epetra_Export> get_element_exporter () const { return element_exporter; };   
  
  

private:
  // Epetra_Map objects
  Teuchos::RCP<const Epetra_Map> vertex_map; 
  Teuchos::RCP<const Epetra_Map> vertex_overlap_map; 
  
  Teuchos::RCP<const Epetra_Map> face_map; 
  Teuchos::RCP<const Epetra_Map> face_overlap_map;   

  Teuchos::RCP<const Epetra_Map> element_map; 
  Teuchos::RCP<const Epetra_Map> element_overlap_map; 

  // we will also need Epetra_CrsGraph objects
  // ...


  // Epetra_Eporter
  Teuchos::RCP<Epetra_Export> vertex_exporter;
  Teuchos::RCP<Epetra_Export> face_exporter;
  Teuchos::RCP<Epetra_Export> element_exporter;
  
};


#endif
