#ifndef _MESH_MAPS_BASE_H_
#define _MESH_MAPS_BASE_H_

#include <Epetra_Map.h>
#include <Epetra_MpiComm.h>
#include <Epetra_SerialComm.h>

#include "../mesh_data/Entity_kind.hh"

#include <memory>
#include <vector>


enum Element_Category
{
    OWNED = 1,
    GHOST = 2,
    USED  = 3
};


class Mesh_maps_base
{

public:
  
  // Local id interfaces.
  // --------------------
  
  virtual
  void cell_to_faces (unsigned int cell, 
		      std::vector<unsigned int>::iterator begin, 
		      std::vector<unsigned int>::iterator end) {};
  virtual
  void cell_to_faces (unsigned int cell, 
		      unsigned int* begin, unsigned int *end) {};


  virtual
  void cell_to_face_dirs (unsigned int cell, 
			  std::vector<int>::iterator begin, 
			  std::vector<int>::iterator end) {};
  virtual
  void cell_to_face_dirs (unsigned int cell, 
			  int * begin, int * end) {};
  
    

  virtual
  void cell_to_nodes (unsigned int cell, 
		      std::vector<unsigned int>::iterator begin, 
		      std::vector<unsigned int>::iterator end) {};
  virtual
  void cell_to_nodes (unsigned int cell, 
		      unsigned int * begin, unsigned int * end) {};




  virtual
  void face_to_nodes (unsigned int face, 
		      std::vector<unsigned int>::iterator begin, 
		      std::vector<unsigned int>::iterator end) {};
  virtual
  void face_to_nodes (unsigned int face, 
		      unsigned int * begin, unsigned int * end) {};



  virtual
  void node_to_coordinates (unsigned int node, 
			    std::vector<double>::iterator begin, 
			    std::vector<double>::iterator end) {};
  virtual
  void node_to_coordinates (unsigned int node, 
			    double * begin, 
			    double * end) {};

  virtual
  void face_to_coordinates (unsigned int face, 
			    std::vector<double>::iterator begin, 
			    std::vector<double>::iterator end) {};
  virtual
  void face_to_coordinates (unsigned int face, 
			    double * begin, 
			    double * end) {};
   
  virtual 
  void cell_to_coordinates (unsigned int cell, 
			    std::vector<double>::iterator begin,
			    std::vector<double>::iterator end) {};
  virtual 
  void cell_to_coordinates (unsigned int cell, 
			    double * begin,
			    double * end) {};

  virtual
  inline const Epetra_Map& cell_map (bool include_ghost) const {};

  virtual
  inline const Epetra_Map& face_map (bool include_ghost) const {}; 

  virtual
  inline const Epetra_Map& node_map (bool include_ghost) const {};
 
  virtual
  unsigned int count_entities (Mesh_data::Entity_kind kind,
			       Element_Category category) const {};

  virtual
  unsigned int num_sets(Mesh_data::Entity_kind kind) const {};

  virtual
  unsigned int get_set_size (unsigned int set_id, 
			     Mesh_data::Entity_kind kind,
			     Element_Category category) const {};

  // virtual
  // unsigned int get_set_size (const char* name,
  // 			     Mesh_data::Entity_kind kind,
  // 			     Element_Category category) const {};
  
  // Id numbers
  virtual
  void get_set_ids (Mesh_data::Entity_kind kind, 
		    std::vector<unsigned int>::iterator begin, 
		    std::vector<unsigned int>::iterator end) const {};
  virtual
  void get_set_ids (Mesh_data::Entity_kind kind, 
		    unsigned int * begin, 
		    unsigned int * end) const {};

  virtual
  bool valid_set_id (unsigned int id, Mesh_data::Entity_kind kind) const {};

  virtual
  void get_set (unsigned int set_id, Mesh_data::Entity_kind kind, 
		Element_Category category, 
		std::vector<unsigned int>::iterator begin, 
		std::vector<unsigned int>::iterator end) const {};
  virtual
  void get_set (unsigned int set_id, Mesh_data::Entity_kind kind, 
		Element_Category category, 
		unsigned int * begin, 
		unsigned int * end) const {};
  
  // virtual
  // void get_set (const char* name, Mesh_data::Entity_kind kind, 
  // 		Element_Category category,
  // 		std::vector<unsigned int>::iterator begin, 
  // 		std::vector<unsigned int>::iterator end) const {};
  
};

#endif /* _MESH_MAPS_BASE_H_ */
