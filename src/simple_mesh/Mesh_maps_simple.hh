#ifndef _MESH_MAPS_SIMPLE_H_
#define _MESH_MAPS_SIMPLE_H_

#include <Epetra_Map.h>
#include <Epetra_MpiComm.h>
#include <Epetra_SerialComm.h>

#include "../mesh_data/Entity_kind.hh"

#include <memory>
#include <vector>

#include "Teuchos_ParameterList.hpp"
#include "Mesh_maps_base.hh"


class Mesh_maps_simple : public virtual Mesh_maps_base
{

public:
  
  Mesh_maps_simple (double x0, double y0, double z0,
		    double x1, double y1, double z1,
		    int nx, int ny, int nz, Epetra_Comm *communicator);
  
  Mesh_maps_simple ( Teuchos::ParameterList &parameter_list,
		     Epetra_Comm *communicator );


  void update ();
  
  // Local id interfaces.
  // --------------------
  
  void cell_to_faces (unsigned int cell, 
		      std::vector<unsigned int>::iterator begin, 
		      std::vector<unsigned int>::iterator end);
  void cell_to_faces (unsigned int cell, 
		      unsigned int * begin, unsigned int * end);


  void cell_to_nodes (unsigned int cell, 
		      std::vector<unsigned int>::iterator begin, 
		      std::vector<unsigned int>::iterator end);
  void cell_to_nodes (unsigned int cell, 
		      unsigned int * begin, unsigned int * end);


  void face_to_nodes (unsigned int face,
		      std::vector<unsigned int>::iterator begin, 
		      std::vector<unsigned int>::iterator end);
  void face_to_nodes (unsigned int face,
		      unsigned int * begin, unsigned int * end);
  

  void node_to_coordinates (unsigned int node, 
			    std::vector<double>::iterator begin,
			    std::vector<double>::iterator end);
  void node_to_coordinates (unsigned int node, 
			    double * begin,
			    double * end);
  
  void face_to_coordinates (unsigned int face,
			    std::vector<double>::iterator begin, 
			    std::vector<double>::iterator end);
  void face_to_coordinates (unsigned int face,
			    double * begin, 
			    double * end);
  
  void cell_to_coordinates (unsigned int cell, 
			    std::vector<double>::iterator begin,
			    std::vector<double>::iterator end);
  void cell_to_coordinates (unsigned int cell, 
			    double * begin,
			    double * end);
  
  inline const Epetra_Map& cell_map (bool include_ghost) const;

  inline const Epetra_Map& face_map (bool include_ghost) const;

  inline const Epetra_Map& node_map (bool include_ghost) const;
  
  unsigned int count_entities (Mesh_data::Entity_kind kind,
			       Element_Category category) const;

  unsigned int num_sets(Mesh_data::Entity_kind kind) const;
  
  unsigned int get_set_size (unsigned int set_id, 
			     Mesh_data::Entity_kind kind,
			     Element_Category category) const;

  unsigned int get_set_size (const char* name,
			     Mesh_data::Entity_kind kind,
			     Element_Category category) const;
  

  // Id numbers
  void get_set_ids (Mesh_data::Entity_kind kind, 
		    std::vector<unsigned int>::iterator begin, 
		    std::vector<unsigned int>::iterator end) const;
  void get_set_ids (Mesh_data::Entity_kind kind, 
		    unsigned int * begin, 
		    unsigned int * end) const;
 
  bool valid_set_id (unsigned int id, Mesh_data::Entity_kind kind) const;

  void get_set (unsigned int set_id, Mesh_data::Entity_kind kind, 
		Element_Category category,
		std::vector<unsigned int>::iterator begin,
		std::vector<unsigned int>::iterator end) const;
  void get_set (unsigned int set_id, Mesh_data::Entity_kind kind, 
		Element_Category category,
		unsigned int * begin,
		unsigned int * end) const;
  

private:
  void update_internals_();
  void clear_internals_();
  void build_maps_();
 

  Epetra_Map *cell_map_, *face_map_, *node_map_;

  
  std::vector<double> coordinates_;

  inline unsigned int node_index_(int i, int j, int k);
  inline unsigned int xyface_index_(int i, int j, int k);
  inline unsigned int yzface_index_(int i, int j, int k);
  inline unsigned int xzface_index_(int i, int j, int k);
  inline unsigned int cell_index_(int i, int j, int k);

  int nx_, ny_, nz_;  // number of cells in the three coordinate directions
  double x0_, x1_, y0_, y1_, z0_, z1_;  // coordinates of lower left front and upper right back of brick
 
  int num_cells_;
  int num_nodes_;
  int num_faces_;

  // Local-id tables of entities
  std::vector<unsigned int> cell_to_face_;
  std::vector<unsigned int> cell_to_node_;
  std::vector<unsigned int> face_to_node_;
  
  std::vector<std::vector<unsigned int> > side_sets_;
  std::vector<std::vector<unsigned int> > element_blocks_;
  
  Epetra_Comm  *communicator_;

};

// -------------------------
// Template & inline members
// ------------------------

// Inlined

const Epetra_Map& Mesh_maps_simple::cell_map (bool include_ghost) const
{
    return *cell_map_;
}

const Epetra_Map& Mesh_maps_simple::face_map (bool include_ghost) const
{
    return *face_map_;
}

const Epetra_Map& Mesh_maps_simple::node_map (bool include_ghost) const
{
    return *node_map_;
}


unsigned int Mesh_maps_simple::node_index_(int i, int j, int k)
{
  return i + j*(nx_+1) + k*(nx_+1)*(ny_+1);
}

unsigned int Mesh_maps_simple::cell_index_(int i, int j, int k)
{
  return i + j*nx_ + k*nx_*ny_;
}

unsigned int Mesh_maps_simple::xyface_index_(int i, int j, int k)
{
  return i + j*nx_ + k*nx_*ny_;
}

unsigned int Mesh_maps_simple::xzface_index_(int i, int j, int k)
{
  return i + j*nx_ + k*nx_*(ny_+1) +  xyface_index_(0,0,nz_+1);
}

unsigned int Mesh_maps_simple::yzface_index_(int i, int j, int k)
{
  return i + j*(nx_+1) + k*(nx_+1)*ny_ + xzface_index_(0,0,nz_);
}





#endif /* _MESH_MAPS_H_ */
