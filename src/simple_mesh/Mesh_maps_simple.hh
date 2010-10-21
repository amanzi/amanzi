#ifndef _MESH_MAPS_SIMPLE_H_
#define _MESH_MAPS_SIMPLE_H_

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



class Mesh_maps_simple
{

public:
  
  Mesh_maps_simple (double x0, double y0, double z0,
		    double x1, double y1, double z1,
		    int nx, int ny, int nz, Epetra_Comm *communicator);
  
  void update ();
  
  // Local id interfaces.
  // --------------------
  
  template <typename IT>
  void cell_to_faces (unsigned int cell, IT begin, IT end) const;

  template <typename IT>
  void cell_to_nodes (unsigned int cell, IT begin, IT end) const;
  
  template <typename IT>
  void face_to_nodes (unsigned int face, IT begin, IT end) const;
  
  template <typename IT>
  void node_to_coordinates (unsigned int node, IT begin, IT end) const;
  
  template <typename IT>
  void face_to_coordinates (unsigned int face, IT begin, IT end) const;
  
  template <typename IT>
  void cell_to_coordinates (unsigned int cell, IT begin, IT end) const;
  
  inline const Epetra_Map& cell_map (bool include_ghost) const;
  inline const Epetra_Map& face_map (bool include_ghost) const;
  inline const Epetra_Map& node_map (bool include_ghost) const;
  
  unsigned int count_entities (Mesh_data::Entity_kind kind, Element_Category category) const;

  unsigned int num_sets(Mesh_data::Entity_kind kind) const;
  unsigned int get_set_size (unsigned int set_id, 
			     Mesh_data::Entity_kind kind,
			     Element_Category category) const;

  unsigned int get_set_size (const char* name,
			     Mesh_data::Entity_kind kind,
			     Element_Category category) const;
  

  // Id numbers
  template <typename IT>
  void get_set_ids (Mesh_data::Entity_kind kind, IT begin, IT end) const;
 
  bool valid_set_id (unsigned int id, Mesh_data::Entity_kind kind) const;

  template <typename IT>
  void get_set (unsigned int set_id, Mesh_data::Entity_kind kind, 
		Element_Category category, IT begin, IT end) const;
  

  template <typename IT>
  void get_set (const char* name, Mesh_data::Entity_kind kind, Element_Category category,
		IT begin, IT end) const;
  
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




template <typename IT>
void Mesh_maps_simple::cell_to_faces (unsigned int cell, IT destination_begin, IT destination_end) const
{
  // ASSERT ((unsigned int) (destination_end - destination_begin) == 6);
    const unsigned int index = 6*cell;
    std::vector<unsigned int>::const_iterator begin = cell_to_face_.begin () + index;
    std::vector<unsigned int>::const_iterator end = begin + 6;
    std::copy (begin, end, destination_begin);
}


template <typename IT>
void Mesh_maps_simple::cell_to_nodes (unsigned int cell, IT destination_begin, IT destination_end) const
{
  // ASSERT ((unsigned int) (destination_end - destination_begin) == 8);
    const unsigned int index = 8*cell;
    std::vector<unsigned int>::const_iterator begin = cell_to_node_.begin () + index;
    std::vector<unsigned int>::const_iterator end   = begin + 8;
    std::copy (begin, end, destination_begin);
}


template <typename IT>
void Mesh_maps_simple::face_to_nodes (unsigned int face, IT destination_begin, IT destination_end) const
{
  // ASSERT ((unsigned int) (destination_end - destination_begin) == 4);
    const unsigned int index = 4*face;
    std::vector<unsigned int>::const_iterator begin = face_to_node_.begin () + index;
    std::vector<unsigned int>::const_iterator end   = begin + 4;
    std::copy (begin, end, destination_begin);
}


// Cooordinate Getters
// -------------------

template <typename IT>
void Mesh_maps_simple::node_to_coordinates (unsigned int local_node_id, IT destination_begin, IT destination_end) const
{
  //  ASSERT ((unsigned int) (end-begin) == 3);
  const unsigned int index = 3*local_node_id;
  std::vector<double>::const_iterator begin = coordinates_.begin() + index;
  std::vector<double>::const_iterator end   = begin + 3;
  std::copy (begin, end, destination_begin);
}


template <typename IT>
void Mesh_maps_simple::face_to_coordinates (unsigned int local_face_id, IT begin, IT end) const
{
  // ASSERT ((unsigned int) (end-begin) == 12);

    unsigned int node_indices [4];
    face_to_nodes (local_face_id, node_indices, node_indices+4);
    for (int i = 0; i < 4; ++i)
    {
        node_to_coordinates (node_indices [i], begin, begin+3);
        begin+=3;
    }


}

template <typename IT>
void Mesh_maps_simple::cell_to_coordinates (unsigned int local_cell_id, IT begin, IT end) const
{
  // ASSERT ((unsigned int) (end-begin) == 24);

    unsigned int node_indices [8];
    cell_to_nodes (local_cell_id, node_indices, node_indices+8);
    for (int i = 0; i < 8; ++i)
    {
        node_to_coordinates (node_indices [i], begin, begin+3);
        begin+=3;
    }

}



// Set getters
// -----------

template <typename IT>
void Mesh_maps_simple::get_set_ids (Mesh_data::Entity_kind kind, IT begin, IT end) const
{
  std::vector<unsigned int> ids(6);
  
  switch (kind) {
  case Mesh_data::FACE: 
    for (int i=0; i<6; i++) ids[i]=i;
    
    std::copy (ids.begin (), ids.end (), begin);
    break;
  default:
    // we do not have anything for CELL and NODE, yet
    throw std::exception();
  }
}

template <typename IT>
void Mesh_maps_simple::get_set (unsigned int set_id, Mesh_data::Entity_kind kind, 
				Element_Category category, IT begin, IT end) const
{
  // we ignore category since this is a serial implementation
  
  switch (kind) {
  case Mesh_data::FACE:
    std::copy(side_sets_[set_id].begin(), side_sets_[set_id].end(), begin) ;
    break;
  default:
    // we do not have anything for CELL and NODE, yet
    throw std::exception();
  }
}



#endif /* _MESH_MAPS_H_ */
