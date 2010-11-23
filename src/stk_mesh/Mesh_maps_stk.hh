#ifndef _MESH_MAPS_STK_H_
#define _MESH_MAPS_STK_H_

#include "Mesh_maps_base.hh"
#include "Mesh.hh"
#include "Entity_map.hh"
#include "Data_structures.hh"

#include "dbc.hh"

#include <Epetra_Map.h>
#include <Epetra_MpiComm.h>

#include <memory>

namespace STK_mesh
{

  class Mesh_maps_stk : public Mesh_maps_base
  {

  private:

    Mesh_p mesh_;
    const Entity_map& entity_map_;
    Epetra_MpiComm communicator_;

    void update_internals_ ();
    void clear_internals_ ();

    // Maps, Accessors and setters.
    // ----------------------------
    std::auto_ptr<Epetra_Map> maps_ [6];
    const Epetra_Map& map_  (Mesh_data::Entity_kind kind, bool include_ghost) const;
    unsigned int map_index_ (Mesh_data::Entity_kind kind, bool include_ghost) const;
    void assign_map_        (Mesh_data::Entity_kind kind, bool include_ghost, Epetra_Map *map);

    void build_maps_ ();
    void build_tables_ ();

    stk::mesh::EntityRank kind_to_rank_ (Mesh_data::Entity_kind kind) const {
      return entity_map_.kind_to_rank (kind);
    }

    bool valid_entity_kind_ (Mesh_data::Entity_kind kind) const;

    // Local-id tables of entities
    std::vector<unsigned int> cell_to_face_;
    std::vector<unsigned int> cell_to_node_;
    std::vector<unsigned int> face_to_node_;


    // Global to local index maps and associated bookkeeping.
    Index_map global_to_local_maps_ [3];
    static unsigned int num_kinds_;
    static Mesh_data::Entity_kind kinds_ [3];
    unsigned int kind_to_index_ (Mesh_data::Entity_kind type) const;
    const Index_map& kind_to_map_ (Mesh_data::Entity_kind kind) const;
    Mesh_data::Entity_kind index_to_kind_ (unsigned int index) const;

    unsigned int global_to_local_ (unsigned int global_id, Mesh_data::Entity_kind kind) const;

    // Builds the global->local maps.
    template <typename F, typename D, typename M>
    void add_global_ids_ (F from, F to, D destination, M& inverse);

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

    template <typename IT>
    void get_set_ids (Mesh_data::Entity_kind kind, IT begin, IT end) const;

    template <typename IT>
    void get_set (unsigned int set_id, Mesh_data::Entity_kind kind, Element_Category category,
                  IT begin, IT end) const;

    template <typename IT>
    void get_set (const char* name, Mesh_data::Entity_kind kind, Element_Category category,
                  IT begin, IT end) const;

  public:

    explicit Mesh_maps_stk (Mesh_p mesh);

    void update ();

    // Local id interfaces
    // --------------------


    void cell_to_faces (unsigned int cell, 
                        std::vector<unsigned int>::iterator begin, 
                        std::vector<unsigned int>::iterator end);

    void cell_to_faces (unsigned int cell, 
                        unsigned int* begin, unsigned int *end);
    

    void cell_to_face_dirs (unsigned int cell, 
                            std::vector<int>::iterator begin, 
                            std::vector<int>::iterator end);

    void cell_to_face_dirs (unsigned int cell, 
                            int * begin, int * end);
  

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

    unsigned int count_entities (Mesh_data::Entity_kind kind, Element_Category category) const;

    // Entity Sets (cell, side, node)
    // ------------------------------

    // Number and sizes

    unsigned int num_sets () const;
    unsigned int num_sets (Mesh_data::Entity_kind kind) const;

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

    // communicator access
    const Epetra_Comm* get_comm() { return &communicator_; };

    // this should be used with extreme caution:
    // modify coordinates  
    void set_coordinate(unsigned int local_node_id, 
                        double* source_begin, double* source_end);
  };

  // -------------------------
  // Template & inline members
  // ------------------------

  // Inlined

  const Epetra_Map& Mesh_maps_stk::cell_map (bool include_ghost) const
  {
    return map_ (Mesh_data::CELL, include_ghost);
  }

  const Epetra_Map& Mesh_maps_stk::face_map (bool include_ghost) const
  {
    return map_ (Mesh_data::FACE, include_ghost);
  }

  const Epetra_Map& Mesh_maps_stk::node_map (bool include_ghost) const
  {
    return map_ (Mesh_data::NODE, include_ghost);
  }


  // Internal template functions
  // ---------------------------


  template <typename F, typename D, typename M>
  void Mesh_maps_stk::add_global_ids_ (F from, F to, D destination, M& global_to_local_map)
  {
    int local_id = 0;
    for (F it = from; it != to;  ++it)
      {
        const unsigned int global_id = (*it)->identifier ();
        *destination = global_id;
        global_to_local_map.insert (std::make_pair (global_id, local_id));

        destination++;
        local_id++;
      }
  }

}



#endif /* _MESH_MAPS_H_ */
