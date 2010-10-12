#ifndef _MESH_HH_
#define _MESH_HH_

#include "Element_category.hh"
#include "Entity_map.hh"
#include "Data_structures.hh"

#include <Teuchos_RCP.hpp>
#include <Epetra_MpiComm.h>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/fem/FieldTraits.hpp>

#include <map>
#include <memory>

namespace STK_mesh
{

struct Mesh_view;

class Mesh
{

private:

    Epetra_MpiComm communicator_;
    Teuchos::RCP<Entity_map> entity_map_;
    
    int space_dimension_;
    bool consistent_;
    
    std::auto_ptr<stk::mesh::MetaData> meta_data_;
    std::auto_ptr<stk::mesh::BulkData> bulk_data_;
    
    stk::mesh::Selector universal_selector_;
    stk::mesh::Selector owned_selector_;
    stk::mesh::Selector ghost_selector_;
    stk::mesh::Selector used_selector_;
    
    Vector_field_type &coordinate_field_;
    
    const stk::mesh::Selector& selector_ (Element_Category category) const;

    Id_map part_to_set_;
    Id_map set_to_part_;
    
    void update_ ();
    void notify_views_ () {  }
    
    // Internal Validators
    bool element_type_ok_ () const;
    bool dimension_ok_ () const;
    
    // Disable copy and assignment.
    Mesh (const Mesh& rhs);
    Mesh& operator=(const Mesh& rhs);
    
public:
    
    
    // Structors
    // ---------
    
    Mesh (int space_dimension, 
          const Epetra_MpiComm& communicator, 
          Entity_map* entity_map, 
          stk::mesh::MetaData *meta_data, 
          stk::mesh::BulkData *bulk_data,
          const Id_map &part_to_set,
          const Id_map &set_to_part,
          Vector_field_type& coordinate_field);
    
    virtual ~Mesh () { }
    
    
    // Accessors
    // ---------
    
    int space_dimension () const { return space_dimension_;  }
    
    const stk::mesh::MetaData& meta_data    () const { return *meta_data_; }
    const stk::mesh::BulkData& build_data   () const { return *bulk_data_; }
    const Entity_map&          entity_map   () const { return *entity_map_; }
    const Epetra_MpiComm&      communicator () const { return communicator_; }
    bool                       consistent   () const { return consistent_; }
    unsigned int               rank_id      () const { return communicator_.MyPID (); }
    
    unsigned int count_entities (stk::mesh::EntityRank rank, Element_Category category) const;
    
    void get_entities (stk::mesh::EntityRank, Element_Category category, Entity_vector& entities) const;
    
    void element_to_faces (stk::mesh::EntityId element, Entity_Ids& ids) const;
    void element_to_nodes (stk::mesh::EntityId element, Entity_Ids& ids) const;
    void face_to_nodes    (stk::mesh::EntityId element, Entity_Ids& ids) const;
    
    
    double const * coordinates (stk::mesh::EntityId node) const;
    double const * coordinates (stk::mesh::Entity* node)  const;
    
    stk::mesh::Entity* id_to_entity (stk::mesh::EntityRank rank, 
                                     stk::mesh::EntityId id,
                                     Element_Category category) const;


    // Sets
    // ----

    unsigned int num_sets (stk::mesh::EntityRank) const;
    
    bool valid_id (stk::mesh::EntityRank) const;

    template <typename T>
    void get_set_ids (stk::mesh::EntityRank rank, T begin, T end) const;

    template <typename T>
    void get_set_element_ids (stk::mesh::EntityRank, unsigned int set_id, T begin, T end);

    template <typename T>
    void get_set_element_ids (stk::mesh::EntityRank, const char* name, T begin, T end);

    
    

    // Manipulators
    // ------------

    void modify_bulk_data () { bulk_data_->modification_begin ();              consistent_ = false; }
    void freeze_bulk_data () { bulk_data_->modification_end ();   update_ ();  consistent_ = true; }
    void rebalance_mesh (const Entity_map& entity_map);

    void add_view (Mesh_view* view) { }



    // Static information
    // ------------------

    static stk::mesh::EntityRank get_element_type (int space_dimension);
    static stk::mesh::EntityRank get_face_type    (int space_dimension);


    // Validators
    // ----------

    static bool valid_dimension (int space_dimension);
    static bool valid_rank (stk::mesh::EntityRank);

};

}

#endif
