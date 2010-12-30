#include <algorithm>
#include <boost/lambda/lambda.hpp>
namespace bl = boost::lambda;

#include "Mesh_maps_stk.hh"
#include "dbc.hh"
#include "Utils.hh"

#include <Teuchos_RCP.hpp>
#include <Epetra_DataAccess.h>
#include <Epetra_CrsMatrix.h>
#include <stk_mesh/fem/EntityRanks.hpp>
#include <Isorropia_EpetraPartitioner.hpp>

static const int ZERO = 0;

namespace STK_mesh
{


unsigned int Mesh_maps_stk::num_kinds_ = 3;
Mesh_data::Entity_kind Mesh_maps_stk::kinds_ [3] = {Mesh_data::NODE, 
                                                    Mesh_data::FACE, 
                                                    Mesh_data::CELL};


Mesh_maps_stk::Mesh_maps_stk (Mesh_p mesh) : mesh_ (mesh),
                                             entity_map_ (mesh->entity_map ()),
                                             communicator_ ((mesh_->communicator ()).Clone())
{
    build_maps_();
}


// -------------------------------------------------------------
// Mesh_maps_stk::build_maps_
// -------------------------------------------------------------
static void
extract_global_ids(const Entity_vector& entities,
                   std::vector<int>& ids)
{
    ids.clear();
    ids.reserve(entities.size());
    for (Entity_vector::const_iterator i = entities.begin(); i != entities.end(); i++) {
        ids.push_back((*i)->identifier());
    }
    std::for_each(ids.begin(), ids.end(), bl::_1 -= 1);
}

void 
Mesh_maps_stk::build_maps_ ()
{
    map_owned_.clear();
    map_used_.clear();

    // For each of Elements, Faces, Nodes:
    for (int entity_kind_index = 0; entity_kind_index < 3; ++entity_kind_index)
    {

        Mesh_data::Entity_kind kind = index_to_kind_ (entity_kind_index);
        stk::mesh::EntityRank rank = entity_map_.kind_to_rank (kind);

        Entity_vector entities;
        std::vector<int> entity_ids;
        Teuchos::RCP<Epetra_Map> map;

        // Get the collection of "owned" entities
        mesh_->get_entities (rank, OWNED, entities);
        extract_global_ids(entities, entity_ids);

        map.reset(new Epetra_Map(-1, entity_ids.size(), &entity_ids[0], ZERO, *communicator_));
        map_owned_.insert(MapSet::value_type(kind, map));

        // Get the collection of "ghost" entities
        Entity_vector ghost_entities;
        mesh_->get_entities (rank, GHOST, ghost_entities);
        std::vector<int> ghost_entity_ids;
        extract_global_ids(ghost_entities, ghost_entity_ids);
        std::copy(ghost_entity_ids.begin(), ghost_entity_ids.end(), 
                  std::back_inserter(entity_ids));

        map.reset(new Epetra_Map(-1, entity_ids.size(), &entity_ids[0], ZERO, *communicator_));
        map_used_.insert(MapSet::value_type(kind, map));
    }
}

// -------------------------------------------------------------
// Mesh_maps_stk::get_map_
// -------------------------------------------------------------
const Epetra_Map&
Mesh_maps_stk::get_map_(const Mesh_data::Entity_kind& kind, const bool& include_ghost) const
{
    MapSet::const_iterator m;
    if (include_ghost) {
        m = map_used_.find(kind);
        ASSERT(m != map_used_.end());
    } else {
        m = map_owned_.find(kind);
        ASSERT(m != map_owned_.end());
    }
    return *(m->second);
}


// Internal validators
// -------------------

bool Mesh_maps_stk::valid_entity_kind_ (const Mesh_data::Entity_kind kind) const
{
    for (Mesh_data::Entity_kind* kind_it = kinds_; kind_it != kinds_+num_kinds_; ++kind_it)
        if (*kind_it == kind) return true;

    return false;
}


// Bookkeeping for the internal relationship maps
// ----------------------------------------------


unsigned int Mesh_maps_stk::kind_to_index_ (const Mesh_data::Entity_kind kind) const
{
    ASSERT (valid_entity_kind_ (kind));

    if (kind == Mesh_data::NODE) return 0;
    if (kind == Mesh_data::FACE) return 1;
    if (kind == Mesh_data::CELL) return 2;
}

Mesh_data::Entity_kind Mesh_maps_stk::index_to_kind_ (const unsigned int index) const
{
    ASSERT (index >= 0 && index < num_kinds_);
    return kinds_ [index];
}

// Bookkeeping for the collection of Epetra_maps
// ---------------------------------------------




// Public Accesor Functions
// ------------------------


unsigned int Mesh_maps_stk::count_entities (Mesh_data::Entity_kind kind, Element_Category category) const
{
    return mesh_->count_entities (kind_to_rank_ (kind), category);
}

bool Mesh_maps_stk::valid_set_id (unsigned int set_id, Mesh_data::Entity_kind kind) const
{
    return mesh_->valid_id (set_id, kind_to_rank_ (kind));
}

unsigned int Mesh_maps_stk::num_sets () const
{
    return mesh_->num_sets ();
}

unsigned int Mesh_maps_stk::num_sets (Mesh_data::Entity_kind kind) const
{
    return mesh_->num_sets (kind_to_rank_ (kind));
}

unsigned int Mesh_maps_stk::get_set_size (unsigned int set_id,
                                          Mesh_data::Entity_kind kind,
                                          Element_Category category) const
{
    ASSERT (valid_set_id (set_id, kind));
    stk::mesh::Part* part = mesh_->get_set (set_id, kind_to_rank_ (kind));
    return mesh_->count_entities (*part, category);
}



// -------------------------------------------------------------
// Mesh_maps_stk::cell_to_faces
// -------------------------------------------------------------
template <typename IT>
void Mesh_maps_stk::cell_to_faces (unsigned int cell, 
                                   IT destination_begin, IT destination_end) const
{
    stk::mesh::EntityId global_cell_id = 
        get_map_(Mesh_data::CELL, true).GID(cell);
    global_cell_id += 1;        // need 1-based for stk::mesh

    Entity_Ids face_ids;
    mesh_->element_to_faces(global_cell_id, face_ids);
    std::for_each(face_ids.begin(), face_ids.end(), bl::_1 -= 1);

    IT i = destination_begin;
    for (Entity_Ids::iterator f = face_ids.begin(); f != face_ids.end(); f++, i++) {
        stk::mesh::EntityId global_face_id(*f);
        if (!get_map_(Mesh_data::FACE, true).MyGID(global_face_id)) {
            std::cerr << "Process " << communicator_->MyPID() << ": "
                      << "cell " << cell << " (" << global_cell_id << ") has bogus face id "
                      << global_face_id << std::endl;
        }
        ASSERT(get_map_(Mesh_data::FACE, true).MyGID(global_face_id));
        stk::mesh::EntityId local_face_id = 
            get_map_(Mesh_data::FACE, true).LID(global_face_id);
        *i = local_face_id;
    }
}

void 
Mesh_maps_stk::cell_to_faces (unsigned int cell, 
                              std::vector<unsigned int>::iterator begin, 
                              std::vector<unsigned int>::iterator end)
{
  cell_to_faces< std::vector<unsigned int>::iterator >(cell, begin, end);
}

void
Mesh_maps_stk::cell_to_faces (unsigned int cell, 
                              unsigned int* begin, unsigned int *end) 
{
  cell_to_faces< unsigned int * >(cell, begin, end);
}

// -------------------------------------------------------------
// Mesh_maps_stk::cell_to_face_dirs
// -------------------------------------------------------------
template <typename IT>
void Mesh_maps_stk::cell_to_face_dirs (unsigned int cell, 
                                   IT destination_begin, IT destination_end) const
{
    stk::mesh::EntityId global_cell_id = 
        get_map_(Mesh_data::CELL, true).GID(cell);
    global_cell_id += 1;        // need 1-based for stk::mesh

    std::vector<int> dirs;
    mesh_->element_to_face_dirs(global_cell_id, dirs);

    std::copy(dirs.begin(), dirs.end(), destination_begin);
}


void
Mesh_maps_stk::cell_to_face_dirs(unsigned int cell, 
                                 std::vector<int>::iterator begin, 
                                 std::vector<int>::iterator end)
{
    cell_to_face_dirs<std::vector<int>::iterator>(cell, begin, end);
}
                                   
void
Mesh_maps_stk::cell_to_face_dirs(unsigned int cell, 
                                 int * begin, int * end)
{
    cell_to_face_dirs<int *>(cell, begin, end);
}  

// -------------------------------------------------------------
// Mesh_maps_stk::cell_to_nodes
// -------------------------------------------------------------
template <typename IT>
void 
Mesh_maps_stk::cell_to_nodes (unsigned int cell, 
                              IT destination_begin, IT destination_end) const
{
    stk::mesh::EntityId global_cell_id = 
        get_map_(Mesh_data::CELL, true).GID(cell);
    global_cell_id += 1;        // need 1-based for stk::mesh

    Entity_Ids node_ids;
    mesh_->element_to_nodes(global_cell_id, node_ids);
    std::for_each(node_ids.begin(), node_ids.end(), bl::_1 -= 1); // 0-based for Epetra_Map

    IT i = destination_begin;
    for (Entity_Ids::iterator n = node_ids.begin(); n != node_ids.end(); n++, i++) {
        stk::mesh::EntityId global_node_id(*n);
        if (!get_map_(Mesh_data::NODE, true).MyGID(global_node_id)) {
            std::cerr << "Process " << communicator_->MyPID() << ": "
                      << "cell " << cell << " (" << global_cell_id << ") has bogus node id "
                      << global_node_id << std::endl;
        }
        ASSERT(get_map_(Mesh_data::NODE, true).MyGID(global_node_id));
        stk::mesh::EntityId local_node_id = 
            get_map_(Mesh_data::NODE, true).LID(global_node_id);
        *i = local_node_id;
    }
}

void 
Mesh_maps_stk::cell_to_nodes (unsigned int cell, 
                              std::vector<unsigned int>::iterator begin, 
                              std::vector<unsigned int>::iterator end) 
{
  cell_to_nodes< std::vector<unsigned int>::iterator > (cell, begin, end);
}

void 
Mesh_maps_stk::cell_to_nodes (unsigned int cell, 
                              unsigned int * begin, unsigned int * end)
{
  cell_to_nodes< unsigned int * > (cell, begin, end);
}


// -------------------------------------------------------------
// Mesh_maps_stk::face_to_nodes
// -------------------------------------------------------------
template <typename IT>
void 
Mesh_maps_stk::face_to_nodes (unsigned int face, 
                                 IT destination_begin, IT destination_end) const
{
    stk::mesh::EntityId global_face_id = 
        get_map_(Mesh_data::FACE, true).GID(face);
    global_face_id += 1;        // need 1-based for stk::mesh

    Entity_Ids node_ids;
    mesh_->face_to_nodes(global_face_id, node_ids);
    std::for_each(node_ids.begin(), node_ids.end(), bl::_1 -= 1);

    IT i = destination_begin;
    for (Entity_Ids::iterator n = node_ids.begin(); n != node_ids.end(); n++, i++) {
        stk::mesh::EntityId global_node_id(*n);
        // ASSERT(get_map_(Mesh_data::NODE, true).MyGID(global_node_id));
        if (!get_map_(Mesh_data::NODE, true).MyGID(global_node_id)) {
            std::cerr << "Process " << communicator_->MyPID() << ": "
                      << "face " << face << " (" << global_face_id << ") has bogus node id "
                      << global_node_id << std::endl;
        }
        stk::mesh::EntityId local_node_id = 
            get_map_(Mesh_data::NODE, true).LID(global_node_id);
        *i = local_node_id;
    }
}

void 
Mesh_maps_stk::face_to_nodes (unsigned int face, 
                    std::vector<unsigned int>::iterator begin, 
                    std::vector<unsigned int>::iterator end) 
{
  face_to_nodes< std::vector<unsigned int>::iterator >(face, begin, end);
}

void 
Mesh_maps_stk::face_to_nodes (unsigned int face, 
                              unsigned int * begin, unsigned int * end) 
{
  face_to_nodes< unsigned int * > (face, begin, end);
}

// Cooordinate Getters
// -------------------


// -------------------------------------------------------------
// Mesh_maps_stk::node_to_coordinates
// -------------------------------------------------------------
template <typename IT>
void Mesh_maps_stk::node_to_coordinates (unsigned int local_node_id, IT begin, IT end) const
{
    // Convert local node to global node Id.
    stk::mesh::EntityId global_node_id = 
        get_map_(Mesh_data::NODE, true).GID(local_node_id);
    global_node_id += 1;        // 1-based for stk::mesh

    const double * coordinates = mesh_->coordinates (global_node_id);
    std::copy (coordinates, coordinates+3, begin);
}

void 
Mesh_maps_stk::node_to_coordinates (unsigned int node, 
                          std::vector<double>::iterator begin, 
                          std::vector<double>::iterator end) 
{
  node_to_coordinates< std::vector<double>::iterator >(node, begin, end);
}

void 
Mesh_maps_stk::node_to_coordinates (unsigned int node, 
                                    double * begin, 
                                    double * end) 
{
  node_to_coordinates< double * >(node, begin, end);
}


// -------------------------------------------------------------
// face_to_coordinates
// -------------------------------------------------------------
template <typename IT>
void Mesh_maps_stk::face_to_coordinates (unsigned int local_face_id, IT begin, IT end) const
{
    stk::mesh::EntityId global_face_id = 
        get_map_(Mesh_data::FACE, true).GID(local_face_id);
    global_face_id += 1;        // need 1-based for stk::mesh

    Entity_Ids node_ids;
    mesh_->face_to_nodes(global_face_id, node_ids);

    IT i = begin;
    const int ndim(mesh_->space_dimension());
    for (Entity_Ids::iterator n = node_ids.begin(); n != node_ids.end(); n++) {
        const double *coord = mesh_->coordinates(*n);
        for (int d = 0; d < ndim; d++, i++) 
            *i = coord[d];
    }
}

void 
Mesh_maps_stk::face_to_coordinates (unsigned int face, 
                                    std::vector<double>::iterator begin, 
                                    std::vector<double>::iterator end) 
{
  face_to_coordinates< std::vector<double>::iterator >(face, begin, end);
}

void 
Mesh_maps_stk::face_to_coordinates (unsigned int face, 
                                    double * begin, 
                                    double * end) 
{
  face_to_coordinates< double * >(face, begin, end);
}
   


// -------------------------------------------------------------
// Mesh_maps_stk::cell_to_coordinates
// -------------------------------------------------------------
template <typename IT>
void Mesh_maps_stk::cell_to_coordinates (unsigned int local_cell_id, IT begin, IT end) const
{
    stk::mesh::EntityId global_cell_id = 
        get_map_(Mesh_data::CELL, true).GID(local_cell_id);
    global_cell_id += 1;        // need 1-based for stk::mesh

    Entity_Ids node_ids;
    mesh_->element_to_nodes(global_cell_id, node_ids);

    IT i = begin;
    const int ndim(mesh_->space_dimension());
    for (Entity_Ids::iterator n = node_ids.begin(); n != node_ids.end(); n++) {
        const double *coord = mesh_->coordinates(*n);
        for (int d = 0; d < ndim; d++, i++) 
            *i = coord[d];
    }
}

void 
Mesh_maps_stk::cell_to_coordinates (unsigned int cell, 
                                    std::vector<double>::iterator begin,
                                    std::vector<double>::iterator end)
{
  cell_to_coordinates< std::vector<double>::iterator > (cell, begin, end);
}

void 
Mesh_maps_stk::cell_to_coordinates (unsigned int cell, 
                              double * begin,
                              double * end) 
{
  cell_to_coordinates< double * > (cell, begin, end);
}



// Set getters
// -----------


// -------------------------------------------------------------
// Mesh_maps_stk::get_set_ids
// -------------------------------------------------------------
template <typename IT>
void 
Mesh_maps_stk::get_set_ids (Mesh_data::Entity_kind kind, IT begin, IT end) const
{
  std::vector<unsigned int> ids;
  mesh_->get_set_ids (kind_to_rank_ (kind), ids);
  ASSERT (ids.size () == num_sets (kind));
  
  IT last = std::copy (ids.begin (), ids.end (), begin);
  
  ASSERT (last == end);
}

void
Mesh_maps_stk::get_set_ids(Mesh_data::Entity_kind kind, 
                           std::vector<unsigned int>::iterator begin, 
                           std::vector<unsigned int>::iterator end) const
{
  get_set_ids< std::vector<unsigned int>::iterator >(kind, begin, end);
}

void
Mesh_maps_stk::get_set_ids (Mesh_data::Entity_kind kind, 
                            unsigned int * begin, 
                            unsigned int * end) const
{
  get_set_ids< unsigned int * >(kind, begin, end);
}


// Connectivity accessors
// ----------------------



// -------------------------------------------------------------
// Mesh_data::get_set
// -------------------------------------------------------------
template <typename IT>
void Mesh_maps_stk::get_set (stk::mesh::Part& set_part, Mesh_data::Entity_kind kind, 
                             Element_Category category,
                             IT begin, IT end) const
{
    ASSERT(category != GHOST);
    Entity_vector entities;
    mesh_->get_entities (set_part, category, entities);
    const Epetra_Map& themap(get_map_(kind, (category == USED)));

    // Convert to local ids.
    for (Entity_vector::const_iterator it = entities.begin (); 
         it != entities.end (); ++it) {
        unsigned int global_id((*it)->identifier());
        global_id -= 1;         // 0-based for Epetra_Map
        ASSERT(themap.MyGID(global_id));
        unsigned int local_id(themap.LID(global_id));
        *begin = local_id;
        begin++;
    }
}

template <typename IT>
void Mesh_maps_stk::get_set (unsigned int set_id, Mesh_data::Entity_kind kind, 
                             Element_Category category,
                             IT begin, IT end) const
{
    stk::mesh::EntityRank rank(kind_to_rank_(kind));
    stk::mesh::Part* set_part = mesh_->get_set (set_id, rank);
    get_set(*set_part, kind, category, begin, end);
}

void 
Mesh_maps_stk::get_set (unsigned int set_id, Mesh_data::Entity_kind kind, 
                        Element_Category category, 
                        std::vector<unsigned int>::iterator begin, 
                        std::vector<unsigned int>::iterator end) const
{
    get_set<std::vector<unsigned int>::iterator>(set_id, kind, category, begin, end);
}

void 
Mesh_maps_stk::get_set (unsigned int set_id, Mesh_data::Entity_kind kind, 
                        Element_Category category, 
                        unsigned int * begin, 
                        unsigned int * end) const
{
    get_set<unsigned int *>(set_id, kind, category, begin, end);
}


template <typename IT>
void Mesh_maps_stk::get_set (const char* name, 
                             Mesh_data::Entity_kind kind, Element_Category category,
                             IT begin, IT end) const
{
    stk::mesh::Part* set_part = mesh_->get_set (name, kind_to_rank_ (kind));
    stk::mesh::EntityRank rank(kind_to_rank_(kind));
    get_set(*set_part, kind, category, begin, end);
}

// -------------------------------------------------------------
// Mesh_maps_stk::set_coordinate
// -------------------------------------------------------------
void
Mesh_maps_stk::set_coordinate(unsigned int local_node_id, 
                              double* source_begin, double* source_end)
{
    // FIXME: not implemented
    ASSERT(false);
}

// -------------------------------------------------------------
// Mesh_maps_stk::cellgraph
// -------------------------------------------------------------
/** 
 * Make as cell-to-cell graph for the mesh.  Cell indexes are 0-based.
 * 
 * 
 * @return smart pointer to a graph instance
 */
Teuchos::RCP<Epetra_CrsGraph> 
Mesh_maps_stk::cellgraph() const
{
  const Epetra_Map& cmap(this->cell_map(false)); // 0-based
  Teuchos::RCP<Epetra_CrsGraph> result;
  result.reset(new Epetra_CrsGraph(Copy, cmap, 6));

  for (int local_idx = 0; local_idx < cmap.NumMyElements(); local_idx++) {

    const int global_idx(cmap.GID(local_idx));
      
    Entity_Ids faceids;         // 1-based
      
    mesh_->element_to_faces(global_idx + 1, faceids);

    for (Entity_Ids::iterator f = faceids.begin(); f != faceids.end(); f++) {
      Entity_Ids nbrids;        // 1-based
      mesh_->face_to_elements(*f, nbrids);
      
      int junk;
      junk = nbrids.front() - 1;
      result->InsertGlobalIndices(global_idx, 1, &junk);
      junk = nbrids.back() - 1;
      result->InsertGlobalIndices(global_idx, 1, &junk);
    }
  }  
  result->FillComplete();
  return result;
}

// // -------------------------------------------------------------
// // Mesh_maps_stk::cellgraph
// // -------------------------------------------------------------
// /** 
//  * Make as cell-to-cell graph for the mesh.  Cell indexes in the graph
//  * are 0-based.
//  * 
//  * 
//  * @return smart pointer to a graph instance
//  */
// Teuchos::RCP<Epetra_CrsGraph> 
// Mesh_maps_stk::cellgraph() const
// {
//   static const double ONE(1.0);
//   const Epetra_Map& cmap(this->cell_map(false)); // 0-based
//   Epetra_CrsMatrix cmat(Copy, cmap, 0, false);   // 0-based

//   int junk;

//   for (int local_idx = 0; local_idx < cmap.NumMyElements(); local_idx++) {

//     std::set<int> nbrset;       // 0-based
//     const int global_idx(cmap.GID(local_idx));
//     junk = global_idx;
//     cmat.InsertGlobalValues(global_idx, 1, &ONE, &junk);
      
//     Entity_Ids faceids;         // 1-based

//     mesh_->element_to_faces(global_idx + 1, faceids);
    
//     for (Entity_Ids::iterator f = faceids.begin(); f != faceids.end(); f++) {
//       Entity_Ids nbrids;        // 1-based
//       mesh_->face_to_elements(*f, nbrids);
//       if (nbrids.front() != global_idx) {
//         junk = nbrids.front() - 1;
//         cmat.InsertGlobalValues(global_idx, 1, &ONE, &junk);
//       }
//       if (nbrids.back() != global_idx) {
//         junk = nbrids.back() - 1;
//         cmat.InsertGlobalValues(global_idx, 1, &ONE, &junk);
//       }
//     }
//   }
//   cmat.FillComplete();

//   Teuchos::RCP<Epetra_CrsGraph> result(new Epetra_CrsGraph(cmat.Graph()));

//   return result;
// }



// -------------------------------------------------------------
// Mesh_maps_stk::redistribute
// -------------------------------------------------------------
/** 
 * This routine redistributes cells amongst the processors according
 * to the specified @c cellmap0.  The indexes in @c cellmap0 are
 * 0-based.
 * 
 * @param cellmap 0-based map of desired cell ownership
 */
void 
Mesh_maps_stk::redistribute(const Epetra_Map& cellmap0)
{
  // change the 0-based map to a 1-based one

  
  std::vector<int> gids(cellmap0.NumMyElements());
  std::copy(cellmap0.MyGlobalElements(), 
            cellmap0.MyGlobalElements() + cellmap0.NumMyElements(), 
            gids.begin());
  std::for_each(gids.begin(), gids.end(), bl::_1 += 1);
  Teuchos::RCP< Epetra_Map > cellmap1(new Epetra_Map(cellmap0.NumGlobalElements(), 
                                                     gids.size(), &gids[0], 1, *communicator_));

  mesh_->redistribute(*cellmap1);

  this->build_maps_();
}

/** 
 * This routine uses the Isorropia::Epetra::Partitioner, with the
 * specified parameter list, to repartition and redistribute the mesh.
 * 
 * @param paramlist parameters for Isorropia::Epetra::Partitioner
 */
void 
Mesh_maps_stk::redistribute(const Teuchos::ParameterList &paramlist)
{
  Teuchos::RCP<const Epetra_CrsGraph> cgraph(this->cellgraph());

  Isorropia::Epetra::Partitioner partitioner(cgraph, paramlist, false);
  partitioner.partition();
  Teuchos::RCP< Epetra_Map > newcmap(partitioner.createNewMap()); // 0-based
  this->redistribute(*newcmap);
}


} // close namespace STK_mesh
