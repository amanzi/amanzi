// -------------------------------------------------------------
// file: Mesh_STK.cc
// -------------------------------------------------------------
/**
 * @file   Mesh_STK.cc
 * @author William A. Perkins
 * @date Mon Aug  8 12:16:12 2011
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// Created May  2, 2011 by William A. Perkins
// Last Change: Mon Aug  8 12:16:12 2011 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

#include <algorithm>
#include <boost/lambda/lambda.hpp>
namespace bl = boost::lambda;

#include <Shards_CellTopology.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>
#include <Isorropia_EpetraPartitioner.hpp>

#include "Entity_map.hh"
#include "Mesh_STK_Impl.hh"
#include "Mesh_STK.hh"
#include "stk_mesh_error.hh"

static const int ZERO = 0;

namespace Amanzi {
namespace AmanziMesh {

// -------------------------------------------------------------
//  class Mesh_STK
// -------------------------------------------------------------

// -------------------------------------------------------------
// Mesh_STK static members
// -------------------------------------------------------------

const Entity_kind Mesh_STK::kinds_[] = { NODE, FACE, CELL };
const unsigned int Mesh_STK::num_kinds_ = 
    sizeof(Mesh_STK::kinds_)/sizeof(Entity_kind);


// -------------------------------------------------------------
// Mesh_STK:: constructors / destructor
// -------------------------------------------------------------
Mesh_STK::Mesh_STK(STK::Mesh_STK_Impl_p mesh)
    : communicator_(new Epetra_MpiComm(mesh->communicator())),
      mesh_(mesh), 
      entity_map_(3),                   // FIXME: can be 2; take from mesh
      map_owned_(), map_used_()
{
  Mesh::set_comm(communicator_->GetMpiComm());
  build_maps_();
}


Mesh_STK::~Mesh_STK(void)
{
  // empty
}

// -------------------------------------------------------------
// Mesh_STK::entity_get_ptype
// -------------------------------------------------------------
// FIXME: can be implemented in Mesh
/** 
 * 
 * 
 * @param kind 
 * @param entid the @em local identifier
 * 
 * @return either @c OWNED or @c GHOST
 */
Parallel_type 
Mesh_STK::entity_get_ptype(const Entity_kind kind, 
                           const Entity_ID entid) const
{
  ASSERT(entity_valid_kind(kind));

  bool owned, used;
  switch (kind) {
    case CELL:
      owned = this->cell_epetra_map(false).MyLID(entid);
      used = this->cell_epetra_map(true).MyLID(entid);
      break;
    case FACE:
      owned = this->face_epetra_map(false).MyLID(entid);
      used = this->face_epetra_map(true).MyLID(entid);
      break;
    case NODE:
      owned = this->node_epetra_map(false).MyLID(entid);
      used = this->node_epetra_map(true).MyLID(entid);
      break;
    default:
      Exceptions::amanzi_throw( STK::Error ("Unknown Entity Kind") );
  }
  if (owned) {
    return OWNED;
  } else if (used) {
    return GHOST;
  } else {
    Exceptions::amanzi_throw( STK::Error ("Invalid local identifier") );
  }
}


// -------------------------------------------------------------
// Mesh_STK::cell_get_type
// -------------------------------------------------------------
Cell_type 
Mesh_STK::cell_get_type(const Entity_ID cellid) const
{
  stk::mesh::EntityId global_cell_id = this->GID(cellid, CELL);
  global_cell_id += 1;        // need 1-based for stk::mesh

  stk::mesh::EntityRank rank(entity_map_.kind_to_rank(CELL));
  stk::mesh::Entity* cell(mesh_->id_to_entity(rank, global_cell_id));

  // FIXME: Throw instead?
  ASSERT(cell != NULL);

  const CellTopologyData* topo = stk::mesh::get_cell_topology (*cell);

  // FIXME: Polyhedral, 2D not yet supported

  Cell_type result(UNKNOWN);
  ASSERT(topo != NULL);
  switch (topo->node_count) {
    case (8):
      result = HEX;
      break;
    case (6):
      result = PRISM;
      break;
    case (5):
      result = PYRAMID;
      break;
    case (4): 
      result = TET;
      break;
    default:
      Exceptions::amanzi_throw( STK::Error ("Unsupported cell type") );
  }
  return result;
}

// -------------------------------------------------------------
// Mesh_STK::num_entities
// -------------------------------------------------------------
unsigned int 
Mesh_STK::num_entities (const Entity_kind kind,
                        const Parallel_type ptype) const
{
  ASSERT(entity_valid_kind(kind));
  stk::mesh::EntityRank rank(entity_map_.kind_to_rank(kind));
  return mesh_->count_entities(rank, ptype);
}

// -------------------------------------------------------------
// Mesh_STK::GID
// -------------------------------------------------------------
unsigned int 
Mesh_STK::GID(const Entity_ID lid, const Entity_kind kind) const
{
  ASSERT(entity_valid_kind(kind));
  unsigned int result;

  switch (kind) {
    case CELL:
      ASSERT(this->cell_epetra_map(true).MyLID(lid));
      result = this->cell_epetra_map(true).GID(lid);
      break;
    case FACE:
      ASSERT(this->face_epetra_map(true).MyLID(lid));
      result = this->face_epetra_map(true).GID(lid);
      break;
    case NODE:
      ASSERT(this->node_epetra_map(true).MyLID(lid));
      result = this->node_epetra_map(true).GID(lid);
      break;
    default:
      Exceptions::amanzi_throw( STK::Error ("Unknown Entity Kind") );
  }
  
  return result;
}

// -------------------------------------------------------------
// Mesh_STK::LID
// -------------------------------------------------------------
Entity_ID 
Mesh_STK::LID(const Entity_ID& gid, const Entity_kind& kind) const
{
  ASSERT (entity_valid_kind(kind));
  unsigned int result;

  switch (kind) {
    case CELL:
      ASSERT(this->cell_epetra_map(true).MyGID(gid));
      result = this->cell_epetra_map(true).LID(gid);
      break;
    case FACE:
      ASSERT(this->face_epetra_map(true).MyGID(gid));
      result = this->face_epetra_map(true).LID(gid);
      break;
    case NODE:
      ASSERT(this->node_epetra_map(true).MyGID(gid));
      result = this->node_epetra_map(true).LID(gid);
      break;
    default:
      Exceptions::amanzi_throw( STK::Error ("Unknown Entity Kind") );
  }
  return result;
}

// -------------------------------------------------------------
// Mesh_STK::cell_get_faces
// -------------------------------------------------------------
void 
Mesh_STK::cell_get_faces (const Entity_ID cellid, 
                          std::vector<Entity_ID> *outfaceids) const
{
  stk::mesh::EntityId global_cell_id = this->GID(cellid, CELL);
  global_cell_id += 1;                  // need 1-based for stk::mesh

  STK::Entity_Ids stk_face_ids;
  mesh_->element_to_faces(global_cell_id, stk_face_ids);
  // 0-based for Epetra_Map
  std::for_each(stk_face_ids.begin(), stk_face_ids.end(), bl::_1 -= 1);

  outfaceids->clear();
  for (STK::Entity_Ids::iterator f = stk_face_ids.begin(); 
       f != stk_face_ids.end(); f++) {
    stk::mesh::EntityId global_face_id(*f);
    ASSERT(this->face_epetra_map(true).MyGID(global_face_id));
    stk::mesh::EntityId local_face_id = 
        this->face_epetra_map(true).LID(global_face_id);
    outfaceids->push_back(local_face_id);
  }
}  

// -------------------------------------------------------------
// Mesh_STK::cell_get_face_dirs
// -------------------------------------------------------------
void 
Mesh_STK::cell_get_face_dirs (const Entity_ID cellid, 
                              std::vector<int> *face_dirs) const
{
  stk::mesh::EntityId global_cell_id = this->GID(cellid, CELL);
  global_cell_id += 1;        // need 1-based for stk::mesh
  
  face_dirs->clear();
  mesh_->element_to_face_dirs(global_cell_id, *face_dirs);
}


// -------------------------------------------------------------
// Mesh_STK::cell_get_nodes
// -------------------------------------------------------------
void 
Mesh_STK::cell_get_nodes (const Entity_ID cellid, 
                          Entity_ID_List *outnodeids) const
{
  stk::mesh::EntityId global_cell_id = this->GID(cellid, CELL);
  global_cell_id += 1;        // need 1-based for stk::mesh

  STK::Entity_Ids node_ids;
  mesh_->element_to_nodes(global_cell_id, node_ids);
  std::for_each(node_ids.begin(), node_ids.end(), bl::_1 -= 1); // 0-based for Epetra_Map

  outnodeids->clear();
  for (STK::Entity_Ids::iterator n = node_ids.begin(); n != node_ids.end(); n++) {
    stk::mesh::EntityId global_node_id(*n);
    ASSERT(this->node_epetra_map(true).MyGID(global_node_id));
    stk::mesh::EntityId local_node_id = 
        this->node_epetra_map(true).LID(global_node_id);
    outnodeids->push_back(local_node_id);
  }
}

// -------------------------------------------------------------
// Mesh_STK::face_get_nodes
// -------------------------------------------------------------
void 
Mesh_STK::face_get_nodes (const Entity_ID faceid, 
                          Entity_ID_List *outnodeids) const
{
  stk::mesh::EntityId global_face_id = this->GID(faceid, FACE);
  global_face_id += 1;        // need 1-based for stk::mesh

  STK::Entity_Ids node_ids;
  mesh_->face_to_nodes(global_face_id, node_ids);
  std::for_each(node_ids.begin(), node_ids.end(), bl::_1 -= 1); // 0-based for Epetra_Map

  outnodeids->clear();
  for (STK::Entity_Ids::iterator n = node_ids.begin(); n != node_ids.end(); n++) {
    stk::mesh::EntityId global_node_id(*n);
    ASSERT(this->node_epetra_map(true).MyGID(global_node_id));
    stk::mesh::EntityId local_node_id = 
        this->node_epetra_map(true).LID(global_node_id);
    outnodeids->push_back(local_node_id);
  }
  ASSERT(!outnodeids->empty());
}

// -------------------------------------------------------------
// Mesh_STK::node_get_cells
// -------------------------------------------------------------
void 
Mesh_STK::node_get_cells(const Entity_ID nodeid, 
                         const Parallel_type ptype,
                         Entity_ID_List *outcellids) const
{
  stk::mesh::EntityId global_node_id = this->GID(nodeid, NODE);
  global_node_id += 1;        // need 1-based for stk::mesh
  STK::Entity_Ids cell_ids;
  mesh_->node_to_elements(global_node_id, cell_ids);
  std::for_each(cell_ids.begin(), cell_ids.end(), bl::_1 -= 1); // 0-based for Epetra_Map

  outcellids->clear();
  for (STK::Entity_Ids::iterator i = cell_ids.begin(); i != cell_ids.end(); i++) {
    Entity_ID local_cell_id(this->cell_epetra_map(true).LID(*i));
    Parallel_type theptype(this->entity_get_ptype(CELL, local_cell_id));
    if (theptype == OWNED && (ptype == OWNED || ptype == USED)) {
      outcellids->push_back(local_cell_id);
    } else if (theptype == GHOST && (ptype == GHOST || ptype == USED)) {
      outcellids->push_back(local_cell_id);
    }
  }
}

// -------------------------------------------------------------
// Mesh_STK::node_get_faces
// -------------------------------------------------------------
void 
Mesh_STK::node_get_faces(const Entity_ID nodeid, 
                         const Parallel_type ptype,
                         Entity_ID_List *outfaceids) const
{
  stk::mesh::EntityId global_node_id = this->GID(nodeid, NODE);
  global_node_id += 1;        // need 1-based for stk::mesh
  STK::Entity_Ids face_ids;
  mesh_->node_to_faces(global_node_id, face_ids);
  std::for_each(face_ids.begin(), face_ids.end(), bl::_1 -= 1); // 0-based for Epetra_Map
  outfaceids->clear();
  for (STK::Entity_Ids::iterator i = face_ids.begin(); i != face_ids.end(); i++) {
    Entity_ID local_face_id(this->cell_epetra_map(true).LID(*i));
    Parallel_type theptype(this->entity_get_ptype(FACE, local_face_id));
    if (theptype == OWNED && (ptype == OWNED || ptype == USED)) {
      outfaceids->push_back(local_face_id);
    } else if (theptype == GHOST && (ptype == GHOST || ptype == USED)) {
      outfaceids->push_back(local_face_id);
    }
  }
}

// -------------------------------------------------------------
// Mesh_STK::node_get_cell_faces
// FIXME: should be implemented like this generally by Mesh
// -------------------------------------------------------------
void 
Mesh_STK::node_get_cell_faces(const Entity_ID nodeid, 
                              const Entity_ID cellid,
                              const Parallel_type ptype,
                              Entity_ID_List *outfaceids) const
{
  Entity_ID_List node_faces;
  Entity_ID_List cell_faces;
  Entity_ID_List outids;

  this->node_get_faces(nodeid, ptype, &node_faces);
  std::sort(node_faces.begin(), node_faces.end());
  this->node_get_cells(cellid, ptype, &cell_faces);
  std::sort(cell_faces.begin(), cell_faces.end());

  std::set_intersection(node_faces.begin(), node_faces.end(),
                        cell_faces.begin(), cell_faces.end(),
                        std::back_inserter(*outfaceids));

  ASSERT(!outfaceids->empty());
}

// -------------------------------------------------------------
// Mesh_STK::face_get_cells
// -------------------------------------------------------------
void
Mesh_STK::face_get_cells(const Entity_ID faceid, 
                         const Parallel_type ptype,
                         Entity_ID_List *outcellids) const
{
  stk::mesh::EntityId global_face_id = 
      this->face_epetra_map(true).GID(faceid);
  
  STK::Entity_Ids cell_ids;
  mesh_->face_to_elements(global_face_id, cell_ids);
  std::for_each(cell_ids.begin(), cell_ids.end(), bl::_1 -= 1); // 0-based for Epetra_Map
  for (STK::Entity_Ids::iterator i = cell_ids.begin(); i != cell_ids.end(); i++) {
    Entity_ID local_cell_id(this->cell_epetra_map(true).LID(*i));
    Parallel_type theptype(this->entity_get_ptype(FACE, local_cell_id));
    if (theptype == OWNED && (ptype == OWNED || ptype == USED)) {
      outcellids->push_back(local_cell_id);
    } else if (theptype == GHOST && (ptype == GHOST || ptype == USED)) {
      outcellids->push_back(local_cell_id);
    }
  }
}

// -------------------------------------------------------------
// Mesh_STK::cell_get_face_adj_cells
// FIXME: not Mesh_STK specific
// -------------------------------------------------------------
void 
Mesh_STK::cell_get_face_adj_cells(const Entity_ID cellid,
                                  const Parallel_type ptype,
                                  Entity_ID_List *fadj_cellids) const
{
  Entity_ID_List faces;
  cell_get_faces(cellid, &faces);

  fadj_cellids->clear();
  for (Entity_ID_List::iterator f = faces.begin(); 
       f != faces.end(); f++) {
    Entity_ID_List fcells;
    face_get_cells(*f, ptype, &fcells);
    if (!fcells.empty()) {
      if (fcells.front() != cellid) {
        fadj_cellids->push_back(fcells.front());
      } else if (fcells.back() != cellid) {
        fadj_cellids->push_back(fcells.back());
      }
    }
  }
}

// -------------------------------------------------------------
// Mesh_STK::cell_get_node_adj_cells
// FIXME: can be implemented in Mesh
// -------------------------------------------------------------
void Mesh_STK::cell_get_node_adj_cells(const Entity_ID cellid,
                                       const Parallel_type ptype,
                                       Entity_ID_List *nadj_cellids) const
{
  Entity_ID_List nodes;
  cell_get_nodes(cellid, &nodes);

  nadj_cellids->clear();
  for (Entity_ID_List::iterator n = nodes.begin(); 
       n != nodes.end(); n++) {
    Entity_ID_List ncells;
    node_get_cells(*n, ptype, &ncells);
    for (Entity_ID_List::iterator c = ncells.begin(); 
         c != ncells.end(); c++) {
      if (*c != cellid) {
        nadj_cellids->push_back(*c);
      }
    }
  }
}


// -------------------------------------------------------------
// Mesh_STK::cell_get_type_4viz
// -------------------------------------------------------------
Cell_type 
Mesh_STK::cell_get_type_4viz(const Entity_ID cellid) const
{
  // FIXME: degenerate cells
  return cell_get_type(cellid);
}


// -------------------------------------------------------------
// Mesh_STK::cell_get_nodes_4viz
// -------------------------------------------------------------
void 
Mesh_STK::cell_get_nodes_4viz (const Entity_ID cellid, 
                               Entity_ID_List *nodeids) const
{
  // FIXME: degenerate cells
  cell_get_nodes(cellid, nodeids);
}


// -------------------------------------------------------------
// Mesh_STK::node_get_coordinates
// -------------------------------------------------------------
void 
Mesh_STK::node_get_coordinates (const Entity_ID nodeid, 
                                AmanziGeometry::Point *x) const
{
  stk::mesh::EntityId gid(this->GID(nodeid, NODE));
  gid++;                                // need 1-based for stk::mesh
  const double *c = mesh_->coordinates(gid);
  switch (space_dimension()) {
    case 2:
      x->init(2);
      x->set(c[0], c[1]);
      break;
    case 3:
      x->init(3);
      x->set(c[0], c[1], c[2]);
      break;
    default:
      Exceptions::amanzi_throw( STK::Error ("Unknown dimension") );      
  }
}

// -------------------------------------------------------------
// Mesh_STK::face_get_coordinates
// -------------------------------------------------------------
void 
Mesh_STK::face_get_coordinates (const Entity_ID faceid, 
                                std::vector<AmanziGeometry::Point> *fcoords) const
{
  stk::mesh::EntityId gid(this->GID(faceid, FACE));
  gid++;                                // need 1-based for stk::mesh
  
  STK::Entity_Ids nodes;
  mesh_->face_to_nodes(gid, nodes);
  fcoords->clear();
  for (STK::Entity_Ids::iterator n = nodes.begin(); n != nodes.end(); n++) {
    const double *c(mesh_->coordinates(*n));
    AmanziGeometry::Point p;
    switch (space_dimension()) {
      case (2):
        p.init(2);
        p.set(c[0], c[1]);
        break;
      case (3):
        p.init(3);
        p.set(c[0], c[1], c[2]);
        break;
      default:
        Exceptions::amanzi_throw( STK::Error ("Unknown dimension") );      
    }
    fcoords->push_back(p);
  }
}

// -------------------------------------------------------------
// Mesh_STK::cell_get_coordinates
// -------------------------------------------------------------
void 
Mesh_STK::cell_get_coordinates (const Entity_ID cellid, 
                                std::vector<AmanziGeometry::Point> *ccoords) const
{
  stk::mesh::EntityId gid(this->GID(cellid, CELL));
  gid++;                                // need 1-based for stk::mesh
  
  STK::Entity_Ids nodes;
  mesh_->element_to_nodes(gid, nodes);
  ccoords->clear();
  for (STK::Entity_Ids::iterator n = nodes.begin(); n != nodes.end(); n++) {
    const double *c(mesh_->coordinates(*n));
    AmanziGeometry::Point p;
    switch (space_dimension()) {
      case (2):
        p.init(2);
        p.set(c[0], c[1]);
        break;
      case (3):
        p.init(3);
        p.set(c[0], c[1], c[2]);
        break;
      default:
        Exceptions::amanzi_throw( STK::Error ("Unknown dimension") );      
    }
    ccoords->push_back(p);
  }
}
    


// -------------------------------------------------------------
// Mesh_STK::cell_epetra_map
// -------------------------------------------------------------
const Epetra_Map& 
Mesh_STK::cell_epetra_map(bool include_ghost) const
{
  return get_map_(CELL, include_ghost);
}

// -------------------------------------------------------------
// Mesh_STK::cell_epetra_map
// -------------------------------------------------------------
const Epetra_Map& 
Mesh_STK::face_epetra_map(bool include_ghost) const
{
  return get_map_(FACE, include_ghost);
}

// -------------------------------------------------------------
// Mesh_STK::cell_epetra_map
// -------------------------------------------------------------
const Epetra_Map& 
Mesh_STK::node_epetra_map (bool include_ghost) const
{
  return get_map_(NODE, include_ghost);
}


// -------------------------------------------------------------
// Mesh_STK::get_set_size (by int setid)
// -------------------------------------------------------------
unsigned int 
Mesh_STK::get_set_size (const Set_ID setid, const Entity_kind kind, 
                        const Parallel_type ptype) const
{
  Entity_ID_List setents;

  get_set_entities(setid, kind, ptype, &setents);

  return setents.size();
}


// -------------------------------------------------------------
// Mesh_STK::get_set_size (by char *setname)
// -------------------------------------------------------------
unsigned int 
Mesh_STK::get_set_size (const char *setname, const Entity_kind kind, 
                        const Parallel_type ptype) const
{
  Entity_ID_List setents;

  get_set_entities(setname, kind, ptype, &setents);

  return setents.size();
}


// -------------------------------------------------------------
// Mesh_STK::get_set_size (by std::string setname)
// -------------------------------------------------------------
unsigned int 
Mesh_STK::get_set_size (const std::string setname, const Entity_kind kind, 
                        const Parallel_type ptype) const
{
  Entity_ID_List setents;

  get_set_entities(setname, kind, ptype, &setents);

  return setents.size();
}


// -------------------------------------------------------------
// Mesh_STK::get_set_entities  (by int setid)
// -------------------------------------------------------------
void 
Mesh_STK::get_set_entities (const Set_ID setid,
                            const Entity_kind kind, 
                            const Parallel_type ptype, 
                            Entity_ID_List *entids) const
{  
  AmanziGeometry::GeometricModelPtr gm = Mesh::geometric_model();
  AmanziGeometry::RegionPtr rgn = gm->FindRegion(setid);

  std::cerr << "DEPRECATED METHOD!" << std::endl;
  std::cerr << "Call get_set_entities with setname instead of setid" << std::endl;

  if (!rgn)
    {
      std::cerr << "Mesh_STK::get_set_entities: No region with id" << setid << std::endl;
      std::cerr << "Cannot construct requested set" << std::endl;
      throw std::exception();
    }

  
  get_set_entities(rgn->name(),kind,ptype,entids);
}


// -------------------------------------------------------------
// Mesh_STK::get_set_entities  (by char *setname)
// -------------------------------------------------------------
void 
Mesh_STK::get_set_entities (const char *setname,
                            const Entity_kind kind, 
                            const Parallel_type ptype, 
                            Entity_ID_List *entids) const
{
  std::string setname1(setname);

  get_set_entities(setname1, kind, ptype, entids);
}

// -------------------------------------------------------------
// Mesh_STK::get_set_entities  (by std::string setname)
// -------------------------------------------------------------
void 
Mesh_STK::get_set_entities (const std::string setname,
                            const Entity_kind kind, 
                            const Parallel_type ptype, 
                            Entity_ID_List *entids) const
{

  ASSERT (entity_valid_ptype(ptype));
  stk::mesh::EntityRank rank(entity_map_.kind_to_rank(kind));

  stk::mesh::Part *part;

  int celldim = cell_dimension();
  int spacedim = space_dimension();

  AmanziGeometry::GeometricModelPtr gm = geometric_model();
  AmanziGeometry::RegionPtr rgn = gm->FindRegion(setname);

  // Did not find the region
  
  if (rgn == NULL) 
    {
      std::cerr << "Geometric model has no region named " << setname << std::endl;
      std::cerr << "Cannot construct set by this name" << std::endl;
      throw std::exception();
    }
  
  if (rgn->type() == AmanziGeometry::LABELEDSET)
    {

      // Region is of type labeled set and a mesh set should have been
      // initialized from the input file
  
      AmanziGeometry::LabeledSetRegionPtr lsrgn = dynamic_cast<AmanziGeometry::LabeledSetRegionPtr> (rgn);
      std::string label = lsrgn->label();
      std::string entity_type = lsrgn->entity_str();
      
      
      if ((kind == CELL && entity_type != "CELL") ||
          (kind == FACE && entity_type != "FACE") ||
          (kind == NODE && entity_type != "NODE"))
        {
          std::cerr << "Found labeled set region named " << setname << " but it"<< std::endl;
          std::cerr << "contains entities of type " << entity_type << ", not the requested type" << std::endl;
          
          throw std::exception();
        } 
      
      std::stringstream internal_name;

      if (kind == CELL)
        internal_name << "element block " << label;
      else if (kind == FACE)
        internal_name << "side set " << label;
      else if (kind == NODE)
        internal_name << "node set " << label;

 
      part = mesh_->get_set(internal_name.str(), rank);

      if (part == NULL)
        {
          std::cerr << "Mesh set " << setname << " should have been read in" << std::endl;
          std::cerr << "as set " << label << " from the mesh file. Its absence" << std::endl;
          std::cerr << "indicates an error in the input file or in the mesh file" << std::endl;
          
          throw std::exception();
        }
    }
  else
    {
      part = mesh_->get_set(setname, rank);
    }



  if (part == NULL) 
    {
      // This mesh entity set did not exist in the input mesh file
      // Create it based on the region defintion

      Entity_ID_List entids2;

      switch (kind) 
        {
        case CELL:           // cellsets

          if (rgn->dimension() != celldim) 
            {
              std::cerr << "No region of dimension " << celldim << " defined in geometric model" << std::endl;
              std::cerr << "Cannot construct cell set from region " << setname << std::endl;
            }

          if (rgn->type() == AmanziGeometry::BOX) 
            {
              if (celldim == 3)  // 3D Mesh
                {
                  int ncell = num_entities(CELL,USED);

                  for (int icell = 0; icell < ncell; icell++) 
                    {
                      if (rgn->inside(cell_centroid(icell))) 
                        {
                          entids2.push_back(icell);
                        }
                    }
                }
              else 
                {
                  throw std::exception();   // FIXME: Not Implemented
                }
            }
          else
            {
              std::cerr << "Region type not applicable/supported for cell sets" << std::endl;
              throw std::exception();
            }
      
          break;

        case FACE:      // sidesets
          if (rgn->dimension() != celldim-1) 
            {
              std::cerr << "No region of dimension " << celldim-1 << " defined in geometric model" << std::endl;
              std::cerr << "Cannot construct cell set from region " << setname << std::endl;
            }

          if (rgn->type() == AmanziGeometry::BOX) 
            {
              if (celldim == 3)   // 3D meshes, 2D sidesets
                {
                  int nface = num_entities(FACE,USED);

                  for (int iface = 0; iface < nface; iface++)
                    {
                      if (rgn->inside(face_centroid(iface))) 
                        {
                          entids2.push_back(iface);
                        }
                    }
                }
              else if (celldim == 2)  // 2D meshes, 1D sidesets
                {
                  throw std::exception();    // FIXME: Not implemented
                }
            }
          else if (rgn->type() == AmanziGeometry::PLANE) 
            {
              if (celldim == 3)
                {
                  int nface = num_entities(FACE,USED);

                  for (int iface = 0; iface < nface; iface++) 
                    {
                      std::vector<AmanziGeometry::Point> fcoords;
                      face_get_coordinates(iface, &fcoords);

                      bool on_plane = true;
                      for (int j = 0; j < fcoords.size(); j++)
                        {
                          if (!rgn->inside(fcoords[j])) 
                            {                          
                              on_plane = false;
                              break;
                            }
                        }

                      if (on_plane)
                        entids2.push_back(iface);
                    }

                }
              else if (celldim == 2)  // 2D Meshes, 1D sidesets
                {
                  throw std::exception();   // FIXME: Not implemented
                }
            }
          else 
            {
              std::cerr << "Region type not applicable/supported for face sets" << std::endl;
              throw std::exception();
            }

          break;

        case NODE:          // Nodesets

          if (rgn->type() == AmanziGeometry::BOX ||
              rgn->type() == AmanziGeometry::PLANE) 
            {
              int nnode = num_entities(NODE,USED);

              for (int inode = 0; inode < nnode; inode++) 
                {
                  AmanziGeometry::Point vcoord;
                  node_get_coordinates(inode, &vcoord);

                  if (rgn->inside(vcoord)) 
                    {                          
                      entids2.push_back(inode);
                    }
                }
            }
          else 
            {
              std::cerr << "Region type not applicable/supported for node sets" << std::endl;
              throw std::exception();
            }

          break;
        }

      part = mesh_->add_set(setname,rgn->id(),rank,entids2);
    }


  ASSERT(part != NULL);

  STK::Entity_vector entities;
  mesh_->get_entities(*part, ptype, entities);

  entids->clear();
  for (STK::Entity_vector::iterator e = entities.begin(); 
       e != entities.end(); e++) {
    Entity_ID gid((*e)->identifier());
    gid--;                              // 0-based for Epetra_Map
    Entity_ID lid(this->LID(gid, kind));
    entids->push_back(lid);
  }
}
    


// -------------------------------------------------------------
// Mesh_STK::build_maps_
// -------------------------------------------------------------
static void
extract_global_ids(const STK::Entity_vector& entities,
                   std::vector<int>& ids)
{
  ids.clear();
  ids.reserve(entities.size());
  for (STK::Entity_vector::const_iterator i = entities.begin(); i != entities.end(); i++) {
    ids.push_back((*i)->identifier());
  }
  std::for_each(ids.begin(), ids.end(), bl::_1 -= 1);
}

void 
Mesh_STK::build_maps_ ()
{
  map_owned_.clear();
  map_used_.clear();

  // For each of Elements, Faces, Nodes:

  for (unsigned int i = 0; i < num_kinds_; i++) {
    Entity_kind kind(kinds_[i]);
    stk::mesh::EntityRank rank = entity_map_.kind_to_rank (kind);

    STK::Entity_vector entities;
    std::vector<int> entity_ids;
    Teuchos::RCP<Epetra_Map> map;

    // Get the collection of "owned" entities
    mesh_->get_entities (rank, OWNED, entities);
    extract_global_ids(entities, entity_ids);

    map.reset(new Epetra_Map(-1, entity_ids.size(), &entity_ids[0], ZERO, *communicator_));
    map_owned_.insert(MapSet::value_type(kind, map));

    // Get the collection of "ghost" entities
    STK::Entity_vector ghost_entities;
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
Mesh_STK::get_map_(const Entity_kind& kind, const bool& include_ghost) const
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

// -------------------------------------------------------------
// Mesh_STK::cellgraph
// -------------------------------------------------------------
/** 
 * Make as cell-to-cell graph for the mesh.  Cell indexes are 0-based.
 * 
 * 
 * @return smart pointer to a graph instance
 */
Teuchos::RCP<Epetra_CrsGraph> 
Mesh_STK::cellgraph() const
{
  const Epetra_Map& cmap(this->cell_map(false)); // 0-based
  Teuchos::RCP<Epetra_CrsGraph> result;
  result.reset(new Epetra_CrsGraph(Copy, cmap, 6));

  for (int local_idx = 0; local_idx < cmap.NumMyElements(); local_idx++) {

    const int global_idx(cmap.GID(local_idx));
      
    STK::Entity_Ids faceids;         // 1-based
      
    mesh_->element_to_faces(global_idx + 1, faceids);

    for (STK::Entity_Ids::iterator f = faceids.begin(); f != faceids.end(); f++) {
      STK::Entity_Ids nbrids;        // 1-based
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

// -------------------------------------------------------------
// Mesh_STK::redistribute
// -------------------------------------------------------------
/** 
 * This routine redistributes cells amongst the processors according
 * to the specified @c cellmap0.  The indexes in @c cellmap0 are
 * 0-based.
 * 
 * @param cellmap 0-based map of desired cell ownership
 */
void 
Mesh_STK::redistribute(const Epetra_Map& cellmap0)
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
Mesh_STK::redistribute(const Teuchos::ParameterList &paramlist)
{
  Teuchos::RCP<const Epetra_CrsGraph> cgraph(this->cellgraph());

  Isorropia::Epetra::Partitioner partitioner(cgraph, paramlist, false);
  partitioner.partition();
  Teuchos::RCP< Epetra_Map > newcmap(partitioner.createNewMap()); // 0-based
  this->redistribute(*newcmap);
}

} // namespace AmanziMesh
} // namespace Amanzi
