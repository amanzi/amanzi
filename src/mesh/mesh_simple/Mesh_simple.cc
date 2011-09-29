/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

#include <algorithm>

#include "Mesh_simple.hh"
#include "GenerationSpec.hh"
#include <Teuchos_RCP.hpp>

namespace Amanzi
{

namespace AmanziMesh
{

Mesh_simple::Mesh_simple (double x0, double y0, double z0,
                          double x1, double y1, double z1,
                          int nx, int ny, int nz,
                          Epetra_Comm *communicator):
    nx_(nx), ny_(ny), nz_(nz),
    x0_(x0), x1_(x1),
    y0_(y0), y1_(y1),
    z0_(z0), z1_(z1),
    communicator_(communicator)

{
  update();

  Mesh::set_comm(communicator);
}

Mesh_simple::Mesh_simple (const GenerationSpec& gspec,
                          Epetra_Comm *communicator ) :
    communicator_(communicator)
{

  // read the parameters from the specification

  nx_ = gspec.xcells();
  ny_ = gspec.ycells();
  nz_ = gspec.zcells();
  
  x0_ = gspec.domain().point0().x();
  x1_ = gspec.domain().point1().x();

  y0_ = gspec.domain().point0().y();
  y1_ = gspec.domain().point1().y();

  z0_ = gspec.domain().point0().z();
  z1_ = gspec.domain().point1().z();


  std::copy(gspec.block_begin(), gspec.block_end(), 
            std::back_inserter(mesh_blocks_));
  
  update();

  Mesh::set_comm(communicator);
}


Mesh_simple::~Mesh_simple()
{
  delete cell_map_;
  delete face_map_;
  delete node_map_;
}



void Mesh_simple::update ()
{
  clear_internals_ ();
  update_internals_ ();
}

void Mesh_simple::clear_internals_ ()
{

  coordinates_.resize(0);

  cell_to_face_.resize(0);
  cell_to_node_.resize(0);
  face_to_node_.resize(0);

  side_sets_.resize(0);

}


void Mesh_simple::update_internals_()
{

  num_cells_ = nx_ * ny_ * nz_;
  num_nodes_ = (nx_+1)*(ny_+1)*(nz_+1);
  num_faces_ = (nx_+1)*(ny_)*(nz_) + (nx_)*(ny_+1)*(nz_) + (nx_)*(ny_)*(nz_+1);

  coordinates_.resize(3*num_nodes_);

  double hx = (x1_ - x0_)/nx_;
  double hy = (y1_ - y0_)/ny_;
  double hz = (z1_ - z0_)/nz_;

  for (int iz=0; iz<=nz_; iz++)
    for (int iy=0; iy<=ny_; iy++)
      for (int ix=0; ix<=nx_; ix++)
      {
        int istart = 3*node_index_(ix,iy,iz);

        coordinates_[ istart ]     = x0_ + ix*hx;
        coordinates_[ istart + 1 ] = y0_ + iy*hy;
        coordinates_[ istart + 2 ] = z0_ + iz*hz;
      }

  cell_to_face_.resize(6*num_cells_);
  cell_to_face_dirs_.resize(6*num_cells_);
  cell_to_node_.resize(8*num_cells_);
  face_to_node_.resize(4*num_faces_);
  node_to_face_.resize(13*num_nodes_); // 1 extra for num faces
  node_to_cell_.resize(9*num_nodes_); // 1 extra for num cells
  face_to_cell_.resize(2*num_faces_); 
  face_to_cell_.assign(2*num_faces_,-1); 
                                          

  // loop over cells and initialize cell_to_node_
  for (int iz=0; iz<nz_; iz++)
    for (int iy=0; iy<ny_; iy++)
      for (int ix=0; ix<nx_; ix++)
      {
        int jstart=0;
        int istart = 8 * cell_index_(ix,iy,iz);
        int ncell=0;

        cell_to_node_[istart]   = node_index_(ix,iy,iz);
        cell_to_node_[istart+1] = node_index_(ix+1,iy,iz);
        cell_to_node_[istart+2] = node_index_(ix+1,iy+1,iz);
        cell_to_node_[istart+3] = node_index_(ix,iy+1,iz);
        cell_to_node_[istart+4] = node_index_(ix,iy,iz+1);
        cell_to_node_[istart+5] = node_index_(ix+1,iy,iz+1);
        cell_to_node_[istart+6] = node_index_(ix+1,iy+1,iz+1);
        cell_to_node_[istart+7] = node_index_(ix,iy+1,iz+1);

        jstart = 9 * node_index_(ix,iy,iz);
        ncell = node_to_cell_[jstart];
        node_to_cell_[jstart+1+ncell] = cell_index_(ix,iy,iz);
        (node_to_cell_[jstart])++;

        jstart = 9 * node_index_(ix+1,iy,iz); 
        ncell = node_to_cell_[jstart];
        node_to_cell_[jstart+1+ncell] = cell_index_(ix,iy,iz);
        (node_to_cell_[jstart])++;

        jstart = 9 * node_index_(ix+1,iy+1,iz); 
        ncell = node_to_cell_[jstart];
        node_to_cell_[jstart+1+ncell] = cell_index_(ix,iy,iz);
        (node_to_cell_[jstart])++;

        jstart = 9 * node_index_(ix,iy+1,iz); 
        ncell = node_to_cell_[jstart];
        node_to_cell_[jstart+1+ncell] = cell_index_(ix,iy,iz);
        (node_to_cell_[jstart])++;

        jstart = 9 * node_index_(ix,iy,iz+1);
        ncell = node_to_cell_[jstart];
        node_to_cell_[jstart+1+ncell] = cell_index_(ix,iy,iz);
        node_to_cell_[jstart]++;

        jstart = 9 * node_index_(ix+1,iy,iz+1); // 1 extra for num cells
        ncell = node_to_cell_[jstart];
        node_to_cell_[jstart+1+ncell] = cell_index_(ix,iy,iz);
        (node_to_cell_[jstart])++;

        jstart = 9 * node_index_(ix+1,iy+1,iz+1); // 1 extra for num cells
        ncell = node_to_cell_[jstart];
        node_to_cell_[jstart+1+ncell] = cell_index_(ix,iy,iz);
        (node_to_cell_[jstart])++;

        jstart = 9 * node_index_(ix,iy+1,iz+1); // 1 extra for num cells
        ncell = node_to_cell_[jstart];
        node_to_cell_[jstart+1+ncell] = cell_index_(ix,iy,iz);
        (node_to_cell_[jstart])++;

      }


  // loop over cells and initialize cell_to_face_
  for (int iz=0; iz<nz_; iz++)
    for (int iy=0; iy<ny_; iy++)
      for (int ix=0; ix<nx_; ix++)
      {
        int istart = 6 * cell_index_(ix,iy,iz);
        int jstart = 0;

        cell_to_face_[istart]    = xzface_index_(ix,iy,iz);
        cell_to_face_[istart+1]  = yzface_index_(ix+1,iy,iz);
        cell_to_face_[istart+2]  = xzface_index_(ix,iy+1,iz);
        cell_to_face_[istart+3]  = yzface_index_(ix,iy,iz);
        cell_to_face_[istart+4]  = xyface_index_(ix,iy,iz);
        cell_to_face_[istart+5]  = xyface_index_(ix,iy,iz+1);

        cell_to_face_dirs_[istart]   = 1;
        cell_to_face_dirs_[istart+1] = 1;
        cell_to_face_dirs_[istart+2] = -1;
        cell_to_face_dirs_[istart+3] = -1;
        cell_to_face_dirs_[istart+4] = -1;
        cell_to_face_dirs_[istart+5] = 1;

        jstart = 2*xzface_index_(ix,iy,iz);
        face_to_cell_[jstart+1] = cell_index_(ix,iy,iz);

        jstart = 2*yzface_index_(ix+1,iy,iz);
        face_to_cell_[jstart+1] = cell_index_(ix,iy,iz);

        jstart = 2*xzface_index_(ix,iy+1,iz);
        face_to_cell_[jstart] = cell_index_(ix,iy,iz);

        jstart = 2*yzface_index_(ix,iy,iz);
        face_to_cell_[jstart] = cell_index_(ix,iy,iz);

        jstart = 2*xyface_index_(ix,iy,iz);
        face_to_cell_[jstart] = cell_index_(ix,iy,iz);

        jstart = 2*xyface_index_(ix,iy,iz+1);
        face_to_cell_[jstart+1] = cell_index_(ix,iy,iz);
      }


  // loop over faces and initialize face_to_node_
  // first we do the xy faces
  for (int iz=0; iz<=nz_; iz++)
    for (int iy=0; iy<ny_; iy++)
      for (int ix=0; ix<nx_; ix++)  
      {
        int istart = 4 * xyface_index_(ix,iy,iz);
        int jstart = 0;
        int nfaces = 0;

        face_to_node_[istart]   = node_index_(ix,iy,iz);
        face_to_node_[istart+1] = node_index_(ix+1,iy,iz);
        face_to_node_[istart+2] = node_index_(ix+1,iy+1,iz);
        face_to_node_[istart+3] = node_index_(ix,iy+1,iz);

        jstart = 13*node_index_(ix,iy,iz);
        nfaces = node_to_face_[jstart];
        node_to_face_[jstart+1+nfaces] = xyface_index_(ix,iy,iz);
        (node_to_face_[jstart])++;

        jstart = 13*node_index_(ix+1,iy,iz);
        nfaces = node_to_face_[jstart];
        node_to_face_[jstart+1+nfaces] = xyface_index_(ix,iy,iz);
        (node_to_face_[jstart])++;

        jstart = 13*node_index_(ix+1,iy+1,iz);
        nfaces = node_to_face_[jstart];
        node_to_face_[jstart+1+nfaces] = xyface_index_(ix,iy,iz);
        (node_to_face_[jstart])++;

        jstart = 13*node_index_(ix,iy+1,iz);
        nfaces = node_to_face_[jstart];
        node_to_face_[jstart+1+nfaces] = xyface_index_(ix,iy,iz);
        (node_to_face_[jstart])++;

      }
  // then we do the xz faces
  for (int iz=0; iz<nz_; iz++)
    for (int iy=0; iy<=ny_; iy++)
      for (int ix=0; ix<nx_; ix++)  
      {
        int istart = 4 * xzface_index_(ix,iy,iz);
        int jstart = 0;
        int nfaces = 0;

        face_to_node_[istart]   = node_index_(ix,iy,iz);
        face_to_node_[istart+1] = node_index_(ix+1,iy,iz);
        face_to_node_[istart+2] = node_index_(ix+1,iy,iz+1);
        face_to_node_[istart+3] = node_index_(ix,iy,iz+1);

        jstart = 13*node_index_(ix,iy,iz);
        nfaces = node_to_face_[jstart];
        node_to_face_[jstart+1+nfaces] = xyface_index_(ix,iy,iz);
        (node_to_face_[jstart])++;

        jstart = 13*node_index_(ix+1,iy,iz);
        nfaces = node_to_face_[jstart];
        node_to_face_[jstart+1+nfaces] = xyface_index_(ix,iy,iz);
        (node_to_face_[jstart])++;

        jstart = 13*node_index_(ix+1,iy,iz+1);
        nfaces = node_to_face_[jstart];
        node_to_face_[jstart+1+nfaces] = xyface_index_(ix,iy,iz);
        (node_to_face_[jstart])++;

        jstart = 13*node_index_(ix,iy,iz+1);
        nfaces = node_to_face_[jstart];
        node_to_face_[jstart+1+nfaces] = xyface_index_(ix,iy,iz);
        (node_to_face_[jstart])++;

      }
  // finally we do the yz faces
  for (int iz=0; iz<nz_; iz++)
    for (int iy=0; iy<ny_; iy++)
      for (int ix=0; ix<=nx_; ix++)  
      {
        int istart = 4 * yzface_index_(ix,iy,iz);
        int jstart = 0;
        int nfaces = 0;

        face_to_node_[istart]   = node_index_(ix,iy,iz);
        face_to_node_[istart+1] = node_index_(ix,iy+1,iz);
        face_to_node_[istart+2] = node_index_(ix,iy+1,iz+1);
        face_to_node_[istart+3] = node_index_(ix,iy,iz+1);

        jstart = 13*node_index_(ix,iy,iz);
        nfaces = node_to_face_[jstart];
        node_to_face_[jstart+1+nfaces] = xyface_index_(ix,iy,iz);
        (node_to_face_[jstart])++;

        jstart = 13*node_index_(ix,iy+1,iz);
        nfaces = node_to_face_[jstart];
        node_to_face_[jstart+1+nfaces] = xyface_index_(ix,iy,iz);
        (node_to_face_[jstart])++;

        jstart = 13*node_index_(ix,iy+1,iz+1);
        nfaces = node_to_face_[jstart];
        node_to_face_[jstart+1+nfaces] = xyface_index_(ix,iy,iz);
        (node_to_face_[jstart])++;

        jstart = 13*node_index_(ix,iy,iz+1);
        nfaces = node_to_face_[jstart];
        node_to_face_[jstart+1+nfaces] = xyface_index_(ix,iy,iz);
        (node_to_face_[jstart])++;

      }


  // we only have 6 side sets

  side_sets_.resize(6);
   
  side_sets_[0].resize( nx_ * nz_ );
  side_sets_[1].resize( ny_ * nz_ );
  side_sets_[2].resize( nx_ * nz_ );
  side_sets_[3].resize( ny_ * nz_ );
  side_sets_[4].resize( nx_ * ny_ );
  side_sets_[5].resize( nx_ * ny_ );

  int count = 0;
  for (int ix=0; ix<nx_; ix++)
    for (int iz=0; iz<nz_; iz++) 
    {
      side_sets_[0][count] = xzface_index_(ix,0,iz);
      side_sets_[2][count] = xzface_index_(ix,ny_,iz);
      count ++;
    }
   
  count = 0;
  for (int iy=0; iy<ny_; iy++)
    for (int iz=0; iz<nz_; iz++) 
    {
      side_sets_[1][count] = yzface_index_(nx_,iy,iz);
      side_sets_[3][count] = yzface_index_(0,iy,iz);
      count++;
    }
   
  count = 0;
  for (int ix=0; ix<nx_; ix++)
    for (int iy=0; iy<ny_; iy++) 
    {
      side_sets_[4][count] = xyface_index_(ix,iy,0);
      side_sets_[5][count] = xyface_index_(ix,iy,nz_);
      count++;
    }
      


  // create element blocks

   
  if (mesh_blocks_.empty()) 
  {
    // we only have one element block
       
    element_blocks_.resize(1);
    element_blocks_[0].resize(nx_*ny_*nz_);
       
    count=0;
    for (int iz=0; iz<nz_; iz++)
      for (int iy=0; iy<ny_; iy++)
        for (int ix=0; ix<nx_; ix++) 
        {
          element_blocks_[0][count] = cell_index_(ix,iy,iz);
          count++;
        }
    
  } 
  else
  {
    // make 1 extra block (w/ id = 0) for those cells that don't
    // fall into the declared blocks

    element_blocks_.resize(mesh_blocks_.size() + 1);

    // The first pass just counts

    std::vector<int> count(element_blocks_.size(), 0);
    for (int pass = 0; pass < 2; ++pass) {
      for (int ic=0; ic<num_cells_; ic++) {
        std::vector<AmanziGeometry::Point> coords(8);
        cell_get_coordinates(ic,&coords);
        AmanziGeometry::Point centroid(0.0, 0.0, 0.0);
        for (std::vector<AmanziGeometry::Point>::iterator c = coords.begin();
             c != coords.end(); ++c) {
          centroid += *c;
        }
        centroid /= 8.0;

        int nb(0);
        AmanziGeometry::RegionVector::const_iterator r;
        for (r = mesh_blocks_.begin(); r != mesh_blocks_.end(); ++r, nb++) {
          bool inside = (*r)->inside(centroid);
          if (inside) break;
        }

        if (nb >= mesh_blocks_.size()) {
          nb = 0;
        } else {
          nb = nb + 1;
        }

        switch (pass) {
          case 0:
            count[nb] += 1;
            break;
          case 1:
            element_blocks_[nb].push_back(ic);
            break;
        }
      }

      switch (pass) {
        case 0:
          for (int nb = 0; nb < element_blocks_.size(); ++nb) {
            element_blocks_[nb].reserve(count[nb]);
          }
          break;
        default:
          // do nothing
          break;
      }
    }

    // check that all cells have been added to an element block
    int ncb = 0;
    for (int nb=0; nb< element_blocks_.size(); nb++) 
    {
      ncb += element_blocks_[nb].size();
       
      std::cout << "Mesh_simple: mesh block " << nb << " has " << element_blocks_[nb].size() << " cells" << std::endl;
    }


    if (ncb < num_cells_) 
    {
      std::cout << "Mesh_simple: WARNING... not all cells are in a mesh block\n";
      std::cout << "  number of cells in mesh blocks = " << ncb << std::endl;
      std::cout << "  number of cells in total       = " << num_cells_ << std::endl;
    }
    if (ncb > num_cells_)
      std::cout << "Mesh_simple: WARNING... overlapping mesh blocks" << std::endl;
       
       
  }
  build_maps_ ();
}


void Mesh_simple::build_maps_ ()
{
  std::vector<int> cells( num_cells_ );
  for (int i=0; i< num_cells_; i++) cells[i] = i;
  
  std::vector<int> nodes( num_nodes_ );
  for (int i=0; i< num_nodes_; i++) nodes[i] = i;

  std::vector<int> faces( num_faces_ );
  for (int i=0; i< num_faces_; i++) faces[i] = i;

  
  cell_map_ = new Epetra_Map(-1, num_cells_, &cells[0], 0, *communicator_ );
  face_map_ = new Epetra_Map(-1, num_faces_, &faces[0], 0, *communicator_ );
  node_map_ = new Epetra_Map(-1, num_nodes_, &nodes[0], 0, *communicator_ );

}





Parallel_type Mesh_simple::entity_get_ptype(const Entity_kind kind, 
					    const Entity_ID entid) const 
{
  return OWNED; // Its a serial code
}




// Get cell type
    
AmanziMesh::Cell_type Mesh_simple::cell_get_type(const AmanziMesh::Entity_ID cellid) const 
{
  return HEX;
}
        
    
unsigned int Mesh_simple::GID(const AmanziMesh::Entity_ID lid, 
			      const AmanziMesh::Entity_kind kind) const
{
  return lid;  // Its a serial code
}



unsigned int Mesh_simple::num_entities (AmanziMesh::Entity_kind kind, 
					AmanziMesh::Parallel_type ptype) const
{
  switch (kind) {
    case AmanziMesh::FACE: 
      return num_faces_;
      break;
    case AmanziMesh::NODE:
      return num_nodes_;
      break;
    case AmanziMesh::CELL:
      return num_cells_;
      break;
    default:
      throw std::exception();
      break;
  }
}


unsigned int Mesh_simple::num_sets(AmanziMesh::Entity_kind kind) const
{
  switch (kind) {
    case AmanziMesh::FACE:
      return side_sets_.size();
      break;
    case AmanziMesh::CELL:
      return element_blocks_.size();
      break;
    default:
      return 0;
      break;
  }
}

unsigned int Mesh_simple::get_set_size (AmanziMesh::Set_ID set_id, 
					AmanziMesh::Entity_kind kind,
					AmanziMesh::Parallel_type ptype) const
{
  // we ignore ptype, since this is a serial implementation

  switch (kind) {
    case AmanziMesh::FACE:
      return side_sets_[set_id].size();
      break;
    case AmanziMesh::CELL:
      return element_blocks_[set_id].size();
      break;
    default:
      throw std::exception();
      break;
  }
}


void Mesh_simple::cell_get_faces (AmanziMesh::Entity_ID cell, 
				  AmanziMesh::Entity_ID_List *faceids)  const 
{
  unsigned int offset = 6*cell;

  faceids->clear();

  for (int i = 0; i < 6; i++) {
    faceids->push_back(*(cell_to_face_.begin()+offset));
    offset++;
  }
}



void Mesh_simple::cell_get_nodes (AmanziMesh::Entity_ID cell, 
				  AmanziMesh::Entity_ID_List *nodeids) const
{
  unsigned int offset = 8*cell;

  nodeids->clear();

  for (int i = 0; i < 8; i++) {
    nodeids->push_back(*(cell_to_node_.begin()+offset));
    offset++;
  }
}


void Mesh_simple::face_get_nodes (AmanziMesh::Entity_ID face, 
				  AmanziMesh::Entity_ID_List *nodeids) const
{
  unsigned int offset = 4*face;

  nodeids->clear();

  for (int i = 0; i < 4; i++) {
    nodeids->push_back(*(face_to_node_.begin()+offset));
    offset++;
  }
}


// Cooordinate Getters
// -------------------

void Mesh_simple::node_get_coordinates (const AmanziMesh::Entity_ID local_node_id, 
					AmanziGeometry::Point *ncoords) const
{
  unsigned int offset = 3*local_node_id;

  ncoords->init(3);
  
  ncoords->set( &(coordinates_[offset]) );
}


void Mesh_simple::face_get_coordinates (AmanziMesh::Entity_ID local_face_id, 
					std::vector<AmanziGeometry::Point> *fcoords) const
{
  Entity_ID_List node_indices(4);

  face_get_nodes (local_face_id, &node_indices);

  fcoords->clear();

  AmanziGeometry::Point xyz(3);
  for (std::vector<unsigned int>::iterator it = node_indices.begin(); 
       it != node_indices.end(); ++it)
  {
    node_get_coordinates ( *it, &xyz);
    fcoords->push_back(xyz);
  }

}

void Mesh_simple::cell_get_coordinates (AmanziMesh::Entity_ID local_cell_id, 
					std::vector<AmanziGeometry::Point> *ccoords) const
{  
  std::vector<unsigned int> node_indices(8);

  cell_get_nodes (local_cell_id, &node_indices);

  ccoords->clear();

  AmanziGeometry::Point xyz(3);
  for (std::vector<unsigned int>::iterator it = node_indices.begin(); 
       it != node_indices.end(); ++it)
  {      
    node_get_coordinates ( *it, &xyz);
    ccoords->push_back(xyz);
  }
}



// Set getters
// -----------

void Mesh_simple::get_set_ids (AmanziMesh::Entity_kind kind, 
			       AmanziMesh::Set_ID_List *setids) const
{
  
  setids->clear();

  switch (kind) {
    case AmanziMesh::FACE: 
      {
        for (int i=0; i<6; i++)
          setids->push_back(i);	
        break;
      }
    case AmanziMesh::CELL:
      {
        for (int i = 0; i < element_blocks_.size(); ++i)
          setids->push_back(i);
        break;
      }      
    default:
      // we do not have anything for NODE, yet
      return;
  }
}


void Mesh_simple::get_set_entities (AmanziMesh::Set_ID set_id, 
				    AmanziMesh::Entity_kind kind, 
				    AmanziMesh::Parallel_type ptype, 
				    AmanziMesh::Set_ID_List *setents) const
{
  // we ignore ptype since this is a serial implementation

  setents->clear();
  
  switch (kind) {
    case AmanziMesh::FACE:
      *setents = side_sets_[set_id];
      break;
    case AmanziMesh::CELL:
      *setents = element_blocks_[set_id];
      break;
    default:
      // we do not have anything for NODE, yet
      throw std::exception();
  }
}


bool Mesh_simple::valid_set_id (AmanziMesh::Entity_ID id, 
				AmanziMesh::Entity_kind kind) const 
{
  switch (kind) {
    case AmanziMesh::FACE:
      return (id<6) ;
      break;
    case AmanziMesh::CELL:
      return (id < element_blocks_.size());
      break;
    default:
      // we do not have anything for NODE, yet
      return false;
  }
    
}

void Mesh_simple::cell_get_face_dirs (AmanziMesh::Entity_ID cell, 
				      std::vector<int> *cfacedirs) const
{
  unsigned int offset = 6*cell;

  cfacedirs->clear();
  for (int i = 0; i < 6; i++) {
    cfacedirs->push_back(*(cell_to_face_dirs_.begin()+offset));
    offset++;
  }
}


void Mesh_simple::set_coordinate(AmanziMesh::Entity_ID local_node_id, 
				 double *source_begin, 
				 double *source_end)
{
  unsigned int offset = 3*local_node_id;
  
  std::vector<double>::iterator destination_begin = coordinates_.begin() + offset;
  std::copy(source_begin, source_end, destination_begin);

}



void Mesh_simple::node_get_cells (const AmanziMesh::Entity_ID nodeid, 
                                  const AmanziMesh::Parallel_type ptype,
                                  AmanziMesh::Entity_ID_List *cellids) const 
{
  unsigned int offset = 9*nodeid;
  unsigned int ncells = node_to_cell_[offset];

  cellids->clear();

  for (int i = 0; i < ncells; i++) 
    cellids->push_back(node_to_cell_[offset+i]);
}
    


// Faces of type 'ptype' connected to a node
    
void Mesh_simple::node_get_faces (const AmanziMesh::Entity_ID nodeid, 
                                  const AmanziMesh::Parallel_type ptype,
                                  AmanziMesh::Entity_ID_List *faceids) const
{
  unsigned int offset = 13*nodeid;
  unsigned int nfaces = node_to_face_[offset];

  faceids->clear();

  for (int i = 0; i < nfaces; i++) 
    faceids->push_back(node_to_face_[offset+i]);
}   


 
// Get faces of ptype of a particular cell that are connected to the
// given node

void Mesh_simple::node_get_cell_faces (const AmanziMesh::Entity_ID nodeid, 
                                       const AmanziMesh::Entity_ID cellid,
                                       const AmanziMesh::Parallel_type ptype,
                                       AmanziMesh::Entity_ID_List *faceids) const
{
  unsigned int offset = 6*cellid;

  faceids->clear();

  for (int i = 0; i < 6; i++) {
    unsigned int cellfaceid = face_to_cell_[offset+i];

    unsigned int offset2 = 4*cellfaceid;
    
    for (int j = 0; j < 4; j++) {
      if (face_to_node_[offset2+j] == nodeid) {
	faceids->push_back(cellfaceid);
	break;
      }
    }
  }

}
    
// Cells connected to a face
    
void Mesh_simple::face_get_cells (const AmanziMesh::Entity_ID faceid, 
                                  const AmanziMesh::Parallel_type ptype,
                                  AmanziMesh::Entity_ID_List *cellids) const
{
  unsigned int offset = 2*faceid;

  cellids->clear();

  if (face_to_cell_[offset] != -1)
    cellids->push_back(face_to_cell_[offset]);
  if (face_to_cell_[offset+1] != -1)
    cellids->push_back(face_to_cell_[offset+1]);
}
    


// Same level adjacencies
//-----------------------

// Face connected neighboring cells of given cell of a particular ptype
// (e.g. a hex has 6 face neighbors)

// The order in which the cellids are returned cannot be
// guaranteed in general except when ptype = USED, in which case
// the cellids will correcpond to cells across the respective
// faces given by cell_get_faces

void Mesh_simple::cell_get_face_adj_cells(const AmanziMesh::Entity_ID cellid,
                                          const AmanziMesh::Parallel_type ptype,
                                          AmanziMesh::Entity_ID_List *fadj_cellids) const
{
  unsigned int offset = 6*cellid;

  fadj_cellids->clear();

  for (int i = 0; i < 6; i++) {    
    unsigned int faceid = cell_to_face_[offset];

    unsigned int foffset = 2*faceid;

    unsigned int adjcell0 = face_to_cell_[foffset];
    if (adjcell0 != -1 && adjcell0 != cellid)
      fadj_cellids->push_back(adjcell0);    
    else {
      unsigned int adjcell1 = face_to_cell_[foffset+1];
      if (adjcell1 != -1 && adjcell1 != cellid)
	fadj_cellids->push_back(adjcell1);
    }

    offset++;
  }
}

// Node connected neighboring cells of given cell
// (a hex in a structured mesh has 26 node connected neighbors)
// The cells are returned in no particular order

void Mesh_simple::cell_get_node_adj_cells(const AmanziMesh::Entity_ID cellid,
                                          const AmanziMesh::Parallel_type ptype,
                                          AmanziMesh::Entity_ID_List *nadj_cellids) const
{
  unsigned int offset = 8*cellid;

  nadj_cellids->clear();
  
  for (int i = 0; i < 8; i++) {
    unsigned int nodeid = cell_to_node_[offset+i];

    unsigned int offset2 = 9*nodeid;
    unsigned int ncell = node_to_cell_[offset2];

    for (int j = 0; j < 8; j++) {
      unsigned int nodecell = node_to_cell_[offset2+j];
      
      unsigned int found = 0;
      unsigned int size = nadj_cellids->size();
      for (int k = 0; k < size; k++) {
	if ((*nadj_cellids)[k] == nodecell) {
	  found = 1;
	  break;
	}
      }

      if (!found)
	nadj_cellids->push_back(nodecell);
    }
  }

}


    
//
// Mesh Topology for viz  
//----------------------
//
// We need a special function because certain types of degenerate
// hexes will not be recognized as any standard element type (hex,
// pyramid, prism or tet). The original topology of this element 
// without any collapsed nodes will be returned by this call.


// Original cell type 

AmanziMesh::Cell_type Mesh_simple::cell_get_type_4viz(const AmanziMesh::Entity_ID cellid) const 
{
  return AmanziMesh::HEX;
}
    
    
// See cell_get_nodes for details on node ordering
    
void Mesh_simple::cell_get_nodes_4viz (const AmanziMesh::Entity_ID cellid, 
                                       AmanziMesh::Entity_ID_List *nodeids) const 
{
  cell_get_nodes(cellid, nodeids);
}
    
    
const Epetra_Map& Mesh_simple::cell_epetra_map (bool include_ghost) const
{
  return *cell_map_;
}

const Epetra_Map& Mesh_simple::face_epetra_map (bool include_ghost) const
{
  return *face_map_;
}

const Epetra_Map& Mesh_simple::node_epetra_map (bool include_ghost) const
{
  return *node_map_;
}

} // close namespace AmanziMesh
} // close namespace Amanzi
