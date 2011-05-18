#include "Mesh_simple.hh"
#include <Teuchos_RCP.hpp>


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
  number_of_mesh_blocks=0;
  update();
}

Mesh_simple::Mesh_simple ( Teuchos::ParameterList &parameter_list,
				     Epetra_Comm *communicator ) :
  communicator_(communicator)
{

  // read the parameters from the parameter list

  nx_ = parameter_list.get<int>("Numer of Cells in X");
  ny_ = parameter_list.get<int>("Numer of Cells in Y");
  nz_ = parameter_list.get<int>("Numer of Cells in Z");
  
  x0_ = parameter_list.get<double>("X_Min");
  x1_ = parameter_list.get<double>("X_Max");

  y0_ = parameter_list.get<double>("Y_Min");
  y1_ = parameter_list.get<double>("Y_Max");

  z0_ = parameter_list.get<double>("Z_Min");
  z1_ = parameter_list.get<double>("Z_Max");


  number_of_mesh_blocks = parameter_list.get<int>("Number of mesh blocks",0);
  
  
  if (number_of_mesh_blocks > 0) {
    mesh_block_z0 = new double [number_of_mesh_blocks];
    mesh_block_z1 = new double [number_of_mesh_blocks];
    
    for (int nb=1; nb<=number_of_mesh_blocks; nb++) {
      std::stringstream s; 
      s << "Mesh block " << nb;

      Teuchos::ParameterList sublist = parameter_list.sublist(s.str());
      
      mesh_block_z0[nb-1] = sublist.get<double>("Z0");
      mesh_block_z1[nb-1] = sublist.get<double>("Z1");
    }
  }

  update();
}







Mesh_simple::~Mesh_simple()
{
  delete cell_map_;
  delete face_map_;
  delete node_map_;
  if (!mesh_block_z0)  delete [] mesh_block_z0;
  if (!mesh_block_z1)  delete [] mesh_block_z1;
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

   // loop over cells and initialize cell_to_node_
   for (int iz=0; iz<nz_; iz++)
     for (int iy=0; iy<ny_; iy++)
       for (int ix=0; ix<nx_; ix++)
   	 {
   	   int istart = 8 * cell_index_(ix,iy,iz);

   	   cell_to_node_[istart]   = node_index_(ix,iy,iz);
   	   cell_to_node_[istart+1] = node_index_(ix+1,iy,iz);
   	   cell_to_node_[istart+2] = node_index_(ix+1,iy+1,iz);
   	   cell_to_node_[istart+3] = node_index_(ix,iy+1,iz);
   	   cell_to_node_[istart+4] = node_index_(ix,iy,iz+1);
   	   cell_to_node_[istart+5] = node_index_(ix+1,iy,iz+1);
   	   cell_to_node_[istart+6] = node_index_(ix+1,iy+1,iz+1);
   	   cell_to_node_[istart+7] = node_index_(ix,iy+1,iz+1);
   	 }


   // loop over cells and initialize cell_to_face_
   for (int iz=0; iz<nz_; iz++)
     for (int iy=0; iy<ny_; iy++)
       for (int ix=0; ix<nx_; ix++)
   	 {
   	   int istart = 6 * cell_index_(ix,iy,iz);

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
	   
   	 }


   // loop over faces and initialize face_to_node_
   // first we do the xy faces
   for (int iz=0; iz<=nz_; iz++)
     for (int iy=0; iy<ny_; iy++)
       for (int ix=0; ix<nx_; ix++)  
   	 {
   	   int istart = 4 * xyface_index_(ix,iy,iz);

   	   face_to_node_[istart]   = node_index_(ix,iy,iz);
   	   face_to_node_[istart+1] = node_index_(ix+1,iy,iz);
   	   face_to_node_[istart+2] = node_index_(ix+1,iy+1,iz);
   	   face_to_node_[istart+3] = node_index_(ix,iy+1,iz);
   	 }
   // then we do the xz faces
   for (int iz=0; iz<nz_; iz++)
     for (int iy=0; iy<=ny_; iy++)
       for (int ix=0; ix<nx_; ix++)  
   	 {
   	   int istart = 4 * xzface_index_(ix,iy,iz);

   	   face_to_node_[istart]   = node_index_(ix,iy,iz);
   	   face_to_node_[istart+1] = node_index_(ix+1,iy,iz);
   	   face_to_node_[istart+2] = node_index_(ix+1,iy,iz+1);
   	   face_to_node_[istart+3] = node_index_(ix,iy,iz+1);
   	 }
   // finally we do the yz faces
   for (int iz=0; iz<nz_; iz++)
     for (int iy=0; iy<ny_; iy++)
       for (int ix=0; ix<=nx_; ix++)  
   	 {
   	   int istart = 4 * yzface_index_(ix,iy,iz);

   	   face_to_node_[istart]   = node_index_(ix,iy,iz);
   	   face_to_node_[istart+1] = node_index_(ix,iy+1,iz);
   	   face_to_node_[istart+2] = node_index_(ix,iy+1,iz+1);
   	   face_to_node_[istart+3] = node_index_(ix,iy,iz+1);
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

   
   if (number_of_mesh_blocks == 0) 
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
       
       number_of_mesh_blocks = 1;

     } 
   else
     {
       element_blocks_.resize(number_of_mesh_blocks);

       for (int nb=0; nb< number_of_mesh_blocks; nb++) 
	 {
	   
	   // count the nunber of cells in mesh block nb
	   count = 0;
	   for (int ic=0; ic<num_cells_; ic++) 
	     {
	       std::vector<Point> coords(8);
	       cell_get_coordinates(ic,&coords);
	       
	       // first check if all z-coordinates are larger than
	       // mesh_block_z0[nb]
	       bool inside = true;
	       for (int i=4; i<8; i++) 
		 if ( coords[i][2] <= mesh_block_z0[nb] ) inside = false;
	       for (int i=0; i<8; i++) 
		 if ( coords[i][2] > mesh_block_z1[nb] ) inside = false;

	       if (inside) count++;
	       
	     }
	   element_blocks_[nb].resize(count);
	   
	   // add cells to mesh block nb
	   count = 0;
	   for (int ic=0; ic<num_cells_; ic++) 
	     {
	       std::vector<Point> coords;
	       cell_get_coordinates(ic,&coords);
	       
	       // first check if all z-coordinates are larger than
	       // mesh_block_z0[nb]
	       bool inside = true;
	       for (int i=4; i<8; i++) 
		 if ( coords[i][2] <= mesh_block_z0[nb] ) inside = false;
	       for (int i=0; i<8; i++) 
		 if ( coords[3][2] > mesh_block_z1[nb] ) inside = false;
	       
	       if (inside) {
		 element_blocks_[nb][count] = ic;
		 count++;
	       }
	       
	     }	   

	   
	 }
       
       // check that all cells have been added to an element block
       int ncb = 0;
       for (int nb=0; nb< number_of_mesh_blocks; nb++) 
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
    
Amanzi::AmanziMesh::Cell_type Mesh_simple::cell_get_type(const Amanzi::AmanziMesh::Entity_ID cellid) const 
{
  return HEX;
}
        
    
unsigned int Mesh_simple::GID(const Amanzi::AmanziMesh::Entity_ID lid, 
			      const Amanzi::AmanziMesh::Entity_kind kind) const
{
  return lid;  // Its a serial code
}



unsigned int Mesh_simple::num_entities (Amanzi::AmanziMesh::Entity_kind kind, 
					Amanzi::AmanziMesh::Parallel_type ptype) const
{
  switch (kind) {
  case Amanzi::AmanziMesh::FACE: 
    return num_faces_;
    break;
  case Amanzi::AmanziMesh::NODE:
    return num_nodes_;
    break;
  case Amanzi::AmanziMesh::CELL:
    return num_cells_;
    break;
  default:
    throw std::exception();
    break;
  }
}


unsigned int Mesh_simple::num_sets(Amanzi::AmanziMesh::Entity_kind kind) const
{
  switch (kind) {
  case Amanzi::AmanziMesh::FACE:
    return 6;
    break;
  case Amanzi::AmanziMesh::CELL:
    return 1;
    break;
  default:
    // nothing yet for NODE or CELL
    throw std::exception();
    break;
  }
}

unsigned int Mesh_simple::get_set_size (Amanzi::AmanziMesh::Set_ID set_id, 
					Amanzi::AmanziMesh::Entity_kind kind,
					Amanzi::AmanziMesh::Parallel_type ptype) const
{
  // we ignore ptype, since this is a serial implementation

  switch (kind) {
  case Amanzi::AmanziMesh::FACE:
    return side_sets_[set_id].size();
    break;
  case Amanzi::AmanziMesh::CELL:
    return element_blocks_[set_id].size();
    break;
  default:
    throw std::exception();
    break;
  }
}


void Mesh_simple::cell_get_faces (Amanzi::AmanziMesh::Entity_ID cell, 
				  Amanzi::AmanziMesh::Entity_ID_List *faceids)  const 
{
  unsigned int index = 6*cell;

  for (int i = 0; i < 6; i++) {
    faceids->push_back(*(cell_to_face_.begin()+index));
    index++;
  }
}



void Mesh_simple::cell_get_nodes (Amanzi::AmanziMesh::Entity_ID cell, 
				  Amanzi::AmanziMesh::Entity_ID_List *nodeids) const
{
  unsigned int index = 8*cell;

  for (int i = 0; i < 8; i++) {
    nodeids->push_back(*(cell_to_node_.begin()+index));
    index++;
  }
}


void Mesh_simple::face_get_nodes (Amanzi::AmanziMesh::Entity_ID face, 
				  Amanzi::AmanziMesh::Entity_ID_List *nodeids) const
{
  unsigned int index = 4*face;

  for (int i = 0; i < 4; i++) {
    nodeids->push_back(*(face_to_node_.begin()+index));
    index++;
  }
}


// Cooordinate Getters
// -------------------

void Mesh_simple::node_get_coordinates (Amanzi::AmanziMesh::Entity_ID local_node_id, 
					Amanzi::AmanziGeometry::Point *ncoords) const
{
  unsigned int index = 3*local_node_id;
  std::vector<double>::iterator begin = coordinates_.begin() + index;

  ncoords->init(3);
  
  ncoords->set( *(begin), *(begin+1), *(begin+2));
}


void Mesh_simple::face_get_coordinates (Amanzi::AmanziMesh::Entity_ID local_face_id, 
					std::vector<Amanzi::AmanziGeometry::Point> *fcoords) const
{
  Entity_ID_List node_indices(4);

  face_get_nodes (local_face_id, &node_indices);

  Point xyz(3);
  for (std::vector<unsigned int>::iterator it = node_indices.begin(); 
       it != node_indices.end(); ++it)
    {
      node_get_coordinates ( *it, &xyz);
      fcoords->push_back(xyz);
    }

}

void Mesh_simple::cell_get_coordinates (Amanzi::AmanziMesh::Entity_ID local_cell_id, 
					std::vector<Amanzi::AmanziGeometry::Point> *ccoords) const
{  
  std::vector<unsigned int> node_indices(8);
  cell_get_nodes (local_cell_id, &node_indices);

  Point xyz(3);
  for (std::vector<unsigned int>::iterator it = node_indices.begin(); 
       it != node_indices.end(); ++it)
    {      
      node_get_coordinates ( *it, &xyz);
      ccoords->push_back(xyz);
    }
}



// Set getters
// -----------

void Mesh_simple::get_set_ids (Amanzi::AmanziMesh::Entity_kind kind, 
			       Amanzi::AmanziMesh::Set_ID_List *setids) const
{
  
  std::vector<int> ids(0);

  switch (kind) {
  case Amanzi::AmanziMesh::FACE: 
    {
      for (int i=0; i<6; i++)
	setids->push_back(i);	
      break;
    }
  case Amanzi::AmanziMesh::CELL:
    {
      setids->push_back(0);
      break;
    }      
  default:
    // we do not have anything for NODE, yet
    throw std::exception();
  }
}


void Mesh_simple::get_set_entities (Amanzi::AmanziMesh::Set_ID set_id, 
				    Amanzi::AmanziMesh::Entity_kind kind, 
				    Amanzi::AmanziMesh::Parallel_type ptype, 
				    Amanzi::AmanziMesh::Set_ID_List *setents) const
{
  // we ignore ptype since this is a serial implementation
  
  switch (kind) {
  case Amanzi::AmanziMesh::FACE:
    *setents = side_sets_[set_id];
    break;
  case Amanzi::AmanziMesh::CELL:
    *setents = element_blocks_[set_id];
    break;
  default:
    // we do not have anything for NODE, yet
    throw std::exception();
  }
}


bool Mesh_simple::valid_set_id (Amanzi::AmanziMesh::Entity_ID id, 
				Amanzi::AmanziMesh::Entity_kind kind) const 
{
  switch (kind) {
  case Amanzi::AmanziMesh::FACE:
    return (id<6) ;
    break;
  case Amanzi::AmanziMesh::CELL:
    return (id < number_of_mesh_blocks);
    break;
  default:
    // we do not have anything for NODE, yet
    throw std::exception();  
  }
    
}

void Mesh_simple::cell_get_face_dirs (Amanzi::AmanziMesh::Entity_ID cell, 
				      std::vector<int> *cfacedirs) const
{
  unsigned int index = 6*cell;
  for (int i = 0; i < 6; i++) {
    cfacedirs->push_back(*(cell_to_face_dirs_.begin()+index));
    index++;
  }
}


void Mesh_simple::set_coordinate(Amanzi::AmanziMesh::Entity_ID local_node_id, 
				 double *source_begin, 
				 double *source_end) const
{
  unsigned int index = 3*local_node_id;
  
  std::vector<double>::iterator destination_begin = coordinates_.begin() + index;
  std::copy(source_begin, source_end, destination_begin);

}



void Mesh_simple::node_get_cells (const Amanzi::AmanziMesh::Entity_ID nodeid, 
			   const Amanzi::AmanziMesh::Parallel_type ptype,
			   Amanzi::AmanziMesh::Entity_ID_List *cellids) const 
{
  throw std::exception();
}
    
// Faces of type 'ptype' connected to a node
    
void Mesh_simple::node_get_faces (const Amanzi::AmanziMesh::Entity_ID nodeid, 
		     const Amanzi::AmanziMesh::Parallel_type ptype,
		     Amanzi::AmanziMesh::Entity_ID_List *faceids) const
{
  throw std::exception();
}   
 
// Get faces of ptype of a particular cell that are connected to the
// given node

void Mesh_simple::node_get_cell_faces (const Amanzi::AmanziMesh::Entity_ID nodeid, 
			  const Amanzi::AmanziMesh::Entity_ID cellid,
			  const Amanzi::AmanziMesh::Parallel_type ptype,
			  Amanzi::AmanziMesh::Entity_ID_List *faceids) const
{
  throw std::exception();
}
    
// Cells connected to a face
    
void Mesh_simple::face_get_cells (const Amanzi::AmanziMesh::Entity_ID faceid, 
		     const Amanzi::AmanziMesh::Parallel_type ptype,
		     Amanzi::AmanziMesh::Entity_ID_List *cellids) const
{
  throw std::exception();
}
    


// Same level adjacencies
//-----------------------

// Face connected neighboring cells of given cell of a particular ptype
// (e.g. a hex has 6 face neighbors)

// The order in which the cellids are returned cannot be
// guaranteed in general except when ptype = USED, in which case
// the cellids will correcpond to cells across the respective
// faces given by cell_get_faces

void Mesh_simple::cell_get_face_adj_cells(const Amanzi::AmanziMesh::Entity_ID cellid,
			     const Amanzi::AmanziMesh::Parallel_type ptype,
			     Amanzi::AmanziMesh::Entity_ID_List *fadj_cellids) const
{
  throw std::exception();
}

// Node connected neighboring cells of given cell
// (a hex in a structured mesh has 26 node connected neighbors)
// The cells are returned in no particular order

void Mesh_simple::cell_get_node_adj_cells(const Amanzi::AmanziMesh::Entity_ID cellid,
			     const Amanzi::AmanziMesh::Parallel_type ptype,
			     Amanzi::AmanziMesh::Entity_ID_List *nadj_cellids) const
{
  throw std::exception();
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

Amanzi::AmanziMesh::Cell_type Mesh_simple::cell_get_type_4viz(const Amanzi::AmanziMesh::Entity_ID cellid) const 
{
  return Amanzi::AmanziMesh::HEX;
}
    
    
// See cell_get_nodes for details on node ordering
    
void Mesh_simple::cell_get_nodes_4viz (const Amanzi::AmanziMesh::Entity_ID cellid, 
				Amanzi::AmanziMesh::Entity_ID_List *nodeids) const 
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

