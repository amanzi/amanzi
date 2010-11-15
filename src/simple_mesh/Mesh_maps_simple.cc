#include "Mesh_maps_simple.hh"
#include <Teuchos_RCP.hpp>


Mesh_maps_simple::Mesh_maps_simple (double x0, double y0, double z0,
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

Mesh_maps_simple::Mesh_maps_simple ( Teuchos::ParameterList &parameter_list,
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







Mesh_maps_simple::~Mesh_maps_simple()
{
  delete cell_map_;
  delete face_map_;
  delete node_map_;
  if (!mesh_block_z0)  delete [] mesh_block_z0;
  if (!mesh_block_z1)  delete [] mesh_block_z1;
}



 void Mesh_maps_simple::update ()
 {
     clear_internals_ ();
     update_internals_ ();
 }

 void Mesh_maps_simple::clear_internals_ ()
 {

   coordinates_.resize(0);

   cell_to_face_.resize(0);
   cell_to_node_.resize(0);
   face_to_node_.resize(0);

   side_sets_.resize(0);

 }


 void Mesh_maps_simple::update_internals_()
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
	       std::vector<double> coords(24);
	       cell_to_coordinates(ic,coords.begin(),coords.end());
	       
	       // first check if all z-coordinates are larger than
	       // mesh_block_z0[nb]
	       bool inside = true;
	       for (int i=4; i<8; i++) 
		 if ( coords[3*i+2] <= mesh_block_z0[nb] ) inside = false;
	       for (int i=0; i<8; i++) 
		 if ( coords[3*i+2] > mesh_block_z1[nb] ) inside = false;

	       if (inside) count++;
	       
	     }
	   element_blocks_[nb].resize(count);
	   
	   // add cells to mesh block nb
	   count = 0;
	   for (int ic=0; ic<num_cells_; ic++) 
	     {
	       std::vector<double> coords(24);
	       cell_to_coordinates(ic,coords.begin(),coords.end());
	       
	       // first check if all z-coordinates are larger than
	       // mesh_block_z0[nb]
	       bool inside = true;
	       for (int i=4; i<8; i++) 
		 if ( coords[3*i+2] <= mesh_block_z0[nb] ) inside = false;
	       for (int i=0; i<8; i++) 
		 if ( coords[3*i+2] > mesh_block_z1[nb] ) inside = false;
	       
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
       
	   std::cout << "Mesh_maps_simple: mesh block " << nb << " has " << element_blocks_[nb].size() << " cells" << std::endl;
	 }


       if (ncb < num_cells_) 
	 {
	   std::cout << "Mesh_maps_simple: WARNING... not all cells are in a mesh block\n";
	   std::cout << "  number of cells in mesh blocks = " << ncb << std::endl;
	   std::cout << "  number of cells in total       = " << num_cells_ << std::endl;
	 }
       if (ncb > num_cells_)
	 std::cout << "Mesh_maps_simple: WARNING... overlapping mesh blocks" << std::endl;
       
       
     }
   
   
   





   build_maps_ ();
}


void Mesh_maps_simple::build_maps_ ()
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


unsigned int Mesh_maps_simple::count_entities (Mesh_data::Entity_kind kind, Element_Category category) const
{
  switch (kind) {
  case Mesh_data::FACE: 
    return num_faces_;
    break;
  case Mesh_data::NODE:
    return num_nodes_;
    break;
  case Mesh_data::CELL:
    return num_cells_;
    break;
  default:
    throw std::exception();
    break;
  }
}


unsigned int Mesh_maps_simple::num_sets(Mesh_data::Entity_kind kind) const
{
  switch (kind) {
  case Mesh_data::FACE:
    return 6;
    break;
  case Mesh_data::CELL:
    return 1;
    break;
  default:
    // nothing yet for NODE or CELL
    throw std::exception();
    break;
  }
}

unsigned int Mesh_maps_simple::get_set_size (unsigned int set_id, 
					     Mesh_data::Entity_kind kind,
					     Element_Category category) const
{
  // we ignore category, since this is a serial implementation

  switch (kind) {
  case Mesh_data::FACE:
    return side_sets_[set_id].size();
    break;
  case Mesh_data::CELL:
    return element_blocks_[set_id].size();
    break;
  default:
    throw std::exception();
    break;
  }
}


void Mesh_maps_simple::cell_to_faces (unsigned int cell, 
				      std::vector<unsigned int>::iterator destination_begin, 
				      std::vector<unsigned int>::iterator destination_end)
{
  // ASSERT ((unsigned int) (destination_end - destination_begin) == 6);
    const unsigned int index = 6*cell;
    std::vector<unsigned int>::iterator begin = cell_to_face_.begin() + index;
    std::vector<unsigned int>::iterator end = begin + 6;
    std::copy (begin, end, destination_begin);
}


void Mesh_maps_simple::cell_to_faces (unsigned int cell, 
				      unsigned int * destination_begin, 
				      unsigned int * destination_end)
{
  // ASSERT ((unsigned int) (destination_end - destination_begin) == 6);
    const unsigned int index = 6*cell;
    std::vector<unsigned int>::iterator begin = cell_to_face_.begin() + index;
    std::vector<unsigned int>::iterator end = begin + 6;
    std::copy (begin, end, destination_begin);
}


void Mesh_maps_simple::cell_to_nodes (unsigned int cell, 
				      std::vector<unsigned int>::iterator destination_begin, 
				      std::vector<unsigned int>::iterator destination_end) 
{
  // ASSERT ((unsigned int) (destination_end - destination_begin) == 8);
    const unsigned int index = 8*cell;
    std::vector<unsigned int>::iterator begin = cell_to_node_.begin () + index;
    std::vector<unsigned int>::iterator end   = begin + 8;
    std::copy (begin, end, destination_begin);
}

void Mesh_maps_simple::cell_to_nodes (unsigned int cell, 
				      unsigned int * destination_begin, 
				      unsigned int * destination_end) 
{
  // ASSERT ((unsigned int) (destination_end - destination_begin) == 8);
    const unsigned int index = 8*cell;
    std::vector<unsigned int>::iterator begin = cell_to_node_.begin () + index;
    std::vector<unsigned int>::iterator end   = begin + 8;
    std::copy (begin, end, destination_begin);
}


void Mesh_maps_simple::face_to_nodes (unsigned int face, 
				      std::vector<unsigned int>::iterator destination_begin,
				      std::vector<unsigned int>::iterator destination_end) 
{
  // ASSERT ((unsigned int) (destination_end - destination_begin) == 4);
    const unsigned int index = 4*face;
    std::vector<unsigned int>::iterator begin = face_to_node_.begin () + index;
    std::vector<unsigned int>::iterator end   = begin + 4;
    
    for (std::vector<unsigned int>::iterator it=begin; it !=end; it++) {
      *destination_begin = *it;
      destination_begin++;
    }

      //std::copy (begin, end, destination_begin);
}

void Mesh_maps_simple::face_to_nodes (unsigned int face, 
				      unsigned int * destination_begin,
				      unsigned int * destination_end) 
{
  // ASSERT ((unsigned int) (destination_end - destination_begin) == 4);
    const unsigned int index = 4*face;
    std::vector<unsigned int>::iterator begin = face_to_node_.begin () + index;
    std::vector<unsigned int>::iterator end   = begin + 4;
    
    for (std::vector<unsigned int>::iterator it=begin; it !=end; it++) {
      *destination_begin = *it;
      destination_begin++;
    }

      //std::copy (begin, end, destination_begin);
}


// Cooordinate Getters
// -------------------

void Mesh_maps_simple::node_to_coordinates (unsigned int local_node_id, 
					    std::vector<double>::iterator destination_begin, 
					    std::vector<double>::iterator destination_end) 
{
  //  ASSERT ((unsigned int) (end-begin) == 3);
  const unsigned int index = 3*local_node_id;
  std::vector<double>::iterator begin = coordinates_.begin() + index;
  std::vector<double>::iterator end   = begin + 3;
  std::copy (begin, end, destination_begin);
}
void Mesh_maps_simple::node_to_coordinates (unsigned int local_node_id, 
					    double * destination_begin, 
					    double * destination_end) 
{
  //  ASSERT ((unsigned int) (end-begin) == 3);
  const unsigned int index = 3*local_node_id;
  std::vector<double>::iterator begin = coordinates_.begin() + index;
  std::vector<double>::iterator end   = begin + 3;
  std::copy (begin, end, destination_begin);
}


void Mesh_maps_simple::face_to_coordinates (unsigned int local_face_id, 
					    std::vector<double>::iterator begin, 
					    std::vector<double>::iterator end)
{
  // ASSERT ((unsigned int) (end-begin) == 12);

  std::vector<unsigned int> node_indices(4);
  face_to_nodes (local_face_id, node_indices.begin(), node_indices.end());
  for (std::vector<unsigned int>::iterator it = node_indices.begin(); 
       it != node_indices.end(); ++it)
    {
      node_to_coordinates ( *it, begin, begin+3);
      begin+=3;
     }


}
void Mesh_maps_simple::face_to_coordinates (unsigned int local_face_id, 
					    double * begin, 
					    double * end)
{
  // ASSERT ((unsigned int) (end-begin) == 12);

  std::vector<unsigned int> node_indices(4);
  face_to_nodes (local_face_id, node_indices.begin(), node_indices.end());
  for (std::vector<unsigned int>::iterator it = node_indices.begin(); 
       it != node_indices.end(); ++it)
    {
      node_to_coordinates ( *it, begin, begin+3);
      begin+=3;
     }


}
void Mesh_maps_simple::cell_to_coordinates (unsigned int local_cell_id, 
					    std::vector<double>::iterator begin, 
					    std::vector<double>::iterator end)
{
  // ASSERT ((unsigned int) (end-begin) == 24);
  
  std::vector<unsigned int> node_indices(8);
  cell_to_nodes (local_cell_id, node_indices.begin(), node_indices.end());
  for (std::vector<unsigned int>::iterator it = node_indices.begin(); 
       it != node_indices.end(); ++it)
    {
        node_to_coordinates ( *it, begin, begin+3);
        begin+=3;
    }

}

void Mesh_maps_simple::cell_to_coordinates (unsigned int local_cell_id, 
					    double * begin, 
					    double * end)
{
  // ASSERT ((unsigned int) (end-begin) == 24);
  
  std::vector<unsigned int> node_indices(8);
  cell_to_nodes (local_cell_id, node_indices.begin(), node_indices.end());
  for (std::vector<unsigned int>::iterator it = node_indices.begin(); 
       it != node_indices.end(); ++it)
    {
        node_to_coordinates ( *it, begin, begin+3);
        begin+=3;
    }

}



// Set getters
// -----------

void Mesh_maps_simple::get_set_ids (Mesh_data::Entity_kind kind, 
				    std::vector<unsigned int>::iterator begin, 
				    std::vector<unsigned int>::iterator end) const
{
  
  std::vector<int> ids(0);

  switch (kind) {
  case Mesh_data::FACE: 
    {
      ids.resize(6);
      for (int i=0; i<6; i++) ids[i]=i;
      
      std::copy (ids.begin(), ids.end(), begin);
      break;
    }
  case Mesh_data::CELL:
    {
      ids.resize(1);
      ids[0] = 0;
      
      std::copy (ids.begin(), ids.end(), begin);
      break;
    }      
  default:
    // we do not have anything for CELL and NODE, yet
    throw std::exception();
  }
}

void Mesh_maps_simple::get_set_ids (Mesh_data::Entity_kind kind, 
				    unsigned int * begin, 
				    unsigned int * end) const
{
  
  std::vector<int> ids(0);

  switch (kind) {
  case Mesh_data::FACE: 
    {
      ids.resize(6);
      for (int i=0; i<6; i++) ids[i]=i;
      
      std::copy (ids.begin(), ids.end(), begin);
      break;
    }
  case Mesh_data::CELL:
    {
      ids.resize(1);
      ids[0] = 0;
      
      std::copy (ids.begin(), ids.end(), begin);
      break;
    }      
  default:
    // we do not have anything for CELL and NODE, yet
    throw std::exception();
  }
}

void Mesh_maps_simple::get_set (unsigned int set_id, Mesh_data::Entity_kind kind, 
				Element_Category category, 
				std::vector<unsigned int>::iterator begin, 
				std::vector<unsigned int>::iterator end) const
{
  // we ignore category since this is a serial implementation
  
  switch (kind) {
  case Mesh_data::FACE:
    std::copy(side_sets_[set_id].begin(), side_sets_[set_id].end(), begin) ;
    break;
  case Mesh_data::CELL:
    std::copy(element_blocks_[set_id].begin(), element_blocks_[set_id].end(), begin) ;
    break;
  default:
    // we do not have anything for CELL and NODE, yet
    throw std::exception();
  }
}
void Mesh_maps_simple::get_set (unsigned int set_id, Mesh_data::Entity_kind kind, 
				Element_Category category, 
				unsigned int * begin, 
				unsigned int * end) const
{
  // we ignore category since this is a serial implementation
  
  switch (kind) {
  case Mesh_data::FACE:
    std::copy(side_sets_[set_id].begin(), side_sets_[set_id].end(), begin) ;
    break;
  case Mesh_data::CELL:
    std::copy(element_blocks_[set_id].begin(), element_blocks_[set_id].end(), begin) ;
    break;
  default:
    // we do not have anything for CELL and NODE, yet
    throw std::exception();
  }
}


bool Mesh_maps_simple::valid_set_id (unsigned int id, Mesh_data::Entity_kind kind) const 
{
  switch (kind) {
  case Mesh_data::FACE:
    return (id<6) ;
    break;
  case Mesh_data::CELL:
    return (id < number_of_mesh_blocks);
    break;
  default:
    // we do not have anything for CELL and NODE, yet
    throw std::exception();  
  }
    
}

void Mesh_maps_simple::cell_to_face_dirs (unsigned int cell, 
					  std::vector<int>::iterator destination_begin, 
					  std::vector<int>::iterator destination_end)
{
  const unsigned int index = 6*cell;
  std::vector<int>::iterator begin = cell_to_face_dirs_.begin() + index;
  std::vector<int>::iterator end = begin + 6;
  std::copy (begin, end, destination_begin);  
}


void Mesh_maps_simple::cell_to_face_dirs (unsigned int cell, 
					  int * destination_begin, int * destination_end)
{
  const unsigned int index = 6*cell;
  std::vector<int>::iterator begin = cell_to_face_dirs_.begin() + index;
  std::vector<int>::iterator end = begin + 6;
  std::copy (begin, end, destination_begin);  
}


void Mesh_maps_simple::set_coordinate(unsigned int local_node_id, 
				      double *source_begin, 
				      double *source_end)
{
  unsigned int index = 3*local_node_id;
  
  std::vector<double>::iterator destination_begin = coordinates_.begin() + index;
  std::copy(source_begin, source_end, destination_begin);

}
