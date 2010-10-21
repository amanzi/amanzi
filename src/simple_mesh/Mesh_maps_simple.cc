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
 
  update();
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

   // now we need to fix the orientation of boundary xy faces for iz=0 
   for (int iy=0; iy<ny_; iy++)
     for (int ix=0; ix<nx_; ix++)
       {
   	 int istart = 4 * xyface_index_(ix,iy,0);

   	 double dummy = face_to_node_[istart];
   	 face_to_node_[istart] = face_to_node_[istart+2];
   	 face_to_node_[istart+2] = dummy;
       }

   // now we need to fix the orientation of boundary xz faces for iy=0 
   for (int iz=0; iz<nz_; iz++)
     for (int ix=0; ix<nx_; ix++)
       {
   	 int istart = 4 * xzface_index_(ix,0,iz);

   	 double dummy = face_to_node_[istart];
   	 face_to_node_[istart] = face_to_node_[istart+2];
   	 face_to_node_[istart+2] = dummy;
       }

   // now we need to fix the orientation of boundary yz faces for ix=0 
   for (int iz=0; iz<nz_; iz++)
     for (int iy=0; iy<ny_; iy++)
       {
   	 int istart = 4 * yzface_index_(0,iy,iz);

   	 double dummy = face_to_node_[istart];
   	 face_to_node_[istart] = face_to_node_[istart+2];
   	 face_to_node_[istart+2] = dummy;
       }


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
  default:
    throw std::exception();
    break;
  }
}
