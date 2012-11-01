#include "deform.hh"

DeformMesh::DeformMesh( Teuchos::ParameterList& plist,
			const Teuchos::RCP<AmanziMesh::Mesh>& mesh ) : 
  plist_(plist), mesh_(mesh) {}

void DeformMesh::print_goodbye() {
  std::cout << "goodbye!" << std::endl;
}

void DeformMesh::check_mesh_nodes() {
  
  // mesh dimension
  int dim = mesh_->space_dimension();
  
  // number of vertices
  int nV = mesh_->num_entities(Amanzi::AmanziMesh::NODE, Amanzi::AmanziMesh::OWNED);

  // get and print the coordinates of the mesh nodes
  for (int iV=0; iV<nV; iV++) {
    Amanzi::AmanziGeometry::Point coords;
    coords.init(dim);
    mesh_->node_get_coordinates(iV,&coords);
    // output
    printf("iV=%3i (",iV) ;
    for ( int s=0; s<dim; ++s ) {
      printf("%14.7e",coords[s]) ;
      if ( s!=dim-1 ) printf(", ") ;
    }
    printf(")\n");
  }

}

void DeformMesh::print_VTK_unstructured_mesh( string fname ) {
  // mesh dimension
  int dim = mesh_->space_dimension();

  // number of vertices
  int nV = mesh_->num_entities(Amanzi::AmanziMesh::NODE, Amanzi::AmanziMesh::OWNED);
  
  // number of faces
  int nF = mesh_->num_entities(Amanzi::AmanziMesh::FACE, Amanzi::AmanziMesh::OWNED);

  // output file name
  string output_filename = fname+string(".vtk");

  // write VTK header
  ofstream vtk_out(output_filename.c_str());
  vtk_out << "# vtk DataFile Version 2.0" << endl;
  vtk_out << "deform example"             << endl;
  vtk_out << "ASCII"                      << endl;
  vtk_out << "DATASET POLYDATA"           << endl;
  vtk_out << "POINTS " << nV << "  float" << endl;
  
  // get and print the coordinates of the mesh nodes
  for (int iV=0; iV<nV; iV++) {
    Amanzi::AmanziGeometry::Point coords;
    coords.init(dim);
    mesh_->node_get_coordinates(iV,&coords);
    if ( dim==2 ) {
      vtk_out << coords[0] << "  " 
	      << coords[1] << "  0.0" << endl; 
    } else {
      vtk_out << coords[0] << "  " 
	      << coords[1] << "  " 
	      << coords[2] << endl;
    }
  }
  // count the number of polygonal entries (size)
  int size = 0 ;
  for (int iF=0; iF<nF; iF++) {
    Entity_ID_List nodeids;
    mesh_->face_get_nodes(iF,&nodeids);
    size += 1+nodeids.size();
  }

  // mesh dimension (dim=2 LINES, dim=3 POLYGONS)
  if ( dim == 2 ) {
    vtk_out << "LINES " << nF << " " << size << endl;
  } else {
    vtk_out << "POLYGONS " << nF << " " << size << endl;
  }

  // print the face node ids
  for (int iF=0; iF<nF; iF++) {
    Entity_ID_List nodeids;
    mesh_->face_get_nodes(iF,&nodeids);
    int nFV = nodeids.size();
    vtk_out << nFV << "  ";
    for (int ilV=0; ilV<nFV; ilV++) {
      vtk_out << nodeids[ilV] << "  ";
    }
    vtk_out << endl;
  }
  vtk_out.close();
}


/*
int Mesh::deform (const Entity_ID_List nodeids,
                  const AmanziGeometry::Point_List new_positions,
                  const bool keep_valid,
                  AmanziGeometry::Point_List *final_positions) {
*/

void DeformMesh::move_a_single_node() {

  // get space dimensions
  int dim = mesh_->space_dimension();

  // set the list of the new position coords
  AmanziGeometry::Point_List newpos, finpos;
  
  // move the mid point on the mesh top
  Entity_ID_List nodeids;

  // number of vertices
  int nV = mesh_->num_entities(Amanzi::AmanziMesh::NODE, Amanzi::AmanziMesh::OWNED);

  // search the id of the mid point on the top
  for (int iV=0; iV<nV; iV++) {
    Amanzi::AmanziGeometry::Point coords;
    coords.init(dim);

    mesh_->node_get_coordinates(iV,&coords);
    bool ok_node = true;
    for ( int s=0; s<dim-1; ++s ) {
      ok_node &= abs( coords[s]-0.5 )<1e-10; 
    }
    ok_node &= abs( coords[dim-1]-1. )<1e-10; 

    if ( ok_node ) {
      nodeids.push_back(iV);

      Amanzi::AmanziGeometry::Point new_coords;
      new_coords.init(dim);

      for ( int s=0; s<dim-1; ++s ) {
	new_coords[s] = coords[s];
      }
      new_coords[dim-1] = stretch( coords[dim-1] );
      newpos.push_back( new_coords );
    }
  }

  // print the nodes that will be changed
  for ( int i=0; i<nodeids.size(); ++i ) {
    printf("found --> nodeids[%2i]=%3i\n",i,nodeids[i]) ;
  }
  LINE(--) ;
  
  // keep_valid --> check mesh consistency
  bool keep_valid = true;

  // compute the deformed mesh
  mesh_->deform( nodeids, newpos, keep_valid, &finpos);
}

// stretching: y(s) = a s^2 + b s
// y(1)   = 1/2
// y(1/2) = 1/3
double DeformMesh::stretch( double s ) {
  return -1./3.*s*s+5./6.*s;
}

void DeformMesh::move_a_node_column() {

  // get space dimensions
  int dim = mesh_->space_dimension();

  // set the list of the new position coords
  AmanziGeometry::Point_List newpos, finpos;

  // move the mid point on the mesh top
  Entity_ID_List nodeids;

  // number of vertices
  int nV = mesh_->num_entities(Amanzi::AmanziMesh::NODE, Amanzi::AmanziMesh::OWNED);

  // get and print the coordinates of the old mesh nodes
  for (int iV=0; iV<nV; iV++) {
    Amanzi::AmanziGeometry::Point coords;
    coords.init(dim);
    mesh_->node_get_coordinates(iV,&coords);

    // stretch the last coordinate
    bool ok_column = true;
    for ( int s=0; s<dim-1; ++s ) {
      ok_column &= abs(coords[s]-0.5)<1.e-10;
    }
    if ( ok_column ) {
      nodeids.push_back( iV );
      Amanzi::AmanziGeometry::Point new_coords;
      new_coords.init(dim);
      for ( int s=0; s<dim-1; ++s ) {
	new_coords[s] = coords[s];
      }
      new_coords[dim-1] = stretch( coords[dim-1] );
      newpos.push_back( new_coords );
    }
  }

  // print the nodes that will be changed
  for ( int i=0; i<nodeids.size(); ++i ) {
    printf("found --> nodeids[%2i]=%3i\n",i,nodeids[i]) ;
  }
  LINE(--) ;

  // keep_valid --> check mesh consistency
  bool keep_valid = true;

  // compute the deformed mesh
  mesh_->deform( nodeids, newpos, keep_valid, &finpos);
}

void DeformMesh::parabolic_profile( double ss ) {
  bool verbose(false);

  if ( verbose ) { cout << "step ss = " << ss << endl; }
  
  // get space dimensions
  int dim = mesh_->space_dimension();

  Amanzi::AmanziGeometry::Point P0;
  P0.init(dim);
  if ( dim==2 ) { P0[0]=0.5; P0[1]=1.; }
  else { P0[0]=0.5; P0[1]=0.5; P0[2]=1.; }

    // set the list of the new position coords
  AmanziGeometry::Point_List newpos, finpos;
  
  // move the mid point on the mesh top
  Entity_ID_List nodeids;

  // number of vertices
  int nV = mesh_->num_entities(Amanzi::AmanziMesh::NODE, Amanzi::AmanziMesh::OWNED);

  Amanzi::AmanziGeometry::Point coords, new_coords;
  coords.init(dim);
  new_coords.init(dim);

  // search the id of the mid point on the top
  for (int iV=0; iV<nV; iV++) {
    
    // get the coords of the node
    mesh_->node_get_coordinates(iV,&coords);

    // check if coords must be changed
    double dist = 0.;
    for ( int s=0; s<dim-1; ++s ) {
      dist += pow( coords[s]-P0[s], 2 );
    }
    dist = sqrt(dist);

    double alpha  = 1./20.;
    double radius = 0.25;

    // if yes, change the coords of the node
    bool ok_node = dist<0.5;
    if ( ok_node ) {
      
      for ( int s=0; s<dim-1; ++s ) {
	new_coords[s] = coords[s];
      }
      double fac = 1.-(1.-ss)*alpha*exp(-pow(dist/radius,2)) ;
      new_coords[dim-1] = coords[dim-1] * fac;

      // puch back for deform method
      nodeids.push_back(iV);
      newpos.push_back( new_coords );

      if ( verbose ) {
	printf("found: iV=%3i,  dist=%14.7e,  fac=%14.7e\n",iV,dist,fac);
	printf("coords    =(%14.7e, %14.7e)\n",coords[0],coords[1]);
	printf("new_coords=(%14.7e, %14.7e)\n",new_coords[0],new_coords[1]);
	LINE(-);
      }
    }
  }

  // print the nodes that will be changed
  if ( verbose ) {
    for ( int i=0; i<nodeids.size(); ++i ) {
      printf("found: nodeids[%2i]=%3i  coords=%14.7e --> new_coords=%14.7e\n",i,nodeids[i],coords[dim-1],new_coords[dim-1]);
    }
    LINE(--);
    LINE(--);
  }

  // keep_valid --> check mesh consistency
  bool keep_valid = true;

  // compute the deformed mesh
  mesh_->deform( nodeids, newpos, keep_valid, &finpos);
}

/*
 Mesh_MSTK(const double x0, const double y0,
	    const double x1, const double y1,
	    const int nx, const int ny, 
	    const Epetra_MpiComm *comm,
	    const AmanziGeometry::GeometricModelPtr& gm = 
	    (AmanziGeometry::GeometricModelPtr) NULL);
*/

void DeformMesh::loop_monitor( int i, int imax ){
  if ( imax>=100 ) {
    int ni_10 = imax/10 ;
    int ni_02 = imax/50 ;
    if      ( i       == 0 ) { cout << "0% " ; }
    else if ( i%ni_10 == 0 ) { cout << 10*i/ni_10 << "% " ; }
    else if ( i%ni_02 == 0 ) { cout << "." ; }
    cout << flush;
  } else {
    if ( i%10 == 0 ) { cout << "." ; }
  }
}

void DeformMesh::parabolic_profile(){
  string fname("");
  int kmax = 20;
  for ( int k=0; k<=kmax; ++k ) {
    //loop_monitor(k,kmax) ;
    double ss = double(k)/double(kmax);
    ostringstream oss;
    oss << k;
    if      (  0<=k && k<10  ) { fname = string("mesh_0") + oss.str(); }
    else if ( 10<=k && k<100 ) { fname = string("mesh_")  + oss.str(); }
    cout << k << "  " << fname << endl ;
    //else if ( 100<=k && k<1000 ) { fname = string("mesh_")   + oss.str(); }
    print_VTK_unstructured_mesh( fname );
    //---
    parabolic_profile(ss);
  }

  // print final mesh
  ostringstream oss;
  oss << kmax+1;
  fname = string("mesh_") + oss.str();
  print_VTK_unstructured_mesh( fname );
}
