#include "deform.hh"

namespace Amanzi {
namespace Deform {


// -------------------------------------------------------------
// CONSTRUCTORS
// -------------------------------------------------------------
DeformMesh::DeformMesh( Teuchos::ParameterList& plist,
			const Teuchos::RCP<AmanziMesh::Mesh>& mesh0,
			const Teuchos::RCP<AmanziMesh::Mesh>& mesh1,
			const Teuchos::RCP<AmanziMesh::Mesh>& mesh ) : 
  plist_(plist), mesh0_(mesh0), mesh1_(mesh1), mesh_(mesh) {}

DeformMesh::DeformMesh( Teuchos::ParameterList& plist,
			const Teuchos::RCP<AmanziMesh::Mesh>& mesh ) : 
  plist_(plist), mesh0_(mesh), mesh1_(mesh), mesh_(mesh) {}

// -------------------------------------------------------------
// MAIN METHODS
// -------------------------------------------------------------

// -------------------------------------------------------------
// Setup data
// -------------------------------------------------------------
void DeformMesh::setup(const Teuchos::Ptr<State>& S) {
};

// -------------------------------------------------------------
// Initialize PK
// -------------------------------------------------------------
void DeformMesh::initialize(const Teuchos::Ptr<State>& S) {
};

// -------------------------------------------------------------
// Update any secondary (dependent) variables given a solution.
// -------------------------------------------------------------
void DeformMesh::commit_state(double dt, const Teuchos::RCP<State>& S) {
}

// -------------------------------------------------------------
// Update any diagnostic variables prior to vis
// -------------------------------------------------------------
void DeformMesh::calculate_diagnostics(const Teuchos::RCP<State>& S) {
}

// -------------------------------------------------------------
// AUX METHODS
// -------------------------------------------------------------

// say goodbye
void DeformMesh::print_goodbye() {
  std::cout << "goodbye!" << std::endl;
}

// monitor iterations
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

// -------------------------------------------------------------
// METHODS TO OUTPUT DATA IN VTK FORMAT
// -------------------------------------------------------------

// write mesh_ in VTK format  
void DeformMesh::print_VTK_unstructured_mesh( string fname ) {
  cout << "--->mesh VTK output: " << fname << endl ;

  // mesh dimension
  int dim = mesh_->space_dimension();

  // number of vertices
  int nV = mesh_->num_entities(Amanzi::AmanziMesh::NODE,
				Amanzi::AmanziMesh::OWNED);
  
  // number of faces
  int nF = mesh_->num_entities(Amanzi::AmanziMesh::FACE, 
				Amanzi::AmanziMesh::OWNED);

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
    mesh_->face_get_nodes( iF, &nodeids );
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

// write the boundary faces of mesh_ in VTK format
void DeformMesh::print_VTK_domain_boundary( string fname ) {
  cout << "--->mesh VTK output: " << fname << endl ;

  // mesh dimension
  int dim = mesh_->space_dimension();

  // number of vertices
  int nV = mesh_->num_entities(Amanzi::AmanziMesh::NODE, 
				Amanzi::AmanziMesh::OWNED);
  
  // number of faces
  int nF = mesh_->num_entities(Amanzi::AmanziMesh::FACE, 
				Amanzi::AmanziMesh::OWNED);

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
  int nbF  = 0 ;
  Entity_ID_List bnd_face_ids;
  for (int iF=0; iF<nF; iF++) {

    Entity_ID_List cellids;
    mesh_->face_get_cells ( iF, Amanzi::AmanziMesh::OWNED, &cellids );

    if ( cellids.size()==1 ) {
      bnd_face_ids.push_back( iF );
      ++nbF;
      Entity_ID_List nodeids;
      mesh_->face_get_nodes(iF,&nodeids);
      size += 1+nodeids.size();
    }
  }
  assert( nbF==bnd_face_ids.size() );

  // mesh dimension (dim=2 LINES, dim=3 POLYGONS)
  if ( dim == 2 ) {
    vtk_out << "LINES " << nbF << " " << size << endl;
  } else {
    vtk_out << "POLYGONS " << nbF << " " << size << endl;
  }

  // print the face node ids
  for (int ilF=0; ilF<nbF; ilF++) {
    int iF = bnd_face_ids[ilF];
    Entity_ID_List nodeids;
    mesh_->face_get_nodes(iF,&nodeids);
    int nFV = nodeids.size();
    vtk_out << nFV << "  ";
    for (int ilV=0; ilV<nFV; ilV++) {
      vtk_out << nodeids[ilV] << "  ";
    }
  }
  vtk_out.close();
}

// write a submesh of mesh_ in VTK format
void DeformMesh::print_VTK_submesh( string fname ) {
  cout << "--->mesh VTK output: " << fname << endl ;

  // mesh dimension
  int dim = mesh_->space_dimension();

  // number of vertices
  int nV = mesh_->num_entities(Amanzi::AmanziMesh::NODE, 
			       Amanzi::AmanziMesh::OWNED);
  
  // flag & renumber the nodes above the threshold (z direction) 
  double threshold = 3.8;
  int sub_nV = 0;
  Entity_ID_List sub_nodeids(nV);
  for (int iV=0; iV<nV; iV++) {
    Amanzi::AmanziGeometry::Point coords;
    coords.init(dim);
    mesh_->node_get_coordinates(iV,&coords);
    if ( coords[dim-1]>threshold ) {
      sub_nodeids[iV]=sub_nV;
      ++sub_nV;
    } else {
      sub_nodeids[iV]=-1;
    }
  }
  
  // output file name
  string output_filename = fname+string(".vtk");

  // write VTK header
  ofstream vtk_out(output_filename.c_str());
  vtk_out << "# vtk DataFile Version 2.0" << endl;
  vtk_out << "deform example"             << endl;
  vtk_out << "ASCII"                      << endl;
  vtk_out << "DATASET POLYDATA"           << endl;
  vtk_out << "POINTS " << sub_nV << "  float" << endl;
  
  // get and print the coordinates of the mesh nodes
  for (int iV=0; iV<nV; iV++) {
    if ( sub_nodeids[iV]!=-1 ) { // check if the node is flagged
      Amanzi::AmanziGeometry::Point coords;
      coords.init(dim);
      mesh_->node_get_coordinates(iV,&coords);
      vtk_out << coords[0] << "  " 
	      << coords[1] << "  " 
	      << coords[2] << endl;
    }
  }

  // number of faces
  int nF = mesh_->num_entities(Amanzi::AmanziMesh::FACE, 
			       Amanzi::AmanziMesh::OWNED);

  // count the number of polygonal entries (size)
  int size = 0 ;
  int nbF  = 0 ;
  Entity_ID_List sub_face_ids;
  for (int iF=0; iF<nF; iF++) {

    Entity_ID_List nodeids;
    mesh_->face_get_nodes ( iF, &nodeids );

    int nFV=nodeids.size();
    bool ok_face=true;
    for ( int ilV=0; ilV<nFV; ++ilV ) {
      Amanzi::AmanziGeometry::Point coords;
      coords.init(dim);
      int iV = nodeids[ilV];
      mesh_->node_get_coordinates(iV,&coords);
      ok_face &= coords[dim-1]>threshold;
    }

    if ( ok_face ) {
      sub_face_ids.push_back( iF );
      ++nbF;
      size += 1+nodeids.size();
    }
  }
  assert( nbF==sub_face_ids.size() );

  // if dim=3 the write POLYGONS
  vtk_out << "POLYGONS " << nbF << " " << size << endl;

  // print the face node ids
  for (int ilF=0; ilF<nbF; ilF++) {
    int iF = sub_face_ids[ilF];
    Entity_ID_List nodeids;
    mesh_->face_get_nodes(iF,&nodeids);
    int nFV = nodeids.size();
    vtk_out << nFV << "  ";
    for (int ilV=0; ilV<nFV; ilV++) {
      int iV = nodeids[ilV];
      int sub_iV = sub_nodeids[iV];
      assert( sub_iV!=-1 );
      vtk_out << sub_iV << "  ";
    }
    vtk_out << endl ;
  }

  // add a colour to each vertex proportional to its height
  vtk_out << "POINT_DATA " << sub_nV << endl;
  vtk_out << "SCALARS scalars float 1" << endl;
  vtk_out << "LOOKUP_TABLE default" << endl;
  
  // get and print the coordinates of the mesh nodes
  int nc=0;
  for (int iV=0; iV<nV; iV++) {
    if ( sub_nodeids[iV]!=-1 ) { // check if the node is flagged
      if ( nc++%20==0 && nc>0 ) { vtk_out<<endl; }
      Amanzi::AmanziGeometry::Point coords;
      coords.init(dim);
      mesh_->node_get_coordinates(iV,&coords);
      vtk_out << coords[2] << "  ";
    }
  }

  vtk_out.close();
}

// -------------------------------------------------------------
// METHODS TO MOVE THE MESH NODES
// -------------------------------------------------------------

// move nodes following a bell_shaped profile
void DeformMesh::bell_shaped_profile( double ss ) {
  bool verbose(false);

  if ( verbose ) { cout << "step ss = " << ss << endl; }
  
  // get space dimensions
  int dim = mesh0_->space_dimension();

  Amanzi::AmanziGeometry::Point P0;
  P0.init(dim);
  if ( dim==2 ) { P0[0]=0.5; P0[1]=1.; }
  else { P0[0]=0.5; P0[1]=0.5; P0[2]=1.; }

    // set the list of the new position coords
  AmanziGeometry::Point_List newpos, finpos;
  
  // move the mid point on the mesh top
  Entity_ID_List nodeids;

  // number of vertices
  int nV = mesh0_->num_entities(Amanzi::AmanziMesh::NODE, 
				Amanzi::AmanziMesh::OWNED);

  Amanzi::AmanziGeometry::Point coords, new_coords;
  coords.init(dim);
  new_coords.init(dim);

  // search the id of the mid point on the top
  for (int iV=0; iV<nV; iV++) {
    
    // get the coords of the node
    mesh0_->node_get_coordinates(iV,&coords);

    // check if coords must be changed
    double dist = 0.;
    for ( int s=0; s<dim-1; ++s ) {
      dist += pow( coords[s]-P0[s], 2 );
    }
    dist = sqrt(dist);

    double alpha  = 1./20.;
    double radius = 0.25;

    // if ok_node is true then change the coords of the node
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
      printf("found: nodeids[%2i]=%3lli  coords=%14.7e --> new_coords=%14.7e\n",
	     i,nodeids[i],coords[dim-1],new_coords[dim-1]);
    }
    LINE(--);
    LINE(--);
  }

  // keep_valid --> check mesh consistency
  bool keep_valid = true;

  // compute the deformed mesh
  mesh0_->deform( nodeids, newpos, keep_valid, &finpos);
}

// move the nodes following a layer profile
void DeformMesh::layer_profile( double ss ) {
  bool verbose(false);

  double alpha = 0.9 ;
  double beta  = 2.  ;

  if ( verbose ) { cout << "step ss = " << ss << endl; }
  
  // get space dimensions
  int dim = mesh0_->space_dimension();

    // set the list of the new position coords
  AmanziGeometry::Point_List newpos, finpos;
  
  // move the mid point on the mesh top
  Entity_ID_List nodeids;

  // number of vertices
  int nV = mesh0_->num_entities(Amanzi::AmanziMesh::NODE, 
				Amanzi::AmanziMesh::OWNED);

  Amanzi::AmanziGeometry::Point coords, new_coords;
  coords.init(dim);
  new_coords.init(dim);

  // search the id of the mid point on the top
  for (int iV=0; iV<nV; iV++) {
    
    // get the coords of the node
    mesh0_->node_get_coordinates(iV,&coords);

    for ( int s=0; s<dim-1; ++s ) {
      new_coords[s] = coords[s];
    }
    double fac = alpha*exp(-beta*ss ) ;
    new_coords[dim-1] = coords[dim-1] * fac;

    // puch back for deform method
    nodeids.push_back(iV);
    newpos.push_back( new_coords );
    
    if ( verbose ) {
      printf("found: iV=%3i,  fac=%14.7e\n",iV,fac);
      printf("coords    =(%14.7e, %14.7e)\n",coords[0],coords[1]);
      printf("new_coords=(%14.7e, %14.7e)\n",new_coords[0],new_coords[1]);
      LINE(-);
    }
  }

  // print the nodes that will be changed
  if ( verbose ) {
    for ( int i=0; i<nodeids.size(); ++i ) {
      printf("found: nodeids[%2i]=%3lli  coords=%14.7e --> new_coords=%14.7e\n",
	     i,nodeids[i],coords[dim-1],new_coords[dim-1]);
    }
    LINE(--);
    LINE(--);
  }

  // keep_valid --> check mesh consistency
  bool keep_valid = true;

  // compute the deformed mesh
  mesh0_->deform( nodeids, newpos, keep_valid, &finpos);
}

void DeformMesh::bell_shaped_profile( double ss, 
				    Amanzi::AmanziGeometry::Point &P0 ) {
  bool verbose(false);

  if ( verbose ) { cout << "step ss = " << ss << endl; }
  
  // get space dimensions
  int dim = mesh0_->space_dimension();
  assert( dim==3 ) ;

  // set the list of the new position coords
  AmanziGeometry::Point_List newpos, finpos;
  
  // move the mid point on the mesh top
  Entity_ID_List nodeids;

  // number of vertices
  int nV = mesh0_->num_entities(Amanzi::AmanziMesh::NODE, 
				Amanzi::AmanziMesh::OWNED);

  Amanzi::AmanziGeometry::Point coords, new_coords;
  coords.init(dim);
  new_coords.init(dim);

  // search the id of the mid point on the top
  cout << "--->modify the coordinates" << endl ;
  for (int iV=0; iV<nV; iV++) {
    
    // get the coords of the node
    mesh0_->node_get_coordinates(iV,&coords);

    // check if coords must be changed
    double dist = 0.;
    for ( int s=0; s<dim-1; ++s ) {
      dist += pow( coords[s]-P0[s], 2 );
    }
    dist = sqrt(dist);

    double alpha  = 1./20.;
    double radius = 12. ;

    // if yes, change the coords of the node
    bool ok_node = dist<10.;
    if ( ok_node ) {
      
      for ( int s=0; s<dim-1; ++s ) {
	new_coords[s] = coords[s];
      }
      double fac = 1.-(1.-ss)*alpha*exp(-pow(dist/radius,2)) ;
      new_coords[dim-1] = coords[dim-1] * fac;

      // push back for deform method
      nodeids.push_back(iV);
      newpos.push_back( new_coords );

      if ( verbose ) {
	printf("found: iV=%3i,  dist=%14.7e,  fac=%14.7e\n",iV,dist,fac);
	printf("coords    =(%14.7e, %14.7e, %14.7e)\n",
	       coords[0],coords[1],coords[2]);
	printf("new_coords=(%14.7e, %14.7e, %14.7e)\n",
	       new_coords[0],new_coords[1],new_coords[2]);
	LINE(-);
      }
    }
  }

  // print the nodes that will be changed
  if ( verbose ) {
    for ( int i=0; i<nodeids.size(); ++i ) {
      printf("found: nodeids[%2i]=%3lli  coords=%14.7e --> new_coords=%14.7e\n",
      	     i,nodeids[i],coords[dim-1],new_coords[dim-1]);
    }
    LINE(--);
    LINE(--);
  }

  // keep_valid --> check mesh consistency
  bool keep_valid = true;

  // compute the deformed mesh
  cout << "--->call deform" << endl ;
  mesh0_->deform( nodeids, newpos, keep_valid, &finpos);
}

// driver routine (to be called by main in UnitTest)
void DeformMesh::bell_shaped_profile(){
  // get space dimensions
  int dim = mesh0_->space_dimension();
  assert( dim==3 ) ;

  // build a bell_shaped profile centered at P0
  Amanzi::AmanziGeometry::Point P0;
  P0.init(dim);
  
  P0[0]=1065; 
  P0[1]= 810; 
  P0[2]=4.56;

  string fname("");
  int kmax = 1;
  for ( int k=0; k<=kmax; ++k ) {
    //loop_monitor(k,kmax) ;
    double ss = double(k)/double(kmax+1);
    ostringstream oss;
    oss << k;
    if      (  0<=k && k<10  ) { fname = string("mesh_0") + oss.str(); }
    else if ( 10<=k && k<100 ) { fname = string("mesh_")  + oss.str(); }
    cout << k << "  " << fname << endl ;
    print_VTK_domain_boundary( fname ) ;
    //print_VTK_unstructured_mesh( fname );
    //---
    bell_shaped_profile(ss,P0);
  }

  // print final mesh
  ostringstream oss;
  oss << kmax+1;
  if      (  0<=kmax+1 && kmax+1<10  ) { fname = string("mesh_0") + oss.str(); }
  else if ( 10<=kmax+1 && kmax+1<100 ) { fname = string("mesh_")  + oss.str(); }
  print_VTK_domain_boundary( fname ) ;
  //print_VTK_unstructured_mesh( fname );
}

// driver routine (to be called by main in UnitTest)
void DeformMesh::layer_profile(){
  // get space dimensions
  int dim = mesh0_->space_dimension();
  assert( dim==3 ) ;

  string fname("");
  int kmax = 20;
  for ( int k=0; k<=kmax; ++k ) {
    //loop_monitor(k,kmax) ;
    double ss = double(k)/double(kmax+1);
    ostringstream oss;
    oss << k;
    if      (  0<=k && k<10  ) { fname = string("mesh_0") + oss.str(); }
    else if ( 10<=k && k<100 ) { fname = string("mesh_")  + oss.str(); }
    cout << k << "  " << fname << endl ;
    //print_VTK_unstructured_mesh( fname );
    print_VTK_domain_boundary( fname ) ;
    //---
    layer_profile(ss);
  }

  // print final mesh
  ostringstream oss;
  oss << kmax+1;
  if      (  0<=kmax+1 && kmax+1<10  ) { fname = string("mesh_0") + oss.str(); }
  else if ( 10<=kmax+1 && kmax+1<100 ) { fname = string("mesh_")  + oss.str(); }
  print_VTK_domain_boundary( fname ) ;
  //print_VTK_unstructured_mesh( fname );
}

// this routine builds the initial mesh from the final one
void DeformMesh::build_the_starting_mesh( Entity_ID_List & newnod ) {
  // get the number of space dimensions
  int dim = mesh0_->space_dimension();
  
  // list of nodes on the top
  Entity_ID_List top_nodeids;

  // number of faces & vertices
  int nF = mesh0_->num_entities(Amanzi::AmanziMesh::FACE, 
				Amanzi::AmanziMesh::OWNED);

  int nV = mesh0_->num_entities(Amanzi::AmanziMesh::NODE, 
				Amanzi::AmanziMesh::OWNED);
  
  bool recompute_flag(false);

  // loop on the mesh faces and get the nodes on the top
  // (nodes may be repeated but this is only to determine the
  // max height)
  for (int iF=0; iF<nF; iF++) {
    
    Entity_ID_List cellids;
    mesh0_->face_get_cells ( iF, Amanzi::AmanziMesh::OWNED, &cellids );
    
    // check if this is a boundary face
    if ( cellids.size()==1 ) {

      // get the normal to the face
      Point norF = mesh0_->face_normal( iF, recompute_flag, cellids[0] );
      
      // check if this face is on the top of the domain
      if ( norF[dim-1]>0 ) {
	Entity_ID_List face_nodeids;
	mesh0_->face_get_nodes ( iF, &face_nodeids ); 
	int nFV = face_nodeids.size();
	for ( int ilV=0; ilV<nFV; ++ilV ) {
	  top_nodeids.push_back( face_nodeids[ilV] );
	}
      }
      
    } // end of --> if ( cellids.size()==1 ) {...
  } // end of --> for (int iF=0; iF<nF; iF++) {...

  int ntV = top_nodeids.size();
  printf("ntV = %i\n",ntV);
  double zmax = 0. ;
  for ( int ilV=0; ilV<ntV; ++ilV ) {

    int iV = top_nodeids[ilV];

    Point ncoord;
    mesh0_->node_get_coordinates ( iV, &ncoord );

    double xV = ncoord[0];
    double yV = ncoord[1];
    double zV = ncoord[2];

    zmax = max( zmax, zV );

    //printf("top_nodeids[%i]=%i, (%14.7e,%14.7e,%14.7e)\n",ilV,iV,xV,yV,zV);
  }
  printf("--------------------\n");
  printf("zmax=%14.7e\n",zmax);

  Amanzi::AmanziGeometry::Point new_coords;
  new_coords.init(dim);

  // set the list of the new position coords
  AmanziGeometry::Point_List newpos, finpos;
  
  // move the mid point on the mesh top
  for ( int ilV=0; ilV<ntV; ++ilV ) {
    
    int iV = top_nodeids[ilV];
    
    Point ncoord;
    mesh0_->node_get_coordinates ( iV, &ncoord );
    
    for ( int s=0; s<dim-1; ++s ) {
      new_coords[s] = ncoord[s];
    }
    new_coords[dim-1] = ncoord[dim-1] + 2.*( zmax-ncoord[dim-1]);

    // push back for deform method
    newnod.push_back(iV);
    newpos.push_back( new_coords );
  }

  // keep_valid --> check mesh consistency
  bool keep_valid = true;

  // compute the deformed mesh
  cout << "--->call deform" << endl ;
  mesh0_->deform( newnod, newpos, keep_valid, &finpos);
}

void DeformMesh::mesh_deformation() {
  bool verbose(false);

  // get the number of space dimensions
  int dim = mesh0_->space_dimension();

  // build the initial mesh and get the list of the reflected nodes
  cout << "--->call build_the_starting_mesh(..)" << endl ;
  Entity_ID_List nodids;

  int nV = mesh0_->num_entities(Amanzi::AmanziMesh::NODE, 
				Amanzi::AmanziMesh::OWNED);

  analyze_final_mesh();

  // set the list of the new position coords
  AmanziGeometry::Point_List newpos, finpos;

  // node coordinates
  Point new_coords, xV0, xV1;
  xV0.init(dim);
  xV1.init(dim);
  new_coords.init(dim);

  // set the number of intermediate steps and the step parameter
  int kmax = 20;
  double dt = 1./double(kmax);

  LINE(---);
  LINE(---);

  // loop on the intermediate steps
  for ( int k=0; k<=kmax; ++k ) {

    cout << endl;
    cout << "k = " << k << endl;

    // set time
    double time = k*dt;

    // reset work arrays
    nodids.resize(0);
    newpos.resize(0);
    finpos.resize(0);
    
    // set the coordinates of the top nodes
    for ( int iV=0; iV<nV; ++iV ) {

      loop_monitor( iV, nV ) ;
      
      mesh0_->node_get_coordinates ( iV, &xV0 );
      mesh1_->node_get_coordinates ( iV, &xV1 );
      
      for ( int s=0; s<dim-1; ++s ) {
	new_coords[s] = xV0[s];
      }
      new_coords[dim-1] = xV0[dim-1]*(1.-time) + xV1[dim-1]*time;
      
      // push back for deform method
      nodids.push_back( iV ) ;
      newpos.push_back( new_coords );

      mesh_->node_set_coordinates(iV,new_coords);
    }

    // keep_valid --> check mesh consistency
    bool keep_valid = true;
    
    // compute the deformed mesh
    cout << endl;
    cout << "--->call deform" << endl ;
    mesh_->deform( nodids, newpos, keep_valid, &finpos);

    // VTK output
    string fname("");
    set_mesh_name(k,fname);
    print_VTK_domain_boundary(string("Mesh-Big/bnd")+fname);
    print_VTK_submesh(string("Mesh-Big/sub")+fname);
    print_VTK_unstructured_mesh(string("Mesh-Big/")+fname);
  }
}

// move the nodes at the top layer of the mesh
void DeformMesh::mesh_deformation_top_nodes() {
  // get the number of space dimensions
  int dim = mesh0_->space_dimension();

  // build the initial mesh and get the list of the reflected nodes
  cout << "--->call build_the_starting_mesh(..)" << endl ;
  Entity_ID_List nodeids;
  build_the_starting_mesh( nodeids );
  int ntV = nodeids.size();

  // set the list of the new position coords
  AmanziGeometry::Point_List newpos, finpos;

  // node coordinates
  Point new_coords, xV0, xV1;
  xV0.init(dim);
  xV1.init(dim);
  new_coords.init(dim);

  // set the number of intermediate steps and the step parameter
  int kmax = 1;
  double dt = 1./double(kmax);

  // loop on the intermediate steps
  for ( int k=0; k<=kmax; ++k ) {

    cout << endl;
    cout << "k = " << k << endl;

    // set time
    double time = k*dt;

    // reset work arrays
    newpos.resize(0);
    finpos.resize(0);
    
    // set the coordinates of the top nodes
    for ( int ilV=0; ilV<ntV; ++ilV ) {
      int iV = nodeids[ilV];
      
      mesh0_->node_get_coordinates ( iV, &xV0 );
      mesh1_->node_get_coordinates ( iV, &xV1 );
      
      for ( int s=0; s<dim-1; ++s ) {
	new_coords[s] = xV0[s];
      }
      new_coords[dim-1] = xV0[dim-1]*(1.-time) + xV1[dim-1]*time;
  
      printf("iV=%i  zV0=%14.7e  zV1=%14.7e\n ",iV,xV0[dim-1],xV1[dim-1]);
    
      // puch back for deform method
      newpos.push_back( new_coords );
    }

    // keep_valid --> check mesh consistency
    bool keep_valid = true;
    
    // compute the deformed mesh
    cout << "--->call deform" << endl ;
    mesh_->deform( nodeids, newpos, keep_valid, &finpos);

    // VTK output
    string fname("");
    set_mesh_name(k,fname);
    print_VTK_submesh(fname);
  }
}

void DeformMesh::set_mesh_name( int k, string &fname ) {
  ostringstream oss;
  oss << k;
  if      (  0<=k && k<10  ) { fname = string("mesh_0") + oss.str(); }
  else if ( 10<=k && k<100 ) { fname = string("mesh_")  + oss.str(); }
}

// get the vertical columns of nodes of the mesh
void DeformMesh::analyze_final_mesh( vector<PCol> & pcol ) {
  bool verbose(false);

  // get the number of space dimensions
  int dim = mesh1_->space_dimension();
  
  // list of nodes on the top
  Entity_ID_List top_nodeids;
  
  // number of vertices
  int nV = mesh1_->num_entities(Amanzi::AmanziMesh::NODE, 
				Amanzi::AmanziMesh::OWNED);
  
  // get the bounding box
  double xmin(+1.e+20), xmax(-1.e+20), ymin(+1.e+20), ymax(-1.e+20);
  for ( int iV=0; iV<nV; ++iV ) {
    Point coords_V;
    mesh1_->node_get_coordinates ( iV, &coords_V );
    xmin = min(xmin,coords_V[0]);
    xmax = max(xmax,coords_V[0]);
    ymin = min(ymin,coords_V[1]);
    ymax = max(ymax,coords_V[1]);
  }
  printf("bbox: (%14.7e,%14.7e)-->(%14.7e,%14.7e)\n",xmin,ymin,xmax,ymax);

  // dx, dy
  double dx=0.25;
  double dy=0.25;

  int nx = int((xmax-xmin)/dx);
  int ny = int((ymax-ymin)/dx);
  
  int size = (nx+1)*(ny+1);
  pcol.resize( size );

  for ( int iV=0; iV<nV; ++iV ) {
    Point coords_V;
    mesh1_->node_get_coordinates ( iV, &coords_V );
    double xV = coords_V[0];
    double yV = coords_V[1];
    int i = int((xV-xmin)/dx);
    int j = int((yV-ymin)/dy);
    int k = j*(nx+1)+i;
    pcol[k].insert(iV,coords_V[2]);
  }

  // print pcol
  if ( verbose ) {
    for ( int ip=0; ip<size; ++ip ) {
      int nk=pcol[ip].get_size();
      LINE(---);
      cout << "ip = " << ip << endl ;
      for ( int k=0; k<nk; ++k ) {
	int    iV = pcol[ip].get_iV(k);
	double zV = pcol[ip].get_zV(k);
	printf("%i %i %14.7e--->\n",k,iV,zV);
      }
    }
  }
}

// analyze the structure of the mesh along the vertical columns
void DeformMesh::analyze_final_mesh() {
  bool verbose(false);
	
  // remap the mesh on columns
  vector<PCol> pcol;
  analyze_final_mesh(pcol);
  
  // get zmax
  double zmax(-1e+20);
  int np = pcol.size();
  for ( int ip=0; ip<np; ++ip ) {
    int nk=pcol[ip].get_size();
    zmax = max( zmax, pcol[ip].get_zV(nk-1) );
  }
  printf("zmax=%14.7e",zmax);
  
  // get the number of space dimensions
  int dim = mesh0_->space_dimension();
  
  // loop 
  int nV = mesh0_->num_entities(Amanzi::AmanziMesh::NODE, 
				Amanzi::AmanziMesh::OWNED);

  // new_coords 
  Amanzi::AmanziGeometry::Point new_coords;
  new_coords.init(dim);
  
  // list of the new position coords
  AmanziGeometry::Point_List newpos, finpos;

  // list of nodes to be moved
  Entity_ID_List newnod;

  // loop on the columns
  for ( int ip=0; ip<np; ++ip ) {
    int nk=pcol[ip].get_size();
    
    // determine the scaling factor
    int    iV_top = pcol[ip].get_iV(nk-1);
    double zV_top = pcol[ip].get_zV(nk-1);
    
    int    iV_bot = pcol[ip].get_iV(0);
    double zV_bot = pcol[ip].get_zV(0);
    
    double new_zV_top = zV_top + 3.*( zmax-zV_top );
    double alpha = (new_zV_top-zV_bot)/(zV_top-zV_bot);

    if ( verbose ) {
      cout << endl;
      LINE(---);
      printf("ip=%i alpha=%14.7e\n",ip,alpha);
      printf("zV_top=%14.7e  new_zV_top=%14.7e\n",zV_top,new_zV_top);
    }

    bool done=false;

    // loop along the column from bottom to top
    double prec_zV=pcol[ip].get_zV(0);
    for ( int k=1; k<nk; ++k ) {
      int iV = pcol[ip].get_iV(k);
      Point coords_V;
      mesh0_->node_get_coordinates ( iV, &coords_V );
      for ( int s=0; s<dim-1; ++s ) {
	new_coords[s] = coords_V[s];
      }
      
      double zV  = pcol[ip].get_zV(k);
      double zVm = pcol[ip].get_zV(k-1);
      double dz  = zV-zVm;

      double new_zV = prec_zV + alpha*dz ;
      new_coords[dim-1] = new_zV;
      
      // push back for deform method
      newnod.push_back( iV );
      newpos.push_back( new_coords );

      mesh0_->node_set_coordinates(iV,new_coords);
      
      if ( verbose ) {
	printf("   k=%i iV=%i\n",k,iV);
	printf("   zV =%14.7e  zVm=%14.7e  dz      =%14.7e\n",zV,zVm,dz);
	printf("   pzV=%14.7e  nzV=%14.7e  alpha*dz=%14.7e\n",prec_zV,new_zV,alpha*dz);
	printf("   alpha=%14.7e\n",alpha);
	cout << endl;
      }

      // reset for the next node
      prec_zV = new_zV; 
    }
  }

  // keep_valid --> check mesh consistency
  bool keep_valid = true;
  
  // check if the mesh deformation is admissible
  // (disabled to be faster)
  bool check_mesh_deformation = false;
  if ( check_mesh_deformation ) {
    cout << "--->call deform" << endl ;
    mesh0_->deform( newnod, newpos, keep_valid, &finpos);
  }
}

}
}
