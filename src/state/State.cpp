#include "State.hpp"
#include "Epetra_Vector.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Mesh.hh"
#include "cell_geometry.hh"
#include "Point.hh"
#include "Geometry.hh"
#include "linear-function.hh"

State::State( int number_of_components_, 
	      Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh_maps_):
  number_of_components(number_of_components_),
  mesh_maps(mesh_maps_)
{
  // create the Eptera_Vector objects
  create_storage();
  create_default_compnames(number_of_components);
};



State::State( Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh_maps_):
  mesh_maps(mesh_maps_)
{
  // this constructor is going to be used in restarts, where we
  // read the number of components from a file before creating
  // storage
};

State::State( Teuchos::ParameterList &parameter_list_,
	      Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh_maps_):
  mesh_maps(mesh_maps_), 
  parameter_list(parameter_list_)
{

  // get the number of component concentrations from the 
  // parameter list
  number_of_components = parameter_list.get<int>("Number of component concentrations");
  
  // create the Eptera_Vector objects

  create_storage();
  initialize_from_parameter_list();
  create_default_compnames(number_of_components);
};


State::~State()
{
  //delete [] (*gravity);
}


void State::create_default_compnames(int n)
{
  compnames.resize(n);
  for (int i=0; i<n; i++)
    {
      std::stringstream ss;
      ss << "component " << i;
      compnames[i] = ss.str();
    }
}



void State::initialize_from_parameter_list()
{
  int num_blocks = parameter_list.get<int>("Number of mesh blocks");

  double u[3];
  u[0] = parameter_list.get<double>("Gravity x");		   
  u[1] = parameter_list.get<double>("Gravity y");		   
  u[2] = parameter_list.get<double>("Gravity z");		   
  set_gravity(u);

  set_zero_total_component_concentration();
  set_water_density(parameter_list.get<double>("Constant water density"));
  // set_water_saturation(parameter_list.get<double>("Constant water saturation"));
  set_viscosity(parameter_list.get<double>("Constant viscosity"));
  
  for (int nb=1; nb<=num_blocks; nb++) {
    
    std::stringstream pname;
    pname << "Mesh block " << nb;

    Teuchos::ParameterList sublist = parameter_list.sublist(pname.str());

    std::string region = sublist.get<std::string>("Region");


    // initialize the arrays with some constants from the input file
    set_porosity(sublist.get<double>("Constant porosity"),region);
    set_permeability(sublist.get<double>("Constant permeability"),region);
    
    u[0] = sublist.get<double>("Constant Darcy flux x",0.0);
    u[1] = sublist.get<double>("Constant Darcy flux y",0.0);
    u[2] = sublist.get<double>("Constant Darcy flux z",0.0);
    set_darcy_flux(u, region);
    
    // set the pressure
    if (sublist.isSublist("uniform pressure"))
      {
	const Teuchos::ParameterList&  unif_p_list = sublist.sublist("uniform pressure");
	set_uniform_pressure( unif_p_list, region );
      }
    else if (sublist.isSublist("linear pressure"))
      {
	const Teuchos::ParameterList&  lin_p_list = sublist.sublist("linear pressure");
	set_linear_pressure( lin_p_list, region );	
      }
    else
      {
	// maybe throw an exception
      }

    // set the saturation
    // set the pressure
    if (sublist.isSublist("uniform saturation"))
      {
	const Teuchos::ParameterList&  unif_s_list = sublist.sublist("uniform saturation");
	set_uniform_saturation( unif_s_list, region );
      }
    else if (sublist.isSublist("linear saturation"))
      {
	const Teuchos::ParameterList&  lin_s_list = sublist.sublist("linear saturation");
	set_uniform_saturation( lin_s_list, region );	
      }
    else
      {
	// maybe throw an exception
      }    

    // read the component names if they are spelled out
    Teuchos::Array<std::string> comp_names;
    if (parameter_list.isParameter("Component Names")) 
      {
	comp_names = parameter_list.get<Teuchos::Array<std::string> >("Component Solutes");
      }
    else
      {
	comp_names.resize(number_of_components);
	for (int i=0; i<number_of_components; i++)
	  {
	    std::stringstream ss;
	    ss << "comp " << i;
	    comp_names[i] = ss.str();
	  }
      }
    // now create the map
    for (int i=0; i<comp_names.size(); i++)
      {
	comp_no[comp_names[i]] = i;
      }


    // read the component concentrations from the xml file
    // and initialize them in mesh block mesh_block_ID
    double tcc_const[number_of_components];
    for (int nc=0; nc<number_of_components; nc++) {
      std::stringstream s; 
      s << "Constant component concentration " << nc;

      tcc_const[nc] = sublist.get<double>(s.str());
    }
    set_total_component_concentration(tcc_const, region);
  }
}




void State::create_storage ()
{
  // create the Eptera_Vector objects

  water_density =    Teuchos::rcp( new Epetra_Vector( mesh_maps->cell_map(false) ) );
  pressure =         Teuchos::rcp( new Epetra_Vector( mesh_maps->cell_map(false) ) );
  darcy_flux =       Teuchos::rcp( new Epetra_Vector( mesh_maps->face_map(false) ) );
  porosity =         Teuchos::rcp( new Epetra_Vector( mesh_maps->cell_map(false) ) );
  water_saturation = Teuchos::rcp( new Epetra_Vector( mesh_maps->cell_map(false) ) ); 
  permeability     = Teuchos::rcp( new Epetra_Vector( mesh_maps->cell_map(false) ) ); 
  total_component_concentration 
    = Teuchos::rcp( new Epetra_MultiVector( mesh_maps->cell_map(false), number_of_components ) );  
  darcy_velocity   = Teuchos::rcp( new Epetra_MultiVector( mesh_maps->cell_map(false), 3));

  density =   Teuchos::rcp(new double);
  viscosity = Teuchos::rcp(new double);
  gravity =   Teuchos::rcp(new double*);
  *gravity = new double[3];
}



void State::set_time ( double new_time ) {

  time = new_time;

}

void State::set_cycle ( int new_cycle ) {

  cycle = new_cycle;

}

void State::set_number_of_components ( const int n )
{
  number_of_components = n;
}


void State::update_total_component_concentration(Teuchos::RCP<Epetra_MultiVector> new_tcc) 
{
  *total_component_concentration = *new_tcc;

}

void State::update_total_component_concentration(const Epetra_MultiVector& new_tcc) 
{
  *total_component_concentration = new_tcc;

}

void State::update_darcy_flux(const Epetra_Vector &new_darcy_flux)
{
  *darcy_flux = new_darcy_flux;
}

void State::update_pressure(const Epetra_Vector &new_pressure)
{
  *pressure = new_pressure;
}



void State::advance_time(double dT)
{
  time = time + dT;
}



void State::set_cell_value_in_region(const double& value, Epetra_Vector& v, 
				     const std::string& region)
{
  if (!mesh_maps->valid_set_name(region,Amanzi::AmanziMesh::CELL)) {
    throw std::exception();
  }
  
  unsigned int mesh_block_size = mesh_maps->get_set_size(region,
							 Amanzi::AmanziMesh::CELL,
							 Amanzi::AmanziMesh::OWNED);

  std::vector<unsigned int> cell_ids(mesh_block_size);
  
  //mesh_maps->get_set(mesh_block_id, Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::OWNED,
  //cell_ids.begin(),cell_ids.end());
  
  mesh_maps->get_set_entities(region, Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::OWNED,
			      &cell_ids);
  

  for( std::vector<unsigned int>::iterator c = cell_ids.begin(); 
       c != cell_ids.end();  c++) {
    v[*c] = value;  
  } 

}


void State::set_cell_value_in_region(const Amanzi::Function& fun, Epetra_Vector& v, 
				     const std::string& region)
{
  if (!mesh_maps->valid_set_name(region,Amanzi::AmanziMesh::CELL)) {
    throw std::exception();
  }
  
  unsigned int mesh_block_size = mesh_maps->get_set_size(region,
							 Amanzi::AmanziMesh::CELL,
							 Amanzi::AmanziMesh::OWNED);
  std::vector<unsigned int> cell_ids(mesh_block_size);
  mesh_maps->get_set_entities(region, Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::OWNED,
			      &cell_ids);
  
  for( std::vector<unsigned int>::iterator c = cell_ids.begin(); c != cell_ids.end();  c++) {
    // get location of centroid for cell *c
    Amanzi::AmanziGeometry::Point p( mesh_maps->cell_centroid(*c) );    
    v[*c] = fun( &p[0] );  
  } 

}


void State::set_cell_value_in_mesh_block(double value, Epetra_Vector &v, 
				    int mesh_block_id)
{
  if (!mesh_maps->valid_set_id(mesh_block_id,Amanzi::AmanziMesh::CELL)) {
    throw std::exception();
  }
  
  unsigned int mesh_block_size = mesh_maps->get_set_size(mesh_block_id,
							 Amanzi::AmanziMesh::CELL,
							 Amanzi::AmanziMesh::OWNED);

  std::vector<unsigned int> cell_ids(mesh_block_size);
  
  mesh_maps->get_set(mesh_block_id, Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::OWNED,
		     cell_ids.begin(),cell_ids.end());
  
  for( std::vector<unsigned int>::iterator c = cell_ids.begin(); 
       c != cell_ids.end();  c++) {
    v[*c] = value;  
  } 

}

void State::set_darcy_flux( const double* u, const int mesh_block_id )
{
  int  i, f;
  double x[4][3], normal[3];

  // Epetra_Map face_map = mesh_maps->face_map(false);

  if (!mesh_maps->valid_set_id(mesh_block_id,Amanzi::AmanziMesh::CELL)) {
    throw std::exception();
  }
  
  unsigned int mesh_block_size = mesh_maps->get_set_size(mesh_block_id,
							 Amanzi::AmanziMesh::CELL,
							 Amanzi::AmanziMesh::OWNED);
  
  std::vector<unsigned int> cell_ids(mesh_block_size);
  


  mesh_maps->get_set(mesh_block_id, Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::OWNED,
		     cell_ids.begin(),cell_ids.end());

  
  for( std::vector<unsigned int>::iterator c = cell_ids.begin(); 
       c != cell_ids.end();  c++) {
    
    std::vector<unsigned int> cface(6);
    mesh_maps->cell_to_faces(*c, cface.begin(), cface.end());

    for (std::vector<unsigned int>::iterator f = cface.begin();
	 f != cface.end(); f++) {
     
      if (mesh_maps->face_map(false).MyLID(*f) ) {
	
	mesh_maps->face_to_coordinates( *f, (double*) x, (double*) x+12 );
	
	cell_geometry::quad_face_normal(x[0], x[1], x[2], x[3], normal);
	
	(*darcy_flux)[*f] = u[0] * normal[0] + u[1] * normal[1] + u[2] * normal[2];
      }
      
    }
  }
}

void State::set_darcy_flux( const double* u, const std::string region )
{
  int  i, f;
  double x[4][3], normal[3];

  if (!mesh_maps->valid_set_name(region,Amanzi::AmanziMesh::CELL)) {
    throw std::exception();
  }
  
  unsigned int mesh_block_size = mesh_maps->get_set_size(region,
							 Amanzi::AmanziMesh::CELL,
							 Amanzi::AmanziMesh::OWNED);
  
  std::vector<unsigned int> cell_ids(mesh_block_size);
  


  mesh_maps->get_set_entities(region, Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::OWNED,
			    &cell_ids);

  
  for( std::vector<unsigned int>::iterator c = cell_ids.begin(); 
       c != cell_ids.end();  c++) {
    
    std::vector<unsigned int> cface(6);
    mesh_maps->cell_to_faces(*c, cface.begin(), cface.end());

    for (std::vector<unsigned int>::iterator f = cface.begin();
	 f != cface.end(); f++) {
      
      if (mesh_maps->face_map(false).MyLID(*f) ) {
	
	mesh_maps->face_to_coordinates( *f, (double*) x, (double*) x+12 );
	
	cell_geometry::quad_face_normal(x[0], x[1], x[2], x[3], normal);
	
	(*darcy_flux)[*f] = u[0] * normal[0] + u[1] * normal[1] + u[2] * normal[2];
      }
      
    }
  }
}


void State::set_water_density( const double wd )
{
  water_density->PutScalar(wd);
  *density = wd;
}


void State::set_water_saturation( const double ws )
{
  water_saturation->PutScalar(ws); 
}


void State::set_porosity( const double phi )
{
  porosity->PutScalar(phi);  
}

void State::set_porosity( const double phi, const int mesh_block_id )
{
  set_cell_value_in_mesh_block(phi,*porosity,mesh_block_id);
}

void State::set_porosity( const double phi, const std::string region )
{
  set_cell_value_in_region(phi,*porosity,region);
}

void State::set_zero_total_component_concentration()
{
  total_component_concentration->PutScalar(0.0);
}


void State::set_total_component_concentration( const double* conc, const int mesh_block_id )
{
  for (int nc=0; nc<number_of_components; nc++) {
    set_cell_value_in_mesh_block(conc[nc], *(*total_component_concentration)(nc),mesh_block_id);
  }
}

void State::set_total_component_concentration( const double* conc, const std::string region )
{
  for (int nc=0; nc<number_of_components; nc++) {
    set_cell_value_in_region(conc[nc], *(*total_component_concentration)(nc),region);
  }
}


void State::set_permeability( const double kappa )
{
  permeability->PutScalar(kappa);
}


void State::set_permeability( const double kappa, const int mesh_block_id)
{
  set_cell_value_in_mesh_block(kappa,*permeability,mesh_block_id);
}

void State::set_permeability( const double kappa, const std::string region)
{
  set_cell_value_in_region(kappa,*permeability,region);
}

void State::set_viscosity(const double mu)
{
  *viscosity = mu;
}


void State::set_gravity(const double *g)
{
  (*gravity)[0] = g[0];
  (*gravity)[1] = g[1];
  (*gravity)[2] = g[2];
}


// void State::init_restart( )
// {
//   int rank = mesh_maps->get_comm()->MyPID();  

//   unsigned int num_nodes = mesh_maps->num_entities(Amanzi::AmanziMesh::NODE,Amanzi::AmanziMesh::OWNED);
//   unsigned int num_cells = mesh_maps->num_entities(Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::OWNED);  
//   unsigned int num_faces = mesh_maps->num_entities(Amanzi::AmanziMesh::FACE,Amanzi::AmanziMesh::OWNED); 

//   int nums[3];
//   int dummy[3];
//   dummy[0] = num_nodes;
//   dummy[1] = num_cells;  
//   dummy[2] = num_faces;

//   mesh_maps->get_comm()->SumAll(dummy,nums,3);  


//   // make the all to one map
//   if (rank == 0) {
//     int *gids = new int[nums[1]];
//     for (int i=0; i<nums[1]; i++) gids[i] = i;
//     all_to_one_cell_map = Teuchos::rcp(new Epetra_Map(nums[1],nums[1],gids,0, * mesh_maps->get_comm() ));
//     delete [] gids;
    
//     gids = new int[nums[0]];
//     for (int i=0; i<nums[0]; i++) gids[i] = i;    
//     all_to_one_node_map = Teuchos::rcp(new Epetra_Map(nums[0],nums[0],gids,0, * mesh_maps->get_comm() ));
//     delete [] gids;

//     // gids = new int[nums[2]];
//     // for (int i=0; i<nums[2]; i++) gids[i] = i;
//     // all_to_one_face_map = Teuchos::rcp(new Epetra_Map(nums[2],nums[2],gids,0, * mesh_maps->get_comm() ));
//     // delete [] gids;  

//     int max_gid = mesh_maps->face_map(false).MaxAllGID();
//     int min_gid = mesh_maps->face_map(false).MinAllGID();
//     gids = new int [max_gid-min_gid+1];
//     for (int i=0; i<max_gid-min_gid+1; i++) gids[i] = min_gid+i;
//     all_to_one_face_map = Teuchos::rcp(new Epetra_Map(max_gid+1,max_gid+1,gids,0, * mesh_maps->get_comm() ));
    

//   } else {
//     int *gids;
//     int max_gid = mesh_maps->face_map(false).MaxAllGID();
//     all_to_one_cell_map = Teuchos::rcp(new Epetra_Map(nums[1],0,gids,0, * mesh_maps->get_comm() ) );
//     all_to_one_node_map = Teuchos::rcp(new Epetra_Map(nums[0],0,gids,0, * mesh_maps->get_comm() ) );
//     all_to_one_face_map = Teuchos::rcp(new Epetra_Map(max_gid+1,0,gids,0, * mesh_maps->get_comm() ) );
//   }

//   // make the all to one exporters 
//   all_to_one_cell_export = Teuchos::rcp(new Epetra_Export(mesh_maps->cell_map(false), *all_to_one_cell_map) );
//   all_to_one_node_export = Teuchos::rcp(new Epetra_Export(mesh_maps->node_map(false), *all_to_one_node_map) );  
//   all_to_one_face_export = Teuchos::rcp(new Epetra_Export(mesh_maps->face_map(false), *all_to_one_face_map) );

// }


// void State::write_restart ( std::string filename )
// {
//   int rank = mesh_maps->get_comm()->MyPID();

//   Epetra_Vector PE0(*all_to_one_cell_map);
//   Epetra_Vector PEF(*all_to_one_face_map);
//   Amanzi::HDF5 restart_output;

//   if (rank == 0) {
//     restart_output.setTrackXdmf(false);
//     restart_output.createDataFile(filename);
//   }
       
//   PE0.Export( *water_density, *all_to_one_cell_export, Insert);
//   if (rank == 0) {
//     restart_output.writeCellData(PE0, "water_density");
//   }

//   PE0.Export( *pressure, *all_to_one_cell_export, Insert);   
//   if (rank == 0) {
//     restart_output.writeCellData(PE0, "pressure");
//   }

//   //darcy_flux->Print(std::cout);
//   PEF.Export( *darcy_flux, *all_to_one_face_export, Insert);
//   //PEF.Print(std::cout);
//   if (rank == 0) {
//     restart_output.writeCellData(PEF, "darcy_flux");
//   }

//   PE0.Export( *porosity, *all_to_one_cell_export, Insert);
//   if (rank == 0) {
//     restart_output.writeCellData(PE0, "porosity");
//   }  

//   PE0.Export( *water_saturation, *all_to_one_cell_export, Insert);
//   if (rank == 0) {
//     restart_output.writeCellData(PE0, "water_saturation");
//   }   
  
//   PE0.Export( *permeability, *all_to_one_cell_export, Insert);
//   if (rank == 0) {
//     restart_output.writeCellData(PE0, "permeability");
//   }   

  
//   for (int i=0; i<darcy_velocity->NumVectors(); i++)
//     {
//       std::stringstream fnss;
//       fnss << "darcy_velocity-" << i;
      
//       PE0.Export( *(*darcy_velocity)(0), *all_to_one_cell_export, Insert);
//       if (rank == 0) {
// 	restart_output.writeCellData(PE0, fnss.str() );
//       }   
//     }
  
//   // for (int i=0; i<total_component_concentration->NumVectors(); i++)
//   //   {
//   //     std::stringstream fnss;
//   //     fnss << "total_component_concentration-" << i;
      
//   //     PE0.Export( *(*total_component_concentration)(0), *all_to_one_cell_export, Insert);
//   //     if (rank == 0) {
//   // 	restart_output.writeCellData(PE0, fnss.str() );
//   //     }   
//   //   }
// }

// void State::read_restart ( std::string filename )
// {
//   int rank = mesh_maps->get_comm()->MyPID();
//   Epetra_Vector PE0(*all_to_one_cell_map); 
//   Epetra_Vector PEF(*all_to_one_face_map); 
  
//   Amanzi::HDF5 *restart_output = new Amanzi::HDF5();

//   if (rank == 0) {
//     restart_output->setTrackXdmf(false);
//     restart_output->setH5DataFilename(filename);
//   }

//   if (rank == 0) {
//     restart_output->readData(PE0, "water_density");
//   }
//   water_density->Import( PE0, *all_to_one_cell_export, Insert);

//   if (rank == 0) {
//     restart_output->readData(PE0, "pressure");
//   }
//   pressure->Import( PE0, *all_to_one_cell_export, Insert);   

//   if (rank == 0) {
//     restart_output->readData(PEF, "darcy_flux");
//   }
//   darcy_flux->Import( PEF, *all_to_one_face_export, Insert);

//   if (rank == 0) {
//     restart_output->readData(PE0, "porosity");
//   }  
//   porosity->Import( PE0, *all_to_one_cell_export, Insert);

//   if (rank == 0) {
//     restart_output->readData(PE0, "water_saturation");
//   }   
//   water_saturation->Import( PE0, *all_to_one_cell_export, Insert);
  
//   if (rank == 0) {
//     restart_output->readData(PE0, "permeability");
//   }   
//   permeability->Import( PE0, *all_to_one_cell_export, Insert);

  
//   for (int i=0; i<darcy_velocity->NumVectors(); i++)
//     {
//       std::stringstream fnss;
//       fnss << "darcy_velocity-" << i;
      
//       if (rank == 0) {
// 	restart_output->readData(PE0, fnss.str() );
//       }   
//       (*darcy_velocity)(i)->Import( PE0, *all_to_one_cell_export, Insert);
//     }
  
//   // for (int i=0; i<total_component_concentration->NumVectors(); i++)
//   //   {
//   //     std::stringstream fnss;
//   //     fnss << "total_component_concentration-" << i;
      
//   //     if (rank == 0) {
//   // 	restart_output->readData(PE0, fnss.str() );
//   //     }   
//   //     PE0.Import( *(*total_component_concentration)(0), *all_to_one_cell_export, Insert);
//   //   }

// }


double State::water_mass()
{
  // compute the total mass of water in the domain
  
  Epetra_Vector wm ( *water_saturation );
  
  wm.Multiply(1.0,*water_density,wm,0.0);
  wm.Multiply(1.0,*porosity,wm,0.0);
  
  Epetra_Vector cell_volume( mesh_maps->cell_map(false) );
  
  for (int i=0; i<(mesh_maps->cell_map(false)).NumMyElements(); i++)
    {
      cell_volume[i] = mesh_maps->cell_volume(i);
    }
	 
  wm.Multiply(1.0,cell_volume,wm,0.0);
       
  double mass;
  wm.Norm1(&mass);
       
  return mass;
}


double State::point_value(const std::string& point_region, const std::string& name)
{
  if (!mesh_maps->valid_set_name(point_region,Amanzi::AmanziMesh::CELL)) 
    {
      // throw
    }
  

  unsigned int mesh_block_size = mesh_maps->get_set_size(point_region,
							 Amanzi::AmanziMesh::CELL,
							 Amanzi::AmanziMesh::OWNED);
  
  if (mesh_block_size > 1)
    {
      // throw
    }
  
  double value(0.0);

  if (mesh_block_size == 1)
    {
      std::vector<unsigned int> cell_ids(mesh_block_size);
      
      mesh_maps->get_set_entities(point_region, Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::OWNED,
				  &cell_ids);
      
      // extract the value if it is a component
      if ( comp_no.find(name) != comp_no.end() )
	{
	  value =   (*(*total_component_concentration)( comp_no[name] )) [cell_ids[0]];
	}

      // extract other point information
      // ...
      
    }  


  // syncronize the result across processors
  
  double result;
  mesh_maps->get_comm()->SumAll(&value,&result,1);


  return result;
}



void State::set_darcy_flux ( const Epetra_Vector& darcy_flux_ )
{
  *darcy_flux = darcy_flux_;
}

void State::set_water_saturation ( const Epetra_Vector& water_saturation_ ) 
{
  *water_saturation = water_saturation_;
};

void State::set_porosity ( const Epetra_Vector& porosity_ ) 
{
  *porosity = porosity_;
}; 

void State::set_permeability ( const Epetra_Vector& permeability_ ) 
{
  *permeability = permeability_;
};

void State::set_pressure ( const Epetra_Vector& pressure_ ) 
{
  *pressure = pressure_;
};

void State::set_water_density ( const Epetra_Vector& water_density_ ) 
{
  *water_density = water_density_;
};

void State::set_darcy_velocity ( const Epetra_MultiVector& darcy_velocity_ ) 
{
  *darcy_velocity = darcy_velocity_;
};

void State::set_total_component_concentration ( const Epetra_MultiVector& total_component_concentration_ ) 
{
  *total_component_concentration = total_component_concentration_;
};

void State::set_uniform_pressure ( const Teuchos::ParameterList& unif_p_list, const std::string& region )
{
  // get value from paramter list
  const double value = unif_p_list.get<double>("value");
  set_cell_value_in_region( value, *pressure, region );
};

void State::set_linear_pressure ( const Teuchos::ParameterList& lin_p_list, const std::string& region )
{
  // get parameters from parameter list
  const double ref_value = lin_p_list.get<double>("reference value");
  const Teuchos::Array<double>& ref_coord = lin_p_list.get<Teuchos::Array<double> >("reference coordinate");
  const Teuchos::Array<double>& gradient = lin_p_list.get<Teuchos::Array<double> >("gradient");

  // create function
  Amanzi::LinearFunction lin_p(ref_value, gradient.toVector(), ref_coord.toVector());

  set_cell_value_in_region ( lin_p, *pressure, region );

};

void State::set_uniform_saturation ( const Teuchos::ParameterList& unif_s_list, const std::string& region )
{
  // get value from paramter list
  const double value = unif_s_list.get<double>("value");
  set_cell_value_in_region( value, *water_saturation, region );
};

void State::set_linear_saturation ( const Teuchos::ParameterList& lin_s_list, const std::string& region )
{
  // get parameters from parameter list
  const double ref_value = lin_s_list.get<double>("reference value");
  const Teuchos::Array<double>& ref_coord = lin_s_list.get<Teuchos::Array<double> >("reference coordinate");
  const Teuchos::Array<double>& gradient = lin_s_list.get<Teuchos::Array<double> >("gradient");

  // create function
  Amanzi::LinearFunction lin_p(ref_value, gradient.toVector(), ref_coord.toVector());

  set_cell_value_in_region ( lin_p, *water_saturation, region );
  


};


void State::write_vis(Amanzi::Vis& vis)
{
  if (vis.dump_requested(get_cycle()) && !vis.is_disabled() )
    {	  
      // create the new time step...
      vis.create_timestep(get_time(),get_cycle());
	  
      // dump all the state vectors into the file
      vis.write_vector(*get_pressure(), "pressure");
      vis.write_vector(*get_porosity(),"porosity");
      vis.write_vector(*get_water_saturation(),"water saturation");
      vis.write_vector(*get_water_density(),"water density");
      vis.write_vector(*get_permeability(),"permeability");

      std::vector<std::string> names(3);
      names[0] = "darcy velocity x";
      names[1] = "darcy velocity y";
      names[1] = "darcy velocity z";
      vis.write_vector(*get_darcy_velocity(), names);

      // write component data
      vis.write_vector( *get_total_component_concentration(), compnames);
    }
}


      
void State::write_vis(Amanzi::Vis& vis, Epetra_MultiVector *auxdata, std::vector<std::string>& auxnames)
{
  write_vis(vis);
  
  if (vis.dump_requested(get_cycle()) && !vis.is_disabled() )
    {
      // write auxillary data
      if (auxdata != NULL) 
	{
	  vis.write_vector( *auxdata , auxnames);
	}
      
      vis.finalize_timestep(); 
    }
}


void State::set_compnames(std::vector<std::string>& compnames_)
{
  compnames = compnames_;
}
