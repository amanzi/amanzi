#include "State.hpp"
#include "Epetra_Vector.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Teuchos_MPISession.hpp"
#include "Mesh.hh"
#include "cell_geometry.hh"
#include "hdf5_mesh.hh"
extern "C" {
#include "gmvwrite.h"
}


State::State( int number_of_components_,
	      Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh_maps_):
  number_of_components(number_of_components_),
  mesh_maps(mesh_maps_)
{
  // create the Eptera_Vector objects

  create_storage();

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
  
  init_restart();

};


State::~State()
{
  //delete [] (*gravity);
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
  set_water_saturation(parameter_list.get<double>("Constant water saturation"));
  set_viscosity(parameter_list.get<double>("Constant viscosity"));
  
  for (int nb=1; nb<=num_blocks; nb++) {
    
    std::stringstream pname;
    pname << "Mesh block " << nb;

    Teuchos::ParameterList sublist = parameter_list.sublist(pname.str());

    int mesh_block_ID = sublist.get<int>("Mesh block ID");

    if (!mesh_maps->valid_set_id(mesh_block_ID,Amanzi::AmanziMesh::CELL)) {
      // there is an inconsistency in the xml input file, report and die
      
      int myrank = Teuchos::MPISession::getRank();

      if (myrank == 0) {
	std::cerr << "State::initialize_from_parameter_list... the mesh block with ID ";
	std::cerr << mesh_block_ID << " does not exist in the mesh" << std::endl;

	// get the mesh block IDs 
	int num_blks = mesh_maps->num_sets(Amanzi::AmanziMesh::CELL);
	std::vector<unsigned int> setids(num_blks);
	mesh_maps->get_set_ids(Amanzi::AmanziMesh::CELL,setids.begin(),setids.end());
	std::cerr << "valid mesh block IDs are: ";
	for (int i=0; i<num_blks; i++) std::cerr << setids[i] << " ";
	std::cerr << std::endl;


	throw std::exception();
      }
    }
	

    // initialize the arrays with some constants from the input file
    set_porosity(sublist.get<double>("Constant porosity"),mesh_block_ID);
    set_permeability(sublist.get<double>("Constant permeability"),mesh_block_ID);
    
    u[0] = sublist.get<double>("Constant Darcy flux x",0.0);
    u[1] = sublist.get<double>("Constant Darcy flux y",0.0);
    u[2] = sublist.get<double>("Constant Darcy flux z",0.0);
    set_darcy_flux(u, mesh_block_ID);
    
    // read the component concentrations from the xml file
    // and initialize them in mesh block mesh_block_ID
    double tcc_const[number_of_components];
    for (int nc=0; nc<number_of_components; nc++) {
      std::stringstream s; 
      s << "Constant component concentration " << nc;

      tcc_const[nc] = sublist.get<double>(s.str());
    }
    set_total_component_concentration(tcc_const, mesh_block_ID);
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


void State::set_permeability( const double kappa )
{
  permeability->PutScalar(kappa);
}


void State::set_permeability( const double kappa, const int mesh_block_id)
{
  set_cell_value_in_mesh_block(kappa,*permeability,mesh_block_id);
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


void State::write_gmv ( std::string filename )
{
  using namespace std;

  gmvwrite_openfile_ir_ascii( (char*) filename.c_str(), 4, 8 );
  
  // first write the node data
  unsigned int num_nodes = mesh_maps->count_entities(Amanzi::AmanziMesh::NODE,Amanzi::AmanziMesh::OWNED);

  double *x = new double [num_nodes];
  double *y = new double [num_nodes];
  double *z = new double [num_nodes];

  std::vector<double> xc(3);
  for (unsigned int i=0; i < num_nodes; i++) {
    mesh_maps->node_to_coordinates(i,xc.begin(),xc.end());

    x[i] = xc[0];
    y[i] = xc[1];
    z[i] = xc[2];
  }
  gmvwrite_node_data(&num_nodes, x, y, z);

  delete x;
  delete y;
  delete z;
  

  // write the cell data
  unsigned int num_cells = mesh_maps->count_entities(Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::OWNED);
  
  gmvwrite_cell_header(&num_cells);
  
  std::vector<unsigned int> xh(8); 
  unsigned int xh_[8];
  for (int i=0; i<num_cells; i++) {
    mesh_maps->cell_to_nodes(i,xh.begin(),xh.end());
    for (int j=0; j<8; j++) xh_[j] = xh[j]+1;

    char cell_type [] = "phex8";
    
    
    gmvwrite_cell_type(cell_type,8,xh_);
  }
  

  
  // write the side sets, we only support 100 side sets
  // with ids less than 100.

  // how many side sets are there?
  unsigned int num_side_sets = mesh_maps->num_sets(Amanzi::AmanziMesh::FACE);
  if (num_side_sets > 0) {
    vector<unsigned int> ssids(num_side_sets);
    mesh_maps->get_set_ids(Amanzi::AmanziMesh::FACE, ssids.begin(), ssids.end());
    
    vector<unsigned int> side_set;
    int total_num_surfaces=0;

    // figure out how many total faces are in side sets
    for (vector<unsigned int>::const_iterator it = ssids.begin(); it != ssids.end(); it++) {
      unsigned int set_size = mesh_maps->get_set_size(*it, Amanzi::AmanziMesh::FACE, Amanzi::AmanziMesh::OWNED);
      total_num_surfaces += set_size;
    }

    if (total_num_surfaces > 0) {
      // now start writing side set information
      gmvwrite_surface_header(&total_num_surfaces);
      
      // loop over the side sets and write all faces belonging to side sets
      for (vector<unsigned int>::const_iterator it = ssids.begin(); it != ssids.end(); it++) {  
	unsigned int set_size = mesh_maps->get_set_size(*it, Amanzi::AmanziMesh::FACE, Amanzi::AmanziMesh::OWNED);
	side_set.resize(set_size);
	
	mesh_maps->get_set(*it,Amanzi::AmanziMesh::FACE,Amanzi::AmanziMesh::OWNED,side_set.begin(),side_set.end());
	
	std::vector<unsigned int> nodes(4);
	unsigned int nodes_[4];
	for (vector<unsigned int>::const_iterator sit = side_set.begin(); sit != side_set.end(); sit++) {
	  mesh_maps->face_to_nodes(*sit,nodes.begin(),nodes.end());
	  for (int j=0; j<4; j++) nodes_[j] = nodes[j]+1;
	  gmvwrite_surface_data(4,nodes_);
	}
      }
      
      
      // now we generate an array of flags, one flag stands for one side set
      gmvwrite_surfflag_header();
      
      char ssname [] = "SideSets";
      gmvwrite_surfflag_name(ssname,num_side_sets);
      
      char sssubname [] = "SideSet_000";

      for (vector<unsigned int>::const_iterator it=ssids.begin(); 
	   it!=ssids.end(); it++) {
	sssubname[8] = '0' + (*it)/100;
	sssubname[9] = '0' + ((*it)/10)%10;
	sssubname[10] = '0' + (*it)%10;

	gmvwrite_surfflag_subname(sssubname);
      }
      
      
      // this is the array of flags
      int *sflagdata = new int[total_num_surfaces];
      unsigned int is=0;
      unsigned int tag=1;
      
      // loop over the side sets
      for (vector<unsigned int>::const_iterator it = ssids.begin(); it != ssids.end(); it++) {  
	unsigned int set_size = mesh_maps->get_set_size(*it, Amanzi::AmanziMesh::FACE, Amanzi::AmanziMesh::OWNED);
	side_set.resize(set_size);
	
	// loop over the faces in the current side set and
	// store the flag for the current side set 
 	for (int i=0; i<set_size; i++) {
	  sflagdata[is] = tag;
	  is++;
	}
	tag++;
      }
      // write the flag array
      gmvwrite_surfflag_data(sflagdata);
      
      delete sflagdata;
      gmvwrite_surfflag_endflag();
    }
    
  }
  

  // write element blocks

  unsigned int num_element_blocks = mesh_maps->num_sets(Amanzi::AmanziMesh::CELL);
  if (num_element_blocks > 0) {
    vector<unsigned int> ebids(num_element_blocks);
    mesh_maps->get_set_ids(Amanzi::AmanziMesh::CELL, ebids.begin(), ebids.end());
    
    vector<unsigned int> element_block;
    int total_num_cells=0;

    // figure out how many total cells are in element blocks
    for (vector<unsigned int>::const_iterator it = ebids.begin(); it != ebids.end(); it++) {
      unsigned int set_size = mesh_maps->get_set_size(*it, Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::OWNED);
      total_num_cells += set_size;
    }

    if (total_num_cells > 0) {
      // now start writing side set information
      gmvwrite_flag_header();

      char flagname [] = "ElementBlocks";
      gmvwrite_flag_name(flagname,num_element_blocks,0);


      // now we generate an array of flags, one flag stands for one side set
      char ebsubname [] = "ElementBlock_000";

      for (vector<unsigned int>::const_iterator it=ebids.begin(); 
	   it!=ebids.end(); it++) {
	ebsubname[13] = '0' + (*it)/100;
	ebsubname[14] = '0' + ((*it)/10)%10;
	ebsubname[15] = '0' + (*it)%10;

	gmvwrite_flag_subname(ebsubname);
      }
      
      // loop over the element blocks and write all cells belonging to element blocks

      // this is the array of flags
      int *ebflagdata = new int[total_num_cells];
      unsigned int is=0;
      unsigned int tag=1;
      
      // loop over the selement blocks sets
      for (vector<unsigned int>::const_iterator it = ebids.begin(); it != ebids.end(); it++) {  
	unsigned int set_size = mesh_maps->get_set_size(*it, Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::OWNED);
	
	// loop over the cells in the current element block set and
	// store the flag for the current side set 
 	for (int i=0; i<set_size; i++) {
	  ebflagdata[is] = tag;
	  is++;
	}
	tag++;
      }
      // write the flag array
      gmvwrite_flag_data(0,ebflagdata);

      delete ebflagdata;

      gmvwrite_flag_endflag();
    }
    
  }  
 

  // done
  gmvwrite_closefile();
}


void State::init_restart( )
{
  int rank = mesh_maps->get_comm()->MyPID();  

  unsigned int num_nodes = mesh_maps->num_entities(Amanzi::AmanziMesh::NODE,Amanzi::AmanziMesh::OWNED);
  unsigned int num_cells = mesh_maps->num_entities(Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::OWNED);  
  unsigned int num_faces = mesh_maps->num_entities(Amanzi::AmanziMesh::FACE,Amanzi::AmanziMesh::OWNED); 

  int nums[3];
  int dummy[3];
  dummy[0] = num_nodes;
  dummy[1] = num_cells;  
  dummy[2] = num_faces;

  mesh_maps->get_comm()->SumAll(dummy,nums,3);  


  // make the all to one map
  if (rank == 0) {
    int *gids = new int[nums[1]];
    for (int i=0; i<nums[1]; i++) gids[i] = i;
    all_to_one_cell_map = Teuchos::rcp(new Epetra_Map(nums[1],nums[1],gids,0, * mesh_maps->get_comm() ));
    delete [] gids;
    
    gids = new int[nums[0]];
    for (int i=0; i<nums[0]; i++) gids[i] = i;    
    all_to_one_node_map = Teuchos::rcp(new Epetra_Map(nums[0],nums[0],gids,0, * mesh_maps->get_comm() ));
    delete [] gids;

    // gids = new int[nums[2]];
    // for (int i=0; i<nums[2]; i++) gids[i] = i;
    // all_to_one_face_map = Teuchos::rcp(new Epetra_Map(nums[2],nums[2],gids,0, * mesh_maps->get_comm() ));
    // delete [] gids;  

    int max_gid = mesh_maps->face_map(false).MaxAllGID();
    int min_gid = mesh_maps->face_map(false).MinAllGID();
    gids = new int [max_gid-min_gid+1];
    for (int i=0; i<max_gid-min_gid+1; i++) gids[i] = min_gid+i;
    all_to_one_face_map = Teuchos::rcp(new Epetra_Map(max_gid+1,max_gid+1,gids,0, * mesh_maps->get_comm() ));
    

  } else {
    int *gids;
    int max_gid = mesh_maps->face_map(false).MaxAllGID();
    all_to_one_cell_map = Teuchos::rcp(new Epetra_Map(nums[1],0,gids,0, * mesh_maps->get_comm() ) );
    all_to_one_node_map = Teuchos::rcp(new Epetra_Map(nums[0],0,gids,0, * mesh_maps->get_comm() ) );
    all_to_one_face_map = Teuchos::rcp(new Epetra_Map(max_gid+1,0,gids,0, * mesh_maps->get_comm() ) );
  }

  // make the all to one exporters 
  all_to_one_cell_export = Teuchos::rcp(new Epetra_Export(mesh_maps->cell_map(false), *all_to_one_cell_map) );
  all_to_one_node_export = Teuchos::rcp(new Epetra_Export(mesh_maps->node_map(false), *all_to_one_node_map) );  
  all_to_one_face_export = Teuchos::rcp(new Epetra_Export(mesh_maps->face_map(false), *all_to_one_face_map) );

}


void State::write_restart ( std::string filename )
{
  int rank = mesh_maps->get_comm()->MyPID();

  Epetra_Vector PE0(*all_to_one_cell_map);
  Epetra_Vector PEF(*all_to_one_face_map);
  Amanzi::HDF5 restart_output;

  if (rank == 0) {
    restart_output.setTrackXdmf(false);
    restart_output.createDataFile(filename);
  }
       
  PE0.Export( *water_density, *all_to_one_cell_export, Insert);
  if (rank == 0) {
    restart_output.writeCellData(PE0, "water_density");
  }

  PE0.Export( *pressure, *all_to_one_cell_export, Insert);   
  if (rank == 0) {
    restart_output.writeCellData(PE0, "pressure");
  }

  //darcy_flux->Print(std::cout);
  PEF.Export( *darcy_flux, *all_to_one_face_export, Insert);
  //PEF.Print(std::cout);
  if (rank == 0) {
    restart_output.writeCellData(PEF, "darcy_flux");
  }

  PE0.Export( *porosity, *all_to_one_cell_export, Insert);
  if (rank == 0) {
    restart_output.writeCellData(PE0, "porosity");
  }  

  PE0.Export( *water_saturation, *all_to_one_cell_export, Insert);
  if (rank == 0) {
    restart_output.writeCellData(PE0, "water_saturation");
  }   
  
  PE0.Export( *permeability, *all_to_one_cell_export, Insert);
  if (rank == 0) {
    restart_output.writeCellData(PE0, "permeability");
  }   

  
  for (int i=0; i<darcy_velocity->NumVectors(); i++)
    {
      std::stringstream fnss;
      fnss << "darcy_velocity-" << i;
      
      PE0.Export( *(*darcy_velocity)(0), *all_to_one_cell_export, Insert);
      if (rank == 0) {
	restart_output.writeCellData(PE0, fnss.str() );
      }   
    }
  
  // for (int i=0; i<total_component_concentration->NumVectors(); i++)
  //   {
  //     std::stringstream fnss;
  //     fnss << "total_component_concentration-" << i;
      
  //     PE0.Export( *(*total_component_concentration)(0), *all_to_one_cell_export, Insert);
  //     if (rank == 0) {
  // 	restart_output.writeCellData(PE0, fnss.str() );
  //     }   
  //   }
}

void State::read_restart ( std::string filename )
{
  int rank = mesh_maps->get_comm()->MyPID();
  Epetra_Vector PE0(*all_to_one_cell_map); 
  Epetra_Vector PEF(*all_to_one_face_map); 
  
  Amanzi::HDF5 *restart_output = new Amanzi::HDF5();

  if (rank == 0) {
    restart_output->setTrackXdmf(false);
    restart_output->setH5DataFilename(filename);
  }

  if (rank == 0) {
    restart_output->readData(PE0, "water_density");
  }
  water_density->Import( PE0, *all_to_one_cell_export, Insert);

  if (rank == 0) {
    restart_output->readData(PE0, "pressure");
  }
  pressure->Import( PE0, *all_to_one_cell_export, Insert);   

  if (rank == 0) {
    restart_output->readData(PEF, "darcy_flux");
  }
  darcy_flux->Import( PEF, *all_to_one_face_export, Insert);

  if (rank == 0) {
    restart_output->readData(PE0, "porosity");
  }  
  porosity->Import( PE0, *all_to_one_cell_export, Insert);

  if (rank == 0) {
    restart_output->readData(PE0, "water_saturation");
  }   
  water_saturation->Import( PE0, *all_to_one_cell_export, Insert);
  
  if (rank == 0) {
    restart_output->readData(PE0, "permeability");
  }   
  permeability->Import( PE0, *all_to_one_cell_export, Insert);

  
  for (int i=0; i<darcy_velocity->NumVectors(); i++)
    {
      std::stringstream fnss;
      fnss << "darcy_velocity-" << i;
      
      if (rank == 0) {
	restart_output->readData(PE0, fnss.str() );
      }   
      (*darcy_velocity)(i)->Import( PE0, *all_to_one_cell_export, Insert);
    }
  
  // for (int i=0; i<total_component_concentration->NumVectors(); i++)
  //   {
  //     std::stringstream fnss;
  //     fnss << "total_component_concentration-" << i;
      
  //     if (rank == 0) {
  // 	restart_output->readData(PE0, fnss.str() );
  //     }   
  //     PE0.Import( *(*total_component_concentration)(0), *all_to_one_cell_export, Insert);
  //   }

}


double State::water_mass()
{
  return 0.0;
}
