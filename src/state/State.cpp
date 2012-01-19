#include "State.hpp"
#include "Epetra_Vector.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Mesh.hh"
#include "cell_geometry.hh"
#include "Point.hh"
#include "Geometry.hh"
#include "linear-function.hh"

#include "dbc.hh"
#include "errors.hh"
#include "exceptions.hh"


State::State(int number_of_components_,
             Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh_maps_):
    number_of_components(number_of_components_), mesh_maps(mesh_maps_)  
{
  init_verbosity(parameter_list);

  // create the Eptera_Vector objects
  create_storage();
  create_default_compnames(number_of_components);
};



State::State( Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh_maps_):
    mesh_maps(mesh_maps_) {
  // this constructor is going to be used in restarts, where we
  // read the number of components from a file before creating
  // storage
  init_verbosity(parameter_list);
};

State::State( Teuchos::ParameterList &parameter_list_,
              Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh_maps_):
    mesh_maps(mesh_maps_),
    parameter_list(parameter_list_)
{
  init_verbosity(parameter_list);
  
  // get the number of component concentrations from the
  // parameter list
  number_of_components = parameter_list.get<int>("Number of component concentrations");

  // create the Eptera_Vector objects

  create_storage();
  initialize_from_parameter_list();
  create_default_compnames(number_of_components);
};


void State::init_verbosity (Teuchos::ParameterList &parameter_list_) {
  // set the line prefix for output
  this->setLinePrefix("Amanzi::State       ");
  // make sure that the line prefix is printed
  this->getOStream()->setShowLinePrefix(true);
  
  // Read the sublist for verbosity settings.
  Teuchos::readVerboseObjectSublist(&parameter_list,this);    

  using Teuchos::OSTab;
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  OSTab tab = this->getOSTab(); // This sets the line prefix and adds one tab
}


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
  using Teuchos::OSTab;
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  OSTab tab = this->getOSTab(); // This sets the line prefix and adds one tab  


  // read the component names if they are spelled out
  Teuchos::Array<std::string> comp_names;
  if (parameter_list.isParameter("Component Solutes"))
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

  double u[3];
  u[0] = parameter_list.get<double>("Gravity x");
  u[1] = parameter_list.get<double>("Gravity y");
  u[2] = parameter_list.get<double>("Gravity z");
  set_gravity(u);

  set_zero_total_component_concentration();
  set_water_density(parameter_list.get<double>("Constant water density"));
  // set_water_saturation(parameter_list.get<double>("Constant water saturation"));
  set_viscosity(parameter_list.get<double>("Constant viscosity"));


  // initialize the material ID array (this is only used for visualization)
  if (parameter_list.isParameter("Region Name to Material ID Map (Material IDs)")) {
    if (parameter_list.isParameter("Region Name to Material ID Map (Region Names)")) {
      // read the region to material id map
      Teuchos::Array<int> matids = parameter_list.get<Teuchos::Array<int> >("Region Name to Material ID Map (Material IDs)");
      Teuchos::Array<std::string> regnames = parameter_list.get<Teuchos::Array<std::string> >("Region Name to Material ID Map (Region Names)");

      // stop if there is a lenght mismatch between the two arrays
      if (matids.size() != regnames.size()) {
        Exceptions::amanzi_throw(Errors::Message("State.cpp: The number of material IDs does not match the number of region names"));
      }
      
      *out << "Region name ---> Material ID:" << std::endl;
      for (int ii=0; ii<regnames.size(); ii++) {
        double value = static_cast<double>(matids[ii]);
        set_cell_value_in_region(value, *get_material_ids(), regnames[ii]);
        
        *out << regnames[ii] << " ---> " << matids[ii] << std::endl;
      }
    }  
  } else {
    material_ids->PutScalar(0.0);
  }




  for (Teuchos::ParameterList::ConstIterator it = parameter_list.begin(); it != parameter_list.end(); it++) {

    if (parameter_list.isSublist(it->first) && (it->first != "VerboseObject")  ) {

      Teuchos::ParameterList& sublist = parameter_list.sublist(it->first);

      std::string region = sublist.get<std::string>("Region");


      // initialize the arrays with some constants from the input file
      set_porosity(sublist.get<double>("Constant porosity"),region);

      if (sublist.isParameter("Constant permeability")) {
        set_permeability(sublist.get<double>("Constant permeability"),region);
      } else {
        set_vertical_permeability(sublist.get<double>("Constant vertical permeability"),region);
        set_horizontal_permeability(sublist.get<double>("Constant horizontal permeability"),region);
      }

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
}




void State::create_storage ()
{
  // create the Eptera_Vector objects

  water_density =    Teuchos::rcp( new Epetra_Vector( mesh_maps->cell_map(false) ) );
  pressure =         Teuchos::rcp( new Epetra_Vector( mesh_maps->cell_map(false) ) );
  darcy_flux =       Teuchos::rcp( new Epetra_Vector( mesh_maps->face_map(false) ) );
  porosity =         Teuchos::rcp( new Epetra_Vector( mesh_maps->cell_map(false) ) );
  water_saturation = Teuchos::rcp( new Epetra_Vector( mesh_maps->cell_map(false) ) );
  vertical_permeability       = Teuchos::rcp( new Epetra_Vector( mesh_maps->cell_map(false) ) );
  horizontal_permeability     = Teuchos::rcp( new Epetra_Vector( mesh_maps->cell_map(false) ) );
  total_component_concentration
      = Teuchos::rcp( new Epetra_MultiVector( mesh_maps->cell_map(false), number_of_components ) );
  darcy_velocity   = Teuchos::rcp( new Epetra_MultiVector( mesh_maps->cell_map(false), 3));
  material_ids =     Teuchos::rcp( new Epetra_Vector( mesh_maps->cell_map(false) ) );

  density =   Teuchos::rcp(new double);
  viscosity = Teuchos::rcp(new double);
  gravity =   Teuchos::rcp(new double*);
  *gravity = new double[3];
}



void State::set_time ( double new_time ) {

  last_time = new_time;
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
  last_time = time;
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
  vertical_permeability->PutScalar(kappa);
  horizontal_permeability->PutScalar(kappa);
}

void State::set_permeability( const double kappa, const int mesh_block_id)
{
  set_cell_value_in_mesh_block(kappa,*vertical_permeability,mesh_block_id);
  set_cell_value_in_mesh_block(kappa,*horizontal_permeability,mesh_block_id);
}

void State::set_permeability( const double kappa, const std::string region)
{
  set_cell_value_in_region(kappa,*vertical_permeability,region);
  set_cell_value_in_region(kappa,*horizontal_permeability,region);
}

void State::set_vertical_permeability( const double kappa )
{
  vertical_permeability->PutScalar(kappa);
}

void State::set_vertical_permeability( const double kappa, const int mesh_block_id)
{
  set_cell_value_in_mesh_block(kappa,*vertical_permeability,mesh_block_id);
}

void State::set_vertical_permeability( const double kappa, const std::string region)
{
  set_cell_value_in_region(kappa,*vertical_permeability,region);
}

void State::set_horizontal_permeability( const double kappa )
{
  horizontal_permeability->PutScalar(kappa);
}

void State::set_horizontal_permeability( const double kappa, const int mesh_block_id)
{
  set_cell_value_in_mesh_block(kappa,*horizontal_permeability,mesh_block_id);
}

void State::set_horizontal_permeability( const double kappa, const std::string region)
{
  set_cell_value_in_region(kappa,*horizontal_permeability,region);
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

  double value(0.0);
  double volume(0.0);

  std::vector<unsigned int> cell_ids(mesh_block_size);

  mesh_maps->get_set_entities(point_region, Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::OWNED,
                              &cell_ids);


  // extract the value if it is a component
  if ( comp_no.find(name) != comp_no.end() )
  {
    value = 0.0;
    volume = 0.0;
    for (int i=0; i<mesh_block_size; i++)
    {
      int ic = cell_ids[i];
      value += (*(*total_component_concentration)( comp_no[name] ))[ic] *  mesh_maps->cell_volume(ic);

      volume += mesh_maps->cell_volume(ic);
    }
  }
  else if ( name == "Water" )
  {
    value = 0.0;
    \
    volume = 0.0;
    for (int i=0; i<mesh_block_size; i++)
    {
      int ic = cell_ids[i];
      value += (*water_density)[ic] * (*porosity)[ic] * (*water_saturation)[ic] * mesh_maps->cell_volume(ic);
      volume += mesh_maps->cell_volume(ic);
    }
  }
  else
  {
    std::stringstream ss;
    ss << "State::point_value: cannot make an observation for variable " << name;
    Errors::Message m(ss.str().c_str());
    Exceptions::amanzi_throw(m);
  }


  // syncronize the result across processors

  double result;
  mesh_maps->get_comm()->SumAll(&value,&result,1);

  double vresult;
  mesh_maps->get_comm()->SumAll(&volume,&vresult,1);

  return result/vresult;
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
  *vertical_permeability = permeability_;
  *horizontal_permeability = permeability_;
};

void State::set_vertical_permeability ( const Epetra_Vector& permeability_ )
{
  *vertical_permeability = permeability_;
};

void State::set_horizontal_permeability ( const Epetra_Vector& permeability_ )
{
  *horizontal_permeability = permeability_;
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

void State::set_material_ids ( const Epetra_Vector& material_ids_ )
{
  *material_ids = material_ids_;
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
    vis.write_vector(*get_vertical_permeability(),"vertical permeability");
    vis.write_vector(*get_horizontal_permeability(),"horizontal permeability");
    vis.write_vector(*get_material_ids(),"material IDs");
    

    std::vector<std::string> names(3);
    names[0] = "darcy velocity x";
    names[1] = "darcy velocity y";
    names[1] = "darcy velocity z";
    vis.write_vector(*get_darcy_velocity(), names);

    // write component data
    vis.write_vector( *get_total_component_concentration(), compnames);
    vis.finalize_timestep();

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
