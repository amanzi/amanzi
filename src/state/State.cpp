#include "Epetra_Vector.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Import.h"

#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_LAPACK.hpp"

#include "errors.hh"
#include "exceptions.hh"

#include "Mesh.hh"
#include "Point.hh"
#include "linear-function.hh"

#include "State.hpp"


/* *******************************************************************/
State::State(int number_of_components_,
             int number_of_minerals,
             Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh_maps_)
    : number_of_components(number_of_components_),
      mesh_maps(mesh_maps_),
      number_of_minerals_(number_of_minerals),
      number_of_ion_exchange_sites_(0),
      number_of_sorption_sites_(0),
      using_sorption_(false),
      use_sorption_isotherms_(false) {
  // the default parameter_list is empty, so we can't sanely implement
  // mineralogy or isotherms here because we can't safely allocate
  // memory.... This constructor can't be used for anything with
  // geochemistry...?!
  init_verbosity(parameter_list);

  SetupSoluteNames();
  // can't call verify mineralogy because no way of setting the 
  //VerifyMaterialChemistry();
  // create the Eptera_Vector objects
  create_storage();
  ExtractVolumeFromMesh();
};


State::State(Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh_maps_)
    : mesh_maps(mesh_maps_),
      number_of_ion_exchange_sites_(0),
      number_of_sorption_sites_(0),
      using_sorption_(false),
      use_sorption_isotherms_(false) {
  // this constructor is going to be used in restarts, where we
  // read the number of components from a file before creating
  // storage
  init_verbosity(parameter_list);
}


State::State(Teuchos::ParameterList &parameter_list_,
             Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh_maps_)
    : mesh_maps(mesh_maps_),
      parameter_list(parameter_list_),
      number_of_minerals_(0),
      number_of_ion_exchange_sites_(0),
      number_of_sorption_sites_(0),
      using_sorption_(false),
      use_sorption_isotherms_(false) {
  init_verbosity(parameter_list);

  // need to set up a few things before we can setup the storage and assign values
  SetupSoluteNames();
  VerifyMaterialChemistry();

  // create the Eptera_Vector objects
  create_storage();
  initialize_from_parameter_list();
  ExtractVolumeFromMesh();
};


/* *******************************************************************/
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


/* *******************************************************************/
State::~State()
{
  //delete [] (*gravity);
}


/* *******************************************************************/
void State::create_default_compnames(int n)
{
  compnames.resize(n);
  for (int i=0; i<n; i++)
  {
    std::stringstream ss;
    ss << "Component " << i;
    compnames[i] = ss.str();
  }
}

void State::SetupSoluteNames(void) {
  // get the number of component concentrations from the
  // parameter list
  if (parameter_list.isParameter("Number of component concentrations")) {
    number_of_components = 
        parameter_list.get<int>("Number of component concentrations");
  } else {
    // if the parameter list does not contain this key, then assume we
    // are being called from a state constructor w/o a valid parameter
    // list (vis/restart) and the number_of_components variable was
    // already set to a valid (non-zero?) value....
  }

  // read the component names if they are spelled out
  Teuchos::Array<std::string> comp_names;
  if (parameter_list.isParameter("Component Solutes"))
  {
    comp_names = 
        parameter_list.get<Teuchos::Array<std::string> >("Component Solutes");
  }

  if (comp_names.size()) {
    set_compnames(comp_names);
  } else {
    create_default_compnames(number_of_components);
  }
  // now create the map
  for (int i = 0; i < comp_names.size(); ++i) {
    comp_no[comp_names[i]] = i;
  }
}  // end SetupSoluteNames()

void State::initialize_from_parameter_list()
{
  using Teuchos::OSTab;
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  OSTab tab = this->getOSTab(); // This sets the line prefix and adds one tab  



  double u[3];
  u[0] = parameter_list.get<double>("Gravity x");
  u[1] = parameter_list.get<double>("Gravity y");
  if (mesh_maps->space_dimension() == 3)
    u[2] = parameter_list.get<double>("Gravity z");
  else
    u[2] = 0.0;
  set_gravity(u);

  set_zero_total_component_concentration();
  set_water_density(parameter_list.get<double>("Constant water density"));
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
      
      if (parameter_list.isParameter("Material Names")) {
        Teuchos::Array<std::string> matnames =  parameter_list.get<Teuchos::Array<std::string> >("Material Names");
        
        *out << std::endl << "Material name ---> Material ID:" << std::endl;
        for (int k=0; k<matnames.size(); k++) {
          *out << matnames[k] << " ---> " << k+1 << std::endl;
        }
        *out << std::endl;
      }

      for (int ii=0; ii<regnames.size(); ii++) {
        double value = static_cast<double>(matids[ii]);
        set_cell_value_in_region(value, *get_material_ids(), regnames[ii]);
      }
    }  
  } else {
    material_ids->PutScalar(0.0);
  }

  for (Teuchos::ParameterList::ConstIterator it = parameter_list.begin(); it != parameter_list.end(); it++) {
    if (parameter_list.isSublist(it->first) && (it->first != "VerboseObject")) {
      Teuchos::ParameterList& sublist = parameter_list.sublist(it->first);

      std::string region = sublist.get<std::string>("Region");

      // initialize the arrays with some constants from the input file
      set_porosity(sublist.get<double>("Constant porosity"), region);

      set_cell_value_in_region(sublist.get<double>("Constant particle density",1.0), *particle_density, region);

      if (sublist.isParameter("Constant permeability")) {
        set_permeability(sublist.get<double>("Constant permeability"), region);
      } else {
        set_vertical_permeability(sublist.get<double>("Constant vertical permeability"), region);
        set_horizontal_permeability(sublist.get<double>("Constant horizontal permeability"), region);
      }

      u[0] = sublist.get<double>("Constant velocity x", 0.0);
      u[1] = sublist.get<double>("Constant velocity y", 0.0);
      u[2] = sublist.get<double>("Constant velocity z", 0.0);
      set_darcy_flux(u, region);

      // set the pressure
      if (sublist.isSublist("uniform pressure")) {
        const Teuchos::ParameterList& unif_p_list = sublist.sublist("uniform pressure");
        set_uniform_pressure( unif_p_list, region );
      } else if (sublist.isSublist("linear pressure")) {
        const Teuchos::ParameterList& lin_p_list = sublist.sublist("linear pressure");
        set_linear_pressure(lin_p_list, region);
      } else if (sublist.isSublist("file pressure")) {
	const Teuchos::ParameterList& file_p_list = sublist.sublist("file pressure");
	set_file_pressure(file_p_list, region);
      }

      // set the saturation
      if (sublist.isSublist("uniform saturation")) {
        const Teuchos::ParameterList& unif_s_list = sublist.sublist("uniform saturation");
        set_uniform_saturation(unif_s_list, region);
      } else if (sublist.isSublist("linear saturation")) {
        const Teuchos::ParameterList& lin_s_list = sublist.sublist("linear saturation");
        set_uniform_saturation(lin_s_list, region);
      }  

      // set specific storage Ss or specific yield Sy in the same state variable
      if (sublist.isParameter("Constant specific storage")) {
        set_specific_storage(sublist.get<double>("Constant specific storage"), region);
      } else if (sublist.isParameter("Constant specific yield")) {
        set_specific_storage(sublist.get<double>("Constant specific yield"), region);
      }

      // read the component concentrations from the xml file
      // and initialize them in mesh block mesh_block_ID
      double tcc_const[number_of_components];
      double free_ion_guess[number_of_components];
      for (int nc=0; nc<number_of_components; nc++) {
        std::stringstream s;
        s << "Constant component concentration " << nc;

        tcc_const[nc] = sublist.get<double>(s.str());
        s.clear();
        s.str("");
        s << "Free Ion Guess " << nc;
        free_ion_guess[nc] = sublist.get<double>(s.str());
        s.clear();
        s.str("");
      }
      set_total_component_concentration(tcc_const, region);
      set_free_ion_concentrations(free_ion_guess, region);
      SetRegionMaterialChemistry(region, &sublist);
    }
  }
}  // end initialize_from_parameter_list()



/* *******************************************************************/
void State::create_storage()
{
  // create the Eptera_Vector objects
  water_density =    Teuchos::rcp( new Epetra_Vector( mesh_maps->cell_map(false) ) );
  pressure =         Teuchos::rcp( new Epetra_Vector( mesh_maps->cell_map(false) ) );
  lambda =           Teuchos::rcp( new Epetra_Vector( mesh_maps->face_map(false) ) );
  darcy_flux =       Teuchos::rcp( new Epetra_Vector( mesh_maps->face_map(false) ) );
  porosity =         Teuchos::rcp( new Epetra_Vector( mesh_maps->cell_map(false) ) );
  water_saturation = Teuchos::rcp( new Epetra_Vector( mesh_maps->cell_map(false) ) );
  prev_water_saturation = Teuchos::rcp( new Epetra_Vector( mesh_maps->cell_map(false) ) );
  vertical_permeability       = Teuchos::rcp( new Epetra_Vector( mesh_maps->cell_map(false) ) );
  horizontal_permeability     = Teuchos::rcp( new Epetra_Vector( mesh_maps->cell_map(false) ) );
  total_component_concentration
      = Teuchos::rcp( new Epetra_MultiVector( mesh_maps->cell_map(false), number_of_components ) );
  darcy_velocity   = Teuchos::rcp( new Epetra_MultiVector( mesh_maps->cell_map(false), 3));
  material_ids =     Teuchos::rcp( new Epetra_Vector( mesh_maps->cell_map(false) ) );
  particle_density = Teuchos::rcp( new Epetra_Vector( mesh_maps->cell_map(false) ) );

  density =   Teuchos::rcp(new double);
  viscosity = Teuchos::rcp(new double);
  gravity =   Teuchos::rcp(new double*);
  *gravity = new double[3];

  specific_storage = Teuchos::rcp(new Epetra_Vector(mesh_maps->cell_map(false)));

  volume_ =     Teuchos::rcp( new Epetra_Vector( mesh_maps->cell_map(false) ) );

  // if chemistry in enabled, we'll always need free_ions stored.
  free_ion_concentrations_ = Teuchos::rcp( new Epetra_MultiVector( mesh_maps->cell_map(false), number_of_components ) );

  // TODO(bandre): activity corrections. Need primaries and
  // secondaries, but don't know number of secondaries when this
  // function is called!

  if (number_of_minerals() > 0) {
    mineral_volume_fractions_ = Teuchos::rcp(
        new Epetra_MultiVector( mesh_maps->cell_map(false), number_of_minerals()));
    mineral_specific_surface_area_ = Teuchos::rcp(
        new Epetra_MultiVector( mesh_maps->cell_map(false), number_of_minerals()));
  } else {
    mineral_volume_fractions_ = Teuchos::null;
    mineral_specific_surface_area_ = Teuchos::null;
  }

  if (using_sorption()) {
    total_sorbed_ = Teuchos::rcp(
        new Epetra_MultiVector(mesh_maps->cell_map(false), number_of_components));
  } else {
    total_sorbed_ = Teuchos::null;
  }

  if (number_of_sorption_sites() > 0) {
    // TODO: this will eventually need to be a 3d array: [cell][mineral][site]
    sorption_sites_ = Teuchos::rcp(
        new Epetra_MultiVector(mesh_maps->cell_map(false), 
                               number_of_sorption_sites()));
  } else {
    sorption_sites_ = Teuchos::null;
  }

  if (number_of_ion_exchange_sites() > 0) {
    // TODO: eventually this probably needs to be a 3d array: [cell][mineral][site]
    // for now we assume [cell][site], but site will always be one?
    ion_exchange_sites_ = Teuchos::rcp(
        new Epetra_MultiVector(mesh_maps->cell_map(false), 
                               number_of_ion_exchange_sites()));    
  } else {
    ion_exchange_sites_ = Teuchos::null;
  }

  if (use_sorption_isotherms()) {
    isotherm_kd_ = Teuchos::rcp(
        new Epetra_MultiVector(mesh_maps->cell_map(false), number_of_components));
    isotherm_freundlich_n_ = Teuchos::rcp(
        new Epetra_MultiVector(mesh_maps->cell_map(false), number_of_components));
    isotherm_langmuir_b_ = Teuchos::rcp(
        new Epetra_MultiVector(mesh_maps->cell_map(false), number_of_components));
  } else {
    isotherm_kd_ = Teuchos::null;
    isotherm_freundlich_n_ = Teuchos::null;
    isotherm_langmuir_b_ = Teuchos::null;
  }

}  // end create_storage()



void State::set_time ( double new_time ) {
  last_time = new_time;
  time = new_time;
}


void State::set_cycle(int new_cycle) {
  cycle = new_cycle;
}


void State::set_number_of_components(const int n)
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


/* *******************************************************************/
void State::update_darcy_flux(const Epetra_Vector &new_darcy_flux)
{
  *darcy_flux = new_darcy_flux;
}


void State::update_pressure(const Epetra_Vector &new_pressure)
{
  *pressure = new_pressure;
}


/* *******************************************************************/
void State::advance_time(double dT)
{
  last_time = time;
  time = time + dT;
}


/* *******************************************************************/
void State::set_cell_value_in_region(const double& value, Epetra_Vector& v,
                                     const std::string& region)
{
  if (!mesh_maps->valid_set_name(region,Amanzi::AmanziMesh::CELL)) {
    throw std::exception();
  }

  unsigned int mesh_block_size = mesh_maps->get_set_size(region,
                                                         Amanzi::AmanziMesh::CELL,
                                                         Amanzi::AmanziMesh::OWNED);

  Amanzi::AmanziMesh::Entity_ID_List cell_ids(mesh_block_size);

  //mesh_maps->get_set(mesh_block_id, Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::OWNED,
  //cell_ids.begin(),cell_ids.end());

  mesh_maps->get_set_entities(region, Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::OWNED,
                              &cell_ids);

  for (Amanzi::AmanziMesh::Entity_ID_List::iterator c = cell_ids.begin(); c != cell_ids.end(); c++) {
    v[*c] = value;
  }
}



/* *******************************************************************/
void State::set_cell_value_in_region(const Epetra_Vector& x, Epetra_Vector& v,
                                     const std::string& region) {

  if (!mesh_maps->valid_set_name(region,Amanzi::AmanziMesh::CELL)) {
    throw std::exception();
  }

  unsigned int mesh_block_size = mesh_maps->get_set_size(region,
                                                         Amanzi::AmanziMesh::CELL,
                                                         Amanzi::AmanziMesh::OWNED);

  Amanzi::AmanziMesh::Entity_ID_List cell_ids(mesh_block_size);

  mesh_maps->get_set_entities(region, Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::OWNED,
                              &cell_ids);

  for(Amanzi::AmanziMesh::Entity_ID_List::iterator c = cell_ids.begin(); c != cell_ids.end(); c++) {
    v[*c] = x[*c];
  }
}



/* *******************************************************************/
void State::set_cell_value_in_region(const Amanzi::Function& fun, Epetra_Vector& v,
                                     const std::string& region)
{
  if (!mesh_maps->valid_set_name(region, Amanzi::AmanziMesh::CELL)) {
    throw std::exception();
  }

  unsigned int mesh_block_size = mesh_maps->get_set_size(region,
                                                         Amanzi::AmanziMesh::CELL,
                                                         Amanzi::AmanziMesh::OWNED);
  Amanzi::AmanziMesh::Entity_ID_List cell_ids(mesh_block_size);
  mesh_maps->get_set_entities(region, Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::OWNED,
                              &cell_ids);

  for( Amanzi::AmanziMesh::Entity_ID_List::iterator c = cell_ids.begin(); c != cell_ids.end();  c++) {
    const Amanzi::AmanziGeometry::Point& p = mesh_maps->cell_centroid(*c);
    v[*c] = fun(&p[0]);
  }
}


/* *******************************************************************/
void State::set_cell_value_in_mesh_block(double value, Epetra_Vector &v,
                                         int mesh_block_id)
{
  if (!mesh_maps->valid_set_id(mesh_block_id,Amanzi::AmanziMesh::CELL)) {
    throw std::exception();
  }

  unsigned int mesh_block_size = mesh_maps->get_set_size(mesh_block_id,
                                                         Amanzi::AmanziMesh::CELL,
                                                         Amanzi::AmanziMesh::OWNED);

  Amanzi::AmanziMesh::Entity_ID_List cell_ids(mesh_block_size);

  mesh_maps->get_set_entities(mesh_block_id, Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::OWNED,
                     &cell_ids);

  for( Amanzi::AmanziMesh::Entity_ID_List::iterator c = cell_ids.begin();
       c != cell_ids.end();  c++) {
    v[*c] = value;
  }
}


/* *******************************************************************/
void State::set_darcy_flux(const double* u, const int mesh_block_id)
{
  // Epetra_Map face_map = mesh_maps->face_map(false);
  if (!mesh_maps->valid_set_id(mesh_block_id,Amanzi::AmanziMesh::CELL)) {
    throw std::exception();
  }

  unsigned int mesh_block_size = mesh_maps->get_set_size(mesh_block_id,
                                                         Amanzi::AmanziMesh::CELL,
                                                         Amanzi::AmanziMesh::OWNED);

  int dim = mesh_maps->space_dimension();
  Amanzi::AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;
  Amanzi::AmanziMesh::Entity_ID_List cell_ids(mesh_block_size);

  mesh_maps->get_set_entities(mesh_block_id, Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::OWNED,
                     &cell_ids);

  for (Amanzi::AmanziMesh::Entity_ID_List::iterator c = cell_ids.begin(); c != cell_ids.end(); c++) {
    mesh_maps->cell_get_faces_and_dirs(*c, &faces, &dirs);
    int nfaces = faces.size();

    for (int i = 0; i < nfaces; i++) {
      int f = faces[i];

      if (mesh_maps->face_map(false).MyLID(f)) {
        const Amanzi::AmanziGeometry::Point& normal = mesh_maps->face_normal(f);
        (*darcy_flux)[f] = 0.0;
        for (int i = 0; i < dim; i++) (*darcy_flux)[f] += u[i] * normal[i];
      }
    }
  }
}


/* *******************************************************************/
void State::set_darcy_flux(const double* u, const std::string region)
{
  if (!mesh_maps->valid_set_name(region,Amanzi::AmanziMesh::CELL)) {
    throw std::exception();
  }

  unsigned int mesh_block_size = mesh_maps->get_set_size(region,
                                                         Amanzi::AmanziMesh::CELL,
                                                         Amanzi::AmanziMesh::OWNED);

  int dim = mesh_maps->space_dimension();
  Amanzi::AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;
  Amanzi::AmanziMesh::Entity_ID_List cell_ids(mesh_block_size);

  mesh_maps->get_set_entities(region, Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::OWNED,
                              &cell_ids);

  for( Amanzi::AmanziMesh::Entity_ID_List::iterator c = cell_ids.begin(); c != cell_ids.end(); c++) {
    mesh_maps->cell_get_faces_and_dirs(*c, &faces, &dirs);
    int nfaces = faces.size();

    for (int i = 0; i< nfaces; i++) {
      int f = faces[i];

      if (mesh_maps->face_map(false).MyLID(f)) {
        const Amanzi::AmanziGeometry::Point& normal = mesh_maps->face_normal(f);
        (*darcy_flux)[f] = 0.0;
        for (int i = 0; i < dim; i++) (*darcy_flux)[f] += u[i] * normal[i];
      }
    }
  }
}


/* *******************************************************************/
void State::set_water_density(const double wd)
{
  water_density->PutScalar(wd);
  *density = wd;
}


void State::set_water_saturation(const double ws)
{
  water_saturation->PutScalar(ws);
}


/* *******************************************************************/
void State::set_porosity(const double phi)
{
  porosity->PutScalar(phi);
}


void State::set_porosity(const double phi, const int mesh_block_id)
{
  set_cell_value_in_mesh_block(phi, *porosity, mesh_block_id);
}


void State::set_porosity( const double phi, const std::string region )
{
  set_cell_value_in_region(phi, *porosity, region);
}


/* *******************************************************************/
void State::set_specific_storage(const double ss, const std::string region)
{
  set_cell_value_in_region(ss, *specific_storage, region);
}


/* *******************************************************************/
void State::set_zero_total_component_concentration()
{
  total_component_concentration->PutScalar(0.0);
}


void State::set_total_component_concentration(const double* conc, const int mesh_block_id)
{
  for (int nc=0; nc<number_of_components; nc++) {
    set_cell_value_in_mesh_block(conc[nc], *(*total_component_concentration)(nc),mesh_block_id);
  }
}

void State::set_total_component_concentration(const double* conc, const std::string region)
{
  for (int nc=0; nc<number_of_components; nc++) {
    set_cell_value_in_region(conc[nc], *(*total_component_concentration)(nc), region);
  }
}


void State::set_free_ion_concentrations( const double* conc, const std::string region )
{
  for (int nc=0; nc<number_of_components; nc++) {
    set_cell_value_in_region(conc[nc], *(*free_ion_concentrations_)(nc),region);
  }
}


/* *******************************************************************/
void State::set_permeability(const double kappa)
{
  vertical_permeability->PutScalar(kappa);
  horizontal_permeability->PutScalar(kappa);
}


void State::set_permeability(const double kappa, const int mesh_block_id)
{
  set_cell_value_in_mesh_block(kappa, *vertical_permeability, mesh_block_id);
  set_cell_value_in_mesh_block(kappa, *horizontal_permeability, mesh_block_id);
}


void State::set_permeability(const double kappa, const std::string region)
{
  set_cell_value_in_region(kappa, *vertical_permeability, region);
  set_cell_value_in_region(kappa, *horizontal_permeability, region);
}


void State::set_vertical_permeability(const double kappa)
{
  vertical_permeability->PutScalar(kappa);
}


void State::set_vertical_permeability(const double kappa, const int mesh_block_id)
{
  set_cell_value_in_mesh_block(kappa, *vertical_permeability, mesh_block_id);
}


void State::set_vertical_permeability(const double kappa, const std::string region)
{
  set_cell_value_in_region(kappa, *vertical_permeability, region);
}


void State::set_horizontal_permeability(const double kappa)
{
  horizontal_permeability->PutScalar(kappa);
}


void State::set_horizontal_permeability(const double kappa, const int mesh_block_id)
{
  set_cell_value_in_mesh_block(kappa,*horizontal_permeability,mesh_block_id);
}


void State::set_horizontal_permeability( const double kappa, const std::string region)
{
  set_cell_value_in_region(kappa,*horizontal_permeability,region);
}


/* *******************************************************************/
void State::set_viscosity(const double mu)
{
  *viscosity = mu;
}


/* *******************************************************************/
void State::set_gravity(const double *g)
{
  (*gravity)[0] = g[0];
  (*gravity)[1] = g[1];
  (*gravity)[2] = g[2];
}


/* ********************************************************************
* Computes the total mass of water in the domain.
**********************************************************************/
double State::water_mass()
{
  Epetra_Vector wm ( *water_saturation );

  wm.Multiply(1.0, *water_density, wm, 0.0);
  wm.Multiply(1.0, *porosity, wm, 0.0);

  Epetra_Vector cell_volume( mesh_maps->cell_map(false) );

  for (int i=0; i<(mesh_maps->cell_map(false)).NumMyElements(); i++)
  {
    cell_volume[i] = mesh_maps->cell_volume(i);
  }

  wm.Multiply(1.0, cell_volume, wm, 0.0);

  double mass;
  wm.Norm1(&mass);

  return mass;
}



/* *******************************************************************/
double State::point_value(const std::string& point_region, const std::string& name)
{
  // if (!mesh_maps->valid_set_name(point_region, Amanzi::AmanziMesh::CELL)) {
  //   throw
  // }

  unsigned int mesh_block_size = mesh_maps->get_set_size(point_region,
                                                         Amanzi::AmanziMesh::CELL,
                                                         Amanzi::AmanziMesh::OWNED);

  double value(0.0);
  double volume(0.0);

  Amanzi::AmanziMesh::Entity_ID_List cell_ids(mesh_block_size);

  mesh_maps->get_set_entities(point_region, Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::OWNED,
                              &cell_ids);

  // check if Aqueous concentration was requested
  std::string var;

  int pos = name.find("Aqueous concentration");
  if (pos != string::npos) {
    var = name.substr(0, pos-1);
  } else {
    var = name;
  }
  
  // extract the value if it is a component
  if (comp_no.find(var) != comp_no.end())  {
    value = 0.0;
    volume = 0.0;
    for (int i=0; i<mesh_block_size; i++) {
      int ic = cell_ids[i];
      value += (*(*total_component_concentration)(comp_no[var]))[ic] * mesh_maps->cell_volume(ic);
      
      volume += mesh_maps->cell_volume(ic);
    }
  } else if (var == "Volumetric water content") {
    value = 0.0;
    volume = 0.0;
    
    for (int i=0; i<mesh_block_size; i++) {
      int ic = cell_ids[i];
      value += (*porosity)[ic] * (*water_saturation)[ic] * mesh_maps->cell_volume(ic);
      volume += mesh_maps->cell_volume(ic);
    }
  } else if (var == "Gravimetric water content") {
    value = 0.0;
    volume = 0.0;
    
    for (int i=0; i<mesh_block_size; i++) {
      int ic = cell_ids[i];
      value += (*porosity)[ic] * (*water_saturation)[ic] * (*water_density)[ic] 
	/ ( (*particle_density)[ic] * (1.0 - (*porosity)[ic] ) )  * mesh_maps->cell_volume(ic);
      volume += mesh_maps->cell_volume(ic);
    }    
  } else if (var == "Aqueous pressure") {
    value = 0.0;
    volume = 0.0;
    
    for (int i=0; i<mesh_block_size; i++) {
      int ic = cell_ids[i];
      value += (*pressure)[ic] * mesh_maps->cell_volume(ic);
      volume += mesh_maps->cell_volume(ic);
    }
  } else if (var == "Aqueous saturation") {
    value = 0.0;
    volume = 0.0;
    
    for (int i=0; i<mesh_block_size; i++) {
      int ic = cell_ids[i];
      value += (*water_saturation)[ic] * mesh_maps->cell_volume(ic);
      volume += mesh_maps->cell_volume(ic);
    }    
  // } else if (var == "Hydrostatic Head") {
  //   value = 0.0;
  //   volume = 0.0;

  //   for (int i=0; i<mesh_block_size; ++i) {
  //     int ic = cell_ids[i];
  //     value += (*pressure)[ic]/ ( (*density) * (*gravity)[2]);
  //     value *= mesh_maps->cell_volume(ic);
  //     volume += mesh_maps->cell_volume(ic);
  //   }
  } else {
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


/* *******************************************************************/
void State::set_darcy_flux(const Epetra_Vector& darcy_flux_)
{
  *darcy_flux = darcy_flux_;
}


void State::set_water_saturation(const Epetra_Vector& water_saturation_)
{
  *water_saturation = water_saturation_;
};


void State::set_prev_water_saturation(const Epetra_Vector& prev_water_saturation_)
{
  *prev_water_saturation = prev_water_saturation_;
};


void State::set_porosity(const Epetra_Vector& porosity_)
{
  *porosity = porosity_;
};

void State::set_particle_density(const Epetra_Vector& particle_density_)
{
  *particle_density = particle_density_;
};



/* *******************************************************************/
void State::set_permeability(const Epetra_Vector& permeability_)
{
  *vertical_permeability = permeability_;
  *horizontal_permeability = permeability_;
};


void State::set_vertical_permeability(const Epetra_Vector& permeability_)
{
  *vertical_permeability = permeability_;
};


void State::set_horizontal_permeability( const Epetra_Vector& permeability_)
{
  *horizontal_permeability = permeability_;
};


/* *******************************************************************/
void State::set_pressure(const Epetra_Vector& pressure_)
{
  *pressure = pressure_;
};


/* *******************************************************************/
void State::set_lambda(const Epetra_Vector& lambda_)
{
  *lambda = lambda_;
};


/* *******************************************************************/
void State::set_water_density(const Epetra_Vector& water_density_)
{
  *water_density = water_density_;
};


void State::set_darcy_velocity(const Epetra_MultiVector& darcy_velocity_)
{
  *darcy_velocity = darcy_velocity_;
};


/* *******************************************************************/
void State::set_total_component_concentration(const Epetra_MultiVector& total_component_concentration_)
{
  *total_component_concentration = total_component_concentration_;
};


/* *******************************************************************/
void State::set_material_ids(const Epetra_Vector& material_ids_)
{
  *material_ids = material_ids_;
}


/* *******************************************************************/
void State::set_uniform_pressure(const Teuchos::ParameterList& unif_p_list, const std::string& region)
{
  // get value from paramter list
  const double value = unif_p_list.get<double>("value");
  set_cell_value_in_region(value, *pressure, region);
}


void State::set_linear_pressure(const Teuchos::ParameterList& lin_p_list, const std::string& region)
{
  // get parameters from parameter list
  const double ref_value = lin_p_list.get<double>("reference value");
  const Teuchos::Array<double>& ref_coord = lin_p_list.get<Teuchos::Array<double> >("reference coordinate");
  const Teuchos::Array<double>& gradient = lin_p_list.get<Teuchos::Array<double> >("gradient");

  // create function
  Amanzi::LinearFunction lin_p(ref_value, gradient.toVector(), ref_coord.toVector());

  set_cell_value_in_region(lin_p, *pressure, region);
}


void State::set_file_pressure ( const Teuchos::ParameterList& file_p_list, const std::string& region )
{
  // get parameters from parameter list
  const std::string filename = file_p_list.get<std::string>("file name");
  const std::string label = file_p_list.get<std::string>("label");

  // read the pressure variable from the checkpoint file
  Amanzi::HDF5_MPI *checkpoint_input = new Amanzi::HDF5_MPI(*mesh_maps->get_comm(), filename);   
  
  Epetra_Vector cell_vector(mesh_maps->cell_epetra_map(false));
  checkpoint_input->readData(cell_vector,label);

  set_cell_value_in_region(cell_vector, *pressure, region);
};


/* *******************************************************************/
void State::set_uniform_saturation(const Teuchos::ParameterList& unif_s_list, const std::string& region)
{
  // get value from paramter list
  const double value = unif_s_list.get<double>("value");
  set_cell_value_in_region(value, *water_saturation, region);
  set_cell_value_in_region(value, *prev_water_saturation, region);
};


void State::set_linear_saturation(const Teuchos::ParameterList& lin_s_list, const std::string& region)
{
  // get parameters from parameter list
  const double ref_value = lin_s_list.get<double>("reference value");
  const Teuchos::Array<double>& ref_coord = lin_s_list.get<Teuchos::Array<double> >("reference coordinate");
  const Teuchos::Array<double>& gradient = lin_s_list.get<Teuchos::Array<double> >("gradient");

  // create function
  Amanzi::LinearFunction lin_p(ref_value, gradient.toVector(), ref_coord.toVector());

  set_cell_value_in_region(lin_p, *water_saturation, region);
  set_cell_value_in_region(lin_p, *prev_water_saturation, region);
  
}


// return component number, -1 if the component does not exist
int State::get_component_number(const std::string component_name) {
  std::map<std::string, int>::const_iterator it = comp_no.find(component_name);
  if (it != comp_no.end()) {
    return it->second;
  } else {
    return -1;
  }
}


// return component name, empty string if number does not exist
std::string State::get_component_name(const int component_number) {
  return compnames[component_number];

  // if ( component_number < 0  || component_number >= compnames.size() ) { 
  //   return compnames[component_number];
  // } else {
  //   return std::string("");
  // }
}



/* *******************************************************************/
void State::write_vis(Amanzi::Vis& vis, bool chemistry_enabled, bool force) {
  if (!vis.is_disabled()) {
    if ( vis.dump_requested(get_cycle(), get_time()) || force )  {
      using Teuchos::OSTab;
      Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
      Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
      OSTab tab = this->getOSTab(); // This sets the line prefix and adds one tab
      
      if (out.get() && includesVerbLevel(verbLevel, Teuchos::VERB_LOW, true)) {
        *out << "Writing visualization dump, cycle = " << get_cycle() << std::endl;
      }
      
      // create the new time step...
      vis.create_timestep(get_time()/(365.25*24*60*60),get_cycle());
      
      // dump all the state vectors into the file
      vis.write_vector(*get_pressure(), "pressure");
      vis.write_vector(*get_porosity(),"porosity");
      vis.write_vector(*get_water_saturation(),"water saturation");
      vis.write_vector(*get_water_density(),"water density");
      vis.write_vector(*get_vertical_permeability(),"vertical permeability");
      vis.write_vector(*get_horizontal_permeability(),"horizontal permeability");
      vis.write_vector(*get_material_ids(),"material IDs");
      
      // compute volumetric water content for visualization (porosity*water_saturation)
      Epetra_Vector vol_water( mesh_maps->cell_map(false) );
      vol_water.Multiply(1.0, *water_saturation, *porosity, 0.0);
      vis.write_vector(vol_water,"volumetric water content");
      
      // compute gravimetric water content for visualization
      // MUST have computed volumetric water content before
      vol_water.Multiply(1.0, *water_density, vol_water, 0.0);
      Epetra_Vector bulk_density( mesh_maps->cell_map(false) );
      bulk_density.PutScalar(1.0);
      bulk_density.Update(-1.0,*porosity,1.0);
      bulk_density.Multiply(1.0,*particle_density,bulk_density,0.0);
      vol_water.ReciprocalMultiply(1.0,bulk_density,vol_water,0.0);
      vis.write_vector(vol_water,"gravimetric water content");

      // compute hydrostatic head
      // TO DO

      std::vector<std::string> names(3);
      names[0] = "darcy velocity x";
      names[1] = "darcy velocity y";
      names[2] = "darcy velocity z";
      vis.write_vector(*get_darcy_velocity(), names);
      
      // write component data
      vis.write_vector( *get_total_component_concentration(), compnames);

      // write the geochemistry data
      if (chemistry_enabled) WriteChemistryToVis(&vis);

      vis.finalize_timestep();
    }
  }
}


/* *******************************************************************/
void State::write_vis(Amanzi::Vis& vis, 
                      Teuchos::RCP<Epetra_MultiVector> auxdata, 
                      const std::vector<std::string>& auxnames, 
                      bool chemistry_enabled, bool force)  {
  write_vis(vis, chemistry_enabled, force);
  
  if ( !vis.is_disabled() ) {
    if ( force || vis.dump_requested(get_cycle())) {
      // write auxillary data
      if (!is_null(auxdata))  {
        vis.write_vector( *auxdata , auxnames);
      }
      
      vis.finalize_timestep();
    }
  }
}


/* *******************************************************************/
void State::set_compnames(std::vector<std::string>& compnames_)
{
  compnames = compnames_;
}

void State::set_compnames(Teuchos::Array<std::string>& compnames_)
{
  compnames.clear();
  compnames.resize(compnames_.size());
  for (int i = 0; i < compnames.size(); ++i) {
    compnames.at(i) = compnames_.at(i);
  }
}

void State::ExtractVolumeFromMesh(void) {
  int ncell = mesh_maps->cell_map(false).NumMyElements();
  if (ncell != volume()->MyLength()) {
    Exceptions::amanzi_throw(Errors::Message("State::ExtractVolumeFromMesh() size error."));
  }
  for (int j = 0; j < ncell; ++j) {
    (*volume_)[j] = mesh_maps->cell_volume(j);
  }
}


void State::VerifyMaterialChemistry(void) {
  /*
  ** This function is setting the mineral_names_, mineral_name_id_map_,
  ** number of ion exchange sites, number of sorption sites, etc. It
  ** must be called *before* create_storage()!
  **
  ** loop through each mesh region/mesh block whatever and verify that
  ** the mineralogy and isotherm data, if provided, is correct.
  **
  */
  const std::string block_key("Mesh block ");

  SetupMineralNames();
  SetupSorptionSiteNames();

  // loop through the state parameter list looking for mesh blocks
  for (Teuchos::ParameterList::ConstIterator item = parameter_list.begin();
       item != parameter_list.end(); ++item) {

    std::string item_name = parameter_list.name(item);
    size_t found_block = item_name.find(block_key);
    if (found_block != std::string::npos) {
      std::string block_name = item_name.substr(found_block + block_key.length(), 
                                                item_name.length());
      Teuchos::ParameterList mesh_block_data = parameter_list.sublist(parameter_list.name(item));
      // block_name should match mesh_block_data.region...
      std::string region_name = mesh_block_data.get<std::string>("Region");

      //
      // check for chemistry related data in the block:
      //

      if (mesh_block_data.isSublist("Mineralogy")) {
        VerifyMineralogy(region_name, mesh_block_data.sublist("Mineralogy"));
      }

      if (mesh_block_data.isSublist("Sorption Isotherms")) {
        VerifySorptionIsotherms(region_name,
                                mesh_block_data.sublist("Sorption Isotherms"));
      }

      if (mesh_block_data.isSublist("Surface Complexation Sites")) {
        VerifySorptionSites(region_name,
                            mesh_block_data.sublist("Surface Complexation Sites"));
      }
      
      if (mesh_block_data.isParameter("Cation Exchange Capacity")) {
        // limit to one ion exchange site for now....
        set_using_sorption(true);
        set_number_of_ion_exchange_sites(1);
      }
    }  // end if(mesh_block)
  }  // end for(parameter_list)


  //
  // error checking?
  //

}  // end VerifyMaterialChemistry()

void State::SetupMineralNames() {
  // do we need to worry about minerals?
  mineral_names_.clear();
  Teuchos::Array<std::string> data;
  if (parameter_list.isParameter("Minerals")) {
    data = parameter_list.get<Teuchos::Array<std::string> >("Minerals");
  } 

  // the mineral_names_ list should be the order expected by the chemistry....
  mineral_names_.clear();
  mineral_name_id_map_.clear();
  for (int m = 0; m < data.size(); ++m) {
    mineral_name_id_map_[data.at(m)] = m;
    mineral_names_.push_back(data.at(m));
  }

  if (mineral_names_.size() > 0) {
    // we read some mineral names, so override any value that may have
    // been set in the constructor
    set_number_of_minerals(mineral_names_.size());
  } else if (number_of_minerals() > 0 && mineral_names_.size() == 0) {
    // assume we are called from the constructor w/o a valid parameter
    // list and the mineral names will be set later....
  }
}  // end SetupMineralNames()

void State::SetupSorptionSiteNames() {
  // could almost generalize the SetupMineralNames and
  // SetupSorptionSiteNames into a single function w/ different
  // parameters, but using_sorption needs to be set...

  // do we need to worry about sorption sites?
  sorption_site_names_.clear();
  Teuchos::Array<std::string> data;
  if (parameter_list.isParameter("Sorption Sites")) {
    data = parameter_list.get<Teuchos::Array<std::string> >("Sorption Sites");
  } 

  // the sorption_site_names_ list should be the order expected by the chemistry...
  sorption_site_names_.clear();
  sorption_site_name_id_map_.clear();
  for (int s = 0; s < data.size(); ++s) {
    sorption_site_name_id_map_[data.at(s)] = s;
    sorption_site_names_.push_back(data.at(s));
  }

  if (sorption_site_names_.size() > 0) {
    // we read some sorption site names, so override any value that
    // may have been set in the constructor and set the sorption flag
    // so we allocate the correct amount of memory
    set_number_of_sorption_sites(sorption_site_names_.size());
    set_using_sorption(true);
  } else if (number_of_sorption_sites() > 0 && sorption_site_names_.size() == 0) {
    // assume we are called from the constructor w/o a valid parameter
    // list and the sorption names will be set later....?
  }
}  // end SetupSorptionSiteNames()

void State::VerifyMineralogy(const std::string& region_name,
                             const Teuchos::ParameterList& minerals_list) {
  // loop through each mineral, verify that the mineral name is known
  for (Teuchos::ParameterList::ConstIterator mineral_iter = minerals_list.begin(); 
       mineral_iter != minerals_list.end(); ++mineral_iter) {
    std::string mineral_name = minerals_list.name(mineral_iter);
    if (!mineral_name_id_map_.count(mineral_name)) {
      std::stringstream message;
      message << "Error: State::VerifyMineralogy(): " << mineral_name
              << " was specified in the mineralogy for region "
              << region_name << " but was not listed in the minerals phase list.\n";
      Exceptions::amanzi_throw(Errors::Message(message.str()));            
    }

    // all minerals will have a volume fraction and specific surface
    // area, but sane defaults can be provided, so we don't bother
    // with them here.

  }  // end for(minerals)
}  // end VerifyMineralogy()

void State::VerifySorptionIsotherms(const std::string& region_name,
                                    const Teuchos::ParameterList& isotherms_list) {
  // verify that every species listed is in the component names list.
  // verify that every listed species has a Kd value (no sane default)
  // langmuir and freundlich values are optional (sane defaults)
  set_using_sorption(true);
  set_use_sorption_isotherms(true);

  // loop through each species in the isotherm list
  for (Teuchos::ParameterList::ConstIterator species_iter = isotherms_list.begin(); 
       species_iter != isotherms_list.end(); ++species_iter) {
    std::string species_name = isotherms_list.name(species_iter);

    // verify that the name is a known species
    if (!comp_no.count(species_name)) {
      std::stringstream message;
      message << "Error: State::VerifySorptionIsotherms(): region: "
              << region_name << " contains isotherm data for solute \'" 
              << species_name 
              << "\' but it is not specified in the component solutes list.\n";
      Exceptions::amanzi_throw(Errors::Message(message.str()));      
    }


    // check that this item is a sublist:
    Teuchos::ParameterList species_data;
    if (!isotherms_list.isSublist(species_name)) {
      std::stringstream message;
      message << "Error: State::VerifySorptionIsotherms(): region: "
              << region_name << " ; species : " << species_name
              << " ; must be a named \'ParameterList\' of isotherm data.\n";
      Exceptions::amanzi_throw(Errors::Message(message.str()));            
    } else {
      species_data = isotherms_list.sublist(species_name);
    }

    // verify that the required parameters are present
    if (!species_data.isParameter("Kd")) {
      std::stringstream message;
      message << "Error: State::VerifySorptionIsotherms(): region: "
              << region_name << " ; species name: " << species_name 
              << " ; each isotherm must have a 'Kd' parameter.\n";
      Exceptions::amanzi_throw(Errors::Message(message.str()));      
    }
    // langmuir and freundlich parameters are optional, we'll assign
    // sane defaults.
  }  // end for(species)
}  // end VerifySorptionIsotherms()

void State::VerifySorptionSites(const std::string& region_name,
                                const Teuchos::ParameterList& sorption_site_list) {
  set_using_sorption(true);
  // loop through each sorption site, verify that the site name is known
  for (Teuchos::ParameterList::ConstIterator site_iter = sorption_site_list.begin(); 
       site_iter != sorption_site_list.end(); ++site_iter) {
    std::string site_name = sorption_site_list.name(site_iter);
    if (!sorption_site_name_id_map_.count(site_name)) {
      std::stringstream message;
      message << "Error: State::VerifySorptionSites(): " << site_name
              << " was specified in the 'Surface Complexation Sites' list for region "
              << region_name << " but was not listed in the sorption sites phase list.\n";
      Exceptions::amanzi_throw(Errors::Message(message.str()));            
    }

    // all sorption sites will have a site density
    // but we can default to zero, so don't do any further checking

  }  // end for(sorption_sites)
}  // end VerifySorptionSites()

void State::SetRegionMaterialChemistry(const std::string& region_name,
                                       Teuchos::ParameterList* region_data) {

  // NOTE(bandre): I don't think this is going to ensure all memory is
  // initialized correctly. If a region doesn't have a particular
  // block (e.g. isotherms) in the XML file, it will have
  // uninitialized data!?!?

  if (region_data->isSublist("Mineralogy")) {
    SetRegionMineralogy(region_name, region_data->sublist("Mineralogy"));
  } else {
    // std::cout << "no mineralogy in region '" << region_name 
    //           << "'..." << std::endl;
  }

  if (region_data->isSublist("Sorption Isotherms")) {
    SetRegionSorptionIsotherms(region_name, region_data->sublist("Sorption Isotherms"));
  } else {
    // std::cout << "no sorption isotherms in region '" << region_name 
    //           << "'..." << std::endl;
  }

  if (region_data->isSublist("Surface Complexation Sites")) {
    SetRegionSorptionSites(region_name, region_data->sublist("Surface Complexation Sites"));
  } else {
    // std::cout << "no surface complexation sites in region '" << region_name 
    //           << "'..." << std::endl;
  }

  double cec = region_data->get<double>("Cation Exchange Capacity", 0.0);
  if (number_of_ion_exchange_sites() > 0) {
    set_cell_value_in_region(cec, *(*ion_exchange_sites_)(0), region_name);
  }

}  // end SetRegionMaterialChemistry()

void State::SetRegionMineralogy(const std::string& region_name,
                                const Teuchos::ParameterList& region_mineralogy) {
  /*
  ** Process a "Mineralogy" parameter list for the current region.
  ** 
  ** NOTE: assumes that error checking in the parameter list was done
  ** during VerifyMaterialChemistry()!
  */

  Teuchos::ParameterList::ConstIterator mineral;
  for (mineral = region_mineralogy.begin(); 
       mineral != region_mineralogy.end(); ++mineral) {
    std::string mineral_name = region_mineralogy.name(mineral);

    // find the correct index for this mineral name
    // std::vector<std::string>::iterator mineral_iterator = 
    //     std::find(mineral_names_.begin(), mineral_names_.end(), mineral_name); 
    // int m = std::distance(mineral_names_.begin(), mineral_iterator);
    
    int m = mineral_name_id_map_[mineral_name];

    Teuchos::ParameterList mineral_data = region_mineralogy.sublist(mineral_name);

    double value = mineral_data.get<double>("Volume Fraction", 0.0);
    set_cell_value_in_region(value, *(*mineral_volume_fractions_)(m), region_name);

    value = mineral_data.get<double>("Specific Surface Area", 1.0);
    set_cell_value_in_region(value, *(*mineral_specific_surface_area_)(m), region_name);

  }  // end for(mineral_list)
}  // end SetRegionMineralogy()


void State::SetRegionSorptionIsotherms(const std::string& region_name,
                                       const Teuchos::ParameterList& region_isotherms) {
  /*
  ** Process a "Sorption Isotherm" parameter list for the current region
  **
  ** Note: assume that the species names have already been verified
  */

  Teuchos::ParameterList::ConstIterator isotherm_species;
  for (isotherm_species = region_isotherms.begin(); 
       isotherm_species != region_isotherms.end(); ++isotherm_species) {

    // assume that all sublists are going to be species names with
    // some optional parameters.
    std::string species_name = region_isotherms.name(isotherm_species);

    // determine the species index for this name
    // std::vector<std::string>::iterator species_iterator = 
    //     std::find(compnames.begin(), compnames.end(), species_name); 
    // int s = std::distance(compnames.begin(), species_iterator);

    int s = comp_no[species_name];

    Teuchos::ParameterList species_data = region_isotherms.sublist(species_name);

    // assign per region per species parameter values with sane defaults
    double value = species_data.get<double>("Kd", 0.0);
    set_cell_value_in_region(value, *(*isotherm_kd_)(s), region_name);

    // TODO(bandre): is this a sane default?
    value = species_data.get<double>("Langmuir b", 1.0);
    set_cell_value_in_region(value, *(*isotherm_langmuir_b_)(s), region_name);

    value = species_data.get<double>("Freundlich N", 1.0);
    set_cell_value_in_region(value, *(*isotherm_freundlich_n_)(s), region_name);
  }  // end for(isotherm_species)
}  // end SetRegionSorptionIsotherms()


void State::SetRegionSorptionSites(const std::string& region_name,
                                   const Teuchos::ParameterList& region_sorption_sites) {
  /*
  ** Process a "Surface Complexation Sites" parameter list for the current region
  **
  ** Note: assume that the surface sites names have already been verified
  */

  std::map<std::string, int>::iterator site_data;
  for (site_data = sorption_site_name_id_map_.begin();
       site_data != sorption_site_name_id_map_.end(); ++site_data) {
    std::string site_name = site_data->first;
    int index = site_data->second;
    double site_density = 0.0;
    if (region_sorption_sites.isSublist(site_name)) {
      Teuchos::ParameterList site_list = region_sorption_sites.sublist(site_name);
      site_density = site_list.get<double>("Site Density", 0.0);
    }
    set_cell_value_in_region(site_density, *(*sorption_sites_)(index), region_name);
  }
}  // end SetRegionSorptionSites()


/*******************************************************************************
 **
 **  Chemistry Vis
 **
 ******************************************************************************/
void State::WriteChemistryToVis(Amanzi::Vis* vis) {
  // TODO(bandre): activity corrections....
  WriteFreeIonsToVis(vis);
  WriteTotalSorbedToVis(vis);
  WriteMineralsToVis(vis);
  WriteIsothermsToVis(vis);
  WriteSorptionSitesToVis(vis);
  WriteIonExchangeSitesToVis(vis);
}  // end WriteChemistryToVis()

void State::WriteFreeIonsToVis(Amanzi::Vis* vis) {
  std::string name;
  for (int i = 0; i < get_number_of_components(); ++i) {
    name = compnames.at(i);
    name += "_free";
    vis->write_vector(*(*free_ion_concentrations_)(i), name);
  }
}  // end WriteFreeIonsToVis()

void State::WriteTotalSorbedToVis(Amanzi::Vis* vis) {
  if (using_sorption()) {
    std::string name;
    for (int i = 0; i < get_number_of_components(); ++i) {
      name = compnames.at(i);
      name += "_sorbed";
      vis->write_vector(*(*total_sorbed_)(i), name);
    }
  }
}  // end WriteTotalSorbedToVis()

void State::WriteMineralsToVis(Amanzi::Vis* vis) {
  std::string name;
  for (int m = 0; m < number_of_minerals(); ++m) {
    name = mineral_names().at(m);
    name += "_volume_fraction";
    vis->write_vector(*(*mineral_volume_fractions_)(m), name);
    name = mineral_names().at(m);
    name += "_specific_surface_area";
    vis->write_vector(*(*mineral_specific_surface_area_)(m), name);
  }
}  // end WriteMineralsToVis()

void State::WriteIsothermsToVis(Amanzi::Vis* vis) {
  if (use_sorption_isotherms()) {
    std::string name;
    for (int i = 0; i < get_number_of_components(); ++i) {
      name = compnames.at(i);
      name += "_Kd";
      vis->write_vector(*(*isotherm_kd_)(i), name);
      name = compnames.at(i);
      name += "_freundlich_n";
      vis->write_vector(*(*isotherm_freundlich_n_)(i), name);
      name = compnames.at(i);
      name += "_langmuir_b";
      vis->write_vector(*(*isotherm_langmuir_b_)(i), name);
    }
  }
}  // end WriteIsothermsToVis()

void State::WriteSorptionSitesToVis(Amanzi::Vis* vis) {
  std::string name;
  for (int s = 0; s < number_of_sorption_sites(); ++s) {
    name = sorption_site_names_.at(s);
    vis->write_vector(*(*sorption_sites_)(s), name);
  }
}  // end WriteSorptionSitesToVis()

void State::WriteIonExchangeSitesToVis(Amanzi::Vis* vis) {
  if (number_of_ion_exchange_sites() > 0) {
    // NOTE: Assume that only one ion exchange site!
    std::string name("CEC");
    vis->write_vector(*(*ion_exchange_sites_)(0), name);
  }
}  // end WriteIonExchangeSitesToVis()




void State::DeriveDarcyVelocity() {
  const Epetra_Map& source_fmap = mesh_maps->face_map(false);
  const Epetra_Map& target_fmap = mesh_maps->face_map(true);
  Epetra_Import face_importer(target_fmap, source_fmap);  

#ifdef HAVE_MPI
  Epetra_Vector darcy_flux_wghost(mesh_maps->face_map(true));
  darcy_flux_wghost.Import(*darcy_flux, face_importer, Insert);
#else
  Epetra_Vector& darcy_flux_wghost = *darcy_flux;
#endif

  Teuchos::LAPACK<int, double> lapack;

  int dim = mesh_maps->space_dimension();
  Teuchos::SerialDenseMatrix<int, double> matrix(dim, dim);
  double rhs_cell[dim];

  Amanzi::AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  int ncells_owned = mesh_maps->num_entities(Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::OWNED);
  for (int c = 0; c < ncells_owned; c++) {
    mesh_maps->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    for (int i = 0; i < dim; i++) rhs_cell[i] = 0.0;
    matrix.putScalar(0.0);

    for (int n = 0; n < nfaces; n++) {  // populate least-square matrix
      int f = faces[n];
      const Amanzi::AmanziGeometry::Point& normal = mesh_maps->face_normal(f);
      double area = mesh_maps->face_area(f);

      for (int i = 0; i < dim; i++) {
        rhs_cell[i] += normal[i] * darcy_flux_wghost[f];
        matrix(i, i) += normal[i] * normal[i];
        for (int j = i+1; j < dim; j++) {
          matrix(j, i) = matrix(i, j) += normal[i] * normal[j];
        }
      }
    }

    int info;
    lapack.POSV('U', dim, 1, matrix.values(), dim, rhs_cell, dim, &info);

    for (int i = 0; i < dim; i++) (*darcy_velocity)[i][c] = rhs_cell[i];
  }
}
