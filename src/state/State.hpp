#ifndef __State_hpp__
#define __State_hpp__

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Map.h"
#include "Epetra_Export.h"
#include "Mesh.hh"
#include "Vis.hpp"
#include "function.hh"

typedef enum { COMPLETE, UPDATING } status_type;

class State : public Teuchos::VerboseObject<State> {
 public:
  State( int, int, Teuchos::RCP<Amanzi::AmanziMesh::Mesh> );

  State(Teuchos::RCP<Amanzi::AmanziMesh::Mesh>);

  State(Teuchos::ParameterList &parameter_list, 
	Teuchos::RCP<Amanzi::AmanziMesh::Mesh>);

  ~State();

  // access methods
  Teuchos::RCP<Epetra_Vector>      get_pressure ()         { return pressure; }; 
  Teuchos::RCP<Epetra_Vector>      get_lambda ()           { return lambda; }
  Teuchos::RCP<Epetra_Vector>      get_darcy_flux ()       { return darcy_flux; };
  Teuchos::RCP<Epetra_Vector>      get_porosity ()         { return porosity; };
  Teuchos::RCP<Epetra_Vector>      get_water_saturation () { return water_saturation; };
  Teuchos::RCP<Epetra_Vector>      get_prev_water_saturation () { return prev_water_saturation; };
  Teuchos::RCP<Epetra_Vector>      get_water_density ()    { return water_density; };
  Teuchos::RCP<Epetra_Vector>      get_vertical_permeability () { return vertical_permeability; };
  Teuchos::RCP<Epetra_Vector>      get_horizontal_permeability () { return horizontal_permeability; };
  Teuchos::RCP<Epetra_Vector>      get_material_ids ()     { return material_ids; };
  Teuchos::RCP<Epetra_Vector>      get_particle_density () { return particle_density; };

  Teuchos::RCP<double>             get_density()   { return density; } 
  Teuchos::RCP<double>             get_viscosity() { return viscosity; }
  Teuchos::RCP<double*>            get_gravity()   { return gravity; }

  Teuchos::RCP<Epetra_MultiVector> get_darcy_velocity () { return darcy_velocity; }
  Teuchos::RCP<Epetra_MultiVector> get_total_component_concentration () { return total_component_concentration; };
  Teuchos::RCP<Epetra_Vector>      get_specific_storage() { return specific_storage; }
  Teuchos::RCP<Epetra_Vector>      get_specific_yield() { return specific_yield; }
  
  const Teuchos::RCP<Amanzi::AmanziMesh::Mesh> get_mesh_maps() const { return mesh_maps; };

  double get_time () const { return time; };
  double get_last_time () const { return last_time; }
  int get_cycle () const { return cycle; };

  int get_number_of_components() const { return number_of_components; };

  const Amanzi::AmanziMesh::Mesh& get_mesh() { return *mesh_maps; };

  int get_component_number(const std::string component_name);
  std::string get_component_name(const int component_number);

  // modify methods
  void set_time ( double new_time );
  void set_cycle ( int new_cycle );
  void advance_time(double dT);
  void update_total_component_concentration(Teuchos::RCP<Epetra_MultiVector>);
  void update_total_component_concentration(const Epetra_MultiVector&);
  void update_darcy_flux(const Epetra_Vector&);
  void update_pressure(const Epetra_Vector&);
  void update_lambda(const Epetra_Vector&);

  // status methods
  status_type get_status () const { return status; };
  void set_status ( status_type new_status ) { status = new_status; }


  // debug helpers
  void set_darcy_flux( const double* u, const int mesh_block_id );
  void set_darcy_flux( const double* u, const std::string region );
  void set_water_saturation(const double ws );
  void set_water_density(const double wd );
  void set_zero_total_component_concentration();
  void set_total_component_concentration(const double* conc, const int mesh_block_id); 
  void set_total_component_concentration(const double* conc, const std::string region ); 
  void set_free_ion_concentrations(const double* conc, const std::string region ); 
  void set_porosity( const double phi );
  void set_porosity( const double phi, const int mesh_block_id );
  void set_porosity( const double phi, const std::string region );
  void set_permeability (const double kappa);
  void set_permeability (const double kappa, const int mesh_block_id);
  void set_permeability (const double kappa, const std::string region);
  void set_horizontal_permeability (const double kappa);
  void set_horizontal_permeability (const double kappa, const int mesh_block_id);
  void set_horizontal_permeability (const double kappa, const std::string region);  
  void set_vertical_permeability (const double kappa);
  void set_vertical_permeability (const double kappa, const int mesh_block_id);
  void set_vertical_permeability (const double kappa, const std::string region);
  void set_viscosity(const double mu);
  void set_gravity(const double *g);
  void set_number_of_components(const int n);
  void set_specific_storage(const double ss, const std::string region);
  void set_specific_yield(const double sy, const std::string region);

  // set methods 
  void set_darcy_flux ( const Epetra_Vector& darcy_flux_ );
  void set_water_saturation ( const Epetra_Vector& water_saturation_ );
  void set_prev_water_saturation ( const Epetra_Vector& prev_water_saturation_ );
  void set_water_density ( const Epetra_Vector& water_density_ );
  void set_porosity ( const Epetra_Vector& porosity_ );
  void set_permeability ( const Epetra_Vector& permeability_ );
  void set_vertical_permeability ( const Epetra_Vector& permeability_ );
  void set_horizontal_permeability ( const Epetra_Vector& permeability_ );
  void set_pressure ( const Epetra_Vector& pressure_ );
  void set_lambda ( const Epetra_Vector& lambda_ );
  void set_darcy_velocity ( const Epetra_MultiVector& darcy_velocity_ );
  void set_total_component_concentration ( const Epetra_MultiVector& total_component_concentration_ );
  void set_material_ids ( const Epetra_Vector& material_ids_ );
  void set_particle_density( const Epetra_Vector& particle_density_);


  void set_linear_pressure ( const Teuchos::ParameterList& ic_list, const std::string& region );
  void set_uniform_pressure ( const Teuchos::ParameterList& ic_list, const std::string& region );
  void set_file_pressure ( const Teuchos::ParameterList& ic_list, const std::string& region );
  void set_linear_saturation ( const Teuchos::ParameterList& ic_list, const std::string& region );
  void set_uniform_saturation ( const Teuchos::ParameterList& ic_list, const std::string& region );

  // observation functions
  double water_mass();
  double point_value(const std::string& point_region, const std::string& comp_name);
      
  void DeriveDarcyVelocity();

  void create_storage();
  void CreateStoragePrimarySpecies(void);
  void CreateStorageSecondaryActivityCoeff(const int size);
  void CreateStorageMinerals(void);
  void CreateStorageTotalSorbed(void);
  void CreateStorageIonExchange(void);
  void CreateStorageSurfaceComplexation(void);
  void CreateStorageSorptionIsotherms(void);

  void write_vis (Amanzi::Vis& vis, bool chemistry_enabled,  bool force=false);
  void write_vis (Amanzi::Vis& vis, 
                  Teuchos::RCP<Epetra_MultiVector> auxdata, 
                  const std::vector<std::string>& auxnames, 
                  bool chemistry_enabled, bool force=false);
  void set_compnames(std::vector<std::string>& compnames_);
  void set_compnames(Teuchos::Array<std::string>& compnames_);

  void SetupSoluteNames(void);
  void SetupMineralNames(void);
  void SetupSorptionSiteNames(void);
  void ExtractVolumeFromMesh(void);
  void VerifyMaterialChemistry(void);
  void VerifyMineralogy(const std::string& region_name,
                        const Teuchos::ParameterList& minerals_list);
  void VerifySorptionIsotherms(const std::string& region_name,
                               const Teuchos::ParameterList& isotherms_list);
  void VerifySorptionSites(const std::string& region_name,
                           const Teuchos::ParameterList& sorption_sites_list);
  void SetRegionMaterialChemistry(const std::string& region_name,
                                  Teuchos::ParameterList* region_data);
  void SetRegionMineralogy(const std::string& region_name,
                           const Teuchos::ParameterList& mineralogy_list);
  void SetRegionSorptionIsotherms(const std::string& region_name,
                          const Teuchos::ParameterList& isotherm_list);
  void SetRegionSorptionSites(const std::string& region_name,
                              const Teuchos::ParameterList& sorption_sites_list);

  void WriteChemistryToVis(Amanzi::Vis* vis);
  void WriteFreeIonsToVis(Amanzi::Vis* vis);
  void WriteTotalSorbedToVis(Amanzi::Vis* vis);
  void WriteMineralsToVis(Amanzi::Vis* vis);
  void WriteIsothermsToVis(Amanzi::Vis* vis);
  void WriteSorptionSitesToVis(Amanzi::Vis* vis);
  void WriteIonExchangeSitesToVis(Amanzi::Vis* vis);

  Teuchos::RCP<const Epetra_Vector> volume() const {
    return volume_;
  }

  Teuchos::RCP<Epetra_MultiVector> free_ion_concentrations() const {
    return free_ion_concentrations_;
  }

  void set_free_ion_concentrations(const Epetra_MultiVector& free_ion_conc) {
    *free_ion_concentrations_ = free_ion_conc;
  }

  Teuchos::RCP<Epetra_MultiVector> primary_activity_coeff() const {
    return primary_activity_coeff_;
  }

  Teuchos::RCP<Epetra_MultiVector> secondary_activity_coeff() const {
    return secondary_activity_coeff_;
  }

  int number_of_minerals(void) const {
    return number_of_minerals_;
  }

  void set_number_of_minerals(const int n) {
    number_of_minerals_ = n;
  }

  std::vector<std::string> mineral_names(void) const {
    return mineral_names_;
  }

  void set_mineral_names(std::vector<std::string> mineral_names) {
    mineral_names_ = mineral_names;
  }
  Teuchos::RCP<Epetra_MultiVector> mineral_volume_fractions() const {
    return mineral_volume_fractions_;
  }

  void set_mineral_volume_fractions(const Epetra_MultiVector& mineral_volume_fractions) {
    *mineral_volume_fractions_ = mineral_volume_fractions;
  }

  Teuchos::RCP<Epetra_MultiVector> mineral_specific_surface_area() const {
    return mineral_specific_surface_area_;
  }

  void set_mineral_specific_surface_area(const Epetra_MultiVector& mineral_specific_surface_area) {
    *mineral_specific_surface_area_ = mineral_specific_surface_area;
  }

  bool using_sorption(void) const {
    return using_sorption_;
  }

  void set_using_sorption(const bool value) {
    using_sorption_ = value;
  }

  Teuchos::RCP<Epetra_MultiVector> total_sorbed() const {
    return total_sorbed_;
  }

  void set_total_sorbed(const Epetra_MultiVector& total_sorbed) {
    *total_sorbed_ = total_sorbed;
  }

  int number_of_ion_exchange_sites(void) const {
    return number_of_ion_exchange_sites_;
  }

  void set_number_of_ion_exchange_sites(const int n) {
    number_of_ion_exchange_sites_ = n;
  }

  Teuchos::RCP<Epetra_MultiVector> ion_exchange_sites() const {
    return ion_exchange_sites_;
  }

  void set_ion_exchange_sites(const Epetra_MultiVector& ion_exchange_sites) {
    *ion_exchange_sites_ = ion_exchange_sites;
  }

  Teuchos::RCP<Epetra_MultiVector> ion_exchange_ref_cation_conc() const {
    return ion_exchange_ref_cation_conc_;
  }

  Teuchos::RCP<Epetra_MultiVector> surface_complex_free_site_conc() const {
    return surface_complex_free_site_conc_;
  }  

  int number_of_sorption_sites(void) const {
    return number_of_sorption_sites_;
  }

  void set_number_of_sorption_sites(const int n) {
    number_of_sorption_sites_ = n;
  }

  Teuchos::RCP<Epetra_MultiVector> sorption_sites() const {
    return sorption_sites_;
  }

  void set_sorption_sites(const Epetra_MultiVector& sorption_sites) {
    *sorption_sites_ = sorption_sites;
  }

  bool use_sorption_isotherms(void) const {
    return use_sorption_isotherms_;
  }

  void set_use_sorption_isotherms(const bool value) {
    use_sorption_isotherms_ = value;
  }

  Teuchos::RCP<Epetra_MultiVector> isotherm_kd() const {
    return isotherm_kd_;
  }

  Teuchos::RCP<Epetra_MultiVector> isotherm_freundlich_n() const {
    return isotherm_freundlich_n_;
  }

  Teuchos::RCP<Epetra_MultiVector> isotherm_langmuir_b() const {
    return isotherm_langmuir_b_;
  }

 private:
  void initialize_from_parameter_list();
  void init_verbosity (Teuchos::ParameterList &parameter_list_);
  void create_default_compnames(int n);
  void set_cell_value_in_mesh_block(double value, Epetra_Vector &v,
				    int mesh_block_id);

  void set_cell_value_in_region(const double& value, Epetra_Vector& v,
				const std::string& region);

  void set_cell_value_in_region(const Amanzi::Function& fun, Epetra_Vector &v,
				const std::string& region);

  void set_cell_value_in_region(const Epetra_Vector& x, Epetra_Vector& v,
				const std::string& region);

  void write_vis_ (Amanzi::Vis& vis, bool chemistry_enabled,  bool force=false);
  void write_vis_ (Amanzi::Vis& vis, 
		   Teuchos::RCP<Epetra_MultiVector> auxdata, 
		   const std::vector<std::string>& auxnames, 
		   bool chemistry_enabled, bool force=false);
  // state vectors
  Teuchos::RCP<Epetra_Vector> water_density;  
  Teuchos::RCP<Epetra_Vector> pressure;
  Teuchos::RCP<Epetra_Vector> lambda;
  Teuchos::RCP<Epetra_Vector> darcy_flux;
  Teuchos::RCP<Epetra_Vector> porosity;
  Teuchos::RCP<Epetra_Vector> particle_density;
  Teuchos::RCP<Epetra_MultiVector> total_component_concentration;
  Teuchos::RCP<Epetra_Vector> water_saturation;
  Teuchos::RCP<Epetra_Vector> prev_water_saturation;
  Teuchos::RCP<Epetra_Vector> horizontal_permeability;
  Teuchos::RCP<Epetra_Vector> vertical_permeability;  
  Teuchos::RCP<Epetra_MultiVector> darcy_velocity;
  Teuchos::RCP<Epetra_Vector> material_ids;
  Teuchos::RCP<Epetra_Vector> volume_;

  Teuchos::RCP<Epetra_MultiVector> free_ion_concentrations_; 
  Teuchos::RCP<Epetra_MultiVector> primary_activity_coeff_;
  Teuchos::RCP<Epetra_MultiVector> secondary_activity_coeff_;
  Teuchos::RCP<Epetra_MultiVector> mineral_volume_fractions_; // [cell][mineral]
  Teuchos::RCP<Epetra_MultiVector> mineral_specific_surface_area_; // [cell][mineral]
  Teuchos::RCP<Epetra_MultiVector> total_sorbed_;  // [cell][species]
  Teuchos::RCP<Epetra_MultiVector> sorption_sites_;  // [cell][site], eventually [cell][mineral][site]
  Teuchos::RCP<Epetra_MultiVector> surface_complex_free_site_conc_;
  Teuchos::RCP<Epetra_MultiVector> ion_exchange_sites_; // CEC, [cell][mineral]
  Teuchos::RCP<Epetra_MultiVector> ion_exchange_ref_cation_conc_;
  Teuchos::RCP<Epetra_MultiVector> isotherm_kd_; // [cell][species]
  Teuchos::RCP<Epetra_MultiVector> isotherm_freundlich_n_; // [cell][species]
  Teuchos::RCP<Epetra_MultiVector> isotherm_langmuir_b_; // [cell][species]
  Teuchos::RCP<Epetra_Vector> specific_storage; 
  Teuchos::RCP<Epetra_Vector> specific_yield; 

  Teuchos::RCP<double*> gravity;
  Teuchos::RCP<double> density;
  Teuchos::RCP<double> viscosity;
  
  int number_of_components;
  std::map<std::string,int> comp_no;

  double time, last_time;
  int cycle;
  status_type status;

  // mesh
  const Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh_maps;

  // parameter list
  Teuchos::ParameterList parameter_list;

  // names for components
  std::vector<std::string> compnames;
  int number_of_minerals_;
  std::vector<std::string> mineral_names_;
  std::map<std::string, int> mineral_name_id_map_;
  bool using_sorption_;
  bool use_sorption_isotherms_;
  int number_of_ion_exchange_sites_;
  int number_of_sorption_sites_;
  std::vector<std::string> sorption_site_names_;
  std::map<std::string, int> sorption_site_name_id_map_;
};

#endif
