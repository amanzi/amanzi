/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>
#include <typeinfo>

#include <UnitTest++.h>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Epetra_SerialComm.h"

#include "MeshFactory.hh"
#include "Domain.hh"
#include "GeometricModel.hh"
#include "GenerationSpec.hh"
#include "State.hpp"

#include "chemistry_pk.hh"
#include "chemistry_state.hh"
#include "species.hh"
#include "chemistry_exception.hh"

#include "dbc.hh"
#include "errors.hh"
#include "exceptions.hh"

/*****************************************************************************
 **
 **  Tests for trilinos based chemistry process kernel in chemistry-pk.cc
 **
 *****************************************************************************/

SUITE(GeochemistryTestsChemistryPK) {
  namespace ac = amanzi::chemistry;
  /*****************************************************************************
   **
   **  Common testing code
   **
   *****************************************************************************/
  class ChemistryPKTest {
   public:
    ChemistryPKTest();
    ~ChemistryPKTest();
 
    void RunTest(const std::string name, double* gamma);

   protected:
    ac::Chemistry_PK* cpk_;
    Teuchos::ParameterList chemistry_parameter_list_;
    Teuchos::RCP<ac::Chemistry_State> chemistry_state_;

   private:
    Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh_;
    Teuchos::RCP<State> state_;
  };  // end class SpeciationTest

  ChemistryPKTest::ChemistryPKTest() {
    // assume that no errors or exceptions will occur in the
    // mesh/state related code....

    // get the parameter list from the input file
    std::string xml_input_filename("test_chemistry_pk.xml");
    Teuchos::ParameterList parameter_list;
    Teuchos::updateParametersFromXmlFile(xml_input_filename, &parameter_list);

    Epetra_MpiComm* comm = new Epetra_MpiComm(MPI_COMM_WORLD);

    Amanzi::AmanziGeometry::Domain *simdomain_ptr = new Amanzi::AmanziGeometry::Domain(3);
    Teuchos::ParameterList& reg_params = parameter_list.sublist("Regions");
    Amanzi::AmanziGeometry::GeometricModelPtr 
      geom_model_ptr( new Amanzi::AmanziGeometry::GeometricModel(3, reg_params, comm) );
    simdomain_ptr->Add_Geometric_Model(geom_model_ptr);    

    // create a dummy mesh?
    // Epetra_SerialComm* comm = new Epetra_SerialComm();
    Teuchos::ParameterList mesh_parameter_list =
        parameter_list.sublist("Mesh Parameters");
 
    // Create a mesh factory for the geometric model
    Amanzi::AmanziMesh::MeshFactory factory(comm) ;
    // Prepare to read/create the mesh specification

    // get the Mesh sublist
    Teuchos::ParameterList mesh_params = parameter_list.sublist("Mesh");
    // Make sure the unstructured mesh option was chosen
    bool unstructured_option = mesh_params.isSublist("Unstructured");
    // Read and initialize the unstructured mesh parameters
    Teuchos::ParameterList unstr_mesh_params = mesh_params.sublist("Unstructured");
    Amanzi::AmanziMesh::FrameworkPreference prefs(Amanzi::AmanziMesh::default_preference());
    prefs.clear(); prefs.push_back(Amanzi::AmanziMesh::MSTK);
    factory.preference(prefs);
    Teuchos::ParameterList gen_params = unstr_mesh_params.sublist("Generate Mesh");
    mesh_ = factory.create(gen_params, geom_model_ptr);
    
    

    // get the state parameter list and create the state object
    Teuchos::ParameterList state_parameter_list = parameter_list.sublist("State");
    state_ = Teuchos::rcp(new State(state_parameter_list, mesh_));
 
    // create the chemistry state object from the state
    chemistry_state_ = Teuchos::rcp(new ac::Chemistry_State(state_));
 
    // create the chemistry parameter list
    chemistry_parameter_list_ = parameter_list.sublist("Chemistry");
  }

  ChemistryPKTest::~ChemistryPKTest() {
    delete cpk_;
  }

  void ChemistryPKTest::RunTest(const std::string name, double * gamma) {
  }  // end ChemistryPKTest::RunTest()

  /*****************************************************************************
   **
   **  individual tests
   **
   *****************************************************************************/
  TEST_FIXTURE(ChemistryPKTest, ChemistryPK_constructor) {
    // just make sure that we can have all the pieces together to set
    // up a chemistry process kernel....
    try {
      cpk_ = new ac::Chemistry_PK(chemistry_parameter_list_, chemistry_state_);
    } catch (ac::ChemistryException chem_error) {
      std::cout << chem_error.what() << std::endl;
    } catch (std::exception e) {
      std::cout << e.what() << std::endl;
    }
    // verbosity should be set after the constructor is finished....
    CHECK_EQUAL(cpk_->verbosity(), 1);
  }  // end TEST_FIXTURE()

  TEST_FIXTURE(ChemistryPKTest, ChemistryPK_initialize) {
    // just make sure that we can have all the pieces together to set
    // up a chemistry process kernel....
    try {
      cpk_ = new ac::Chemistry_PK(chemistry_parameter_list_, chemistry_state_);
      cpk_->InitializeChemistry();
    } catch (ac::ChemistryException chem_error) {
      std::cout << chem_error.what() << std::endl;
    } catch (std::exception e) {
      std::cout << e.what() << std::endl;
    }
    // XXXX is set by InitializeChemistry....
    CHECK_EQUAL(0, 0);
  }  // end TEST_FIXTURE()

  TEST_FIXTURE(ChemistryPKTest, ChemistryPK_get_chem_output_names) {
    cpk_ = new ac::Chemistry_PK(chemistry_parameter_list_, chemistry_state_);
    cpk_->InitializeChemistry();
    std::vector<std::string> names;
    cpk_->set_chemistry_output_names(&names);
    CHECK_EQUAL(names.at(0), "pH");
  }  // end TEST_FIXTURE()

  TEST_FIXTURE(ChemistryPKTest, ChemistryPK_set_component_names) {
    cpk_ = new ac::Chemistry_PK(chemistry_parameter_list_, chemistry_state_);
    cpk_->InitializeChemistry();
    std::vector<std::string> names;
    cpk_->set_component_names(&names);
    CHECK_EQUAL(names.at(0), "Al+++");
    CHECK_EQUAL(names.at(1), "H+");
    CHECK_EQUAL(names.at(2), "HPO4--");
    CHECK_EQUAL(names.at(3), "SiO2(aq)");
    CHECK_EQUAL(names.at(4), "UO2++");
  }  // end TEST_FIXTURE()
}  // end SUITE(GeochemistryTestChemistryPK)
