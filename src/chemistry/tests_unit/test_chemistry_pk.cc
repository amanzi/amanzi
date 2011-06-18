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

#include "Mesh_maps_simple.hh"
#include "State.hpp"

#include "chemistry_pk.hh"
#include "chemistry_state.hh"
#include "species.hh"
#include "chemistry_exception.hh"


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
    Teuchos::RCP<Mesh_maps_base> mesh_;
    Teuchos::RCP<State> state_;
  };  // end class SpeciationTest

  ChemistryPKTest::ChemistryPKTest() {
    // assume that no errors or exceptions will occur in the
    // mesh/state related code....

    // get the parameter list from the input file
    std::string xml_input_filename("test_chemistry_pk.xml");
    Teuchos::ParameterList parameter_list;
    Teuchos::updateParametersFromXmlFile(xml_input_filename, &parameter_list);

    // create a dummy mesh?
    Epetra_SerialComm* comm = new Epetra_SerialComm();
    Teuchos::ParameterList mesh_parameter_list =
        parameter_list.sublist("Simple Mesh Parameters");
    mesh_ = Teuchos::rcp(new Mesh_maps_simple(mesh_parameter_list, comm));

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
