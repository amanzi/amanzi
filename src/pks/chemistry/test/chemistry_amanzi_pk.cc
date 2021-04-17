/*
  Chemistry PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.
*/

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>
#include <typeinfo>

#include <UnitTest++.h>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Epetra_SerialComm.h"
#include "XMLParameterListWriter.hh"

#include "dbc.hh"
#include "chemistry_exception.hh"
#include "errors.hh"
#include "exceptions.hh"
#include "GenerationSpec.hh"
#include "MeshFactory.hh"
#include "State.hh"

#include "Amanzi_PK.hh"
#include "Species.hh"


/*****************************************************************************
 **
 **  Tests for trilinos based chemistry process kernel in chemistry-pk.cc
 **
 *****************************************************************************/

SUITE(GeochemistryTestsChemistryPK) {
  namespace ac = Amanzi::AmanziChemistry;
  namespace ag = Amanzi::AmanziGeometry;
  namespace am = Amanzi::AmanziMesh;

  /* **************************************************************************
   * Common testing code
   * *************************************************************************/
  class ChemistryPKTest {
   public:
    ChemistryPKTest();
    ~ChemistryPKTest();
 
    void RunTest(const std::string name, double* gamma);

   protected:
    ac::Amanzi_PK* cpk_;
    Teuchos::ParameterList pk_tree_;
    Teuchos::RCP<Teuchos::ParameterList> glist_;
    Teuchos::RCP<Amanzi::State> state_;
    Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh_;

   private:
    Amanzi::Comm_ptr_type comm_;
    Teuchos::RCP<ag::GeometricModel> gm_;
  };  // end class SpeciationTest


  ChemistryPKTest::ChemistryPKTest() {
    // assume that no errors or exceptions will occur in the
    // mesh/state related code....
    
    // get the parameter list from the input file.
    std::string xml_input_filename("test/chemistry_amanzi_pk.xml");
    glist_ = Teuchos::getParametersFromXmlFile(xml_input_filename);

    // create a test mesh
    comm_ = Amanzi::getCommSelf();
    Teuchos::ParameterList mesh_parameter_list =
      glist_->sublist("mesh").sublist("unstructured").sublist("generate mesh");

    am::GenerationSpec g(mesh_parameter_list);
    
    Teuchos::ParameterList region_parameter_list = glist_->sublist("regions");
    gm_ = Teuchos::rcp(new ag::GeometricModel(3, region_parameter_list, *comm_));
  
    am::MeshFactory meshfactory(comm_, gm_);
    meshfactory.set_preference(am::Preference({am::Framework::SIMPLE}));

    mesh_ = meshfactory.create(0.,0.,0.,1.,1.,1.,1,1,10);

    // get the state parameter list and create the state object
    Teuchos::ParameterList state_parameter_list = glist_->sublist("state");

    state_ = Teuchos::rcp(new Amanzi::State(state_parameter_list));
    state_->RegisterDomainMesh(mesh_);

    // create the chemistry state object
    Teuchos::ParameterList chemistry_parameter_list = glist_->sublist("PKs").sublist("chemistry");
    std::vector<std::string> component_names;
    component_names.push_back("Al+++");
    component_names.push_back("H+");
    component_names.push_back("HP04--");
    component_names.push_back("SiO2(aq)");
    component_names.push_back("UO2++");

    // other input parameters in the constructor
    pk_tree_ = glist_->sublist("PK tree").sublist("chemistry");
  }

  ChemistryPKTest::~ChemistryPKTest() {}

  void ChemistryPKTest::RunTest(const std::string name, double * gamma) {};


  /* ***************************************************************************
   * Individual tests
   * **************************************************************************/
  TEST_FIXTURE(ChemistryPKTest, ChemistryPK_constructor) {
    // just make sure that we can have all the pieces together to set
    // up a chemistry process kernel....
    try {
      cpk_ = new ac::Amanzi_PK(pk_tree_, glist_, state_, Teuchos::null);
    } catch (ac::ChemistryException chem_error) {
      std::cout << "ERROR test1 "<< chem_error.what() << std::endl;
    } catch (std::exception e) {
      std::cout << "ERROR test1a " << e.what() << std::endl;
    }
  }  


  TEST_FIXTURE(ChemistryPKTest, ChemistryPK_initialize) {
    // make sure that we can initialize the pk and internal chemistry
    // object correctly based on the xml input....
    try {
      cpk_ = new ac::Amanzi_PK(pk_tree_, glist_, state_, Teuchos::null);
      cpk_->Setup(state_.ptr());
      state_->Setup();
      state_->InitializeFields();
      state_->InitializeEvaluators();
      cpk_->Initialize(state_.ptr());
    } catch (std::exception e) {
      std::cout << "ERROR test2 "<<e.what() << std::endl;
      throw e;
    }
    // assume all is right with the world if we exited w/o an error
    CHECK_EQUAL(0, 0);
  }


  TEST_FIXTURE(ChemistryPKTest, ChemistryPK_get_chem_output_names) {
    try {
      cpk_ = new ac::Amanzi_PK(pk_tree_, glist_, state_, Teuchos::null);
      cpk_->Setup(state_.ptr());
      state_->Setup();
      state_->InitializeFields();
      state_->InitializeEvaluators();
      cpk_->Initialize(state_.ptr());
    } catch (std::exception e) {
      std::cout << "ERROR test3 "<<e.what() << std::endl;
      throw e;
    }
    std::vector<std::string> names;
    cpk_->set_chemistry_output_names(&names);
    std::cout<<names.at(0)<<"\n";
    CHECK_EQUAL(names.at(0), "pH");
  }
}
