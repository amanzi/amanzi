/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
   Amanzi

   License: see $AMANZI_DIR/COPYRIGHT
   Author: Svetlana Tokareva (tokareva@lanl.gov)

   Helper class for soil unit tests
   ------------------------------------------------------------------------- */

#ifndef SOIL_THERMO_TEST_CLASS_HH_
#define SOIL_THERMO_TEST_CLASS_HH_

#include <string>
#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "State.hh"
#include "soil_thermo_pk.hh"

using namespace Amanzi;

class SoilThermoTest {

public:
  // constuctor
  SoilThermoTest(Teuchos::ParameterList& plist,
                const Teuchos::RCP<AmanziMesh::Mesh>& mesh,
                int num_components);

  // init
  void initialize();

  // data
  Teuchos::ParameterList parameter_list;
  Teuchos::RCP<State> S0, S1;
  Teuchos::RCP<Soil_Thermo_PK> LTPK;
  Teuchos::RCP<AmanziMesh::Mesh> mesh;
  int num_components;
  double constant_value;

  // helper analytic solutions
  virtual double my_f(const AmanziGeometry::Point& x, double t) = 0;
  virtual AmanziGeometry::Point my_u(const AmanziGeometry::Point& x, double t);

  // helper evaluation methods
  void commit_step();
  void initialize_owned();
  void initialize_mass_flux();
  void evaluate_error_temp(double t, double* L1, double* L2);
};

class SoilThermoTestTestOne : public SoilThermoTest {
public:
  SoilThermoTestTestOne(Teuchos::ParameterList& plist,
                   const Teuchos::RCP<AmanziMesh::Mesh>& mesh,
                   int num_components);

  virtual double my_f(const AmanziGeometry::Point& x, double t);
};

class SoilThermoTestStep : public SoilThermoTest {
public:
  SoilThermoTestStep(Teuchos::ParameterList& plist,
                   const Teuchos::RCP<AmanziMesh::Mesh>& mesh,
                   int num_components);

  virtual double my_f(const AmanziGeometry::Point& x, double t);
};

class SoilThermoTestSmooth : public SoilThermoTest {
public:
  SoilThermoTestSmooth(Teuchos::ParameterList& plist,
                   const Teuchos::RCP<AmanziMesh::Mesh>& mesh,
                   int num_components);

  virtual double my_f(const AmanziGeometry::Point& x, double t);
};

class SoilThermoTestCubic : public SoilThermoTest {
public:
  SoilThermoTestCubic(Teuchos::ParameterList& plist,
                   const Teuchos::RCP<AmanziMesh::Mesh>& mesh,
                   int num_components);

  virtual double my_f(const AmanziGeometry::Point& x, double t);
};

class SoilThermoTestTwoDOne : public SoilThermoTest {
public:
  SoilThermoTestTwoDOne(Teuchos::ParameterList& plist,
                   const Teuchos::RCP<AmanziMesh::Mesh>& mesh,
                   int num_components);

  virtual double my_f(const AmanziGeometry::Point& x, double t);
  virtual AmanziGeometry::Point my_u(const AmanziGeometry::Point& x, double t);
};

#endif

