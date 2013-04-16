/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
   Amanzi

   License: see $AMANZI_DIR/COPYRIGHT
   Author: Ethan Coon

   Helper class for flow unit tests
   ------------------------------------------------------------------------- */

#ifndef FLOW_TEST_CLASS_HH_
#define FLOW_TEST_CLASS_HH_

#include <string>
#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "state.hh"
#include "richards.hh"

using namespace Amanzi;

class FlowTest {

public:
  // constuctor
  FlowTest(Teuchos::ParameterList& plist,
                const Teuchos::RCP<AmanziMesh::Mesh>& mesh,
                int num_components);

  // init
  void initialize();

  // data
  Teuchos::ParameterList parameter_list;
  Teuchos::RCP<State> S0, S1;
  Teuchos::RCP<Flow::Richards> FPK;
  Teuchos::RCP<AmanziMesh::Mesh> mesh;
  int num_components;
  double constant_value;

  // helper analytic solutions
  virtual double my_f(const AmanziGeometry::Point& x, double t) = 0;

  // helper evaluation methods
  void commit_step();
  void initialize_owned();
  void evaluate_error_pressure(double t, double* L1, double* L2);
};

class FlowTestOne : public FlowTest {
public:
  FlowTestOne(Teuchos::ParameterList& plist,
                   const Teuchos::RCP<AmanziMesh::Mesh>& mesh,
                   int num_components);

  virtual double my_f(const AmanziGeometry::Point& x, double t);
};

class FlowTestTwoDOne : public FlowTest {
public:
  FlowTestTwoDOne(Teuchos::ParameterList& plist,
                   const Teuchos::RCP<AmanziMesh::Mesh>& mesh,
                   int num_components);

  virtual double my_f(const AmanziGeometry::Point& x, double t);
};

#endif

