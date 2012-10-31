/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
   Amanzi

   License: see $AMANZI_DIR/COPYRIGHT
   Author: Ethan Coon

   Helper class for transport unit tests
   ------------------------------------------------------------------------- */

#ifndef TRANSPORT_TEST_CLASS_HH_
#define TRANSPORT_TEST_CLASS_HH_

#include <string>
#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "state.hh"
#include "passive_tracer.hh"

using namespace Amanzi;

class TransportTest {

public:
  // constuctor
  TransportTest(Teuchos::ParameterList& plist,
                const Teuchos::RCP<AmanziMesh::Mesh>& mesh,
                int num_components);

  // init
  void initialize();

  // data
  Teuchos::ParameterList parameter_list;
  Teuchos::RCP<State> S0, S1;
  Teuchos::RCP<Transport::PassiveTracer> TPK;
  Teuchos::RCP<AmanziMesh::Mesh> mesh;
  int num_components;
  double constant_value;

  // helper analytic solutions
  virtual double my_f(const AmanziGeometry::Point& x, double t) = 0;
  virtual AmanziGeometry::Point my_u(const AmanziGeometry::Point& x, double t);

  // helper evaluation methods
  void commit_step();
  void initialize_tcc();
  void initialize_darcy_flux();
  void evaluate_error_tcc(double t, double* L1, double* L2);
};

class TransportTestOne : public TransportTest {
public:
  TransportTestOne(Teuchos::ParameterList& plist,
                   const Teuchos::RCP<AmanziMesh::Mesh>& mesh,
                   int num_components);

  virtual double my_f(const AmanziGeometry::Point& x, double t);
};

class TransportTestStep : public TransportTest {
public:
  TransportTestStep(Teuchos::ParameterList& plist,
                   const Teuchos::RCP<AmanziMesh::Mesh>& mesh,
                   int num_components);

  virtual double my_f(const AmanziGeometry::Point& x, double t);
};

class TransportTestSmooth : public TransportTest {
public:
  TransportTestSmooth(Teuchos::ParameterList& plist,
                   const Teuchos::RCP<AmanziMesh::Mesh>& mesh,
                   int num_components);

  virtual double my_f(const AmanziGeometry::Point& x, double t);
};

class TransportTestCubic : public TransportTest {
public:
  TransportTestCubic(Teuchos::ParameterList& plist,
                   const Teuchos::RCP<AmanziMesh::Mesh>& mesh,
                   int num_components);

  virtual double my_f(const AmanziGeometry::Point& x, double t);
};

class TransportTestTwoDOne : public TransportTest {
public:
  TransportTestTwoDOne(Teuchos::ParameterList& plist,
                   const Teuchos::RCP<AmanziMesh::Mesh>& mesh,
                   int num_components);

  virtual double my_f(const AmanziGeometry::Point& x, double t);
  virtual AmanziGeometry::Point my_u(const AmanziGeometry::Point& x, double t);
};

#endif

