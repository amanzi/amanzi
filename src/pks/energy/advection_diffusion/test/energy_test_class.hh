/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
   Amanzi

   License: see $AMANZI_DIR/COPYRIGHT
   Author: Ethan Coon

   Helper class for energy unit tests
   ------------------------------------------------------------------------- */

#ifndef ENERGY_TEST_CLASS_HH_
#define ENERGY_TEST_CLASS_HH_

#include <string>
#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "State.hh"
#include "advection_diffusion.hh"

using namespace Amanzi;

class EnergyTest {

public:
  // constuctor
  EnergyTest(Teuchos::RCP<Teuchos::ParameterList> plist,
                const Teuchos::RCP<AmanziMesh::Mesh>& mesh,
                int num_components);

  // init
  void initialize();

  // data
  Teuchos::RCP<Teuchos::ParameterList> parameter_list;
  Teuchos::RCP<State> S0, S1;
  Teuchos::RCP<Energy::AdvectionDiffusion> EPK;
  Teuchos::RCP<AmanziMesh::Mesh> mesh;
  int num_components;
  double constant_value;

  // helper analytic solutions
  virtual double my_f(const AmanziGeometry::Point& x, double t) = 0;
  virtual AmanziGeometry::Point my_u(const AmanziGeometry::Point& x, double t) = 0;
  virtual double my_K() = 0;

  // helper evaluation methods
  void commit_step(double dt);
  void initialize_owned();
  void initialize_mass_flux();
  void evaluate_error_temp(double t, double* L1, double* L2);
};

class EnergyTestOne : public EnergyTest {
public:
  EnergyTestOne(Teuchos::RCP<Teuchos::ParameterList> plist,
                const Teuchos::RCP<AmanziMesh::Mesh>& mesh,
                int num_components) :
    EnergyTest(plist, mesh, num_components) {}      

  virtual double my_f(const AmanziGeometry::Point& x, double t) { return 1.; }
  virtual double my_K() { return 1.; }

  virtual AmanziGeometry::Point my_u(const AmanziGeometry::Point& x,
          double t) {
    return AmanziGeometry::Point(1.0, 0.0, 0.0);
  }
  
};

class EnergyTestStep : public EnergyTest {
public:
  EnergyTestStep(Teuchos::RCP<Teuchos::ParameterList> plist,
                 const Teuchos::RCP<AmanziMesh::Mesh>& mesh,
                 int num_components) :
    EnergyTest(plist, mesh, num_components) {}      


  virtual double my_f(const AmanziGeometry::Point& x, double t) {
    if (x[0] <= t) return 1.0;
    return 0;
  }
    
  virtual double my_K() { return 0.; }
  virtual AmanziGeometry::Point my_u(const AmanziGeometry::Point& x,
          double t) {
    return AmanziGeometry::Point(1.0, 0.0, 0.0);
  }
};

class EnergyTestDiffusedStep : public EnergyTest {
public:
  EnergyTestDiffusedStep(Teuchos::RCP<Teuchos::ParameterList> plist,
                         const Teuchos::RCP<AmanziMesh::Mesh>& mesh,
                         int num_components) :
    EnergyTest(plist, mesh, num_components) {}      

  virtual double my_f(const AmanziGeometry::Point& x, double t) {
    if (x[0] <= t) return 1.0;
    return 0;
  }

  virtual double my_K() { return 1.; }
  virtual AmanziGeometry::Point my_u(const AmanziGeometry::Point& x,
          double t) {
    return AmanziGeometry::Point(0.0, 0.0, 0.0);
  }
};

class EnergyTestAdvDiffusedStep : public EnergyTest {
public:
  EnergyTestAdvDiffusedStep(Teuchos::RCP<Teuchos::ParameterList> plist,
                         const Teuchos::RCP<AmanziMesh::Mesh>& mesh,
                         int num_components) :
    EnergyTest(plist, mesh, num_components) {}      

  virtual double my_f(const AmanziGeometry::Point& x, double t) {
    if (x[0] <= t) return 1.0;
    return 0;
  }
    
  virtual double my_K() { return 1.; }
  virtual AmanziGeometry::Point my_u(const AmanziGeometry::Point& x,
          double t) {
    return AmanziGeometry::Point(1.0, 0.0, 0.0);
  }
};

class EnergyTestSmooth : public EnergyTest {
public:
  EnergyTestSmooth(Teuchos::RCP<Teuchos::ParameterList> plist,
                   const Teuchos::RCP<AmanziMesh::Mesh>& mesh,
                   int num_components) :
    EnergyTest(plist, mesh, num_components) {}      


  virtual double my_f(const AmanziGeometry::Point& x, double t);
};

class EnergyTestCubic : public EnergyTest {
public:
  EnergyTestCubic(Teuchos::RCP<Teuchos::ParameterList> plist,
                  const Teuchos::RCP<AmanziMesh::Mesh>& mesh,
                  int num_components) :
    EnergyTest(plist, mesh, num_components) {}      

  virtual double my_f(const AmanziGeometry::Point& x, double t);
};

class EnergyTestTwoDOne : public EnergyTest {
public:
  EnergyTestTwoDOne(Teuchos::RCP<Teuchos::ParameterList> plist,
                    const Teuchos::RCP<AmanziMesh::Mesh>& mesh,
                    int num_components) :
    EnergyTest(plist, mesh, num_components) {}      

  virtual double my_f(const AmanziGeometry::Point& x, double t);
  virtual AmanziGeometry::Point my_u(const AmanziGeometry::Point& x, double t);
};

#endif

