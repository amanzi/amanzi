/*
  This is the operators component of the Amanzi code.

  License: BSD
  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)

  Discrete source operator.
*/

#include "AnalyticBase.hh"

class Analytic03 : public AnalyticBase {
 public:
  Analytic03(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh) : AnalyticBase(mesh) {
    k1 = 1.0;
    k2 = 20.0;
    a1 = 1.0 / k1;
    a2 = 1.0 / k2;
    b2 = (a1 - a2) / 4;
  }
  ~Analytic03() {};

  Amanzi::WhetStone::Tensor Tensor(const Amanzi::AmanziGeometry::Point& p, double t) {
    double x = p[0];
    Amanzi::WhetStone::Tensor K(2, 1);
    if (x < 0.5) { 
      K(0, 0) = k1;
    } else {
      K(0, 0) = k2;
    }
    return K;
  }

  double pressure_exact(const Amanzi::AmanziGeometry::Point& p, double t) { 
    double x = p[0];
    double y = p[1];
    if (x < 0.5) { 
      return a1 * x * x + y * y;
    } else {
      return a2 * x * x + y * y + b2;
    }
  }

  Amanzi::AmanziGeometry::Point velocity_exact(const Amanzi::AmanziGeometry::Point& p, double t) { 
    double x = p[0];
    double y = p[1];
    Amanzi::AmanziGeometry::Point v(2);
    if (x < 0.5) { 
      v[0] = -2 * a1 * k1 * x;
      v[1] = -2 * k1 * y;
    } else {
      v[0] = -2 * a2 * k2 * x;
      v[1] = -2 * k2 * y;
    }
    return v;
  }
 
  Amanzi::AmanziGeometry::Point gradient_exact(const Amanzi::AmanziGeometry::Point& p, double t) { 
    double x = p[0];
    double y = p[1];
    Amanzi::AmanziGeometry::Point v(2);
    if (x < 0.5) { 
      v[0] = 2 * a1 * x;
      v[1] = 2 * k1 * y;
    } else {
      v[0] = 2 * a2 * x;
      v[1] = 2 * k2 * y;
    }
    return v;
  }

  double source_exact(const Amanzi::AmanziGeometry::Point& p, double t) { 
    double x = p[0];
    if (x < 0.5) { 
      return -2 * a1 * k1 - 2 * k1;
    } else {
      return -2 * a2 * k2 - 2 * k2;
    }
  }

 private:
  double k1, k2;
  double a1, a2, b2;
};

