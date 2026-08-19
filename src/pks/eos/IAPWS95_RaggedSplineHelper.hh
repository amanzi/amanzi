/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_IAPWS95_RAGGED_SPLINE_HELPER_HH_
#define AMANZI_IAPWS95_RAGGED_SPLINE_HELPER_HH_

#include <array> 
#include <vector> 

namespace Amanzi {
namespace AmanziEOS {

class IAPWS95_RaggedSplineHelper {
 public:
  struct BasisData {
    std::array<int, 4> index{};
    std::array<double, 4> value{};
    std::array<double, 4> d1{};
    std::array<double, 4> d2{};
  };

  struct CubicSpanData {
    int first_basis;

    double x0;
    double inv_h;

    // N_r(x) = poly[r][0] + poly[r][1] x + poly[r][2] x^2 + poly[r][3] x^3
    double poly[4][4];
  };

  struct Mesh {
    // Physical coordinate lines used to describe and sample the domain
    std::vector<double> x_lines;  // rho or p
    std::vector<double> y_lines;  // T or h

    // Open/clamped degree-3 knot vectors in delta and tau
    std::vector<double> x_knots;
    std::vector<double> y_knots;

    int nx_basis_ = 0;
    int ny_basis_ = 0;

    int CoefficientIndex(int i, int j) const noexcept { return i * ny_basis_ + j; }

    // cache to avoid repeated calcualtions
    std::vector<CubicSpanData> x_span_data;
    std::vector<CubicSpanData> y_span_data;
  };

  std::vector<double> MakeClampedCubicKnots(const std::vector<double>& coordinate_lines);
  BasisData EvaluateCubicBasis(const std::vector<double>& knots, double x) const;
  BasisData EvaluateCachedCubicBasis(const std::vector<double>& U,
                                     const std::vector<CubicSpanData>& cache,
                                     double x) const;

  std::array<double, 6> Evaluate(const Mesh& mesh,
                                 const std::vector<double>& coefficients,
                                 double x, double y) const;

  double MinSpacing(const std::vector<double>& x) const;
  int LowerCell(const std::vector<double>& x, double value);

 protected:
  std::vector<CubicSpanData> BuildCubicSpanCache(const std::vector<double>& U);
};

}  // namespace AmanziEOS
}  // namespace Amanzi

#endif
