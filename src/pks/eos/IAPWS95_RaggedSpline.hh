/*
  Globally C2 cubic tensor-product B-spline approximation of the residual
  dimensionless Helmholtz energy phi^r(delta,tau) from IAPWS-95.
 
  Public inputs use physical variables: rho [kg/m^3], T [K].

  Internal coordinates are: delta = rho / RHOC, tau = TC / T.
 
  ResidualPart(rho,T) returns, in this order,
    ar, 
    d(ar)/d(delta), d(ar)/d(tau), 
    d2(ar)/d(delta)^2, d2(ar)/(d(delta)d(tau)), d2(ar)/d(tau)^2.
 
  The spline space is defined on a rectangular background knot grid, while
  fitting samples are retained only in a ragged domain T >= F(rho), plus a
  narrow extension band below F. The extension can use direct homogeneous
  IAPWS-95 values inside the saturation dome.
*/

#ifndef AMANZI_IAPWS95_RAGGED_SPLINE_HH_
#define AMANZI_IAPWS95_RAGGED_SPLINE_HH_

#include <array>
#include <cstddef>
#include <functional>
#include <memory>
#include <string>
#include <utility>
#include <vector>

// Amanzi::EOS
#include "IAPWS95.hh"

namespace Amanzi {
namespace AmanziEOS {

class IAPWS95_RaggedSpline : public IAPWS95 {
 public:
  struct Options {
    Options() {};

    // Temperature range used to construct the ragged domain.
    double T_min = 275.0;
    double T_max = 1000.0;

    // Background density range. rho_min must be positive.
    double rho_min = 1.0e-5;
    double rho_max = 1200.0;

    // Initial numbers of intervals. Adaptive refinement can add more.
    int initial_rho_intervals = 40;
    int initial_T_intervals   = 40;

    // Maximum numbers of intervals after adaptive refinement.
    int max_rho_intervals = 120;
    int max_T_intervals   = 120;

    // Number of fitting points per active background cell and direction.
    // Two gives four interior fitting points per active cell.
    unsigned samples_per_cell_direction = 2;

    // Extension width below F(rho), measured in local T-cell widths.
    double extension_cells = 4.0;
    double extension_weight = 0.15;

    // Relative fitting weights for value and delta/tau derivatives.
    std::array<double, 6> fit_weights = {1.0, 0.25, 0.25, 0.05, 0.05, 0.05};

    // Adaptive indicators are normalized errors in the six outputs.
    std::array<double, 6> abs_tolerance = {1e-10, 1e-9, 1e-9, 1e-8, 1e-8, 1e-8};
    std::array<double, 6> rel_tolerance = {1e-8,  1e-7, 1e-7, 1e-6, 1e-6, 1e-6};

    unsigned max_adaptive_passes = 8;
  };

  struct SaturationPoint {
    double T = 0.0;
    double rho_l = 0.0;
    double rho_v = 0.0;
    double p_sat = 0.0;
  };

  struct RaggedColumn {
    double rho = 0.0;
    double boundary_T = 0.0;  // SaturationLine(rho)
    int first_physical_T = 0;
    int first_extended_T = 0;
  };

  struct Mesh {
    // Physical coordinate lines used to describe and sample the domain
    std::vector<double> rho_lines;
    std::vector<double> T_lines;

    // Open/clamped degree-3 knot vectors in delta and tau
    std::vector<double> delta_knots;
    std::vector<double> tau_knots;

    std::vector<SaturationPoint> saturation;
    std::vector<RaggedColumn> columns;
  };

  struct Sample {
    double rho = 0.0;
    double T = 0.0;
    double weight = 1.0;
  };

  struct BasisData {
    std::array<int, 4> index{};
    std::array<double, 4> value{};
    std::array<double, 4> d1{};
    std::array<double, 4> d2{};
  };

  IAPWS95_RaggedSpline(Teuchos::ParameterList& plist,
                       Options options = IAPWS95_RaggedSpline::Options());
  ~IAPWS95_RaggedSpline() {};

  virtual std::array<double, 6> ResidualPart(double rho, double T) override;

  // Construct saturation data, the parabola-like boundary F(rho), the
  // rectangular background grid, and the ragged column offsets.
  //
  // The default boundary construction interpolates the sampled saturation
  // branches and produces a smooth exclusion boundary above the dome. Users
  // can replace BoundaryTemperature() if their F(rho) is prescribed directly.
  void CreateRaggedMesh();

  // Assemble a global overdetermined least-squares system for the active
  // cubic tensor-product basis functions and solve for all spline
  // coefficients. Samples include values and optionally first/second
  // derivatives. Physical samples receive full weight; extension-band
  // samples receive extension_weight.
  void BuildSplineCoefficients();

  // access
  const Mesh& GetMesh() const noexcept { return mesh_; }

 private:
  std::array<double, 6> ResidualPartDimensionless_(double delta, double tau) const;

  void ValidateOptions_() const;
  void BuildInitialCoordinateLines_();
  void BuildRaggedColumns_();
  void AdaptiveRefineCoordinateLines_();
  void BuildKnotVectors_();
  std::vector<double> MakeClampedCubicKnots_(const std::vector<double>& coordinate_lines);

  std::vector<Sample> BuildSamples_() const;
  bool IsPhysical_(double rho, double T) const { return T >= BoundaryTemperature(rho); }
  bool IsExtended_(double rho, double T) const;
  BasisData EvaluateCubicBasis_(const std::vector<double>& knots, double x) const;

  double BoundaryTemperature(double rho) const;
  double LocalTemperatureSpacing(double rho) const;

  int CoefficientIndex(int i, int j) const noexcept { return i * n_tau_basis_ + j; }

 private:
  std::shared_ptr<IAPWS95> eos95_;

  static constexpr int degree_ = 3;

  Options options_;
  Mesh mesh_;

  // Coefficients are stored with tau basis index changing fastest:
  // coefficient(i_delta, j_tau) = coefficients_[i_delta * n_tau_basis + j_tau].
  std::vector<double> coefficients_;
  int n_delta_basis_ = 0;
  int n_tau_basis_ = 0;
  bool built_ = false;
};

} // namespace AmanziEOS
} // namespace Amanzi

#endif
