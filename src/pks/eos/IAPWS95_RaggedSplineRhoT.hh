/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

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

  Cells entirely below the extended boundary contribute no samples 
  (delta_p, tau_p). Cut cells contribute only the interior points that pass 
  IsExtended(rho, T). Additional smaples are dded only where Tb(rho) lies 
  on the flat critical cap. They constrain basis functions whose support
  crosses the flat critical cutoff. The fitting equations are

     sum_ij c_ij B_i(delta_p) B_j(tau_p) = Helmholtz_r(delta_p, tau_p)

  Because only four basis functions are nonzero in each direction, each of 
  these equations contains only 16 nonzero coefficient entries.

  The spline coefficeint minimize

    J(c) = sum_p wp sum_i g_i s_i^2 [L_i B(delta_p, tau_p) - f_pi]^2

  where L_i are six derivative operators of spline B. The operator weights
  are properly normalized:

    g_i = [1, 0.25, 0.25, 0.05, 0.05, 0.05],
    s_i = [1, hd,   ht,   hdd,  hdt,  htt].

  The banded Cholesky function DPBSV is used to solve for the vector of 
  coefficients c which minimizes J(c). Since coefficeints in the saturation
  dome may not be supported by samples, corresponding equations in the 
  normal system are nearly singular. We add a small diagonal regularization. 

  Adding extension samples can modify coefficients whose rectangular supports 
  cross the boundary, and those coefficients can also affect physical points 
  above the boundary over approximately four knot intervals in each coordinate.
  They redistribute error into nearby physical regions.
*/

#ifndef AMANZI_IAPWS95_RAGGED_SPLINE_RHO_T_HH_
#define AMANZI_IAPWS95_RAGGED_SPLINE_RHO_T_HH_

#include <array>
#include <cstddef>
#include <functional>
#include <memory>
#include <string>
#include <utility>
#include <vector>

// Amanzi::EOS
#include "IAPWS95.hh"
#include "IAPWS95_RaggedSplineHelper.hh"

namespace Amanzi {
namespace AmanziEOS {

class IAPWS95_RaggedSplineRhoT : public IAPWS95,
                                 public IAPWS95_RaggedSplineHelper {
 public:
  struct Options {
    Options() {};

    // Rectangle supported by spline
    double rho_min = 1.0e-5;
    double rho_max = 1200.0;

    double T_min = 275.0;
    double T_max = 1000.0;

    // The initial and maximum numbers of intervals.
    int initial_rho_intervals = 40;
    int initial_T_intervals   = 40;

    int max_rho_intervals = 120;
    int max_T_intervals   = 120;

    unsigned max_adaptive_passes = 8;

    // Number of fitting points per active background cell and direction.
    // Two gives four interior fitting points per active cell.
    unsigned samples_per_cell_direction = 2;

    // Extension inside the metastable region.
    double extension_cells = 4.0;  // must be obsolete
    int extension_samples = 4;
    double extension_weight = 0.15;

    double metastable_stability_fraction = 0.25;
    double max_metastable_extension_K = 40.0;

    // initial temperature decrement used to locate the stability margin, [K]
    double metastable_scan_step_K = 0.5;

    // Relative fitting weights for value and delta/tau derivatives
    std::array<double, 6> fit_weights = {1.0, 0.25, 0.25, 0.05, 0.05, 0.05};

    // minimum permitted temperature, [K]
    double minimum_extension_temperature_K = 273.16;

    // Number of bisection iterations after bracketing the margin.
    int metastable_bisection_iterations = 50;

    double critical_cutoff_temperature_K = 637.5;
    double critical_cap_extension_K = 4.0;
    double critical_cap_tolerance_K = 1.0e-8;
  };

  struct SaturationPoint {
    double T = 0.0;
    double rho_l = 0.0;
    double rho_v = 0.0;
    double p = 0.0;
  };

  struct RaggedColumn {
    double rho = 0.0;  // column coordinate
    double boundary_T = 0.0;  // saturation line / boundary of spline domain
    double extension_T = 0.0;  // lower endpoint of metastable extension
    int first_physical_T = 0;
    int first_extended_T = 0;
  };

  struct Sample {
    double rho = 0.0;
    double T = 0.0;
    double weight = 1.0;
  };

  IAPWS95_RaggedSplineRhoT(Teuchos::ParameterList& plist,
                           Options options = IAPWS95_RaggedSplineRhoT::Options());
  ~IAPWS95_RaggedSplineRhoT() {};

  virtual std::array<double, 6> ResidualPart(double rho, double T) override;

  // Construct saturation data bounded by parabola-like function F(rho), the
  // rectangular background grid, and the ragged column offsets.
  //
  // The default boundary construction interpolates the sampled saturation
  // branches and produces a smooth exclusion boundary above the dome. 
  void CreateRaggedMesh();

  // Assemble a global overdetermined least-squares system for the active
  // cubic tensor-product basis functions and solve for all spline coefficients.
  // Samples include values and optionally first/second derivatives. Physical 
  // samples receive full weight; the metastable extension-band samples receive
  // extension_weight.
  void BuildSplineCoefficients(const std::vector<Sample>& samples);

  // Creates the set of (rho, T) points/samples used to fit spline coefficients.
  std::vector<Sample> BuildSamples();

  // access
  const Mesh& GetMesh() const noexcept { return mesh_; }

 private:
  void ValidateOptions_() const;
  void BuildInitialCoordinateLines_();
  void BuildRaggedColumns_();
  void AdaptiveRefineCoordinateLines_();
  void BuildKnotVectors_();

  bool IsPhysical_(double rho, double T) const { return T >= BoundaryTemperature(rho); }
  bool IsExtended_(double rho, double T);

  double BoundaryTemperature(double rho) const;
  double ExtensionTemperature(double rho) const;
  double MinTemperatureSpacing() const;

  double FindMetastableMarginTemperature(double rho, double boundary_T);
  double StabilityFactor(double rho, double T);

 private:
  std::shared_ptr<IAPWS95> eos95_;

  Options options_;
  Mesh mesh_;

  std::vector<SaturationPoint> saturation_;
  std::vector<RaggedColumn> columns_;

  // Coefficients are stored with tau basis index changing fastest:
  // coefficient(i_delta, j_tau) = coefficients_[i_delta * ny_basis + j_tau].
  std::vector<double> coefficients_;
  bool built_ = false;
};

} // namespace AmanziEOS
} // namespace Amanzi

#endif
