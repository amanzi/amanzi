/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Globally C2 cubic tensor-product B-spline approximation of the 
  dimensionless entropy from IAPWS-95.
 
  Public inputs use physical variables: p [MPa], h [kJ / kg].

  Internal coordinates are: pi = p / PC, theta = h / HC.
 
  Derivatives(pi, theta) returns, in this order,
    s, 
    d(s)/d(pi), d(s)/d(theta), 
    d2(s)/d(pi)^2, d2(s)/(d(pi)d(theta)), d2(s)/d(theta)^2.
 
  The spline space is defined on a rectangular background knot grid, while
  fitting samples are retained only in a ragged domain, plus a narrow 
  extension band. The extension can use direct homogeneous IAPWS-95 values 
  inside the saturation dome.

  Cells entirely below the extended boundary contribute no samples 
  (pi_p, theta_p). Cut cells contribute only the interior points that pass 
  IsExtended(rho, T). Additional smaples are dded only where Tb(rho) lies 
  on the flat critical cap. They constrain basis functions whose support
  crosses the flat critical cutoff. The fitting equations are

     sum_ij c_ij B_i(pi_p) B_j(theta_p) = Enthalpy(pi_p, theta_p)

  Because only four basis functions are nonzero in each direction, each of 
  these equations contains only 16 nonzero coefficient entries.

  The spline coefficeint minimize

    J(c) = sum_p wp sum_i g_i s_i^2 [L_i B(pi_p, theta_p) - f_pi]^2

  where L_i are six derivative operators of spline B. The operator weights
  are properly normalized:

    g_i = [1, 0.25, 0.25, 0.05, 0.05, 0.05],
    s_i = [1, hp,   hh,   hpp,  hph,  hhh].

  The banded Cholesky function DPBSV is used to solve for the vector of 
  coefficients c which minimizes J(c). Since coefficeints in the saturation
  dome may not be supported by samples, corresponding equations in the 
  normal system are nearly singular. We add a small diagonal regularization. 

  Adding extension samples can modify coefficients whose rectangular supports 
  cross the boundary, and those coefficients can also affect physical points 
  above the boundary over approximately four knot intervals in each coordinate.
  They redistribute error into nearby physical regions.
*/

#ifndef AMANZI_IAPWS95_RAGGED_SPLINE_P_H_HH_
#define AMANZI_IAPWS95_RAGGED_SPLINE_P_H_HH_

#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <mutex>
#include <string>
#include <utility>
#include <vector>

// Amanzi::EOS
#include "IAPWS95.hh"
#include "IAPWS95_RaggedSplineHelper.hh"

namespace Amanzi {
namespace AmanziEOS {

class IAPWS95_RaggedSplinePH : public IAPWS95,
                               public IAPWS95_RaggedSplineHelper {
 public:
  struct Options {
    Options() {};

    // background density range, rho_min must be positive
    double p_min = 0.10;
    double p_max = 50.0;

    // temperature range used to construct the ragged domain
    double h_min = 150.0;
    double h_max = 3800.0;

    // Initial numbers of intervals. Adaptive refinement can add more.
    int initial_p_intervals = 40;
    int initial_h_intervals = 40;

    // Maximum numbers of intervals after adaptive refinement.
    int max_p_intervals = 120;
    int max_h_intervals = 120;

    // Number of fitting points per active background cell and direction.
    // Two gives four interior fitting points per active cell.
    unsigned samples_per_cell_direction = 2;

    // Extension width below F(rho), measured in local T-cell widths.
    double extension_cells = 4.0;  // must be obsolete
    int extension_samples = 4;
    double extension_weight = 0.15;

    // Relative fitting weights for value and delta/tau derivatives
    std::array<double, 6> fit_weights = {1.0, 0.25, 0.25, 0.05, 0.05, 0.05};

    unsigned max_adaptive_passes = 8;

    // metastable extension fraction, must lie in (0,1)
    double metastable_stability_fraction = 0.25;

    // largest allowed extension below two-phase boundary, [K]
    double max_metastable_extension_K = 40.0;

    // initial density decrement used to locate the stability margin
    double metastable_density_fraction_step = 0.005;

    // minimum permitted temperature, [K]
    double minimum_extension_temperature_K = 273.16;

    // Number of bisection iterations after bracketing the margin.
    int metastable_bisection_iterations = 50;

    double critical_cutoff_temperature_K = 637.5;
    double critical_cap_extension_K = 4.0;
    double critical_cap_tolerance_K = 1.0e-8;

    // adaptive refinement
    double anisotropic_refinement_fraction = 0.20;
  };

  struct ActiveInterval {
    double lower = 0.0;
    double upper = 0.0;
    double weight = 1.0;
  };

  struct RaggedColumn {
    double p = 0.0;
    std::vector<ActiveInterval> physical_intervals;
    std::vector<ActiveInterval> extension_intervals;
  };

  struct Sample {
    double p = 0.0;
    double h = 0.0;
    double rho = 0.0;
    double T = 0.0;
    double weight = 1.0;
    bool is_physical = false;
  };

  struct CellError {
    double error;
    double error_p;
    double error_h;

    double p_mid;
    double h_mid;
  };

  IAPWS95_RaggedSplinePH(Teuchos::ParameterList& plist,
                         Options options = IAPWS95_RaggedSplinePH::Options());
  ~IAPWS95_RaggedSplinePH() {};

  virtual std::array<double, 6> EntropyDerivativesPH(double p, double h) override;

  // initialize shared data 
  std::vector<Sample> InitializeSharedData();

  // Construct rectangular background grid and active column intervals.
  void CreateRaggedMesh();

  // Assemble a global overdetermined least-squares system for the active
  // cubic tensor-product basis functions and solve for all spline
  // coefficients. Samples include values and optionally first/second
  // derivatives. Physical samples receive full weight; extension-band
  // samples receive extension_weight.
  void BuildSplineCoefficients(const std::vector<Sample>& samples);

  // Creates the set of (rho, T) points used to fit spline coefficients
  std::vector<Sample> BuildSamples();

  // point classification
  bool IsPhysical(double p, double h, const SaturationState& sat);
  bool IsPhysical(double p, double h) {
    const auto& sat = eos95_->SaturationLineP(p);
    return IsPhysical(p, h, sat);
  }

  double FindLiquidMetastableMargin(double p, const SaturationState& sat);
  double FindVaporMetastableMargin(double p, const SaturationState& sat);

  double StabilityFactor(double rho, double T);

  // mesh refinement
  void AnisotropicRefinement();

  // access
  const Mesh& GetMesh() const noexcept { return mesh_; }
  const Options& GetOptions() const { return options_; }

 private:
  void ValidateOptions_() const;
  void BuildInitialCoordinateLines_();
  void BuildRaggedColumns_();
  void AdaptiveRefineCoordinateLines_();
  void BuildKnotVectors_();

  // boundaries
  bool IsExtended_(double p, double h, const SaturationState& sat);
  std::pair<double, double> BoundaryEnthalpies(double p) const;

  // support of samples
  struct MetastableWork {
    int index;

    // Normalized distance from the appropriate saturation boundary.
    // Used only to order continuation from saturation into metastability.
    double depth;

    // Saturation state on the appropriate homogeneous branch at sample.p.
    double rho_sat;
    double T_sat;
  };

  void
  BuildMetastableWorkLists_(const std::vector<Sample>& samples,
                            std::vector<MetastableWork>& liquid,
                            std::vector<MetastableWork>& vapor);
  void
  PopulateMetastableSamples_(std::vector<Sample>& samples,
                             const std::vector<MetastableWork>& work);
  bool
  SolveMetastableSample_(Sample& sample,
                         const std::vector<Sample>& solved,
                         double rho_sat,
                         double T_sat);
  int FindNearestSolvedState_(const std::vector<Sample>& solved, double p, double h) const;

 public:
  std::uint64_t residual_calls = 0;  // statistics 

 private:
  std::shared_ptr<IAPWS95> eos95_;

  Options options_;
  std::vector<RaggedColumn> columns_;

  // shared data
  inline static Mesh mesh_;
  inline static std::vector<double> coefficients_;
  inline static std::once_flag init_shared_;
};

} // namespace AmanziEOS
} // namespace Amanzi

#endif
