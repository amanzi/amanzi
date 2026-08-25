/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <sstream>

#include "lapack.hh"

#include "IAPWS95_RaggedSplineRhoT.hh"

namespace Amanzi {
namespace AmanziEOS {

/* ******************************************************************
* Constructor
****************************************************************** */
IAPWS95_RaggedSplineRhoT::IAPWS95_RaggedSplineRhoT(Teuchos::ParameterList& plist, Options options)
  : IAPWS95(plist),
    options_(std::move(options))
{
  ValidateOptions_();
  eos95_ = std::make_shared<IAPWS95>(plist);
}


/* ******************************************************************
* Residual part of the Helmholtz free energy
****************************************************************** */
std::array<double, 6>
IAPWS95_RaggedSplineRhoT::ResidualPart(double rho, double T)
{ 
  AMANZI_ASSERT(built_);

  // AMANZI_ASSERT(rho >= options_.rho_min && rho <= options_.rho_max);
  // AMANZI_ASSERT(T >= options_.T_min && T <= options_.T_max);

  // double bndT = BoundaryTemperature(rho);
  // double extT = options_.extension_cells * min_T_spacing_;
  // AMANZI_ASSERT(T >= bndT - extT);

  return Evaluate(mesh_, coefficients_, rho / RHOC, TC / T);
}


/* ******************************************************************
* Wrapper for initialization of shared data
****************************************************************** */
void
IAPWS95_RaggedSplineRhoT::InitializeSharedData()
{
  CreateRaggedMesh();
  auto samples = BuildSamples();
  BuildSplineCoefficients(samples);

  AnisotropicRefinement();
  samples = BuildSamples();
  BuildSplineCoefficients(samples);
}


/* ******************************************************************
* Create a ragged mesh
****************************************************************** */
void
IAPWS95_RaggedSplineRhoT::CreateRaggedMesh()
{
  mesh_ = Mesh{};
  coefficients_.clear();
  built_ = false;

  BuildInitialCoordinateLines_();
  BuildSaturationData_();
  BuildRaggedColumns_();

  AdaptiveRefineCoordinateLines_();
  BuildKnotVectors_();

  min_T_spacing_ = MinSpacing(mesh_.y_lines);
}


/* ******************************************************************
* Assemble a global overdetermined least-squares system for the active
* cubic tensor-product basis functions and solve for all spline coefficients.
****************************************************************** */
void
IAPWS95_RaggedSplineRhoT::BuildSplineCoefficients(const std::vector<Sample>& samples)
{
  // CreateRaggedMesh() must be called before fitting
  AMANZI_ASSERT(mesh_.x_knots.size() > 0 && mesh_.y_knots.size() > 0);

  int n = mesh_.nx_basis_ * mesh_.ny_basis_;
  AMANZI_ASSERT(n != 0);  // The spline has no coefficients

  // With coefficient(i,j) = i * ny_basis + j, two overlapping cubic
  // tensor-product basis functions differ by at most 3 in each index.
  int kd = 3 * mesh_.ny_basis_ + 3;
  int ldab = kd + 1;

  // Upper-band LAPACK storage for G = A^T A and rhs = A^T b.
  // AB(kd + i - j, j) stores G(i,j), i <= j.
  std::vector<double> ab(ldab * n, 0.0);
  std::vector<double> rhs(n, 0.0);

  auto add_upper = [&](int i, int j, double value) {
    if (i > j) std::swap(i, j);
    const int d = j - i;
    AMANZI_ASSERT(d <= kd);  // Normal-matrix entry lies outside the band
    ab[(kd - d) + j * ldab] += value;
  };

  std::array<int, 16> index{};
  std::array<double, 16> row_value{};

  for (const Sample& sample : samples) {
    double delta = sample.rho / RHOC;
    double tau = TC / sample.T;
    const BasisData Bd = EvaluateCubicBasis(mesh_.x_knots, delta);
    const BasisData Bt = EvaluateCubicBasis(mesh_.y_knots, tau);

    // target values for fitting
    const std::array<double, 6> target = eos95_->ResidualPart(sample.rho, sample.T);

    int id = LowerCell(mesh_.x_lines, sample.rho);
    int it = LowerCell(mesh_.y_lines, sample.T);

    double tau0 = TC / mesh_.y_lines[it];
    double tau1 = TC / mesh_.y_lines[it + 1];

    const double h_delta = (mesh_.x_lines[id + 1] - mesh_.x_lines[id]) / RHOC;
    const double h_tau = std::fabs(tau1 - tau0);

    const std::array<double, 6> derivative_scale = { 1.0,
                                                     h_delta,
                                                     h_tau,
                                                     h_delta * h_delta,
                                                     h_delta * h_tau,
                                                     h_tau * h_tau };

    for (int component = 0; component < 6; ++component) {
      double s = derivative_scale[component];
      double weight = sample.weight * options_.fit_weights[component] * s * s;
      if (!(weight > 0.0)) continue;

      int q = 0;
      for (int a = 0; a < 4; ++a) {
        for (int b = 0; b < 4; ++b) {
          double xd, xt;
          switch (component) {
            case 0:
              xd = Bd.value[a];
              xt = Bt.value[b];
              break;
            case 1:
              xd = Bd.d1[a];
              xt = Bt.value[b];
              break;
            case 2:
              xd = Bd.value[a];
              xt = Bt.d1[b];
              break;
            case 3:
              xd = Bd.d2[a];
              xt = Bt.value[b];
              break;
            case 4:
              xd = Bd.d1[a];
              xt = Bt.d1[b];
              break;
            case 5:
              xd = Bd.value[a];
              xt = Bt.d2[b];
              break;
            default:
              ;  // pass
          }
          index[q] = mesh_.CoefficientIndex(Bd.index[a], Bt.index[b]);
          row_value[q] = xd * xt;
          q++;
        }
      }

      for (int a = 0; a < 16; ++a) {
        rhs[index[a]] += weight * row_value[a] * target[component];
        for (int b = a; b < 16; ++b) {
          add_upper(index[a], index[b], weight * row_value[a] * row_value[b]);
        }
      }
    }
  }

  // Small diagonal shift protects coefficients whose support is mostly in
  // the excluded part of a cut cell.
  double diagonal_sum = 0.0;
  for (int j = 0; j < n; ++j) {
     diagonal_sum += std::abs(ab[kd + j * ldab]);
  }
  double regularization = 1.0e-12 * std::max(1.0, diagonal_sum / n);
  for (int j = 0; j < n; ++j) {
    ab[kd + j * ldab] += regularization;
  }

  int nrhs(1), info(0);
  WhetStone::DPBSV_F77("U", &n, &kd, &nrhs, ab.data(), &ldab, rhs.data(), &n, &info);
  if (info != 0) AMANZI_ASSERT(false);
  // info < 0 : illegal argument 
  // info > 0 : normal matrix is not positive definite at pivot info

  coefficients_ = std::move(rhs);
  built_ = true;
}


/* ******************************************************************
* Build initial mesh
****************************************************************** */
void
IAPWS95_RaggedSplineRhoT::BuildInitialCoordinateLines_()
{
  mesh_.x_lines.resize(options_.initial_rho_intervals + 1);
  mesh_.y_lines.resize(options_.initial_T_intervals + 1);

  // Logarithmic density distribution handles the dilute-vapor scale while
  // remaining monotone and simple. We can replace by another distribution
  // if the liquid region requires stronger clustering.
  double log_min = std::log(options_.rho_min);
  double log_max = std::log(options_.rho_max);
  for (int i = 0; i < mesh_.x_lines.size(); ++i) {
    double s = (double)i / (mesh_.x_lines.size() - 1);
    mesh_.x_lines[i] = std::exp((1.0 - s) * log_min + s * log_max);
  }

  for (int j = 0; j < mesh_.y_lines.size(); ++j) {
    double s = (double)j / (mesh_.y_lines.size() - 1);
    mesh_.y_lines[j] = (1.0 - s) * options_.T_min + s * options_.T_max;
  }
}


/* ******************************************************************
* Should be called every time the coordinate lines change, e.g.
* during adaptive refinement.
****************************************************************** */
void
IAPWS95_RaggedSplineRhoT::BuildRaggedColumns_()
{
  columns_.clear();
  columns_.reserve(mesh_.x_lines.size());

  for (double rho : mesh_.x_lines) {
    double bndT = BoundaryTemperature(rho);

    RaggedColumn c;
    c.rho = rho;
    c.boundary_T = bndT;
    double extT = FindMetastableMarginTemperature(rho, bndT);
    c.extension_T = bndT;

    c.first_physical_T = static_cast<int>(
      std::lower_bound(mesh_.y_lines.begin(), mesh_.y_lines.end(), bndT) - mesh_.y_lines.begin());
    c.first_extended_T = static_cast<int>(
      std::lower_bound(mesh_.y_lines.begin(), mesh_.y_lines.end(), extT) - mesh_.y_lines.begin());

    columns_.push_back(c);
  }
}


/* ******************************************************************
* Validate options
****************************************************************** */
void IAPWS95_RaggedSplineRhoT::ValidateOptions_() const
{
  AMANZI_ASSERT(options_.rho_min > 0.0 && options_.rho_max > options_.rho_min);
  AMANZI_ASSERT(options_.T_max > options_.T_min && options_.T_min > 0.0);

  // At least four intervals per direction are required.
  AMANZI_ASSERT(options_.initial_rho_intervals >= 4 && options_.initial_T_intervals >= 4);

  // Samples_per_cell_direction must be positive.
  AMANZI_ASSERT(options_.samples_per_cell_direction != 0);
}


/* ******************************************************************
*
****************************************************************** */
double
IAPWS95_RaggedSplineRhoT::BoundaryTemperature(double rho) const
{
  // The saturation arrays parameterize two branches rho_v(T) and rho_l(T).
  // For a density between the vapor and liquid branches, return the highest
  // sampled saturation temperature whose dome interval contains rho. This 
  // provides a piecewise-linear approximation to the lower exclusion boundary. 
  double result = options_.T_min;
  bool found = false;

  for (int k = 0; k + 1 < saturation_.size(); ++k) {
    const auto& a = saturation_[k];
    const auto& b = saturation_[k + 1];

    auto branch_intersection = [&](double r0, double r1) {
      double rlo = std::min(r0, r1);
      double rhi = std::max(r0, r1);
      if (rho < rlo || rho > rhi || r0 == r1) return;
      double s = (rho - r0) / (r1 - r0);
      double T = a.T + s * (b.T - a.T);
      result = std::max(result, T);
      found = true;
    };

    branch_intersection(a.rho_v, b.rho_v);
    branch_intersection(a.rho_l, b.rho_l);

    // Densities strictly between both branches at this temperature lie
    // inside the dome. The upper sampled T is a conservative boundary.
    if (rho >= b.rho_v && rho <= b.rho_l) {
      result = std::max(result, b.T);
      found = true;
    }
  }

  return found ? std::min(result, TC) : options_.T_min;
}


/* ******************************************************************
*
****************************************************************** */
double
IAPWS95_RaggedSplineRhoT::ExtensionTemperature(double rho) const
{
  if (columns_.empty()) {
    return BoundaryTemperature(rho);
  }

  if (rho <= columns_.front().rho) {
    return columns_.front().extension_T;
  }

  if (rho >= columns_.back().rho) {
    return columns_.back().extension_T;
  }

  const auto it = std::upper_bound(columns_.begin(),
                                   columns_.end(),
                                   rho,
                                   [](double value, const RaggedColumn& column) {
                                     return value < column.rho;
                                   });

  int i1 = static_cast<int>(std::distance(columns_.begin(), it));
  int i0 = i1 - 1;

  const RaggedColumn& c0 = columns_[i0];
  const RaggedColumn& c1 = columns_[i1];

  double s = (rho - c0.rho) / (c1.rho - c0.rho);

  return (1.0 - s) * c0.extension_T + s * c1.extension_T;
}


/* ******************************************************************
* Metastable-margin function searches downward from two-phase boundary
****************************************************************** */
double
IAPWS95_RaggedSplineRhoT::FindMetastableMarginTemperature(double rho, double boundary_T)
{
  if (std::fabs(boundary_T - options_.critical_cutoff_temperature_K) < options_.critical_cap_tolerance_K) {
    return boundary_T - options_.critical_cap_extension_K;
  }

  double fraction = options_.metastable_stability_fraction;
  double D_boundary = StabilityFactor(rho, boundary_T);

  if (!std::isfinite(D_boundary) || !(D_boundary > 0.0)) return boundary_T;

  double target = fraction * D_boundary;
  double lower_limit = std::max(options_.minimum_extension_temperature_K,
                                boundary_T - options_.max_metastable_extension_K);

  double T_high = boundary_T;
  double D_high = D_boundary;

  double T_low = boundary_T;
  double D_low = D_boundary;

  bool found = false;
  while (T_low > lower_limit) {
    double next_T = std::max(lower_limit, T_low - options_.metastable_scan_step_K);
    double next_D = StabilityFactor(rho, next_T);

    if (!std::isfinite(next_D) || next_D <= target) {
      T_high = T_low;
      D_high = D_low;

      T_low = next_T;
      D_low = next_D;

      found = true;
      break;
    }

    T_low = next_T;
    D_low = next_D;

    if (T_low == lower_limit) break;
  }

  if (!found) return lower_limit;

  for (unsigned itr = 0; itr < options_.metastable_bisection_iterations; ++itr) {
    double T_mid = (T_low + T_high) / 2;
    double D_mid = StabilityFactor(rho, T_mid);

    if (!std::isfinite(D_mid) || D_mid <= target) {
      T_low = T_mid;
    } else {
      T_high = T_mid;
    }
  }

  return T_high;
}


/* ******************************************************************
* Evaluates factor propotional to (dp/drho)_T
****************************************************************** */
double
IAPWS95_RaggedSplineRhoT::StabilityFactor(double rho, double T)
{
  const auto& exact = IAPWS95::ResidualPart(rho, T);

  double delta = rho / RHOC;
  return 1.0 + 2.0 * delta * exact[1] + delta * delta * exact[3];
}


/* ******************************************************************
* Approximation of the saturation dome.
****************************************************************** */
void
IAPWS95_RaggedSplineRhoT::BuildSaturationData_()
{
  saturation_.reserve(mesh_.y_lines.size());
  saturation_.clear();
  for (double T : mesh_.y_lines) {
    if (T >= TC) break;

    double rhol0 = eos95_->DensityLiquid(T);
    double rhov0 = eos95_->DensityVapor(T);

    auto [rhol, rhov, psat] = eos95_->SaturationLineT(T, rhol0, rhov0);
    AMANZI_ASSERT(rhol > rhov && rhov > 0.0 && psat > 0.0);

    saturation_.push_back({T, rhol, rhov, psat});
  }
}


/* ******************************************************************
* Mesh adaptation
****************************************************************** */
void
IAPWS95_RaggedSplineRhoT::AdaptiveRefineCoordinateLines_()
{
  // A complete production version can perform solve-estimate-refine cycles.
  // Here we implement geometry-driven refinement before the first solve:
  // refine cells crossed by the saturation boundary and cells whose endpoint
  // exact derivatives vary more than the configured normalized tolerance.
  for (int pass = 0; pass < options_.max_adaptive_passes; ++pass) {
    bool changed = false;
    std::vector<double> add_rho;
    std::vector<double> add_T;

    if (mesh_.x_lines.size() - 1 < options_.max_rho_intervals) {
      for (int i = 0; i + 1 < mesh_.x_lines.size(); ++i) {
        double r0 = mesh_.x_lines[i];
        double r1 = mesh_.x_lines[i + 1];
        double rm = std::sqrt(r0 * r1);

        double T0 = BoundaryTemperature(r0);
        double T1 = BoundaryTemperature(r1);
        double curvature = std::fabs(T0 - 2 * BoundaryTemperature(rm) + T1);
        // double hT = MinSpacing(mesh_.y_lines);
        double hT = std::fabs(T1 - T0);
        if (hT > 1e-12 && curvature > 0.04 * hT) add_rho.push_back(rm);
      }
    }

    if (mesh_.y_lines.size() - 1 < options_.max_T_intervals) {
      for (int j = 0; j + 1 < mesh_.y_lines.size(); ++j) {
        double T0 = mesh_.y_lines[j];
        double T1 = mesh_.y_lines[j + 1];
        double Tm = (T0 + T1) / 2;

        bool crossed(false);
        for (double rho : mesh_.x_lines) {
          double F = BoundaryTemperature(rho);
          if (F > T0 && F < T1) {
            crossed = true;
            break;
          }
        }
        if (crossed) add_T.push_back(Tm);
      }
    }

    auto insert_unique = [&](std::vector<double>& lines,
                             std::vector<double>& additions,
                             int max_intervals) {
      std::sort(additions.begin(), additions.end());
      additions.erase(std::unique(additions.begin(), additions.end()), additions.end());
      for (double x : additions) {
        if (lines.size() - 1 >= max_intervals) break;
        const auto it = std::lower_bound(lines.begin(), lines.end(), x);
        if (it == lines.end() || std::abs(*it - x) > 1e-14 * std::max(1.0, std::abs(x))) {
          lines.insert(it, x);
          changed = true;
        }
      }
    };

    insert_unique(mesh_.x_lines, add_rho, options_.max_rho_intervals);
    insert_unique(mesh_.y_lines, add_T, options_.max_T_intervals);
    if (!changed) break;

    BuildSaturationData_();
    BuildRaggedColumns_();
  }

  // additional lines
  auto insert_unique = [&](std::vector<double>& lines,
                           const std::vector<double>& additions) {
    for (double x : additions) {
      const auto it = std::lower_bound(lines.begin(), lines.end(), x);
      if (it == lines.begin() || it == lines.end()) continue;
      if (std::abs(*it - x) > 1e-14 * std::max(1.0, std::abs(x))) lines.insert(it, x);
    }
  };

  insert_unique(mesh_.x_lines, { 320.0, 324.0, 1075.0 });
  insert_unique(mesh_.y_lines, { 283, 286, 645.0, 646.0, 647.0, 648.5, 650.5 });
  BuildSaturationData_();
  BuildRaggedColumns_();
}


/* ******************************************************************
* Refine based on error
****************************************************************** */
void
IAPWS95_RaggedSplineRhoT::AnisotropicRefinement()
{
  std::vector<CellError> cells;

  auto Error = [&](double rho, double T) {
    const auto& exact  = eos95_->ResidualPart(rho, T);
    const auto& spline = this->ResidualPart(rho, T);

    int id = LowerCell(mesh_.x_lines, rho);
    int it = LowerCell(mesh_.y_lines, T);

    double tau0 = TC / mesh_.y_lines[it];
    double tau1 = TC / mesh_.y_lines[it + 1];

    const double h_delta = (mesh_.x_lines[id + 1] - mesh_.x_lines[id]) / RHOC;
    const double h_tau = std::fabs(tau1 - tau0);
    const std::array<double, 6> derivative_scale = { 1.0,
                                                     h_delta,
                                                     h_tau,
                                                     h_delta * h_delta,
                                                     h_delta * h_tau,
                                                     h_tau * h_tau };

    double e = 0.0;
    for (int k = 0; k < 6; ++k) {
      double ek = std::sqrt(options_.fit_weights[k]) * derivative_scale[k] * std::fabs(spline[k] - exact[k]);
      e = std::max(e, ek);
    }
    return e;
  };

  for (int i = 0; i + 1 < mesh_.x_lines.size(); ++i) {
    double rho0 = mesh_.x_lines[i];
    double rho1 = mesh_.x_lines[i + 1];
    double rhom = std::sqrt(rho0 * rho1);

    // Quarter points in log(rho).
    double rhoL = std::sqrt(rho0 * rhom);
    double rhoR = std::sqrt(rhom * rho1);

    for (int j = 0; j + 1 < mesh_.y_lines.size(); ++j) {
      double T0 = mesh_.y_lines[j];
      double T1 = mesh_.y_lines[j + 1];
      double Tm = 0.5 * (T0 + T1);
      if (!IsExtended_(rhom, Tm)) continue;

      // Quarter points in T.
      double TL = 0.5 * (T0 + Tm);
      double TR = 0.5 * (Tm + T1);

      double error_rho = Error(rhom, Tm);
      double error_T = error_rho;

      if (IsExtended_(rhoL, Tm)) error_rho = std::max(error_rho, Error(rhoL, Tm));
      if (IsExtended_(rhoR, Tm)) error_rho = std::max(error_rho, Error(rhoR, Tm));

      if (IsExtended_(rhom, TL)) error_T = std::max(error_T, Error(rhom, TL));
      if (IsExtended_(rhom, TR)) error_T = std::max(error_T, Error(rhom, TR));

      cells.push_back({std::max(error_rho, error_T), error_rho, error_T, rhom, Tm});
    }
  }

  if (cells.empty()) return;

  // extract cells with the largest error
  std::vector<double> add_rho;
  std::vector<double> add_T;

  std::sort(cells.begin(),
            cells.end(),
            [](const CellError& a, const CellError& b) { return a.error > b.error; });

  int nbase = std::min(mesh_.x_lines.size(), mesh_.y_lines.size());
  int nrefine = std::max(1, (int)std::ceil(options_.refinement_fraction * nbase));

  double anisotropy_factor = 1.25;

  for (int n = 0; n < nrefine; ++n) {
    const CellError& cell = cells[n];

    if (cell.error_rho > anisotropy_factor * cell.error_T) {
      add_rho.push_back(cell.rho_mid);
    } else if (cell.error_T > anisotropy_factor * cell.error_rho) {
      add_T.push_back(cell.T_mid);
    } else {
      add_rho.push_back(cell.rho_mid);
      add_T.push_back(cell.T_mid);
    }
  }

  // remove duplicates
  std::sort(add_rho.begin(), add_rho.end());
  add_rho.erase(std::unique(add_rho.begin(), add_rho.end()), add_rho.end());

  std::sort(add_T.begin(), add_T.end());
  add_T.erase(std::unique(add_T.begin(), add_T.end()), add_T.end());

  // add to the mesh
  auto insert_unique = [&](std::vector<double>& lines,
                           const std::vector<double>& additions) {
    for (double x : additions) {
      const auto it = std::lower_bound(lines.begin(), lines.end(), x);
      if (it == lines.begin() || it == lines.end()) continue;
      if (std::abs(*it - x) > 1e-14 * std::max(1.0, std::abs(x))) lines.insert(it, x);
    }
  };

  insert_unique(mesh_.x_lines, add_rho);
  insert_unique(mesh_.y_lines, add_T);
  BuildSaturationData_();
  BuildRaggedColumns_();

  BuildKnotVectors_();
  min_T_spacing_ = MinSpacing(mesh_.y_lines);
}


/* ******************************************************************
*
****************************************************************** */
void
IAPWS95_RaggedSplineRhoT::BuildKnotVectors_()
{
  std::vector<double> delta_lines(mesh_.x_lines.size());
  std::transform(mesh_.x_lines.begin(), mesh_.x_lines.end(),
                 delta_lines.begin(),
                 [&](double rho) { return rho / RHOC; });

  // tau decreases as T increases. Build an increasing tau line array.
  std::vector<double> tau_lines(mesh_.y_lines.size());
  std::transform(mesh_.y_lines.begin(), mesh_.y_lines.end(),
                 tau_lines.begin(),
                 [&](double T) { return TC / T; });
  std::reverse(tau_lines.begin(), tau_lines.end());

  mesh_.x_knots = MakeClampedCubicKnots(delta_lines);
  mesh_.y_knots = MakeClampedCubicKnots(tau_lines);

  static constexpr int degree = 3;
  mesh_.nx_basis_ = mesh_.x_knots.size() - degree - 1;
  mesh_.ny_basis_ = mesh_.y_knots.size() - degree - 1;

  mesh_.x_span_data = BuildCubicSpanCache(mesh_.x_knots);
  mesh_.y_span_data = BuildCubicSpanCache(mesh_.y_knots);
}


/* ******************************************************************
* Creates the set of (rho, T) points used to fit spline coefficients
****************************************************************** */
std::vector<IAPWS95_RaggedSplineRhoT::Sample>
IAPWS95_RaggedSplineRhoT::BuildSamples()
{
  std::vector<Sample> samples;
  const unsigned q = options_.samples_per_cell_direction;

  // interior cell samples
  for (int i = 0; i + 1 < mesh_.x_lines.size(); ++i) {
    double r0 = mesh_.x_lines[i];
    double r1 = mesh_.x_lines[i + 1];
    for (int j = 0; j + 1 < mesh_.y_lines.size(); ++j) {
      double T0 = mesh_.y_lines[j];
      double T1 = mesh_.y_lines[j + 1];

      for (unsigned ir = 0; ir < q; ++ir) {
        double sr = (ir + 0.5) / q;
        double rho = (1.0 - sr) * r0 + sr * r1;
        for (unsigned jt = 0; jt < q; ++jt) {
          double st = (jt + 0.5) / q;
          double T = (1.0 - st) * T0 + st * T1;

          if (!IsExtended_(rho, T)) continue;
          double weight = IsPhysical(rho, T) ? 1.0 : options_.extension_weight;
          // double D = StabilityFactor(rho, T);
          // if (D >= 0.0) samples.push_back({rho, T, weight});
          samples.push_back({rho, T, weight});
        }
      }
    }
  }

  // Add all background-grid vertices. This directly constrains clamped
  // endpoint coefficients, especially corners such as (rho_max, T_min).
  for (double rho : mesh_.x_lines) {
    for (double T : mesh_.y_lines) {
      if (!IsExtended_(rho, T)) continue;
      double weight = IsPhysical(rho, T) ? 2.0 : options_.extension_weight;
      samples.push_back({rho, T, weight});
    }
  }

  // Additional smaples only under the flat critical cutoff.
  for (double rho : mesh_.x_lines) {
    double Tb = BoundaryTemperature(rho);
    if (Tb < options_.critical_cutoff_temperature_K) continue;

    double Te = Tb - options_.critical_cap_extension_K;

    for (int i = 0; i < options_.extension_samples; ++i) {
      double s = (i + 0.5) / options_.extension_samples;
      double T = Te + s * (Tb - Te);
      samples.push_back({rho, T, options_.extension_weight});
    }
  }

  return samples;
}


/* ******************************************************************
* Returns true when point is inside the metastable extension.
****************************************************************** */
bool
IAPWS95_RaggedSplineRhoT::IsExtended_(double rho, double T)
{
  double Tb = BoundaryTemperature(rho);
  double Te = FindMetastableMarginTemperature(rho, Tb);
  return T >= Te;
}

} // namespace AmanziEOS
} // namespace Amanzi

