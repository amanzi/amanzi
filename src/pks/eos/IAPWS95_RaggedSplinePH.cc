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

#include "Brent.hh"
#include "lapack.hh"

#include "IAPWS95_RaggedSplinePH.hh"

namespace Amanzi {
namespace AmanziEOS {

/* ******************************************************************
* Constructor
****************************************************************** */
IAPWS95_RaggedSplinePH::IAPWS95_RaggedSplinePH(Teuchos::ParameterList& plist, Options options)
  : IAPWS95(plist),
    options_(std::move(options))
{
  ValidateOptions_();
  eos95_ = std::make_shared<IAPWS95>(plist);
}


/* ******************************************************************
* Six derivatives of entropy 
****************************************************************** */
std::array<double, 6>
IAPWS95_RaggedSplinePH::EntropyDerivativesPH(double p, double h)
{ 
  residual_calls++;

  double x = std::log(p / PC);
  double theta = h / HC;
  auto d = Evaluate(mesh_, coefficients_, x, theta);

  double sp = d[1] / p;
  double sh = d[2] / HC;

  double spp = (d[3] - d[1]) / (p * p);
  double sph = d[4] / (p * HC);
  double shh = d[5] / (HC * HC);

  return { d[0], sp, sh, spp, sph, shh };
}


/* ******************************************************************
* Wrapper for initialization of shared data
****************************************************************** */
std::vector<IAPWS95_RaggedSplinePH::Sample>
IAPWS95_RaggedSplinePH::InitializeSharedData()
{
  std::vector<Sample> samples;

  std::call_once(init_shared_, [this, &samples]() {
    CreateRaggedMesh();
    samples = BuildSamples();
    BuildSplineCoefficients(samples);

    constexpr int lookup_bins = 512;
    mesh_.x_lookup = MakeSpanLookupLeft(mesh_.x_knots, lookup_bins);
    mesh_.y_lookup = MakeSpanLookupLeft(mesh_.y_knots, lookup_bins);

    AnisotropicRefinement();
    samples = BuildSamples();
    BuildSplineCoefficients(samples);

    mesh_.x_lookup = MakeSpanLookupLeft(mesh_.x_knots, lookup_bins);
    mesh_.y_lookup = MakeSpanLookupLeft(mesh_.y_knots, lookup_bins);
  });

  return samples;
}


/* ******************************************************************
* Create a ragged mesh
****************************************************************** */
void
IAPWS95_RaggedSplinePH::CreateRaggedMesh()
{
  mesh_ = Mesh{};
  coefficients_.clear();

  BuildInitialCoordinateLines_();
  BuildRaggedColumns_();
  AdaptiveRefineCoordinateLines_();
  BuildKnotVectors_();
}


/* ******************************************************************
* Assemble a global overdetermined least-squares system for the active
* cubic tensor-product basis functions and solve for all spline coefficients.
****************************************************************** */
void
IAPWS95_RaggedSplinePH::BuildSplineCoefficients(const std::vector<Sample>& samples)
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
    double x = std::log(sample.p / PC);
    double theta = sample.h / HC;
    const BasisData Bp = EvaluateCubicBasis(mesh_.x_knots, x);
    const BasisData Bt = EvaluateCubicBasis(mesh_.y_knots, theta);

    const double p = sample.p;
    const auto& der = eos95_->EntropyDerivativesPHbase(sample.rho, sample.T);
    const std::array<double, 6> target = { der[0],
                                           p * der[1],
                                           HC * der[2],
                                           p * der[1] + p * p * der[3],
                                           p * HC * der[4],
                                           HC * HC * der[5] };

    int id = LowerCell(mesh_.x_lines, sample.p);
    int it = LowerCell(mesh_.y_lines, sample.h);

    const double h_x = std::log(mesh_.x_lines[id + 1] / mesh_.x_lines[id]);
    const double h_theta = (mesh_.y_lines[it + 1] - mesh_.y_lines[it]) / HC;

    const std::array<double, 6> derivative_scale = { 1.0,
                                                     h_x,
                                                     h_theta,
                                                     h_x * h_x,
                                                     h_x * h_theta,
                                                     h_theta * h_theta };

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
              xd = Bp.value[a];
              xt = Bt.value[b];
              break;
            case 1:
              xd = Bp.d1[a];
              xt = Bt.value[b];
              break;
            case 2:
              xd = Bp.value[a];
              xt = Bt.d1[b];
              break;
            case 3:
              xd = Bp.d2[a];
              xt = Bt.value[b];
              break;
            case 4:
              xd = Bp.d1[a];
              xt = Bt.d1[b];
              break;
            case 5:
              xd = Bp.value[a];
              xt = Bt.d2[b];
              break;
            default:
              ;  // pass
          }

          index[q] = mesh_.CoefficientIndex(Bp.index[a], Bt.index[b]);
          row_value[q] = xd * xt;
          ++q;
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
  double diagonal_mean = diagonal_sum / static_cast<double>(n);
  double regularization = 1.0e-12 * std::max(1.0, diagonal_mean);
  for (int j = 0; j < n; ++j) {
    ab[kd + j * ldab] += regularization;
  }

  int nrhs(1), info(0);
  WhetStone::DPBSV_F77("U", &n, &kd, &nrhs, ab.data(), &ldab, rhs.data(), &n, &info);
  if (info != 0) AMANZI_ASSERT(false);
  // info < 0 : illegal argument 
  // info > 0 : normal matrix is not positive definite at pivot info

  coefficients_ = std::move(rhs);

  constexpr int lookup_bins = 512;
  mesh_.x_lookup = MakeSpanLookupLeft(mesh_.x_knots, lookup_bins);
  mesh_.y_lookup = MakeSpanLookupLeft(mesh_.y_knots, lookup_bins);
}


/* ******************************************************************
* Build initial mesh
****************************************************************** */
void
IAPWS95_RaggedSplinePH::BuildInitialCoordinateLines_()
{
  mesh_.x_lines.resize(options_.initial_p_intervals + 1);
  mesh_.y_lines.resize(options_.initial_h_intervals + 1);

  // Uniform initial spacing in the spline coordinate x = log(p / PC).
  // x_lines themselves remain stored as physical pressures.
  double log_min = std::log(options_.p_min);
  double log_max = std::log(options_.p_max);
  for (int i = 0; i < mesh_.x_lines.size(); ++i) {
    double s = (double)i / (mesh_.x_lines.size() - 1);
    mesh_.x_lines[i] = std::exp((1.0 - s) * log_min + s * log_max);
  }

  for (int j = 0; j < mesh_.y_lines.size(); ++j) {
    double s = (double)j / (mesh_.y_lines.size() - 1);
    mesh_.y_lines[j] = (1.0 - s) * options_.h_min + s * options_.h_max;
  }
}


/* ******************************************************************
* Should be called every time the coordinate lines change, e.g.
* during adaptive refinement. Will be used for boundary interpolation.
****************************************************************** */
void
IAPWS95_RaggedSplinePH::BuildRaggedColumns_()
{
  columns_.clear();
  columns_.reserve(mesh_.x_lines.size());

  for (double p : mesh_.x_lines) {
    RaggedColumn c;
    c.p = p;

    if (p < PC) {
      auto sat = eos95_->SaturationLineP(p);

      auto liquid = PopulateProperties(sat.rhol, sat.Tsat);
      auto vapor = PopulateProperties(sat.rhov, sat.Tsat);

      double ext_hl = FindLiquidMetastableMargin(p, sat);
      double ext_hv = FindVaporMetastableMargin(p, sat);

      c.physical_intervals = { {options_.h_min, liquid.h, 1.0},
                               {vapor.h, options_.h_max, 1.0} };
      c.extension_intervals = { {liquid.h, ext_hl, options_.extension_weight},
                                {ext_hv, vapor.h, options_.extension_weight} };
    } else {
      c.physical_intervals = { {options_.h_min, options_.h_max, 1.0} };
      c.extension_intervals = { {options_.h_min, options_.h_max, 1.0} };
    }

    columns_.push_back(c);
  }
}


/* ******************************************************************
* Validate options
****************************************************************** */
void IAPWS95_RaggedSplinePH::ValidateOptions_() const
{
  AMANZI_ASSERT(options_.p_min > 0.0 && options_.p_max > options_.p_min);
  AMANZI_ASSERT(options_.h_max > options_.h_min && options_.h_min > 0.0);

  // At least four intervals per direction are required.
  AMANZI_ASSERT(options_.initial_p_intervals >= 4 && options_.initial_h_intervals >= 4);

  // Samples_per_cell_direction must be positive.
  AMANZI_ASSERT(options_.samples_per_cell_direction != 0);
}


/* ******************************************************************
*
****************************************************************** */
std::pair<double, double>
IAPWS95_RaggedSplinePH::BoundaryEnthalpies(double p) const
{
  if (p > PC) return { options_.h_min, options_.h_min };

  const auto& sat = eos95_->SaturationLineP(p);
  return { sat.hl, sat.hv };
}


/* ******************************************************************
* Metastable-margin function searches downward from two-phase boundary
****************************************************************** */
struct Frho95s {
  Frho95s(double p, double rho, IAPWS95* eos) 
    : p_(p), rho_(rho), Rrho_(IAPWS95::R * rho_ / 1000.0), eos_(eos) {};
  double operator()(double T) const {
    double delta = rho_ / IAPWS95::RHOC;

    double gd = eos_->ResidualPart(rho_, T)[1];
    double po = (1.0 + delta * gd) * Rrho_ * T;
    return po - p_;
  }

  double p_, rho_, Rrho_;
  IAPWS95* eos_;
};


double
IAPWS95_RaggedSplinePH::FindLiquidMetastableMargin(double p,
                                                   const SaturationState& sat)
{
  double T0 = sat.Tsat;
  double rho0 = sat.rhol;
  double rho_mid = (rho0 + sat.rhov) / 2; 

  double D0 = StabilityFactor(rho0, T0);
  double target = options_.metastable_stability_fraction * D0;
  double fraction_step = options_.metastable_density_fraction_step;

  for (;;) {
    double rho1 = rho0 * (1.0 - fraction_step);

    int itrs = 20;
    double tol = 1e-8;
    Frho95s f(p, rho1, eos95_.get());
    auto [Tmin, Tmax] = Utils::bracketRootSymmetric(f, T0, T0 * 0.01, &itrs);
    brent_bracket_itrs += itrs;
    if (itrs < 0) return sat.hl;

    itrs = 20;
    double T1 = Utils::findRootBrent(f, Tmin, Tmax, tol, &itrs);
    AMANZI_ASSERT(itrs >= 0);
    brent_root_itrs += itrs;

    double D1 = StabilityFactor(rho1, T1);
    if (D1 <= 0.0) {
      fraction_step *= 0.5;
    } else if (D1 <= target || rho1 < rho_mid) {
      auto prop = eos95_->PopulateProperties(rho1, T1);
      return prop.h;
    } else {
      T0 = T1;
      rho0 = rho1;
    }
  }

  return sat.hl;
}


double
IAPWS95_RaggedSplinePH::FindVaporMetastableMargin(double p,
                                                  const SaturationState& sat)
{
  double T0 = sat.Tsat;
  double rho0 = sat.rhov;
  double rho_mid = (rho0 + sat.rhol) / 2; 

  double D0 = StabilityFactor(rho0, T0);
  double target = options_.metastable_stability_fraction * D0;
  double fraction_step = 2 * options_.metastable_density_fraction_step;

  for (;;) {
    double rho1 = rho0 * (1.0 + fraction_step);

    int itrs = 20;
    double tol = 1e-8;
    Frho95s f(p, rho1, eos95_.get());
    auto [Tmin, Tmax] = Utils::bracketRootSymmetric(f, T0, T0 * 0.01, &itrs);
    brent_bracket_itrs += itrs;
    if (itrs < 0) return sat.hv;

    itrs = 20;
    double T1 = Utils::findRootBrent(f, Tmin, Tmax, tol, &itrs);
    AMANZI_ASSERT(itrs >= 0);
    brent_root_itrs += itrs;

    double D1 = StabilityFactor(rho1, T1);
    if (D1 < 0.0) {
      fraction_step *= 0.5;
    } else if (D1 <= target || rho1 > rho_mid) {
      auto prop = eos95_->PopulateProperties(rho1, T1);
      return prop.h;
    } else {
      T0 = T1;
      rho0 = rho1;
    }
  }

  return sat.hv;
}


/* ******************************************************************
* Evaluates factor propotional to (dp/drho)_T
****************************************************************** */
double
IAPWS95_RaggedSplinePH::StabilityFactor(double rho, double T)
{
  const auto& exact = IAPWS95::ResidualPart(rho, T);

  double delta = rho / RHOC;
  return 1.0 + 2.0 * delta * exact[1] + delta * delta * exact[3];
}


/* ******************************************************************
* Mesh adaptation
****************************************************************** */
void
IAPWS95_RaggedSplinePH::AdaptiveRefineCoordinateLines_()
{
  // A complete production version can perform solve-estimate-refine cycles.
  // Here we implement geometry-driven refinement before the first solve:
  // refine cells crossed by the saturation boundary and cells whose endpoint
  // exact derivatives vary more than the configured normalized tolerance.
  for (int pass = 0; pass < options_.max_adaptive_passes; ++pass) {
    bool changed = false;
    std::vector<double> add_p;
    std::vector<double> add_h;

    double min_h_spacing = MinSpacing(mesh_.y_lines);

    if (mesh_.x_lines.size() - 1 < options_.max_p_intervals) {
      for (int i = 0; i + 1 < mesh_.x_lines.size(); ++i) {
        double p0 = mesh_.x_lines[i];
        double p1 = mesh_.x_lines[i + 1];
        double pm = std::sqrt(p0 * p1);

        auto b0 = BoundaryEnthalpies(p0);
        auto b1 = BoundaryEnthalpies(p1);
        auto bm = BoundaryEnthalpies(pm);

        double curvature_l = std::fabs(b0.first - 2 * bm.first + b1.first);
        double curvature_v = std::fabs(b0.second - 2 * bm.second + b1.second);
        double hb1 = std::fabs(b1.first - b0.first);
        double hb2 = std::fabs(b1.second - b0.second);
        if (hb1 > 1e-12 && curvature_l > 0.2 * hb1 ||
            hb2 > 1e-12 && curvature_v > 0.2 * hb2) add_p.push_back(pm);
      }
    }

    if (mesh_.y_lines.size() - 1 < options_.max_h_intervals) {
      for (int j = 0; j + 1 < mesh_.y_lines.size(); ++j) {
        double h0 = mesh_.y_lines[j];
        double h1 = mesh_.y_lines[j + 1];
        double hm = 0.5 * (h0 + h1);

        bool crossed = false;
        for (double p : mesh_.x_lines) {
          auto F = BoundaryEnthalpies(p);
          if (F.first > h0 && F.first < h1) {
            crossed = true;
            break;
          }
          if (F.second > h0 && F.second < h1) {
            crossed = true;
            break;
          }
        }
        if (crossed) add_h.push_back(hm);
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

    insert_unique(mesh_.x_lines, add_p, options_.max_p_intervals);
    insert_unique(mesh_.y_lines, add_h, options_.max_h_intervals);
    if (!changed) break;

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

  insert_unique(mesh_.x_lines, { 20.4, 20.8, 21.2 });
  insert_unique(mesh_.y_lines, { 505.0, 2092.0 });
  // insert_unique(mesh_.y_lines, { 505.0 });
  BuildRaggedColumns_();
}


/* ******************************************************************
* Refine based on error
****************************************************************** */
void
IAPWS95_RaggedSplinePH::AnisotropicRefinement()
{
  std::vector<CellError> cells;

  auto Transform = [&](const std::array<double, 6>& d, double p) {
    return std::array<double, 6>{ d[0],
                                  p * d[1],
                                  HC * d[2],
                                  p * d[1] + p * p * d[3],
                                  p * HC * d[4],
                                  HC * HC * d[5] };
  };

  auto Error = [&](double p, double h) {
    const auto exact_ph = eos95_->EntropyDerivativesPH(p, h);
    const auto spline_ph = this->EntropyDerivativesPH(p, h);

    const auto exact = Transform(exact_ph, p);
    const auto spline = Transform(spline_ph, p);

    int id = LowerCell(mesh_.x_lines, p);
    int it = LowerCell(mesh_.y_lines, h);

    const double h_x = std::log(mesh_.x_lines[id + 1] / mesh_.x_lines[id]);
    const double h_theta = (mesh_.y_lines[it + 1] - mesh_.y_lines[it]) / HC;

    const std::array<double, 6> derivative_scale = { 1.0,
                                                     h_x,
                                                     h_theta,
                                                     h_x * h_x,
                                                     h_x * h_theta,
                                                     h_theta * h_theta };

    double e = 0.0;
    for (int k = 0; k < 6; ++k) {
      double ek = std::sqrt(options_.fit_weights[k]) * derivative_scale[k] * std::fabs(spline[k] - exact[k]);
      e = std::max(e, ek);
    }

    return e;
  };

  for (int i = 0; i + 1 < mesh_.x_lines.size(); ++i) {
    double p0 = mesh_.x_lines[i];
    double p1 = mesh_.x_lines[i + 1];
    double pm = std::sqrt(p0 * p1);

    // Quarter points in log(p).
    double pL = std::sqrt(p0 * pm);
    double pR = std::sqrt(pm * p1);

    for (int j = 0; j + 1 < mesh_.y_lines.size(); ++j) {
      double h0 = mesh_.y_lines[j];
      double h1 = mesh_.y_lines[j + 1];
      double hm = 0.5 * (h0 + h1);

      const SaturationState& sat = eos95_->SaturationLineP(pm);
      if (!IsPhysical(pm, hm, sat)) continue;

      // Quarter points in h
      double hL = 0.5 * (h0 + hm);
      double hR = 0.5 * (hm + h1);

      double error_p = Error(pm, hm);
      double error_h = error_p;

      const SaturationState& satL = eos95_->SaturationLineP(pL);
      const SaturationState& satR = eos95_->SaturationLineP(pR);
      if (IsPhysical(pL, hm, satL)) error_p = std::max(error_p, Error(pL, hm));
      if (IsPhysical(pR, hm, satR)) error_p = std::max(error_p, Error(pR, hm));

      if (IsPhysical(pm, hL, sat)) error_h = std::max(error_h, Error(pm, hL));
      if (IsPhysical(pm, hR, sat)) error_h = std::max(error_h, Error(pm, hR));

      cells.push_back({std::max(error_p, error_h), error_p, error_h, pm, hm});
    }
  }

  if (cells.empty()) return;

  // extract cells with the largest error
  std::vector<double> add_p;
  std::vector<double> add_h;

  std::sort(cells.begin(),
            cells.end(),
            [](const CellError& a, const CellError& b) { return a.error > b.error; });

  int nrefine = std::max(1, (int)std::ceil(options_.anisotropic_refinement_fraction * cells.size()));

  double anisotropy_factor = 1.25;

  for (int n = 0; n < nrefine; ++n) {
    const CellError& cell = cells[n];

    if (cell.error_p > anisotropy_factor * cell.error_h) {
      add_p.push_back(cell.p_mid);
    } else if (cell.error_h > anisotropy_factor * cell.error_p) {
      add_h.push_back(cell.h_mid);
    } else {
      add_p.push_back(cell.p_mid);
      add_h.push_back(cell.h_mid);
    }
  }

  // remove duplicates
  std::sort(add_p.begin(), add_p.end());
  add_p.erase(std::unique(add_p.begin(), add_p.end()), add_p.end());

  std::sort(add_h.begin(), add_h.end());
  add_h.erase(std::unique(add_h.begin(), add_h.end()), add_h.end());

  // add to the mesh
  auto insert_unique = [&](std::vector<double>& lines,
                           const std::vector<double>& additions) {
    for (double x : additions) {
      const auto it = std::lower_bound(lines.begin(), lines.end(), x);
      if (it == lines.begin() || it == lines.end()) continue;
      if (std::abs(*it - x) > 1e-14 * std::max(1.0, std::abs(x))) lines.insert(it, x);
    }
  };

  insert_unique(mesh_.x_lines, add_p);
  insert_unique(mesh_.y_lines, add_h);

  BuildRaggedColumns_();
  BuildKnotVectors_();
}


/* ******************************************************************
*
****************************************************************** */
void
IAPWS95_RaggedSplinePH::BuildKnotVectors_()
{
  std::vector<double> x_lines(mesh_.x_lines.size());
  std::transform(mesh_.x_lines.begin(), mesh_.x_lines.end(),
                 x_lines.begin(),
                 [&](double p) { return std::log(p / PC); });

  std::vector<double> theta_lines(mesh_.y_lines.size());
  std::transform(mesh_.y_lines.begin(), mesh_.y_lines.end(),
                 theta_lines.begin(),
                 [&](double h) { return h / HC; });

  mesh_.x_knots = MakeClampedCubicKnots(x_lines);
  mesh_.y_knots = MakeClampedCubicKnots(theta_lines);

  static constexpr int degree = 3;
  mesh_.nx_basis_ = mesh_.x_knots.size() - degree - 1;
  mesh_.ny_basis_ = mesh_.y_knots.size() - degree - 1;

  mesh_.x_span_data = BuildCubicSpanCache(mesh_.x_knots);
  mesh_.y_span_data = BuildCubicSpanCache(mesh_.y_knots);
}


/* ******************************************************************
* Creates the set of (p, h) points used to fit spline coefficients
****************************************************************** */
std::vector<IAPWS95_RaggedSplinePH::Sample>
IAPWS95_RaggedSplinePH::BuildSamples()
{
  std::vector<Sample> samples;
  const unsigned q = options_.samples_per_cell_direction;

  // interior cell samples
  for (int i = 0; i + 1 < mesh_.x_lines.size(); ++i) {
    double p0 = mesh_.x_lines[i];
    double p1 = mesh_.x_lines[i + 1];
    for (int j = 0; j + 1 < mesh_.y_lines.size(); ++j) {
      double h0 = mesh_.y_lines[j];
      double h1 = mesh_.y_lines[j + 1];

      for (unsigned ir = 0; ir < q; ++ir) {
        double sr = (ir + 0.5) / q;
        double p = p0 * std::pow(p1 / p0, sr);
        for (unsigned jt = 0; jt < q; ++jt) {
          double st = (jt + 0.5) / q;
          double h = (1.0 - st) * h0 + st * h1;

          const SaturationState& sat = eos95_->SaturationLineP(p);
          if (!IsExtended_(p, h, sat)) continue;

          bool flag = IsPhysical(p, h, sat);
          double weight = flag ? 1.0 : options_.extension_weight;
          samples.push_back({p, h, 0.0, 0.0, weight, flag});
        }
      }
    }
  }

  // Add all background-grid vertices. This directly constrains clamped
  // endpoint coefficients, especially corners such as (p_max, h_min).
  for (double p : mesh_.x_lines) {
    for (double h : mesh_.y_lines) {
      const SaturationState& sat = eos95_->SaturationLineP(p);
      if (!IsExtended_(p, h, sat)) continue;

      bool flag = IsPhysical(p, h, sat);
      double weight = flag ? 2.0 : options_.extension_weight;
      samples.push_back({p, h, 0.0, 0.0, weight, flag});
    }
  }

  // Additional samples only under the flat critical cutoff.
  // ... pass

  // pass samples in physical region, populate rho/T data
  for (Sample& sample : samples) {
    if (sample.is_physical) {
      auto [prop, liquid, vapor] = eos97_.ThermodynamicsPH(sample.p, sample.h); 
      sample.rho = prop.rho;
      sample.T = prop.T;
    }
  }

  // pass samples in metastable region
  std::vector<MetastableWork> liquid_work;
  std::vector<MetastableWork> vapor_work;

  BuildMetastableWorkLists_(samples, liquid_work, vapor_work);

  PopulateMetastableSamples_(samples, liquid_work);
  PopulateMetastableSamples_(samples, vapor_work);

  return samples;
}


/* ******************************************************************
* Construct liquid and vapor metastable work lists in one pass.
****************************************************************** */
void
IAPWS95_RaggedSplinePH::BuildMetastableWorkLists_(const std::vector<Sample>& samples,
                                                  std::vector<MetastableWork>& liquid,
                                                  std::vector<MetastableWork>& vapor)
{
  liquid.clear();
  vapor.clear();

  liquid.reserve(samples.size());
  vapor.reserve(samples.size());

  for (int k = 0; k < static_cast<int>(samples.size()); ++k) {
    const Sample& sample = samples[k];
    if (sample.is_physical) continue;

    const SaturationState& sat = eos95_->SaturationLineP(sample.p);

    double dl = sample.h - sat.hl;
    double dv = sat.hv - sample.h;
    double latent_h = sat.hv - sat.hl;

    if (dl <= dv) {
      liquid.push_back({k, dl / latent_h, sat.rhol, sat.Tsat});
    } else {
      vapor.push_back({k, dv / latent_h, sat.rhov, sat.Tsat});
    }
  }

  // Solve closest-to-saturation samples first.
  auto compare = [](const MetastableWork& a, const MetastableWork& b) { return a.depth < b.depth; };

  std::sort(liquid.begin(), liquid.end(), compare);
  std::sort(vapor.begin(), vapor.end(), compare);
}


/* ******************************************************************
* Populate one metastable branch.
* The work list already contains the branch-specific saturation state.
****************************************************************** */
void
IAPWS95_RaggedSplinePH::PopulateMetastableSamples_(std::vector<Sample>& samples,
                                                   const std::vector<MetastableWork>& work)
{
  std::vector<Sample> solved;
  solved.reserve(work.size());

  for (const MetastableWork& w : work) {
    Sample& sample = samples[w.index];
    bool success = SolveMetastableSample_(sample, solved, w.rho_sat, w.T_sat);
    AMANZI_ASSERT(success);

    solved.push_back({sample.p, sample.h, sample.rho, sample.T, 1.0, false});
  }
}


/* ******************************************************************
* Solve one metastable (p,h) state. First use the nearest previously 
* converged state on this branch. If fails, restart from the 
* saturation state at the current pressure.
****************************************************************** */
bool
IAPWS95_RaggedSplinePH::SolveMetastableSample_(Sample& sample,
                                               const std::vector<Sample>& solved,
                                               double rho_sat,
                                               double T_sat)
{
  double p = sample.p;
  double h = sample.h;

  constexpr double tol = 1.0e-10;
  FrhoT f(p, h, eos95_.get());
  FrhoT::Vector sol(2), x0(2);

  // Attempt 1: nearest previously converged metastable state.
  int nearest = FindNearestSolvedState_(solved, p, h);

  if (nearest >= 0) {
    x0[0] = solved[nearest].rho;
    x0[1] = solved[nearest].T;

    int itrs = 100;

    sol = PowellHybrid(x0, f, &itrs, tol);
    if (itrs >= 0) {
      double D = StabilityFactor(sol[0], sol[1]);
      if (D > 0.0) {
        sample.rho = sol[0];
        sample.T = sol[1];
        return true;
      }
    }
  }

  // Attempt 2: saturation state at the current pressure.
  x0[0] = rho_sat;
  x0[1] = T_sat;

  int itrs = 100;

  sol = PowellHybrid(x0, f, &itrs, tol);
  if (itrs >= 0) {
    double D = StabilityFactor(sol[0], sol[1]);
    if (D > 0.0) {
      sample.rho = sol[0];
      sample.T = sol[1];
      return true;
    }
  }

  return false;
}


/* ******************************************************************
* Find the nearest previously converged state on the same branch.
****************************************************************** */
int
IAPWS95_RaggedSplinePH::FindNearestSolvedState_(const std::vector<Sample>& solved,
                                                double p,
                                                double h) const
{
  int nearest = -1;
  double dist_min = std::numeric_limits<double>::max();

  for (int k = 0; k < solved.size(); ++k) {
    // Logarithmic pressure distance is more appropriate over the
    // large pressure range of the table.
    double dp = std::log(p / solved[k].p);
    double dh = (h - solved[k].h) / HC;

    double dist = dp * dp + dh * dh;
    if (dist < dist_min) {
      dist_min = dist;
      nearest = k;
    }
  }

  return nearest;
}


/* ******************************************************************
*
****************************************************************** */
bool
IAPWS95_RaggedSplinePH::IsPhysical(double p, double h, const SaturationState& sat)
{ 
  if (p >= eos95_->PC) return true;
  return h <= sat.hl || h >= sat.hv;
}


/* ******************************************************************
*
****************************************************************** */
bool
IAPWS95_RaggedSplinePH::IsExtended_(double p, double h, const SaturationState& sat)
{
  if (p >= eos95_->PC) return true;
  if (h <= sat.hl || h >= sat.hv) return true;

  double ext_hl = FindLiquidMetastableMargin(p, sat);
  double ext_hv = FindVaporMetastableMargin(p, sat);
  return h <= ext_hl || h >= ext_hv;
}

} // namespace AmanziEOS
} // namespace Amanzi

