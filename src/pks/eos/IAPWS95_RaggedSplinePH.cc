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
  AMANZI_ASSERT(built_);

  // check for spline domain
  AMANZI_ASSERT(p >= options_.P_min && p <= options_.P_max);
  AMANZI_ASSERT(h >= options_.H_min && h <= options_.H_max);

  const SaturationState& sat = eos95_->SaturationLineP(p);
  AMANZI_ASSERT(IsPhysical_(p, h, sat));

  return Evaluate(mesh_, coefficients_, p / PC, h / HC);
}


/* ******************************************************************
* Create a ragged mesh
****************************************************************** */
void
IAPWS95_RaggedSplinePH::CreateRaggedMesh()
{
  mesh_ = Mesh{};
  coefficients_.clear();
  built_ = false;

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
    double pi = sample.p / PC;
    double theta = sample.h / HC;
    const BasisData Bp = EvaluateCubicBasis(mesh_.x_knots, pi);
    const BasisData Bt = EvaluateCubicBasis(mesh_.y_knots, theta);
    const std::array<double, 6> target = eos95_->EntropyDerivativesRhoT(sample.rho, sample.T);

    int id = LowerCell(mesh_.x_lines, sample.p);
    int it = LowerCell(mesh_.y_lines, sample.h);

    const double h_pi = (mesh_.x_lines[id + 1] - mesh_.x_lines[id]) / PC;
    const double h_theta = (mesh_.y_lines[it + 1] - mesh_.y_lines[it]) / HC;

    const std::array<double, 3> derivative_scale = { 1.0, h_pi, h_theta };

    for (int component = 0; component < 3; ++component) {
      double s = derivative_scale[component];
      double weight = sample.weight * options_.fit_weights[component] * s * s;
      if (!(weight > 0.0)) continue;

      int q = 0;
      for (int a = 0; a < 4; ++a) {
        const double xd = (component == 0 || component == 2) ? 
                          Bp.value[a] : (component == 1) ? Bp.d1[a] : Bp.d2[a];

        for (int b = 0; b < 4; ++b) {
          const double xt = (component == 0 || component == 1) ? 
                            Bt.value[b] : (component == 2) ? Bt.d1[b] : Bt.d2[b];

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
  built_ = true;
}


/* ******************************************************************
* Build initial mesh
****************************************************************** */
void
IAPWS95_RaggedSplinePH::BuildInitialCoordinateLines_()
{
  mesh_.x_lines.resize(options_.initial_p_intervals + 1);
  mesh_.y_lines.resize(options_.initial_h_intervals + 1);

  // Logarithmic density distribution handles the dilute-vapor scale while
  // remaining monotone and simple. We can replace by another distribution
  // if the liquid region requires stronger clustering.
  double log_min = std::log(options_.P_min);
  double log_max = std::log(options_.P_max);
  for (int i = 0; i < mesh_.x_lines.size(); ++i) {
    double s = (double)i / (mesh_.x_lines.size() - 1);
    mesh_.x_lines[i] = std::exp((1.0 - s) * log_min + s * log_max);
  }

  for (int j = 0; j < mesh_.y_lines.size(); ++j) {
    double s = (double)j / (mesh_.y_lines.size() - 1);
    mesh_.y_lines[j] = (1.0 - s) * options_.H_min + s * options_.H_max;
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

      double ext_hl = FindLiquidMetastableMargin_(p, sat);
      double ext_hv = FindVaporMetastableMargin_(p, sat);

      c.physical_intervals = { {options_.H_min, liquid.h, 1.0},
                               {vapor.h, options_.H_max, 1.0} };
      c.extension_intervals = { {liquid.h, ext_hl, options_.extension_weight},
                                {ext_hv, vapor.h, options_.extension_weight} };
    } else {
      c.physical_intervals = { {options_.H_min, options_.H_max, 1.0} };
      c.extension_intervals = { {options_.H_min, options_.H_max, 1.0} };
    }

    columns_.push_back(c);
  }
}


/* ******************************************************************
* Validate options
****************************************************************** */
void IAPWS95_RaggedSplinePH::ValidateOptions_() const
{
  AMANZI_ASSERT(options_.P_min > 0.0 && options_.P_max > options_.P_min);
  AMANZI_ASSERT(options_.H_max > options_.H_min && options_.H_min > 0.0);

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
  if (p > PC) return { options_.H_min, options_.H_min };

  const auto& sat = eos95_->SaturationLineP(p);
  return { sat.hl, sat.hv };
}


/* ******************************************************************
* Metastable-margin function searches downward from two-phase boundary
****************************************************************** */
struct Frho95s {
  Frho95s(double p, double rho, IAPWS95* eos) : p_(p), rho_(rho), eos_(eos) {};
  double operator()(double T) const {
    double delta = rho_ / eos_->RHOC;

    double gd = eos_->ResidualPart(rho_, T)[1];
    double po = (1.0 + delta * gd) * eos_->R * T * rho_ / 1000.0;
    return po - p_;
  }

  double p_, rho_;
  IAPWS95* eos_;
};


double
IAPWS95_RaggedSplinePH::FindLiquidMetastableMargin_(double p,
                                                    const SaturationState& sat)
{
  double T0 = sat.Tsat;
  double rho0 = sat.rhol;
  double rho_mid = (rho0 + sat.rhov) / 2; 

  double D0 = StabilityFactor(rho0, T0);
  double target = options_.metastable_stability_fraction * D0;

  for (;;) {
    double rho1 = rho0 - options_.metastable_density_fraction_step * rho0;

    itrs_ = 20;
    double tol = 1e-8;
    Frho95s f(p, rho1, eos95_.get());
    auto [Tmin, Tmax] = Utils::bracketRootSymmetric(f, T0, T0 * 0.01, &itrs_);
    if (itrs_ < 0) return sat.hl;

    itrs_ = 20;
    double T1 = Utils::findRootBrent(f, Tmin, Tmax, tol, &itrs_);
    AMANZI_ASSERT(itrs_ >= 0);

    double D1 = StabilityFactor(rho1, T1);
    if (D1 <= target || rho1 < rho_mid) {
      auto prop = eos95_->PopulateProperties(rho1, T1);
      return prop.h;
    }

    T0 = T1;
    rho0 = rho1;
  }

  return sat.hl;
}


double
IAPWS95_RaggedSplinePH::FindVaporMetastableMargin_(double p,
                                                   const SaturationState& sat)
{
  double T0 = sat.Tsat;
  double rho0 = sat.rhov;
  double rho_mid = (rho0 + sat.rhol) / 2; 

  double D0 = StabilityFactor(rho0, T0);
  double target = options_.metastable_stability_fraction * D0;

  for (;;) {
    double rho1 = rho0 + options_.metastable_density_fraction_step * rho0;

    itrs_ = 20;
    double tol = 1e-8;
    Frho95s f(p, rho1, eos95_.get());
    auto [Tmin, Tmax] = Utils::bracketRootSymmetric(f, T0, T0 * 0.01, &itrs_);
    if (itrs_ < 0) return sat.hv;

    itrs_ = 20;
    double T1 = Utils::findRootBrent(f, Tmin, Tmax, tol, &itrs_);
    AMANZI_ASSERT(itrs_ >= 0);

    double D1 = StabilityFactor(rho1, T1);
    if (D1 <= target || rho1 > rho_mid) {
      auto prop = eos95_->PopulateProperties(rho1, T1);
      return prop.h;
    }

    T0 = T1;
    rho0 = rho1;
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
* Minimum dT
****************************************************************** */
double
IAPWS95_RaggedSplinePH::MinEnthalpySpacing(double p) const
{
  double h = std::numeric_limits<double>::max();
  for (int j = 0; j + 1 < mesh_.y_lines.size(); ++j) {
    h = std::min(h, mesh_.y_lines[j + 1] - mesh_.y_lines[j]);
  }
  return h;
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
        double hh = MinEnthalpySpacing(pm);
        if (std::max(curvature_l, curvature_v) > 0.25 * hh) add_p.push_back(pm);
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
}


/* ******************************************************************
*
****************************************************************** */
void
IAPWS95_RaggedSplinePH::BuildKnotVectors_()
{
  std::vector<double> pi_lines(mesh_.x_lines.size());
  std::transform(mesh_.x_lines.begin(), mesh_.x_lines.end(),
                 pi_lines.begin(),
                 [&](double p) { return p / PC; });

  std::vector<double> theta_lines(mesh_.y_lines.size());
  std::transform(mesh_.y_lines.begin(), mesh_.y_lines.end(),
                 theta_lines.begin(),
                 [&](double h) { return h / HC; });

  mesh_.x_knots = MakeClampedCubicKnots(pi_lines);
  mesh_.y_knots = MakeClampedCubicKnots(theta_lines);

  static constexpr int degree = 3;
  mesh_.nx_basis_ = mesh_.x_knots.size() - degree - 1;
  mesh_.ny_basis_ = mesh_.y_knots.size() - degree - 1;
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
        double p = (1.0 - sr) * p0 + sr * p1;
        for (unsigned jt = 0; jt < q; ++jt) {
          double st = (jt + 0.5) / q;
          double h = (1.0 - st) * h0 + st * h1;

          const SaturationState& sat = eos95_->SaturationLineP(p);
          if (!IsExtended_(p, h, sat)) continue;

          bool flag = IsPhysical_(p, h, sat);
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

      bool flag = IsPhysical_(p, h, sat);
      double weight = flag ? 2.0 : options_.extension_weight;
      samples.push_back({p, h, 0.0, 0.0, weight, flag});
    }
  }

  // Additional samples only under the flat critical cutoff.
  // ... pass

  // Recursive build of sample rho/T data
  for (Sample& sample : samples) {
    if (sample.is_physical) {
      auto [prop, liquid, vapor] = eos97_.ThermodynamicsPH(sample.p, sample.h); 
      sample.rho = prop.rho;
      sample.T = prop.T;
    }
  }

  for (Sample& sample : samples) {
    if (!sample.is_physical) {
      double p = sample.p;
      double h = sample.h;
      const SaturationState& sat = eos95_->SaturationLineP(p); 

      int itrs(30);
      double tol(1e-10);
      FrhoT f(p, h, eos95_.get());
      FrhoT::Vector x0(2);
      x0[0] = (h <= HC) ? sat.rhol : sat.rhov;
      x0[1] = sat.Tsat;
      FrhoT::Vector sol = PowellHybrid(x0, f, &itrs, tol);
      // AMANZI_ASSERT(itrs >= 0); // FIXME

      sample.rho = sol[0];
      sample.T = sol[1];
    }
  }

  return samples;
}


/* ******************************************************************
*
****************************************************************** */
bool
IAPWS95_RaggedSplinePH::IsPhysical_(double p, double h, const SaturationState& sat)
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

  double ext_hl = FindLiquidMetastableMargin_(p, sat);
  double ext_hv = FindVaporMetastableMargin_(p, sat);
  return h <= ext_hl || h >= ext_hv;
}

} // namespace AmanziEOS
} // namespace Amanzi

