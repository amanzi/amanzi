#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <sstream>

#include "Brent.hh"
#include "lapack.hh"
#include "dbc.hh"

#include "IAPWS95_RaggedSpline.hh"

namespace Amanzi {
namespace AmanziEOS {

int LowerCell(const std::vector<double>& x, double value)
{
  AMANZI_ASSERT(x.size() >= 2); // Coordinate line has fewer than two points.
  if (value <= x.front()) return 0;
  if (value >= x.back()) return x.size() - 2;
  return std::distance(x.begin(), std::upper_bound(x.begin(), x.end(), value)) - 1;
}


/* ******************************************************************
* Constructor
****************************************************************** */
IAPWS95_RaggedSpline::IAPWS95_RaggedSpline(Teuchos::ParameterList& plist, Options options)
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
IAPWS95_RaggedSpline::ResidualPart(double rho, double T)
{ 
  // check for spline domain
  AMANZI_ASSERT(rho >= options_.rho_min && rho <= options_.rho_max);
  AMANZI_ASSERT(T >= options_.T_min && T <= options_.T_max);

  double bndT = BoundaryTemperature(rho);
  double extT = options_.extension_cells * LocalTemperatureSpacing(bndT);
  AMANZI_ASSERT(T >= bndT - extT);

  return ResidualPartDimensionless_(rho / RHOC, TC / T);
}


/* ******************************************************************
* Create a ragged mesh
****************************************************************** */
void
IAPWS95_RaggedSpline::CreateRaggedMesh()
{
  mesh_ = Mesh{};
  coefficients_.clear();
  built_ = false;

  BuildInitialCoordinateLines_();

  // March along T. The first call uses zero guesses; subsequent calls reuse
  // the previous converged densities for fast and branch-consistent solves.
  mesh_.saturation.reserve(mesh_.T_lines.size());
  for (double T : mesh_.T_lines) {
    if (T >= TC) break;

    double rhol0 = eos95_->DensityLiquid(T);
    double rhov0 = eos95_->DensityVapor(T);

    auto [rhol, rhov, p_sat] = eos95_->SaturationLine(T, rhol0, rhov0);
    AMANZI_ASSERT(rhol > rhov && rhov > 0.0 && p_sat > 0.0);

    mesh_.saturation.push_back({T, rhol, rhov, p_sat});
  }

  BuildRaggedColumns_();
  AdaptiveRefineCoordinateLines_();
  BuildKnotVectors_();
}


/* ******************************************************************
* Assemble a global overdetermined least-squares system for the active
* cubic tensor-product basis functions and solve for all spline coefficients.
****************************************************************** */
void
IAPWS95_RaggedSpline::BuildSplineCoefficients()
{
  // CreateRaggedMesh() must be called before fitting
  AMANZI_ASSERT(mesh_.delta_knots.size() > 0 && mesh_.tau_knots.size() > 0);

  const auto samples = BuildSamples_();
  int n = n_delta_basis_ * n_tau_basis_;
  AMANZI_ASSERT(n != 0);  // The spline has no coefficients

  // With coefficient(i,j) = i * n_tau_basis + j, two overlapping cubic
  // tensor-product basis functions differ by at most 3 in each index.
  int kd = 3 * n_tau_basis_ + 3;
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
    const BasisData Bd = EvaluateCubicBasis_(mesh_.delta_knots, delta);
    const BasisData Bt = EvaluateCubicBasis_(mesh_.tau_knots, tau);
    const std::array<double, 6> target = eos95_->ResidualPart(sample.rho, sample.T);

    int id = LowerCell(mesh_.rho_lines, sample.rho);
    int it = LowerCell(mesh_.T_lines, sample.T);

    double tau0 = TC / mesh_.T_lines[it];
    double tau1 = TC / mesh_.T_lines[it + 1];

    const double h_delta = (mesh_.rho_lines[id + 1] - mesh_.rho_lines[id]) / RHOC;
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
        const double xd =
          (component == 0 || component == 2 || component == 5) ? Bd.value[a] :
          (component == 1 || component == 4) ? Bd.d1[a] : Bd.d2[a];

        for (int b = 0; b < 4; ++b) {
          const double xt =
            (component == 0 || component == 1 || component == 3) ? Bt.value[b] :
            (component == 2 || component == 4) ? Bt.d1[b] : Bt.d2[b];

          index[q] = CoefficientIndex(Bd.index[a], Bt.index[b]);
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
IAPWS95_RaggedSpline::BuildInitialCoordinateLines_()
{
  mesh_.rho_lines.resize(options_.initial_rho_intervals + 1);
  mesh_.T_lines.resize(options_.initial_T_intervals + 1);

  // Logarithmic density distribution handles the dilute-vapor scale while
  // remaining monotone and simple. We can replace by another distribution
  // if the liquid region requires stronger clustering.
  double log_min = std::log(options_.rho_min);
  double log_max = std::log(options_.rho_max);
  for (int i = 0; i < mesh_.rho_lines.size(); ++i) {
    double s = (double)i / (mesh_.rho_lines.size() - 1);
    mesh_.rho_lines[i] = std::exp((1.0 - s) * log_min + s * log_max);
  }

  for (int j = 0; j < mesh_.T_lines.size(); ++j) {
    double s = (double)j / (mesh_.T_lines.size() - 1);
    mesh_.T_lines[j] = (1.0 - s) * options_.T_min + s * options_.T_max;
  }
}


/* ******************************************************************
* Should be called every time the coordinate lines change, e.g.
* during adaptive refinement.
****************************************************************** */
void
IAPWS95_RaggedSpline::BuildRaggedColumns_()
{
  mesh_.columns.clear();
  mesh_.columns.reserve(mesh_.rho_lines.size());

  for (double rho : mesh_.rho_lines) {
    double bndT = BoundaryTemperature(rho);

    RaggedColumn c;
    c.rho = rho;
    c.boundary_T = bndT;
    double extT = bndT - options_.extension_cells * LocalTemperatureSpacing(rho);

    c.first_physical_T = static_cast<int>(
      std::lower_bound(mesh_.T_lines.begin(), mesh_.T_lines.end(), bndT) - mesh_.T_lines.begin());
    c.first_extended_T = static_cast<int>(
      std::lower_bound(mesh_.T_lines.begin(), mesh_.T_lines.end(), extT) - mesh_.T_lines.begin());

    mesh_.columns.push_back(c);
  }
}


/* ******************************************************************
* Validate options
****************************************************************** */
void IAPWS95_RaggedSpline::ValidateOptions_() const
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
IAPWS95_RaggedSpline::BoundaryTemperature(double rho) const
{
  // The saturation arrays parameterize two branches rho_v(T) and rho_l(T).
  // For a density between the vapor and liquid branches, return the highest
  // sampled saturation temperature whose dome interval contains rho.
  // This provides a piecewise-linear approximation to the lower exclusion
  // boundary. A prescribed parabola F(rho) can replace this routine directly.
  double result = options_.T_min;
  bool found = false;

  for (int k = 0; k + 1 < mesh_.saturation.size(); ++k) {
    const auto& a = mesh_.saturation[k];
    const auto& b = mesh_.saturation[k + 1];

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
* Minimum dT
****************************************************************** */
double
IAPWS95_RaggedSpline::LocalTemperatureSpacing(double /*rho*/) const
{
  if (mesh_.T_lines.size() < 2) return 0.0;
  double h = std::numeric_limits<double>::max();
  for (int j = 0; j + 1 < mesh_.T_lines.size(); ++j) {
    h = std::min(h, mesh_.T_lines[j + 1] - mesh_.T_lines[j]);
  }
  return h;
}


/* ******************************************************************
* Mesh adaptation
****************************************************************** */
void
IAPWS95_RaggedSpline::AdaptiveRefineCoordinateLines_()
{
  // A complete production version can perform solve-estimate-refine cycles.
  // Here we implement geometry-driven refinement before the first solve:
  // refine cells crossed by the saturation boundary and cells whose endpoint
  // exact derivatives vary more than the configured normalized tolerance.
  for (int pass = 0; pass < options_.max_adaptive_passes; ++pass) {
    bool changed = false;
    std::vector<double> add_rho;
    std::vector<double> add_T;

    if (mesh_.rho_lines.size() - 1 < options_.max_rho_intervals) {
      for (int i = 0; i + 1 < mesh_.rho_lines.size(); ++i) {
        double r0 = mesh_.rho_lines[i];
        double r1 = mesh_.rho_lines[i + 1];
        double rm = std::sqrt(r0 * r1);
        double curvature = std::fabs(BoundaryTemperature(r0)
                               - 2 * BoundaryTemperature(rm)
                                   + BoundaryTemperature(r1));
        double hT = LocalTemperatureSpacing(rm);
        if (curvature > 0.25 * hT) add_rho.push_back(rm);
      }
    }

    if (mesh_.T_lines.size() - 1 < options_.max_T_intervals) {
      for (int j = 0; j + 1 < mesh_.T_lines.size(); ++j) {
        double T0 = mesh_.T_lines[j];
        double T1 = mesh_.T_lines[j + 1];
        double Tm = 0.5 * (T0 + T1);

        bool crossed = false;
        for (double rho : mesh_.rho_lines) {
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

    insert_unique(mesh_.rho_lines, add_rho, options_.max_rho_intervals);
    insert_unique(mesh_.T_lines, add_T, options_.max_T_intervals);
    if (!changed) break;

    BuildRaggedColumns_();
  }
}


/* ******************************************************************
*
****************************************************************** */
void
IAPWS95_RaggedSpline::BuildKnotVectors_()
{
  std::vector<double> delta_lines(mesh_.rho_lines.size());
  std::transform(mesh_.rho_lines.begin(), mesh_.rho_lines.end(),
                 delta_lines.begin(),
                 [&](double rho) { return rho / RHOC; });

  // tau decreases as T increases. Build an increasing tau line array.
  std::vector<double> tau_lines(mesh_.T_lines.size());
  std::transform(mesh_.T_lines.begin(), mesh_.T_lines.end(),
                 tau_lines.begin(),
                 [&](double T) { return TC / T; });
  std::reverse(tau_lines.begin(), tau_lines.end());

  mesh_.delta_knots = MakeClampedCubicKnots_(delta_lines);
  mesh_.tau_knots = MakeClampedCubicKnots_(tau_lines);
  n_delta_basis_ = mesh_.delta_knots.size() - degree_ - 1;
  n_tau_basis_ = mesh_.tau_knots.size() - degree_ - 1;
}


/* ******************************************************************
* Constructs knot vector for a cubic B-spline from the coordinate lines.
* For cubic splines, the degree is p=3. A clamped knot vector repeats 
* the first and last knot p+1=4 times. This makes the spline interpolate 
* the endpoint behavior in the usual open-knot-vector sense and ensures 
* that the basis spans the full interval.
****************************************************************** */
std::vector<double>
IAPWS95_RaggedSpline::MakeClampedCubicKnots_(const std::vector<double>& lines)
{
  std::vector<double> knots;
  knots.reserve(lines.size() + 6);
  for (int k = 0; k < 4; ++k) knots.push_back(lines.front());
  for (int i = 1; i + 1 < lines.size(); ++i) knots.push_back(lines[i]);
  for (int k = 0; k < 4; ++k) knots.push_back(lines.back());
  return knots;
}


/* ******************************************************************
* Creates the set of (rho, T) points used to fit spline coefficients
****************************************************************** */
std::vector<IAPWS95_RaggedSpline::Sample>
IAPWS95_RaggedSpline::BuildSamples_() const
{
  std::vector<Sample> samples;
  const unsigned q = options_.samples_per_cell_direction;

  // interior cell samples
  for (int i = 0; i + 1 < mesh_.rho_lines.size(); ++i) {
    double r0 = mesh_.rho_lines[i];
    double r1 = mesh_.rho_lines[i + 1];
    for (int j = 0; j + 1 < mesh_.T_lines.size(); ++j) {
      double T0 = mesh_.T_lines[j];
      double T1 = mesh_.T_lines[j + 1];

      for (unsigned ir = 0; ir < q; ++ir) {
        double sr = (ir + 0.5) / q;
        double rho = (1.0 - sr) * r0 + sr * r1;
        for (unsigned jt = 0; jt < q; ++jt) {
          double st = (jt + 0.5) / q;
          double T = (1.0 - st) * T0 + st * T1;

          if (!IsExtended_(rho, T)) continue;
          double weight = IsPhysical_(rho, T) ? 1.0 : options_.extension_weight;
          samples.push_back({rho, T, weight});
        }
      }
    }
  }

  // Add all background-grid vertices. This directly constrains clamped
  // endpoint coefficients, especially corners such as (rho_max, T_min).
  for (double rho : mesh_.rho_lines) {
    for (double T : mesh_.T_lines) {
      if (!IsExtended_(rho, T)) continue;
      double weight = IsPhysical_(rho, T) ? 2.0 : options_.extension_weight;
      samples.push_back({rho, T, weight});
    }
  }

  return samples;
}


/* ******************************************************************
*
****************************************************************** */
IAPWS95_RaggedSpline::BasisData
IAPWS95_RaggedSpline::EvaluateCubicBasis_(const std::vector<double>& U, double x) const
{
  const int p = 3;
  const int n = static_cast<int>(U.size()) - p - 2;
  // Invalid cubic knot vector.
  AMANZI_ASSERT(n >= p);

  x = std::clamp(x, U[p], U[n + 1]);
  int span;
  if (x >= U[n + 1]) {
    span = n;
  } else {
    int low = p, high = n + 1;
    span = (low + high) / 2;
    while (x < U[span] || x >= U[span + 1]) {
        if (x < U[span]) high = span;
        else low = span;
        span = (low + high) / 2;
    }
  }

  // Algorithm A2.3 from The NURBS Book: derivatives of nonzero basis functions.
  double ndu[4][4] = {};
  double left[4] = {}, right[4] = {};
  ndu[0][0] = 1.0;
  for (int j = 1; j <= p; ++j) {
    left[j] = x - U[span + 1 - j];
    right[j] = U[span + j] - x;
    double saved = 0.0;
    for (int r = 0; r < j; ++r) {
      ndu[j][r] = right[r + 1] + left[j - r];
      const double temp = ndu[r][j - 1] / ndu[j][r];
      ndu[r][j] = saved + right[r + 1] * temp;
      saved = left[j - r] * temp;
    }
    ndu[j][j] = saved;
  }

  double ders[3][4] = {};
  for (int j = 0; j <= p; ++j) ders[0][j] = ndu[j][p];

  double a[2][4] = {};
  for (int r = 0; r <= p; ++r) {
    int s1 = 0, s2 = 1;
    a[0][0] = 1.0;
    for (int k = 1; k <= 2; ++k) {
      double d = 0.0;
      const int rk = r - k;
      const int pk = p - k;
      if (r >= k) {
        a[s2][0] = a[s1][0] / ndu[pk + 1][rk];
        d = a[s2][0] * ndu[rk][pk];
      }
      const int j1 = (rk >= -1) ? 1 : -rk;
      const int j2 = (r - 1 <= pk) ? k - 1 : p - r;
      for (int j = j1; j <= j2; ++j) {
        a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][rk + j];
        d += a[s2][j] * ndu[rk + j][pk];
      }
      if (r <= pk) {
        a[s2][k] = -a[s1][k - 1] / ndu[pk + 1][r];
        d += a[s2][k] * ndu[r][pk];
      }
      ders[k][r] = d;
      std::swap(s1, s2);
    }
  }

  int factor = p;
  for (int k = 1; k <= 2; ++k) {
    for (int j = 0; j <= p; ++j) ders[k][j] *= factor;
    factor *= (p - k);
  }

  BasisData out;
  for (int j = 0; j < 4; ++j) {
    out.index[j] = span - p + j;
    out.value[j] = ders[0][j];
    out.d1[j] = ders[1][j];
    out.d2[j] = ders[2][j];
  }
  return out;
}


/* ******************************************************************
* Evaluate residual part uing demensionless input
****************************************************************** */
std::array<double, 6>
IAPWS95_RaggedSpline::ResidualPartDimensionless_(double delta, double tau) const
{
  // Check that spline coefficients have been built
  AMANZI_ASSERT(built_);

  const BasisData Bd = EvaluateCubicBasis_(mesh_.delta_knots, delta);
  const BasisData Bt = EvaluateCubicBasis_(mesh_.tau_knots, tau);
  std::array<double, 6> out{};

  for (int a = 0; a < 4; ++a) {
    int i = Bd.index[a];
    for (int c = 0; c < 4; ++c) {
      int j = Bt.index[c];
      double z = coefficients_[CoefficientIndex(i, j)];
      out[0] += z * Bd.value[a] * Bt.value[c];
      out[1] += z * Bd.d1[a] * Bt.value[c];
      out[2] += z * Bd.value[a] * Bt.d1[c];
      out[3] += z * Bd.d2[a] * Bt.value[c];
      out[4] += z * Bd.d1[a] * Bt.d1[c];
      out[5] += z * Bd.value[a] * Bt.d2[c];
    }
  }
  return out;
}


/* ******************************************************************
*
****************************************************************** */
bool
IAPWS95_RaggedSpline::IsExtended_(double rho, double T) const
{
  double F = BoundaryTemperature(rho);
  double width = options_.extension_cells * LocalTemperatureSpacing(rho);
  return T >= F - width;
}

} // namespace AmanziEOS
} // namespace Amanzi

