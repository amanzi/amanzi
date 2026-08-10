#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <limits>
#include <vector>

#include "UnitTest++.h"
#include "Teuchos_ParameterList.hpp"

#include "dbc.hh"
#include "IAPWS95.hh"
#include "IAPWS95_RaggedSplineRhoT.hh"
#include "IAPWS95_RaggedSplinePH.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziEOS;

void WriteHelmholtzPlotData(IAPWS95_RaggedSplinePH& spline,
                            IAPWS95& eos95,
                            double p_min,
                            double p_max,
                            double h_min,
                            double h_max,
                            int np,
                            int nh,
                            const std::string& filename)
{
  std::ofstream out(filename);
  out << std::setprecision(17);

  // Header
  out << "# p h spline_s exact_s error\n";

  for (int j = 0; j < np; ++j) {
    double p = p_min + (p_max - p_min) * (double)j / (np - 1);

    for (int i = 0; i < nh; ++i) {
      double h = h_min + (h_max - h_min) * (double)i / (nh - 1);

      try {
        const auto& spline_value = spline.EntropyDerivativesPH(p, h);
        const auto& exact_value = eos95.EntropyDerivativesPH(p, h);
        double error = spline_value[0] - exact_value[0];

        out << p << " " << h << " " << spline_value[0] << " " << exact_value[0] << " " << error << "\n";
      }
      catch (const std::exception&) {
        // Point is outside the ragged physical domain.
        // NaN preserves the rectangular plotting structure.
        const double nan = std::numeric_limits<double>::quiet_NaN();

        out << p << " " << h << " " << nan << " " << nan << " " << nan << "\n";
      }
    }

    out << "\n";
  }
  std::cout << "Wrote plotting data to " << filename << '\n';
}


TEST(EOS_IAPWS95_SPLINE_P_H)
{
  IAPWS95_RaggedSplinePH::Options opt;
  opt.P_min = 0.50;
  opt.P_max = 50.0;
  opt.H_min = 500.0;
  opt.H_max = 3600.0;
  opt.initial_p_intervals = 50;
  opt.initial_h_intervals = 50;
  opt.extension_cells = 4.0;
  opt.extension_weight = 0.1;

  Teuchos::ParameterList plist;
  IAPWS95_RaggedSplinePH spline(plist, opt);
  IAPWS95 eos95(plist);

  spline.CreateRaggedMesh();
  const auto samples = spline.BuildSamples();
  spline.BuildSplineCoefficients(samples);

  const auto& mesh = spline.GetMesh();

  WriteHelmholtzPlotData(spline,
                         eos95,
                         mesh.x_lines.front(),
                         mesh.x_lines.back(),
                         mesh.y_lines.front(),
                         mesh.y_lines.back(),
                         301,
                         301,
                         "eos_iapws95_spline_ph.dat");
}

