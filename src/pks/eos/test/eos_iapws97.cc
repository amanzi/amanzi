/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
*/

#include "UnitTest++.h"

#include "IAPWS97.hh"

TEST(EOS_IAPWS97_PT)
{
  using namespace Amanzi::AmanziEOS;

  Teuchos::ParameterList plist;
  IAPWS97 eos(plist);

  double p, T, v;
  T = eos.BoundaryLine23(16.5291643); // after Table 1
  CHECK_CLOSE(623.15, T, 1e-6);

  // region 1
  auto prop = eos.Region1(3.0, 300.0);
  CHECK_CLOSE(0.100215168e-2, prop.v, 1e-11);
  CHECK_CLOSE(0.115331273e+3, prop.h, 1e-6);
  CHECK_CLOSE(0.112324818e+3, prop.u, 1e-6);
  CHECK_CLOSE(0.392294792, prop.s, 1e-9);
  CHECK_CLOSE(0.417301218e+1, prop.cp, 1e-8);
  CHECK_CLOSE(0.150773921e+4, prop.w, 1e-5);

  prop = eos.Region1(80.0, 300.0);
  CHECK_CLOSE(0.971180894e-3, prop.v, 1e-11);
  CHECK_CLOSE(0.184142828e+3, prop.h, 1e-6);
  CHECK_CLOSE(0.106448356e+3, prop.u, 1e-6);
  CHECK_CLOSE(0.368563852, prop.s, 1e-9);
  CHECK_CLOSE(0.401008987e+1, prop.cp, 1e-8);
  CHECK_CLOSE(3.91736606, prop.cv, 1e-8);
  CHECK_CLOSE(0.163469054e+4, prop.w, 1e-5);

  prop = eos.Region1(3.0, 500.0);
  CHECK_CLOSE(0.120241800e-2, prop.v, 1e-11);
  CHECK_CLOSE(0.975542239e+3, prop.h, 1e-6);
  CHECK_CLOSE(0.971934985e+3, prop.u, 1e-6);
  CHECK_CLOSE(0.258041912e+1, prop.s, 1e-8);
  CHECK_CLOSE(0.465580682e+1, prop.cp, 1e-8);
  CHECK_CLOSE(0.124071337e+4, prop.w, 1e-5);
  CHECK_CLOSE(0.112892188e-2, prop.kt, 1e-11);
  CHECK_CLOSE(0.00164118128, prop.av, 1e-11);

  // region 2
  prop = eos.Region2(0.0035, 300.0);
  CHECK_CLOSE(0.394913866e+2, prop.v, 1e-7);
  CHECK_CLOSE(0.254991145e+4, prop.h, 1e-5);
  CHECK_CLOSE(0.241169160e+4, prop.u, 1e-5);
  CHECK_CLOSE(0.852238967e+1, prop.s, 1e-8);
  CHECK_CLOSE(0.191300162e+1, prop.cp, 1e-8);
  CHECK_CLOSE(0.427920172e+3, prop.w, 1e-6);
  CHECK_CLOSE(0.286239651e+3, prop.kt, 1e-6);
  CHECK_CLOSE(0.00337578289, prop.av, 1e-11);

  prop = eos.Region2(0.0035, 700.0);
  CHECK_CLOSE(0.923015898e+2, prop.v, 1e-7);
  CHECK_CLOSE(0.333568375e+4, prop.h, 1e-5);
  CHECK_CLOSE(0.301262819e+4, prop.u, 1e-5);
  CHECK_CLOSE(0.101749996e+2, prop.s, 1e-7);
  CHECK_CLOSE(0.208141274e+1, prop.cp, 1e-8);
  CHECK_CLOSE(0.161978333e+1, prop.cv, 1e-8);
  CHECK_CLOSE(0.644289068e+3, prop.w, 1e-6);

  prop = eos.Region2(30.0, 700.0);
  CHECK_CLOSE(0.542946619e-2, prop.v, 1e-11);
  CHECK_CLOSE(0.263149474e+4, prop.h, 1e-5);
  CHECK_CLOSE(0.246861076e+4, prop.u, 1e-5);
  CHECK_CLOSE(0.517540298e+1, prop.s, 1e-8);
  CHECK_CLOSE(0.103505092e+2, prop.cp, 1e-8);
  CHECK_CLOSE(0.480386523e+3, prop.w, 1e-6);

  // region 3: boundaries
  T = eos.Boundary3Formula1(25.0, IAPWS97::rgn3n_cd);
  CHECK_CLOSE(649.3659208, T, 1e-7);

  T = eos.Boundary3Formula1(23.0, IAPWS97::rgn3n_gh);
  CHECK_CLOSE(649.8873759, T, 1e-7);

  T = eos.Boundary3Formula1(23.0, IAPWS97::rgn3n_ij);
  CHECK_CLOSE(651.5778091, T, 1e-7);

  T = eos.Boundary3Formula1(23.0, IAPWS97::rgn3n_jk);
  CHECK_CLOSE(655.8338344, T, 1e-7);

  T = eos.Boundary3Formula1(22.8, IAPWS97::rgn3n_mn);
  CHECK_CLOSE(649.6054133, T, 1e-7);

  T = eos.Boundary3Formula1(22.0, IAPWS97::rgn3n_qu);
  CHECK_CLOSE(645.6355027, T, 1e-7);

  T = eos.Boundary3Formula1(22.0, IAPWS97::rgn3n_rx);
  CHECK_CLOSE(648.2622754, T, 1e-7);

  T = eos.Boundary3Formula1(22.3, IAPWS97::rgn3n_uv);
  CHECK_CLOSE(647.7996121, T, 1e-7);

  T = eos.Boundary3Formula2(40.0, IAPWS97::rgn3k_ab, IAPWS97::rgn3n_ab);
  CHECK_CLOSE(693.0341408, T, 1e-7);

  T = eos.Boundary3Formula2(22.8, IAPWS97::rgn3k_op, IAPWS97::rgn3n_op);
  CHECK_CLOSE(650.0106943, T, 1e-7);

  T = eos.Boundary3Formula3(40.0);
  CHECK_CLOSE(713.9593992, T, 1e-7);

  // region 3: specific valumes
  v = eos.Region3_SubregionFormula(50.0, 630.0, IAPWS97::ModelId::a);
  CHECK_CLOSE(0.001470853100, v, 1e-11);
  v = eos.Region3_SubregionFormula(80.0, 670.0, IAPWS97::ModelId::a);
  CHECK_CLOSE(0.001503831359, v, 1e-11);

  v = eos.Region3_SubregionFormula(50.0, 710.0, IAPWS97::ModelId::b);
  CHECK_CLOSE(0.002204728587, v, 1e-11);
  v = eos.Region3_SubregionFormula(80.0, 750.0, IAPWS97::ModelId::b);
  CHECK_CLOSE(0.001973692940, v, 1e-11);

  v = eos.Region3_SubregionFormula(20.0, 630.0, IAPWS97::ModelId::c);
  CHECK_CLOSE(0.001761696406, v, 1e-11);
  v = eos.Region3_SubregionFormula(30.0, 650.0, IAPWS97::ModelId::c);
  CHECK_CLOSE(0.001819560617, v, 1e-11);

  v = eos.Region3_SubregionFormula(26.0, 656.0, IAPWS97::ModelId::d);
  CHECK_CLOSE(0.002245587720, v, 1e-11);
  v = eos.Region3_SubregionFormula(30.0, 670.0, IAPWS97::ModelId::d);
  CHECK_CLOSE(0.002506897702, v, 1e-11);

  v = eos.Region3_SubregionFormula(26.0, 661.0, IAPWS97::ModelId::e);
  CHECK_CLOSE(0.002970225962, v, 1e-11);
  v = eos.Region3_SubregionFormula(30.0, 675.0, IAPWS97::ModelId::e);
  CHECK_CLOSE(0.003004627086, v, 1e-11);

  v = eos.Region3_SubregionFormula(26.0, 671.0, IAPWS97::ModelId::f);
  CHECK_CLOSE(0.005019029401, v, 1e-11);
  v = eos.Region3_SubregionFormula(30.0, 690.0, IAPWS97::ModelId::f);
  CHECK_CLOSE(0.004656470142, v, 1e-11);

  v = eos.Region3_SubregionFormula(23.6, 649.0, IAPWS97::ModelId::g);
  CHECK_CLOSE(0.002163198378, v, 1e-11);
  v = eos.Region3_SubregionFormula(24.0, 650.0, IAPWS97::ModelId::g);
  CHECK_CLOSE(0.002166044161, v, 1e-11);

  v = eos.Region3_SubregionFormula(23.6, 652.0, IAPWS97::ModelId::h);
  CHECK_CLOSE(0.002651081407, v, 1e-11);
  v = eos.Region3_SubregionFormula(24.0, 654.0, IAPWS97::ModelId::h);
  CHECK_CLOSE(0.002967802335, v, 1e-11);

  v = eos.Region3_SubregionFormula(23.6, 653.0, IAPWS97::ModelId::i);
  CHECK_CLOSE(0.003273916816, v, 1e-11);
  v = eos.Region3_SubregionFormula(24.0, 655.0, IAPWS97::ModelId::i);
  CHECK_CLOSE(0.003550329864, v, 1e-11);

  v = eos.Region3_SubregionFormula(23.5, 655.0, IAPWS97::ModelId::j);
  CHECK_CLOSE(0.004545001142, v, 1e-11);
  v = eos.Region3_SubregionFormula(24.0, 660.0, IAPWS97::ModelId::j);
  CHECK_CLOSE(0.005100267704, v, 1e-11);

  v = eos.Region3_SubregionFormula(23.0, 660.0, IAPWS97::ModelId::k);
  CHECK_CLOSE(0.006109525997, v, 1e-11);
  v = eos.Region3_SubregionFormula(24.0, 670.0, IAPWS97::ModelId::k);
  CHECK_CLOSE(0.006427325645, v, 1e-11);

  v = eos.Region3_SubregionFormula(22.6, 646.0, IAPWS97::ModelId::l);
  CHECK_CLOSE(0.002117860851, v, 1e-11);
  v = eos.Region3_SubregionFormula(23.0, 646.0, IAPWS97::ModelId::l);
  CHECK_CLOSE(0.002062374674, v, 1e-11);

  v = eos.Region3_SubregionFormula(22.6, 648.6, IAPWS97::ModelId::m);
  CHECK_CLOSE(0.002533063780, v, 1e-11);
  v = eos.Region3_SubregionFormula(22.8, 649.3, IAPWS97::ModelId::m);
  CHECK_CLOSE(0.002572971781, v, 1e-11);

  v = eos.Region3_SubregionFormula(22.6, 649.0, IAPWS97::ModelId::n);
  CHECK_CLOSE(0.002923432711, v, 1e-11);
  v = eos.Region3_SubregionFormula(22.8, 649.7, IAPWS97::ModelId::n);
  CHECK_CLOSE(0.002913311494, v, 1e-11);

  v = eos.Region3_SubregionFormula(22.6, 649.1, IAPWS97::ModelId::o);
  CHECK_CLOSE(0.003131208996, v, 1e-11);
  v = eos.Region3_SubregionFormula(22.8, 649.9, IAPWS97::ModelId::o);
  CHECK_CLOSE(0.003221160278, v, 1e-11);

  v = eos.Region3_SubregionFormula(22.6, 649.4, IAPWS97::ModelId::p);
  CHECK_CLOSE(0.003715596186, v, 1e-11);
  v = eos.Region3_SubregionFormula(22.8, 650.2, IAPWS97::ModelId::p);
  CHECK_CLOSE(0.003664754790, v, 1e-11);

  v = eos.Region3_SubregionFormula(21.1, 640.0, IAPWS97::ModelId::q);
  CHECK_CLOSE(0.001970999272, v, 1e-11);
  v = eos.Region3_SubregionFormula(21.8, 643.0, IAPWS97::ModelId::q);
  CHECK_CLOSE(0.002043919161, v, 1e-11);

  v = eos.Region3_SubregionFormula(21.1, 644.0, IAPWS97::ModelId::r);
  CHECK_CLOSE(0.005251009921, v, 1e-11);
  v = eos.Region3_SubregionFormula(21.8, 648.0, IAPWS97::ModelId::r);
  CHECK_CLOSE(0.005256844741, v, 1e-11);

  v = eos.Region3_SubregionFormula(19.1, 635.0, IAPWS97::ModelId::s);
  CHECK_CLOSE(0.001932829079, v, 1e-11);
  v = eos.Region3_SubregionFormula(20.0, 638.0, IAPWS97::ModelId::s);
  CHECK_CLOSE(0.001985387227, v, 1e-11);

  v = eos.Region3_SubregionFormula(17.0, 626.0, IAPWS97::ModelId::t);
  CHECK_CLOSE(0.008483262001, v, 1e-11);
  v = eos.Region3_SubregionFormula(20.0, 640.0, IAPWS97::ModelId::t);
  CHECK_CLOSE(0.006227528101, v, 1e-11);

  v = eos.Region3_SubregionFormula(21.5, 644.6, IAPWS97::ModelId::u);
  CHECK_CLOSE(0.002268366647, v, 1e-11);
  v = eos.Region3_SubregionFormula(22.0, 646.1, IAPWS97::ModelId::u);
  CHECK_CLOSE(0.002296350553, v, 1e-11);

  v = eos.Region3_SubregionFormula(22.5, 648.6, IAPWS97::ModelId::v);
  CHECK_CLOSE(0.002832373260, v, 1e-11);
  v = eos.Region3_SubregionFormula(22.3, 647.9, IAPWS97::ModelId::v);
  CHECK_CLOSE(0.002811424405, v, 1e-11);

  v = eos.Region3_SubregionFormula(22.15, 647.5, IAPWS97::ModelId::w);
  CHECK_CLOSE(0.003694032281, v, 1e-11);
  v = eos.Region3_SubregionFormula(22.3, 648.1, IAPWS97::ModelId::w);
  CHECK_CLOSE(0.003622226305, v, 1e-11);

  v = eos.Region3_SubregionFormula(22.11, 648.0, IAPWS97::ModelId::x);
  CHECK_CLOSE(0.004528072649, v, 1e-11);
  v = eos.Region3_SubregionFormula(22.3, 649.0, IAPWS97::ModelId::x);
  CHECK_CLOSE(0.004556905799, v, 1e-11);

  v = eos.Region3_SubregionFormula(22.0, 646.84, IAPWS97::ModelId::y);
  CHECK_CLOSE(0.002698354719, v, 1e-11);
  v = eos.Region3_SubregionFormula(22.064, 647.05, IAPWS97::ModelId::y);
  CHECK_CLOSE(0.002717655648, v, 1e-11);

  v = eos.Region3_SubregionFormula(22.0, 646.89, IAPWS97::ModelId::z);
  CHECK_CLOSE(0.003798732962, v, 1e-11);
  v = eos.Region3_SubregionFormula(22.064, 647.15, IAPWS97::ModelId::z);
  CHECK_CLOSE(0.003701940009, v, 1e-11);

  // Region 3
  prop = eos.Region3(500.0, 650.0);
  CHECK_CLOSE(25.5837018, prop.p, 1e-6);
  CHECK_CLOSE(1863.43019, prop.h, 1e-5);
  CHECK_CLOSE(1812.26279, prop.u, 1e-5);
  CHECK_CLOSE(0.405427273e+1, prop.s, 1e-8);
  CHECK_CLOSE(0.138935717e+2, prop.cp, 1e-7);
  CHECK_CLOSE(0.502005554e+3, prop.w, 1e-6);

  prop = eos.Region3(200.0, 650.0);
  CHECK_CLOSE(0.237512401e+4, prop.h, 1e-5);
  CHECK_CLOSE(0.226365868e+4, prop.u, 1e-5);
  CHECK_CLOSE(0.485438792e+1, prop.s, 1e-8);
  CHECK_CLOSE(0.446579342e+2, prop.cp, 1e-7);
  CHECK_CLOSE(4.04118076, prop.cv, 1e-7);
  CHECK_CLOSE(0.383444594e+3, prop.w, 1e-6);

  prop = eos.Region3(500.0, 750.0);
  CHECK_CLOSE(0.225868845e+4, prop.h, 1e-5);
  CHECK_CLOSE(0.210206932e+4, prop.u, 1e-5);
  CHECK_CLOSE(0.446971906e+1, prop.s, 1e-8);
  CHECK_CLOSE(0.634165359e+1, prop.cp, 1e-7);
  CHECK_CLOSE(0.760696041e+3, prop.w, 1e-6);
  CHECK_CLOSE(0.00806710817, prop.kt, 1e-11);
  CHECK_CLOSE(0.00441515098, prop.av, 1e-11);
  eos.Print(prop);

  // Region 4
  T = eos.SaturationLineP(0.1);
  CHECK_CLOSE(372.755919, T, 1e-6); // Table 36

  T = eos.SaturationLineP(1.0);
  CHECK_CLOSE(453.035632, T, 1e-6);

  T = eos.SaturationLineP(10.0);
  CHECK_CLOSE(584.149488, T, 1e-6);

  p = eos.SaturationLineT(500.0);
  CHECK_CLOSE(2.63889776, p, 1e-8);  // Table 35

  p = eos.SaturationLineT(643.15);
  CHECK_CLOSE(IAPWS97::PSAT_643_15, p, 1e-8);

  // region 5
  prop = eos.Region5(0.5, 1500.0);
  CHECK_CLOSE(0.138455090e+1, prop.v, 1e-8);
  CHECK_CLOSE(0.521976855e+4, prop.h, 1e-5);
  CHECK_CLOSE(0.452749310e+4, prop.u, 1e-5);
  CHECK_CLOSE(0.965408875e+1, prop.s, 1e-8);
  CHECK_CLOSE(0.261609445e+1, prop.cp, 1e-8);
  CHECK_CLOSE(0.917068690e+3, prop.w, 1e-6);

  prop = eos.Region5(30.0, 1500.0);
  CHECK_CLOSE(0.230761299e-1, prop.v, 1e-8);
  CHECK_CLOSE(0.516723514e+4, prop.h, 1e-5);
  CHECK_CLOSE(0.447495124e+4, prop.u, 1e-5);
  CHECK_CLOSE(0.772970133e+1, prop.s, 1e-8);
  CHECK_CLOSE(0.272724317e+1, prop.cp, 1e-8);
  CHECK_CLOSE(0.219274829e+1, prop.cv, 1e-8);
  CHECK_CLOSE(0.928548002e+3, prop.w, 1e-6);

  prop = eos.Region5(30.0, 2000.0);
  CHECK_CLOSE(0.311385219e-1, prop.v, 1e-8);
  CHECK_CLOSE(0.657122604e+4, prop.h, 1e-5);
  CHECK_CLOSE(0.563707038e+4, prop.u, 1e-5);
  CHECK_CLOSE(0.853640523e+1, prop.s, 1e-8);
  CHECK_CLOSE(0.288569882e+1, prop.cp, 1e-8);
  CHECK_CLOSE(0.106736948e+4, prop.w, 1e-5);
  CHECK_CLOSE(0.329193892e-1, prop.kt, 1e-10);
  CHECK_CLOSE(0.000508830641, prop.av, 1e-10);

  // full thermodynamics: PT formulation
  prop = eos.ThermodynamicsPT(3.0, 500.0);
  CHECK_CLOSE(0.258041912e+1, prop.s, 1e-8);
  CHECK_EQUAL(prop.rgn, 1);

  prop = eos.ThermodynamicsPT(0.0006112127, 323.15);
  CHECK_CLOSE(2594.66, prop.h, 0.2);
  CHECK_EQUAL(prop.rgn, 2);

  prop = eos.ThermodynamicsPT(30.0, 700.0);
  CHECK_CLOSE(0.103505092e+2, prop.cp, 1e-8);
  CHECK_EQUAL(prop.rgn, 2);

  prop = eos.ThermodynamicsPT(25.5837018, 650.0);
  CHECK_CLOSE(0.502005554e+3, prop.w, 1e-6);
  CHECK_EQUAL(prop.rgn, 3);
  CHECK(eos.get_itrs() < 3);

  prop = eos.ThermodynamicsPT(22.29305999995, 650.0);
  CHECK_CLOSE(44.65762757698, prop.cp, 1e-9);
  CHECK_EQUAL(prop.rgn, 3);
  CHECK(eos.get_itrs() < 3);

  prop = eos.ThermodynamicsPT(78.30956391692, 750.0);
  CHECK_CLOSE(7.606960408824e+02, prop.w, 1e-6);
  CHECK_EQUAL(prop.rgn, 3);
  CHECK(eos.get_itrs() < 3);
}


TEST(EOS_IAPWS97_TABLE_PT)
{
  using namespace Amanzi::AmanziEOS;
  std::cout << "\nTEST: PT table" << std::endl;

  Teuchos::ParameterList plist;
  IAPWS97 eos(plist);

  int nitrs(0), mitrs(0), count(0);
  double dp(5.0e-2), dT(1.0), p, T, rho; 

  for (int i = -50; i < 50; ++i) {
    for (int j = -50; j < 50; ++j) {
      p = eos.PC + i * dp; // Mpa
      T = eos.TC + j * dT;
      rho = (eos.ThermodynamicsPT(p, T)).rho;
      nitrs += eos.get_itrs();
      mitrs = std::max(mitrs, eos.get_itrs());
      count++;
      // std::cout << p << " " << T << " " << rho << std::endl;
    }
  }
  std::cout << "mean itrs = " << double(nitrs) / count << std::endl;
  std::cout << "max itrs = " << mitrs << std::endl;
}


TEST(EOS_IAPWS97_CRITICAL_POINT_PT)
{
  using namespace Amanzi::AmanziEOS;
  std::cout << "\nTEST: critical point" << std::endl;

  Teuchos::ParameterList plist;
  IAPWS97 eos(plist);

  double dp(0.01), p;
  for (int i = -6; i < -2; ++i) {
    p = eos.PC + i * dp; 
    auto prop = eos.ThermodynamicsPT(p, eos.TC);;
    std::cout << "dp = " << p - eos.PC << " cp = " << prop.cp << " cv = " << prop.cv << std::endl;
  }
}


TEST(EOS_IAPWS97_PH)
{
  using namespace Amanzi::AmanziEOS;
  std::cout << "\nTEST: PH-formulation" << std::endl;

  Teuchos::ParameterList plist;
  IAPWS97 eos(plist);

  double p, T, h, v, k, sigma, mu;
  Properties prop;

  p = eos.SaturationLineH(1700.0);
  CHECK_CLOSE(17.24175718, p, 1e-8);
  p = eos.SaturationLineH(2400.0);
  CHECK_CLOSE(20.18090839, p, 1e-8);

  T = eos.Region1_BackwardPH(3.0, 500.0);
  CHECK_CLOSE(391.798509, T, 1e-6);
  T = eos.Region1_BackwardPH(80.0, 1500.0);
  CHECK_CLOSE(611.041229, T, 1e-6);

  T = eos.Region2_SubregionFormulaA(0.001, 3000.0);
  CHECK_CLOSE(534.433241, T, 1e-6);
  T = eos.Region2_SubregionFormulaA(3.0, 4000.0);
  CHECK_CLOSE(1010.77577, T, 1e-5);
  T = eos.Region2_SubregionFormulaB(5.0, 4000.0);
  CHECK_CLOSE(1015.31583, T, 1e-5);
  T = eos.Region2_SubregionFormulaB(25.0, 3500.0);
  CHECK_CLOSE(875.279054, T, 1e-6);
  T = eos.Region2_SubregionFormulaC(40.0, 2700.0);
  CHECK_CLOSE(743.056411, T, 1e-6);
  T = eos.Region2_SubregionFormulaC(60.0, 3200.0);
  CHECK_CLOSE(882.756860, T, 1e-6);

  h = eos.BoundaryLine2b2c(100.0);
  CHECK_CLOSE(3516.004323, h, 1e-6);

  h = eos.BoundaryLine3a3b(25.0);
  CHECK_CLOSE(2095.936454, h, 1e-6);

  v = eos.Region3_BackwardPH_v(20.0, 1700.0); // Table 2.43
  CHECK_CLOSE(0.001749903962, v, 1e-11);
  v = eos.Region3_BackwardPH_v(50.0, 2000.0);
  CHECK_CLOSE(1.908139035e-3, v, 1e-11);
  v = eos.Region3_BackwardPH_v(100.0, 2100.0);
  CHECK_CLOSE(0.001676229776, v, 1e-11);

  v = eos.Region3_BackwardPH_v(20.0, 2500.0); 
  CHECK_CLOSE(6.670547043e-3, v, 1e-11);
  v = eos.Region3_BackwardPH_v(50.0, 2400.0);
  CHECK_CLOSE(2.801244590e-3, v, 1e-11);
  v = eos.Region3_BackwardPH_v(100.0, 2700.0);
  CHECK_CLOSE(2.404234998e-3, v, 1e-11);

  T = eos.Region3_BackwardPH_T(20.0, 1700.0); // Table 2.47
  CHECK_CLOSE(6.293083892e2, T, 1e-7);
  T = eos.Region3_BackwardPH_T(50.0, 2000.0);
  CHECK_CLOSE(6.905718338e2, T, 1e-7);
  T = eos.Region3_BackwardPH_T(100.0, 2100.0);
  CHECK_CLOSE(7.336163014e2, T, 1e-7);

  T = eos.Region3_BackwardPH_T(20.0, 2500.0);
  CHECK_CLOSE(6.418418053e2, T, 1e-7);
  T = eos.Region3_BackwardPH_T(50.0, 2400.0);
  CHECK_CLOSE(7.351848618e2, T, 1e-7);
  T = eos.Region3_BackwardPH_T(100.0, 2700.0);
  CHECK_CLOSE(8.420460876e2, T, 1e-7);

  // other properties
  k = eos.ThermalConductivity(998.0, 298.15, prop); // Table 4
  CHECK_CLOSE(0.607712868, k, 1e-9);
  k = eos.ThermalConductivity(0.0, 873.15, prop);
  CHECK_CLOSE(0.0791034659, k, 1e-10);
  k = eos.ThermalConductivity(1200.0, 298.15, prop);
  CHECK_CLOSE(0.799038144, k, 1e-9);

  sigma = eos.SurfaceTension(300.0);
  CHECK_CLOSE(0.0716859625, sigma, 1e-10);
  sigma = eos.SurfaceTension(450.0);
  CHECK_CLOSE(0.0428914992, sigma, 1e-10);

  mu = eos.Viscosity(998.0, 298.15); // Table 4
  CHECK_CLOSE(8.89735100e-4, mu, 1e-12);
  mu = eos.Viscosity(600.0, 873.15);
  CHECK_CLOSE(7.7430195e-05, mu, 1e-12);
  mu = eos.Viscosity(1.0, 873.15);
  CHECK_CLOSE(3.2619287e-05, mu, 1e-12);
  mu = eos.Viscosity(100.0, 873.15);
  CHECK_CLOSE(3.5802262e-05, mu, 1e-12);
  mu = eos.Viscosity(400.0, 1173.15);
  CHECK_CLOSE(6.4154608e-05, mu, 1e-12);

  {
    auto [prop, liquid, gas] = eos.ThermodynamicsPH(20.0, 2988.8);
    CHECK(prop.rgn == 2);
    CHECK(prop.x == 1.0); // gas phas
  }

  {
    auto [prop, liquid, gas] = eos.ThermodynamicsPH(16.0, 2087.0);
    CHECK(prop.rgn == 4);
    CHECK_CLOSE(prop.h, (1.0 - prop.x) * liquid.h + prop.x * gas.h, 1e-8);
    CHECK_CLOSE(prop.v, (1.0 - prop.x) * liquid.v + prop.x * gas.v, 1e-10);
    CHECK(liquid.mu > gas.mu);
  }

  {
    auto [prop, liquid, gas] = eos.ThermodynamicsPH(40.0, 2200.0);
    CHECK(prop.rgn == 3);
  }

  {
    auto [prop, liquid, gas] = eos.ThermodynamicsPH(20.0, 1750.0);
    CHECK(prop.rgn == 3);
  }
}


TEST(EOS_IAPWS97_TABLE_PH)
{
  using namespace Amanzi::AmanziEOS;
  std::cout << "\nTEST: PH table" << std::endl;

  Teuchos::ParameterList plist;
  IAPWS97 eos(plist);

  int nitrs(0), mitrs(0), count(0);
  double dp(5.0e-2), dh(15.0), p, h, T, Tmin(1000.0), Tmax(100.0); 

  for (int i = -50; i < 50; ++i) {
    for (int j = -50; j < 50; ++j) {
      p = eos.PC + i * dp; // Mpa
      h = eos.HC + j * dh;
      auto [prop, liquid, gas] = eos.ThermodynamicsPH(p, h);
      nitrs += eos.get_itrs();
      mitrs = std::max(mitrs, eos.get_itrs());
      count++;
     
      if (prop.x == 0.0 || prop.x == 1.0) {
        T = prop.T;
        Tmin = std::min(Tmin, T);
        Tmax = std::max(Tmax, T);
      } else {
        T = liquid.T;
        Tmin = std::min(Tmin, T);
        Tmax = std::max(Tmax, T);
      }
      // std::cout << p << " " << h << " " << prop.rho << std::endl;
    }
  }
  std::cout << "mean itrs = " << double(nitrs) / count << std::endl;
  std::cout << "max itrs = " << mitrs << std::endl;
  std::cout << "temperature variation = " << Tmin << " " << Tmax << std::endl;
}


TEST(EOS_IAPWS97_REGION4_CONSTANT_P)
{
  using namespace Amanzi::AmanziEOS;
  std::cout << "\nTEST: region4, p=constant" << std::endl;

  Teuchos::ParameterList plist;
  IAPWS97 eos(plist);

  double h, dh(10.0), p(eos.PC / 2), x0(0.0), x, T, Tmin(1000.0), Tmax(100.0);
  for (int i = -20; i < 20; ++i) {
    h = eos.HC + i * dh; 
    auto [prop, liquid, gas] = eos.ThermodynamicsPH(p, h);

    T = liquid.T;
    Tmin = std::min(Tmin, T);
    Tmax = std::max(Tmax, T);

    // liquid content must grow
    x = liquid.x;
    CHECK(x > x0);
    x0 = x;
    // std::cout << "dh = " << h - eos.HC << " T = " << liquid.T << " x = " << liquid.x << std::endl;
  }
  CHECK_CLOSE(Tmin, Tmax, 1e-6);
  std::cout << "temperature variation along p = constant: " << Tmin << " " << Tmax << std::endl;
}


TEST(EOS_IAPWS97_REGION4_CONSTANT_H)
{
  using namespace Amanzi::AmanziEOS;
  std::cout << "\nTEST: region 4, h=constant" << std::endl;

  Teuchos::ParameterList plist;
  IAPWS97 eos(plist);

  double h(eos.HC - 100.0), p, dp(0.2);
  for (int i = -7; i < 3; ++i) {
    p = eos.PC + i * dp; 
    auto [prop, liquid, gas] = eos.ThermodynamicsPH(p, h);
    if (liquid.T != 0.0) {
      std::cout << "dp=" << p - eos.PC << " T=" << liquid.T << " x=" << liquid.x << " rgn=" << liquid.rgn << std::endl;
    } else {
      std::cout << "dp=" << p - eos.PC << " T=" << gas.T << " x=" << gas.x << " rgn=" << gas.rgn << std::endl;
    }
  }
}


