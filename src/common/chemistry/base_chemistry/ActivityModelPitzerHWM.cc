/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include <cstdlib>
#include <math.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>

#include "Teuchos_XMLParameterListHelpers.hpp"

#include "errors.hh"
#include "exceptions.hh"
#include "VirialCoefficient.hh"

#include "ActivityModelPitzerHWM.hh"

namespace Amanzi {
namespace AmanziChemistry {

const double ActivityModelPitzerHWM::bdh = 1.2;
const double ActivityModelPitzerHWM::cwater = 55.50837;
// -------------------------------------------------------------
// Limiting Debye-Hückel slope to 25º   0.39153  0.392
// -------------------------------------------------------------
const double ActivityModelPitzerHWM::aphi_debye_huckel_slope25 = 0.392;
// -------------------------------------------------------------
// aphi_debye_huckel_slope25=0.39153d0,   & ! Limiting Debye-Hückel slope to 25º   0.39153  0.392
// -------------------------------------------------------------
const double ActivityModelPitzerHWM::c0aphi_debye_huckel_slope = 0.13422;    // Temperature depending coefficients
const double ActivityModelPitzerHWM::c1aphi_debye_huckel_slope = 0.0368329;  // Temperature depending coefficients
const double ActivityModelPitzerHWM::c2aphi_debye_huckel_slope = 14.62718;   // Temperature depending coefficients
const double ActivityModelPitzerHWM::c3aphi_debye_huckel_slope = 1530.1474;  // Temperature depending coefficients
const double ActivityModelPitzerHWM::c4aphi_debye_huckel_slope = 80.40631;   // Temperature depending coefficients
const double ActivityModelPitzerHWM::c5aphi_debye_huckel_slope = 4.1725332;  // Temperature depending coefficients
const double ActivityModelPitzerHWM::c6aphi_debye_huckel_slope = 0.1481291;  // Temperature depending coefficients
const double ActivityModelPitzerHWM::c7aphi_debye_huckel_slope = 1.5188505;  // Temperature depending coefficients
const double ActivityModelPitzerHWM::c8aphi_debye_huckel_slope = 1.8016317;  // Temperature depending coefficients
const double ActivityModelPitzerHWM::c9aphi_debye_huckel_slope = 9.3816144;  // Temperature depending coefficients


/* ******************************************************************
* Create the object
****************************************************************** */
ActivityModelPitzerHWM::ActivityModelPitzerHWM()
  : ActivityModel(),
    aphi_debye_huckel_slope(aphi_debye_huckel_slope25),
    number_b_functions(0),
    number_j_functions(0),
    number_non_zero_beta(0),
    number_non_zero_theta(0),
    number_non_zero_cphi(0),
    number_non_zero_lamda(0),
    number_non_zero_psi(0),
    number_non_zero_q(0),
    index_cl_species(-1),
    index_h2o_species(-1),
    index_k_species(-1),
    number_species(0),
    macinnes_scaled(false),
    jfunction_approach("pitzer1975") {
}


/* ******************************************************************
* Setup the object
****************************************************************** */
void ActivityModelPitzerHWM::Setup(
    const ActivityModelParameters& parameters,
    const std::vector<Species>& primary_species,
    const std::vector<AqueousEquilibriumComplex>& aqueous_complexes)
{
  // initialize storage in the base class:
  ResizeGamma(primary_species.size() + aqueous_complexes.size());

  // Store the J's functions approach name
  jfunction_approach = parameters.pitzer_jfunction;

  // Store constant values for j functions
  // These constants were taken from Table III in Pitzer (1975).
  if (jfunction_approach == "pitzer1975") {
    const_j_functions.push_back(4.118);
    const_j_functions.push_back(7.247);
    const_j_functions.push_back(-4.408);
    const_j_functions.push_back(1.837);
    const_j_functions.push_back(-0.251);
    const_j_functions.push_back(0.0164);
  }

  // collect aquoues species in one array
  for (auto it = primary_species.begin(); it != primary_species.end(); ++it) {
    molality.push_back(0.0);
    charge.push_back((*it).charge());
    name_species.push_back((*it).name());
  }
  for (auto it = aqueous_complexes.begin(); it != aqueous_complexes.end(); ++it) {
    molality.push_back(0.0);
    charge.push_back((*it).charge());
    name_species.push_back((*it).name());
  }
  number_species = primary_species.size() + aqueous_complexes.size();

  auto plist = Teuchos::getParametersFromXmlFile(parameters.database_filename + ".xml");

  ParseVirialCoefficient_(plist, "b0");
  ParseVirialCoefficient_(plist, "b1");
  ParseVirialCoefficient_(plist, "b2");
  ParseVirialCoefficient_(plist, "cfi");
  ParseVirialCoefficient_(plist, "theta");
  ParseVirialCoefficient_(plist, "lamda");
  ParseVirialCoefficient_(plist, "psi");

  // Read Pitzer coefficients database
  AssignIndexBetaFunctions();
  AssignIndexJFunctions();
  number_non_zero_q = number_non_zero_beta + number_non_zero_theta + number_non_zero_lamda;
  PushPrivateVectors();
  Update(273.15, 0.0);

  for (int isp = 0; isp < number_species; isp++) {
    if (name_species.at(isp) == "H2O") index_h2o_species = isp;
    if (name_species.at(isp) == "Cl-") index_cl_species = isp;
    if (name_species.at(isp) == "K+") index_k_species = isp;
  }
  if (index_cl_species > -1 && index_k_species > -1) {
    macinnes_scaled = true;
  }
}


/* ******************************************************************
* Compute the activity coefficient
****************************************************************** */
double ActivityModelPitzerHWM::Evaluate(const Species& species) {
  return 1.0;
}


/* ******************************************************************
* Compute the activity coefficients
****************************************************************** */
void ActivityModelPitzerHWM::EvaluateVector(
    const std::vector<Species>& primary_species,
    const std::vector<AqueousEquilibriumComplex>& aqueous_complexes,
    std::vector<double>* gamma,
    double* water_activity)
{
  int number_species_(primary_species.size() + aqueous_complexes.size());

  // Check the number of species
  if (number_species_ != number_species) {
    std::ostringstream error_stream;
    error_stream << "Error, different number of aqueous species" << "\n";
    Exceptions::amanzi_throw(Errors::Message(error_stream.str()));
  }
  double gcl(1.0);
  double gclm(1.0);
  double osmotic_coefficient(1.0);
  CalculateIonicStrength(primary_species, aqueous_complexes);
  CalculateSumAbsZ(primary_species, aqueous_complexes);
  CalculateSumC(primary_species, aqueous_complexes);

  if (I_ == 0.0 || Z_ == 0.0 || M_ == 0.0) {
    std::ostringstream error_stream;
    error_stream << "Error, zero concentrations" << "\n";
    Exceptions::amanzi_throw(Errors::Message(error_stream.str()));
  }
  int isp(-1);
  for (auto it = primary_species.begin(); it != primary_species.end(); ++it) {
    isp++;
    molality.at(isp) = (*it).molality();
  }
  for (auto it = aqueous_complexes.begin(); it != aqueous_complexes.end(); ++it) {
    isp++;
    molality.at(isp) = (*it).molality();
  }

  ComputeQmatrices();
  ComputeDebyeHuckelTerm(*gamma, osmotic_coefficient, gclm);
  ComputemQmProduct(*gamma, osmotic_coefficient);
  ComputemQlmProduct(osmotic_coefficient);
  ComputemQcmProduct(*gamma, osmotic_coefficient);
  ComputemTmmProduct(*gamma, osmotic_coefficient);

  if (macinnes_scaled && index_cl_species > -1) {
    gcl = gamma->at(index_cl_species);
    for (int i = 0; i < number_species; i++)
      if (i != index_h2o_species) {
        gamma->at(i) *= pow((gcl / gclm), charge.at(i));
      }
  }
  *water_activity = osmotic_coefficient;
}


/* ******************************************************************
* Compute the m(Qc)m product.
****************************************************************** */
void ActivityModelPitzerHWM::ComputemQcmProduct(std::vector<double>& gamma,
                                                double& osmotic_coefficient)
{
  double qc(0.0);
  for (int nz = 0; nz < number_non_zero_cphi; nz++) {
    int i(cphi_virial.at(nz).GetIsp1());
    int j(cphi_virial.at(nz).GetIsp2());
    double mi(molality.at(i));
    double mj(molality.at(j));
    double Cij(cphi_virial.at(nz).GetVirial());
    double mjCij(mj * Cij);
    double miCij(mi * Cij);
    double mimjCij(mi * mj * Cij);
    qc += mimjCij;
    gamma.at(i) += mjCij * Z_;
    gamma.at(j) += miCij * Z_;
  }
  osmotic_coefficient += Z_ * qc;
  for (int i = 0; i < gamma.size(); i++) {
    gamma.at(i) += std::abs(charge.at(i)) * qc;
  }
}


/* ******************************************************************
* Compute the m(Ql)m product.
****************************************************************** */
void ActivityModelPitzerHWM::ComputemQlmProduct(double& osmotic_coefficient)
{
  for (int nz = 0; nz < number_non_zero_lamda; nz++) {
    int i(lamda_virial.at(nz).GetIsp1());
    int j(lamda_virial.at(nz).GetIsp2());
    double Lij(lamda_virial.at(nz).GetVirial());
    double mi(molality.at(i));
    double mj(molality.at(j));
    osmotic_coefficient += mi * mj * Lij;
  }
}


/* ******************************************************************
* Compute mQm product.
****************************************************************** */
void ActivityModelPitzerHWM::ComputemQmProduct(std::vector<double>& gamma,
                                               double& osmotic_coefficient)
{
  double q2phi(0.0), q2prim(0.0);
  for (int nz = 0; nz < number_non_zero_q; nz++) {
    double qij(q_matrix.at(nz));
    double qpriij(qpri_matrix.at(nz));
    double qphiij(qphi_matrix.at(nz));
    int i(index_non_zero_q.at(0).at(nz));
    int j(index_non_zero_q.at(1).at(nz));
    double mi(molality.at(i));
    double mj(molality.at(j));
    q2phi += mi * mj * qphiij;
    q2prim += mi * mj * qpriij;
    gamma.at(i) += 2.0 * qij * mj;
    gamma.at(j) += 2.0 * qij * mi;
  }
  osmotic_coefficient += q2phi;
  for (int i = 0; i < gamma.size(); i++) {
    gamma.at(i) += charge.at(i) * charge.at(i) * q2prim;
  }
}


/* ******************************************************************
* Compute mTmm product.
****************************************************************** */
void ActivityModelPitzerHWM::ComputemTmmProduct(std::vector<double>& gamma,
                                                double& osmotic_coefficient)
{
  double tril(0.0);
  std::vector<double> vector;
  for (int i = 0; i < number_species; i++) {
    vector.push_back(0.0);
  }
  for (int nz = 0; nz < number_non_zero_psi; nz++) {
    double Psiijk(psi_virial.at(nz).GetVirial());
    int i(psi_virial.at(nz).GetIsp1());
    int j(psi_virial.at(nz).GetIsp2());
    int k(psi_virial.at(nz).GetIsp3());
    double mi(molality.at(i));
    double mj(molality.at(j));
    double mk(molality.at(k));
    tril += mi * mj * mk * Psiijk;
    vector.at(i) += Psiijk * mj * mk;
    vector.at(j) += Psiijk * mi * mk;
    vector.at(k) += Psiijk * mj * mi;
  }
  osmotic_coefficient += tril;
  osmotic_coefficient *= 2.0 / (M_ + 1.0e-20);
  osmotic_coefficient += 1.0;
  osmotic_coefficient *= -M_ / cwater;
  osmotic_coefficient = exp(osmotic_coefficient);
  for (int i = 0; i < number_species; i++) {
    if (i == index_h2o_species) {
      gamma.at(index_h2o_species) = osmotic_coefficient;
    } else {
      gamma.at(i) = exp(gamma.at(i) + vector.at(i));
    }
  }
}


/* ******************************************************************
* Compute the Debye-Huckel term.
****************************************************************** */
void ActivityModelPitzerHWM::ComputeDebyeHuckelTerm(std::vector<double>& gamma,
                                                    double& osmotic_coefficient,
                                                    double& gclm)
{
  gclm = 1.0;
  double den(1.0 + bdh * sqrt(I_));
  double dh(-aphi_debye_huckel_slope * sqrt(I_) / den);
  osmotic_coefficient = I_ * dh;
  dh -= (2.0 / bdh) * aphi_debye_huckel_slope * log(den);
  int isp(-1);
  for (std::vector<double>::iterator i = gamma.begin(); i != gamma.end(); i++) {
    isp++;
    (*i) = charge.at(isp) * charge.at(isp) * dh;
  }
  if (macinnes_scaled) {
    gclm = gclm_(dh);
  }
}


/* ******************************************************************
* Compute the activity coefficient for Cl- in a KCl-K system
*  (performed for the MacIness convention).
****************************************************************** */
double ActivityModelPitzerHWM::gclm_(const double& dhterm)
{
  const double mtb0kcl(0.04835);
  const double mtb1kcl(0.2122);
  const double mtc0kcl(-0.00084);
  double x(2.0 * sqrt(I_));
  double xxx(-2.0 * (1.0 - (1.0 + x + 0.5 * x * x)*exp(-x)) / (x * x));
  xxx *= mtb1kcl / I_;
  double yyy(2.0 * (1.0 - (1.0 + x)*exp(-x)) / (x * x));
  yyy *= mtb1kcl;
  yyy += mtb0kcl;
  return exp(dhterm + I_ * I_ * xxx + I_ * (2.0 * yyy + I_ * mtc0kcl) + I_ * I_ * mtc0kcl / 2.0);
}


/* ******************************************************************
* Compute the Q's matrices and store them in sparse storage
****************************************************************** */
void ActivityModelPitzerHWM::ComputeQmatrices()
{
  for (int i = 0; i < number_b_functions; i++) {
    g_function.at(i).at(0) = 0.0;
    g_function.at(i).at(1) = 0.0;
    g_pri_function.at(i).at(0) = 0.0;
    g_pri_function.at(i).at(1) = 0.0;
    f_function.at(i).at(0) = 0.0;
    f_function.at(i).at(1) = 0.0;
  }
  for (int i = 0; i < number_j_functions; i++) {
    j_function.at(i) = 0.0;
    j_pri_function.at(i) = 0.0;
  }
  ComputeBetaFunctions_();
  ComputeJFunctions();
  int nz(-1);
  for (int nz_loc = 0; nz_loc < number_non_zero_beta; nz_loc++) {
    nz++;
    int i(beta0_virial.at(nz_loc).GetIsp1());
    int j(beta0_virial.at(nz_loc).GetIsp2());
    int k(beta0_virial.at(nz_loc).GetIfun1());
    double B0ij(beta0_virial.at(nz_loc).GetVirial());
    double B1ij(beta1_virial.at(nz_loc).GetVirial());
    double B2ij(beta2_virial.at(nz_loc).GetVirial());
    double expo1(f_function.at(k).at(0));
    double expo2(f_function.at(k).at(1));
    double G1(g_function.at(k).at(0));
    double G2(g_function.at(k).at(1));
    double Gp1(g_pri_function.at(k).at(0));
    double Gp2(g_pri_function.at(k).at(1));
    qphi_matrix.at(nz) = B0ij + B1ij * expo1 + B2ij * expo2;
    q_matrix.at(nz) = B0ij + B1ij * G1 + B2ij * G2;
    qpri_matrix.at(nz) = B1ij * (Gp1 / I_) + B2ij * (Gp2 / I_);
    index_non_zero_q.at(0).at(nz) = i;
    index_non_zero_q.at(1).at(nz) = j;
  }

  for (int nz_loc = 0; nz_loc < number_non_zero_theta; nz_loc++) {
    nz++;
    int i(theta_virial.at(nz_loc).GetIsp1());
    int j(theta_virial.at(nz_loc).GetIsp2());
    int funij(theta_virial.at(nz_loc).GetIfun1());
    int funii(theta_virial.at(nz_loc).GetIfun2());
    int funjj(theta_virial.at(nz_loc).GetIfun3());
    double thij(theta_virial.at(nz_loc).GetVirial());
    double jij(j_function.at(funij));
    double jii(j_function.at(funii));
    double jjj(j_function.at(funjj));
    double jpij(j_pri_function.at(funij));
    double jpii(j_pri_function.at(funii));
    double jpjj(j_pri_function.at(funjj));
    double zizj(charge_product.at(funij));
    double zizi(charge_product.at(funii));
    double zjzj(charge_product.at(funjj));
    double xij(2.352 * sqrt(I_)*zizj);
    double xii(2.352 * sqrt(I_)*zizi);
    double xjj(2.352 * sqrt(I_)*zjzj);
    double eth((zizj / (4.0 * I_)) * (jij - 0.5 * jii - 0.5 * jjj));
    q_matrix.at(nz) = thij + eth;
    double eth_i(eth / I_);
    double ethpri(-eth_i + (zizj / (8.0 * I_ * I_)) * (xij * jpij - 0.5 * xii * jpii - 0.5 * xjj * jpjj));
    qpri_matrix.at(nz) = ethpri;
    qphi_matrix.at(nz) = thij + eth + I_ * ethpri;
    index_non_zero_q.at(0).at(nz) = i;
    index_non_zero_q.at(1).at(nz) = j;
  }
  for (int nz_loc = 0; nz_loc < number_non_zero_lamda; nz_loc++) {
    nz++;
    int i(lamda_virial.at(nz_loc).GetIsp1());
    int j(lamda_virial.at(nz_loc).GetIsp2());
    q_matrix.at(nz) = lamda_virial.at(nz_loc).GetVirial();
    index_non_zero_q.at(0).at(nz) = i;
    index_non_zero_q.at(1).at(nz) = j;
  }

  if ((nz + 1) != number_non_zero_q) {
    std::ostringstream error_stream;
    error_stream << "Warning different number non-zero terms in Q's matrices\n";
    Exceptions::amanzi_throw(Errors::Message(error_stream.str()));
  }
}


/* ******************************************************************
* Compute the Beta functions
****************************************************************** */
void ActivityModelPitzerHWM::ComputeBetaFunctions_()
{
  for (int j = 0; j < number_b_functions; j++) {
    double x1(alpha1.at(j)*sqrt(I_));
    double x1q(x1 * x1);
    f_function.at(j).at(0) = exp(-x1);
    g_function.at(j).at(0) = 2.0 * (1.0 - (1.0 + x1) * exp(-x1)) / x1q;
    g_pri_function.at(j).at(0) = -2.0 * (1.0 - (1.0 + x1 + (x1q / 2.0)) * exp(-x1)) / x1q;
    if (alpha2.at(j) != 0.0) {
      double x2(alpha2.at(j)*sqrt(I_));
      double x2q(x2 * x2);
      f_function.at(j).at(1) = exp(-x2);
      g_function.at(j).at(1) = 2.0 * (1.0 - (1.0 + x2) * exp(-x2)) / x2q;
      g_pri_function.at(j).at(1) = -2.0 * (1.0 - (1.0 + x2 + (x2q / 2.0)) * exp(-x2)) / x2q;
    }
  }
}


/* ******************************************************************
* Compute the J's functions.
****************************************************************** */
void ActivityModelPitzerHWM::ComputeJFunctions()
{
  const double e1(4.581), e2(0.7237), e3(0.012), e4(0.528);  // e12(7.8963);

  if (jfunction_approach == "pitzer1975") {
    for (int i = 0; i < number_j_functions; i++) {
      double zizj(charge_product.at(i));
      double x(2.352 * sqrt(I_)*zizj);
      double x2(x * x);
      if (x <= 0.03) {
        double s1(const_j_functions.at(5) / x);
        double s3(6.0 * s1);
        for (int k = 4; k >= 0; k--) {
          s1 = (s1 + const_j_functions.at(k)) / x;
          s3 = (s3 + k * const_j_functions.at(k)) / x;
        }
        s3 = s3 / x;
        double s1q(s1 * s1);
        j_function.at(i) = -(1.0 / 6.0) * x2 * log(x) * exp(-10.0 * x2) + (1.0 / s1);
        j_pri_function.at(i) = ((10.0 * x2 - 1.0) * log(x) - 0.5) * (x / 3.0) * exp(-10.0 * x2) + (s3 / s1q);

      } else {
        double xc4(pow(x, e4));
        double xc2(pow(x, -e2));
        double td1(e1 * xc2 * exp(-e3 * xc4));
        double td(4.0 + td1);
        j_function.at(i) = x / td;
        j_pri_function.at(i) = (j_function.at(i) / x2) * (x + td1 * (e2 + e3 * e4 * xc4) * j_function.at(i));
      }
    }

  } else {
    std::ostringstream error_stream;
    error_stream << "Name for the J's functions approach not recognized" << "\n";
    Exceptions::amanzi_throw(Errors::Message(error_stream.str()));
  }
}


/* ******************************************************************
* Write the attributes
****************************************************************** */
void ActivityModelPitzerHWM::Display() const {
  std::cout << "============================================>" << std::endl;
  std::cout << "Activity model: HWM (Harvie et al., 1984)" << std::endl;
  std::cout << "Species:" << std::endl;

  for (int i = 0; i < number_species; i++) {
    std::cout << name_species.at(i) << " " << charge.at(i) << " " << std::endl;
  }
  std::cout << "--------------------------------------------------------------------" << std::endl;
  std::cout << " Virial coefficients" << std::endl;
  std::cout << "--------------------------------------------------------------------" << std::endl;
  int nvirial(0);

  if (number_non_zero_beta > 0) {
    std::cout << "=================> Beta0 ==============>" << std::endl;
    for (int i = 0; i < number_non_zero_beta; i++) {
      int isp1 = beta0_virial.at(i).GetIsp1();
      int isp2 = beta0_virial.at(i).GetIsp2();
      if (beta0_virial.at(i).GetVirial() != 0.0) {
        nvirial++;
        std::cout << name_species.at(isp1) << "  " << name_species.at(isp2) << "  " << beta0_virial.at(i).GetVirial() << std::endl;
      }
    }

    std::cout << "=================> Beta1 ==============>" << std::endl;
    for (int i = 0; i < number_non_zero_beta; i++) {
      int isp1 = beta1_virial.at(i).GetIsp1();
      int isp2 = beta1_virial.at(i).GetIsp2();
      if (beta1_virial.at(i).GetVirial() != 0.0) {
        nvirial++;
        std::cout << name_species.at(isp1) << "  " << name_species.at(isp2) << "  " << beta1_virial.at(i).GetVirial() << std::endl;
      }
    }

    std::cout << "=================> Beta2 ==============>" << std::endl;
    for (int i = 0; i < number_non_zero_beta; i++) {
      int isp1 = beta2_virial.at(i).GetIsp1();
      int isp2 = beta2_virial.at(i).GetIsp2();
      if (beta2_virial.at(i).GetVirial() != 0.0) {
        nvirial++;
        std::cout << name_species.at(isp1) << "  " << name_species.at(isp2) << "  " << beta2_virial.at(i).GetVirial() << std::endl;
      }
    }
  }

  if (number_non_zero_cphi > 0) {
    std::cout << "=================> Cphi ==============>" << std::endl;
    for (int i = 0; i < number_non_zero_cphi; i++) {
      int isp1 = cphi_virial.at(i).GetIsp1();
      int isp2 = cphi_virial.at(i).GetIsp2();
      if (cphi_virial.at(i).GetVirial() != 0.0) {
        nvirial++;
        std::cout << name_species.at(isp1) << "  " << name_species.at(isp2) << "  " << cphi_virial.at(i).GetVirial() << std::endl;
      }
    }
  }

  if (number_non_zero_theta > 0) {
    std::cout << "=================> Theta ==============>" << std::endl;
    for (int i = 0; i < number_non_zero_theta; i++) {
      int isp1 = theta_virial.at(i).GetIsp1();
      int isp2 = theta_virial.at(i).GetIsp2();;
      if (theta_virial.at(i).GetVirial() != 0.0) {
        nvirial++;
        std::cout << name_species.at(isp1) << "  " << name_species.at(isp2) << "  " << theta_virial.at(i).GetVirial() << std::endl;
      }
    }
  }

  if (number_non_zero_lamda > 0) {
    std::cout << "=================> Lamda ==============>" << std::endl;
    for (int i = 0; i < number_non_zero_lamda; i++) {
      int isp1 = lamda_virial.at(i).GetIsp1();
      int isp2 = lamda_virial.at(i).GetIsp2();
      if (lamda_virial.at(i).GetVirial() != 0.0) {
        nvirial++;
        std::cout << name_species.at(isp1) << "  " << name_species.at(isp2) << "  " << lamda_virial.at(i).GetVirial() << std::endl;
      }
    }
  }

  if (number_non_zero_psi > 0) {
    std::cout << "=================> Psi ==============>" << std::endl;
    for (int i = 0; i < number_non_zero_psi; i++) {
      int isp1 = psi_virial.at(i).GetIsp1();
      int isp2 = psi_virial.at(i).GetIsp2();
      int isp3 = psi_virial.at(i).GetIsp3();
      if (psi_virial.at(i).GetVirial() != 0.0) {
        nvirial++;
        std::cout << name_species.at(isp1) << "  " << name_species.at(isp2) << "  " 
                  << name_species.at(isp3) << "  " << psi_virial.at(i).GetVirial() << std::endl;
      }
    }
  }

  std::cout << "=====================================>" << std::endl;
  std::cout << "Total number of virial coefficients: " << nvirial << std::endl;
  std::cout << "Total number of Beta's functions: " << number_b_functions << std::endl;
  std::cout << "Total number of J's functions: " << number_j_functions << std::endl;
  std::cout << "--------------------------------------" << std::endl;
  for (int i = 0; i < number_j_functions; i++) {
    std::cout << "Zi Zj product:" << charge_product[i] << std::endl;
  }
  std::cout << "=====================================>" << std::endl;
}


/* ******************************************************************
* Set virial coefficients
****************************************************************** */
void ActivityModelPitzerHWM::SetVirialCoefficient(
     const std::vector<double>& virial, const std::string& typevirial,
     const int& isp1, const int& isp2, const int& isp3)
{
  VirialCoefficient vir;
  VirialCoefficient vir0;
  bool not_stored(true);
  if (typevirial == "b0") {
    for (int i = 0; i < number_non_zero_beta && not_stored; i++) {
      int isp1_old = beta0_virial.at(i).GetIsp1();
      int isp2_old = beta0_virial.at(i).GetIsp2();
      if ((isp1 == isp1_old && isp2 == isp2_old) || (isp1 == isp2_old && isp2 == isp1_old)) {
        for (std::vector<double>::const_iterator j = virial.begin(); j != virial.end(); j++) {
          beta0_virial.at(i).SetPol((*j));
        }
        not_stored = false;
      }
    }
    if (not_stored) {
      number_non_zero_beta++;
      for (std::vector<double>::const_iterator j = virial.begin(); j != virial.end(); j++) {
        vir.SetPol((*j));
      }
      vir.SetIsp1(isp1);
      vir.SetIsp2(isp2);
      vir0.SetIsp1(isp1);
      vir0.SetIsp2(isp2);
      vir.SetIfun1(0);
      vir0.SetIfun1(0);
      beta0_virial.push_back(vir);
      beta1_virial.push_back(vir0);
      beta2_virial.push_back(vir0);
    }

  } else if (typevirial == "b1") {

    for (int i = 0; i < number_non_zero_beta && not_stored; i++) {
      int isp1_old = beta1_virial.at(i).GetIsp1();
      int isp2_old = beta1_virial.at(i).GetIsp2();;
      if ((isp1 == isp1_old && isp2 == isp2_old) || (isp1 == isp2_old && isp2 == isp1_old)) {
        for (std::vector<double>::const_iterator j = virial.begin(); j != virial.end(); j++) {
          beta1_virial.at(i).SetPol((*j));
        }
        not_stored = false;
      }
    }
    if (not_stored) {
      number_non_zero_beta++;
      for (std::vector<double>::const_iterator j = virial.begin(); j != virial.end(); j++) {
        vir.SetPol((*j));
      }
      vir.SetIsp1(isp1);
      vir.SetIsp2(isp2);
      vir0.SetIsp1(isp1);
      vir0.SetIsp2(isp2);
      vir.SetIfun1(0);
      vir0.SetIfun1(0);
      beta1_virial.push_back(vir);
      beta0_virial.push_back(vir0);
      beta2_virial.push_back(vir0);
    }

  } else if (typevirial == "b2") {


    for (int i = 0; i < number_non_zero_beta && not_stored; i++) {
      int isp1_old = beta2_virial.at(i).GetIsp1();
      int isp2_old = beta2_virial.at(i).GetIsp1();
      if ((isp1 == isp1_old && isp2 == isp2_old) || (isp1 == isp2_old && isp2 == isp1_old)) {
        for (std::vector<double>::const_iterator j = virial.begin(); j != virial.end(); j++) {
          beta2_virial.at(i).SetPol((*j));
        }
        not_stored = false;
      }
    }
    if (not_stored) {
      number_non_zero_beta++;
      for (std::vector<double>::const_iterator j = virial.begin(); j != virial.end(); j++) {
        vir.SetPol((*j));
      }
      vir.SetIsp1(isp1);
      vir.SetIsp2(isp2);
      vir0.SetIsp1(isp1);
      vir0.SetIsp2(isp2);
      vir.SetIfun1(0);
      vir0.SetIfun1(0);
      beta2_virial.push_back(vir);
      beta0_virial.push_back(vir0);
      beta1_virial.push_back(vir0);
    }
  } else if (typevirial == "cfi") {

    number_non_zero_cphi++;
    for (std::vector<double>::const_iterator j = virial.begin(); j != virial.end(); j++) {
      vir.SetPol((*j));
    }
    vir.SetIsp1(isp1);
    vir.SetIsp2(isp2);
    cphi_virial.push_back(vir);

  } else if (typevirial == "theta") {

    number_non_zero_theta++;
    for (std::vector<double>::const_iterator j = virial.begin(); j != virial.end(); j++) {
      vir.SetPol((*j));
    }
    vir.SetIsp1(isp1);
    vir.SetIsp2(isp2);
    vir.SetIfun1(int(charge.at(isp1)*charge.at(isp2)));
    vir.SetIfun2(int(charge.at(isp1)*charge.at(isp1)));
    vir.SetIfun3(int(charge.at(isp2)*charge.at(isp2)));
    theta_virial.push_back(vir);


  } else if (typevirial == "lamda") {

    number_non_zero_lamda++;
    for (auto j = virial.begin(); j != virial.end(); j++) {
      vir.SetPol((*j));
    }
    vir.SetIsp1(isp1);
    vir.SetIsp2(isp2);
    lamda_virial.push_back(vir);

  } else if (typevirial == "psi") {

    number_non_zero_psi++;
    for (std::vector<double>::const_iterator j = virial.begin(); j != virial.end(); j++) {
      vir.SetPol((*j));
    }
    vir.SetIsp1(isp1);
    vir.SetIsp2(isp2);
    vir.SetIsp3(isp3);
    psi_virial.push_back(vir);

  } else {
    std::ostringstream error_stream;
    error_stream << "Type virial coefficient not defined" << typevirial;
    Exceptions::amanzi_throw(Errors::Message(error_stream.str()));
  }
}


/* ******************************************************************
* Assign Beta's functions
****************************************************************** */
void ActivityModelPitzerHWM::AssignIndexBetaFunctions()
{
  // Local variables and constants
  int l1(0), l2(0), l3(0), n1(0), n2(0), n3(0), nz(0);
  number_b_functions = 0;
  for (nz = 0; nz < number_non_zero_beta; nz++) {
    double z1(charge.at(beta0_virial.at(nz).GetIsp1()));
    double z2(charge.at(beta0_virial.at(nz).GetIsp2()));
    if (std::abs(z1) == 1.0 || std::abs(z2) == 1.0) {
      if (l1 == 0) {
        number_b_functions++;
        beta0_virial.at(nz).SetIfun1(number_b_functions - 1);
        beta1_virial.at(nz).SetIfun1(number_b_functions - 1);
        beta2_virial.at(nz).SetIfun1(number_b_functions - 1);
        n1 = number_b_functions - 1;
        alpha1.push_back(2.0);
        alpha2.push_back(12.0);
        l1 = 1;
      } else {
        beta0_virial.at(nz).SetIfun1(n1);
        beta1_virial.at(nz).SetIfun1(n1);
        beta2_virial.at(nz).SetIfun1(n1);
      }
    } else {
      if (std::abs(z1) != std::abs(z2)) {
        if (l2 == 0) {
          number_b_functions++;
          beta0_virial.at(nz).SetIfun1(number_b_functions - 1);
          beta1_virial.at(nz).SetIfun1(number_b_functions - 1);
          beta2_virial.at(nz).SetIfun1(number_b_functions - 1);
          n2 = number_b_functions - 1;
          alpha1.push_back(2.0);
          alpha2.push_back(50.0);
          l2 = 1;
        } else {
          beta0_virial.at(nz).SetIfun1(n2);
          beta1_virial.at(nz).SetIfun1(n2);
          beta2_virial.at(nz).SetIfun1(n2);
        }
      } else {
        if (l3 == 0) {
          number_b_functions++;
          beta0_virial.at(nz).SetIfun1(number_b_functions - 1);
          beta1_virial.at(nz).SetIfun1(number_b_functions - 1);
          beta2_virial.at(nz).SetIfun1(number_b_functions - 1);
          n3 = number_b_functions - 1;
          alpha1.push_back(1.4);
          alpha2.push_back(12.0);
          l3 = 1;
        } else {
          beta0_virial.at(nz).SetIfun1(n3);
          beta1_virial.at(nz).SetIfun1(n3);
          beta2_virial.at(nz).SetIfun1(n3);
        }
      }
    }
  }
}


/* ******************************************************************
* Assign J's functions
****************************************************************** */
void ActivityModelPitzerHWM::AssignIndexJFunctions()
{
  VirialCoefficient vir;
  for (int i = 0; i < number_species; i++) {
    for (int j = 0; j < number_species; j++) {
      if (i != j && ((charge.at(i) > 0.0 && charge.at(j) > 0.0) || (charge.at(i) < 0.0 && charge.at(j) < 0.0))) {
        bool not_found(true);
        for (int k = 0; k < number_non_zero_theta && not_found; k++) {
          int isp1(theta_virial.at(k).GetIsp1());
          int isp2(theta_virial.at(k).GetIsp2());
          if ((isp1 == i && isp2 == j) || (isp1 == j && isp2 == i)) {
            not_found = false;
          }
        }
        if (not_found) {
          number_non_zero_theta++;
          vir.SetPol(0.0);
          vir.SetIsp1(i);
          vir.SetIsp2(j);
          vir.SetIfun1(int(charge.at(i)*charge.at(j)));
          vir.SetIfun2(int(charge.at(i)*charge.at(i)));
          vir.SetIfun3(int(charge.at(j)*charge.at(j)));
          theta_virial.push_back(vir);
        }
      }
    }
  }

  // Compute number of functions j and save
  if (number_non_zero_theta > 0) {
    number_j_functions = 1;
    charge_product.push_back(double(theta_virial.at(0).GetIfun1()));
    double zz(0.0);
    bool not_found(true);
    for (int nz = 0; nz < number_non_zero_theta; nz++) {
      not_found = true;
      for (int j = 0; j < number_j_functions && not_found; j++) {
        zz = double(theta_virial.at(nz).GetIfun1());
        if (charge_product.at(j) == zz) {
          theta_virial.at(nz).SetIfun1(j);
          not_found = false;
        }
      }
      if (not_found) {
        charge_product.push_back(zz);
        number_j_functions++;
        theta_virial.at(nz).SetIfun1(number_j_functions - 1);
      }
      not_found = true;
      for (int j = 0; j < number_j_functions && not_found; j++) {
        zz = double(theta_virial.at(nz).GetIfun2());
        if (charge_product.at(j) == zz) {
          theta_virial.at(nz).SetIfun2(j);
          not_found = false;
        }
      }
      if (not_found) {
        charge_product.push_back(zz);
        number_j_functions++;
        theta_virial.at(nz).SetIfun2(number_j_functions - 1);
      }
      not_found = true;
      for (int j = 0; j < number_j_functions && not_found; j++) {
        zz = double(theta_virial.at(nz).GetIfun3());
        if (charge_product.at(j) == zz) {
          theta_virial.at(nz).SetIfun3(j);
          not_found = false;
        }
      }
      if (not_found) {
        charge_product.push_back(zz);
        number_j_functions++;
        theta_virial.at(nz).SetIfun3(number_j_functions - 1);
      }
    }
  }
}


/* ******************************************************************
* Push back private vectors.
****************************************************************** */
void ActivityModelPitzerHWM::PushPrivateVectors()
{
  g_function.resize(number_b_functions);
  g_pri_function.resize(number_b_functions);
  g_function.resize(number_b_functions);
  f_function.resize(number_b_functions);
  for (int i = 0; i < number_b_functions; i++) {
    g_function.at(i).push_back(0.0);
    g_function.at(i).push_back(0.0);
    g_pri_function.at(i).push_back(0.0);
    g_pri_function.at(i).push_back(0.0);
    f_function.at(i).push_back(0.0);
    f_function.at(i).push_back(0.0);
  }
  for (int i = 0; i < number_j_functions; i++) {
    j_function.push_back(0.0);
    j_pri_function.push_back(0.0);
  }
  index_non_zero_q.resize(2);
  for (int i = 0; i < number_non_zero_q; i++) {
    index_non_zero_q.at(0).push_back(0);
    index_non_zero_q.at(1).push_back(0);
    q_matrix.push_back(0.0);
    qpri_matrix.push_back(0.0);
    qphi_matrix.push_back(0.0);
  }
}


/* ******************************************************************
* Update virial coefficients with temperature and liquid pressure.
****************************************************************** */
void ActivityModelPitzerHWM::Update(const double& temperature,
                                    const double& pressure)
{
  for (auto i = beta0_virial.begin(); i != beta0_virial.end(); i++) {
    i->UpdateVirial(temperature, pressure);
  }
  for (auto i = beta1_virial.begin(); i != beta1_virial.end(); i++) {
    i->UpdateVirial(temperature, pressure);
  }
  for (auto i = beta2_virial.begin(); i != beta2_virial.end(); i++) {
    i->UpdateVirial(temperature, pressure);
  }
  for (auto i = cphi_virial.begin(); i != cphi_virial.end(); i++) {
    i->UpdateVirial(temperature, pressure);
  }
  for (auto i = theta_virial.begin(); i != theta_virial.end(); i++) {
    i->UpdateVirial(temperature, pressure);
  }
  for (auto i = psi_virial.begin(); i != psi_virial.end(); i++) {
    i->UpdateVirial(temperature, pressure);
  }
  for (auto i = lamda_virial.begin(); i != lamda_virial.end(); i++) {
    i->UpdateVirial(temperature, pressure);
  }
}


/* ******************************************************************
* Search in an vector
****************************************************************** */
void ActivityModelPitzerHWM::ParseVirialCoefficient_(
   Teuchos::RCP<Teuchos::ParameterList> plist, const std::string& block_name)
{
  int isp1, isp2, isp3;
  double virial;
  std::string name;

  const auto& tmp = plist->sublist(block_name);

  for (auto it = tmp.begin(); it != tmp.end(); ++it) {
    std::istringstream iss(tmp.get<std::string>(it->first));
    iss >> name; 
    isp1 = GetIndexSpeciesFromName(name);
    if (isp1 > -1) {
      iss >> name; 
      isp2 = GetIndexSpeciesFromName(name);
    }

    isp3 = -1;
    if (block_name == "psi") {
      if (isp1 > -1 && isp2 > -1) {
        iss >> name; 
        isp3 = GetIndexSpeciesFromName(name);
      }
      if (isp3 == -1) continue;
    }

    if (isp1 > -1 && isp2 > -1) {
      iss >> virial;
      if (block_name == "cfi")
        virial /= 2.0 * sqrt(std::abs(charge.at(isp1) * charge.at(isp2)));
      std::vector<double> virial_vec;
      virial_vec.push_back(virial);
      SetVirialCoefficient(virial_vec, block_name, isp1, isp2, isp3);
    }
  }
}


/* ******************************************************************
* search in an vector
****************************************************************** */
int ActivityModelPitzerHWM::GetIndexSpeciesFromName(const std::string& name)
{
  for (int i = 0; i < name_species.size(); i++) {
    if (name_species.at(i) == name) return i;
  }
  return -1;
}

}  // namespace AmanziChemistry
}  // namespace Amanzi
