/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ben Andre
*/

#ifndef AMANZI_CHEMISTRY_ACTIVITY_MODEL_PITZER_HWM_HH_
#define AMANZI_CHEMISTRY_ACTIVITY_MODEL_PITZER_HWM_HH_

#include <vector>
#include <string>
#include <cstdlib>
#include <math.h>
// Base class for activity calculations
#include "activity_model.hh"
#include <virial_coefficient.hh>


// forward declarations
class Epetra_MultiVector;
class Epetra_Vector;
class Epetra_SerialDenseVector;

namespace Amanzi {
namespace AmanziChemistry {

class Species;

class VirialCoefficient;

class ActivityModelPitzerHWM : public ActivityModel {
 public:
  ActivityModelPitzerHWM();
  ~ActivityModelPitzerHWM() {};

  void Setup(const ActivityModelParameters& parameters,
             const std::vector<Species>& primary_species,
             const std::vector<AqueousEquilibriumComplex>& secondary_species);

  double Evaluate(const Species& species);
  void EvaluateVector(const std::vector<Species>& prim, 
                      const std::vector<AqueousEquilibriumComplex>& sec,
                      std::vector<double>* gamma,
                      double* actw);
  void Display() const;

 private:
  void ReadDataBase(const std::string& database,
		            const std::vector<Species>& primary_species,
		            const std::vector<AqueousEquilibriumComplex>& aqueous_complexes);
  void ParseBeta0VirialCoefficient(const std::string& data);
  void ParseBeta1VirialCoefficient(const std::string& data);
  void ParseBeta2VirialCoefficient(const std::string& data);
  void ParseCfiVirialCoefficient(const std::string& data);
  void ParseThetaVirialCoefficient(const std::string& data);
  void ParseLamdaVirialCoefficient(const std::string& data);
  void ParsePsiVirialCoefficient(const std::string& data);
  void AssignIndexBetaFunctions();
  void AssignIndexJFunctions();
  void ComputemQmProduct(std::vector<double>& gamma, double& osmotic_coefficient);
  void ComputemQlmProduct(double& osmotic_coefficient);
  void ComputemQcmProduct(std::vector<double>& gamma, double& osmotic_coefficient);
  void ComputemTmmProduct(std::vector<double>& gamma, double& osmotic_coefficient);
  void ComputeQmatrices();
  void ComputeBetaFunctions();
  void ComputeJFunctions();
  void ComputeDebyeHuckelTerm(std::vector<double>& gamma, double& osmotic_coefficient, double& gclm);
  double gclm_(const double& dhterm);
  void PushPrivateVectors();
  void Update(const double& temperature, const double& pressure);
  void SetVirialCoefficient(const std::vector<double>& virial, const std::string& typevirial,
    		                const int& isp1, const int& isp2, const int& isp3);
  int GetIndexSpeciesFromName(const std::string& name_species);

  static const double cwater;
  static const double bdh;
  static const double aphi_debye_huckel_slope25;
  //-------------------------------------------------------------
  // Limiting Debye-Hückel slope to 25º   0.39153  0.392
  //-------------------------------------------------------------
  static const double c0aphi_debye_huckel_slope;       // Temperature depending coefficients
  static const double c1aphi_debye_huckel_slope;       // Temperature depending coefficients
  static const double c2aphi_debye_huckel_slope;       // Temperature depending coefficients
  static const double c3aphi_debye_huckel_slope;       // Temperature depending coefficients
  static const double c4aphi_debye_huckel_slope;       // Temperature depending coefficients
  static const double c5aphi_debye_huckel_slope;       // Temperature depending coefficients
  static const double c6aphi_debye_huckel_slope;       // Temperature depending coefficients
  static const double c7aphi_debye_huckel_slope;       // Temperature depending coefficients
  static const double c8aphi_debye_huckel_slope;       // Temperature depending coefficients
  static const double c9aphi_debye_huckel_slope;       // Temperature depending coefficients
  std::vector<VirialCoefficient> beta0_virial;         // Beta0 virial coefficients
  std::vector<VirialCoefficient> beta1_virial;         // Beta1 virial coefficients
  std::vector<VirialCoefficient> beta2_virial;         // Beta2 virial coefficients
  std::vector<VirialCoefficient> theta_virial;         // Theta virial coefficients
  std::vector<VirialCoefficient> lamda_virial;         // Lamda virial coefficients
  std::vector<VirialCoefficient> psi_virial;           // Psi virial coefficients
  std::vector<VirialCoefficient> cphi_virial;          // Cphi virial coefficients
  std::vector<double> charge_product;                  // charge products [number_j_functions]
  std::vector<double> alpha1;                          // [number_b_functions]
  std::vector<double> alpha2;                          // [number_b_functions]
  std::vector<double> const_j_functions;               // constant values for j functions (Table III, Pitzer, 1975)
  double aphi_debye_huckel_slope;                      // Debye-Hückel limiting slope
  int number_b_functions;                              // Number of beta functions
  int number_j_functions;                              // Number of j functions
  int number_non_zero_beta;                            // Number of non cero Beta matrix terms
  int number_non_zero_theta;                           // Number of non cero theta matrix terms
  int number_non_zero_cphi;                            // Number of non cero C_ca matrix terms
  int number_non_zero_lamda;                           // Number of non cero lambda matrix terms
  int number_non_zero_psi;                             // Number of non cero Psi matrix terms
  int number_non_zero_q;                               // Number of non cero q matrix terms (q,q',q'',q-fi and q-fi')
  int index_cl_species;                                // Local indice of Cl species  (usefull for macinnes convention)
  int index_h2o_species;                               // Local indice of water species
  int index_k_species;                                 // Local indice of K species   (usefull for macinnes convention)
  bool macinnes_scaled;                                // macinnes_scaled=true, then activity coefficients will be scaled according macinnes convention
  std::string jfunction_approach;                      // Name of the J's functions approach
  std::vector<std::vector<double> > g_function;
  std::vector<std::vector<double> > g_pri_function;
  std::vector<std::vector<double> > f_function;
  std::vector<double> j_function;
  std::vector<double> j_pri_function;
  std::vector<double> q_matrix;
  std::vector<double> qphi_matrix;
  std::vector<double> qpri_matrix;
  std::vector<std::vector<int> > index_non_zero_q;
  std::vector<double> molality;                        // Molality of aqueous species [number_species]
  std::vector<double> charge;                          // Electric charge of the aqueous species [number_species]
  std::vector<std::string> name_species;               // Name of the aqueous species [number_species]
  int number_species;                                  // Number of aqueous species

};
}  // namespace AmanziChemistry
}  // namespace Amanzi
#endif  // AMANZI_CHEMISTRY_ACTIVITY_MODEL_PITZER_HWM_HH_
