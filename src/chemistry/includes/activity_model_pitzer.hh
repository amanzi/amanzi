/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_ACTIVITY_MODEL_PITZER_HH_
#define AMANZI_CHEMISTRY_ACTIVITY_MODEL_PITZER_HH_

#include <vector>
#include <string>
#include <cstdlib>
#include <math.h>
// Base class for activity calculations
#include "activity_model.hh"


// forward declarations
class Epetra_MultiVector;
class Epetra_Vector;
class Epetra_SerialDenseVector;

namespace amanzi {
namespace chemistry {

class Species;

class ActivityModelPitzer : public ActivityModel {
 public:
	struct theta_i {
	     double virial;
	     int fzizi;
	     int fzjzj;
	     int fzizj;
	     double zizj;
	     double zizi;
	     double zjzj;
	     int i;
	     int j;
	   };
  ActivityModelPitzer(const std::string& database, const std::vector<Species>& prim, const std::vector<AqueousEquilibriumComplex>& sec);
  ~ActivityModelPitzer();

  double Evaluate(const Species& species);

  void EvaluateVector(std::vector<double>& gamma, const std::vector<Species>& prim, const std::vector<AqueousEquilibriumComplex>& sec);

  void Display(void) const;

  void SetVirial (const double& virial, const std::string& typevirial,
  		          const int& isp1, const int& isp2, const int& isp3);

  private:

  //!%-------------------------------------------------------------
  //!% Private services
  //!%-------------------------------------------------------------
  void ReadDataBase (const std::string& database, const std::vector<Species>& prim, const std::vector<AqueousEquilibriumComplex>& sec);
  void ParseB0 (const std::string& data);
  void ParseB1 (const std::string& data);
  void ParseB2 (const std::string& data);
  void ParseCfi (const std::string& data);
  void ParseTheta (const std::string& data);
  void ParseLamda (const std::string& data);
  void ParsePsi (const std::string& data);
  void AssignFbeta ();
  void AssignFj ();
  void ComputeQ(std::vector<double>& gamma, double& osco);
  void ComputeQl(double& osco);
  void ComputeQc(std::vector<double>& gamma, double& osco);
  void ComputeT(std::vector<double>& gamma, double& osco);
  void ComputeQmatrices();
  void ComputeFbeta();
  void ComputeFj();
  void ComputeDH(std::vector<double>& gamma, double& osco, double& gclm);
  double gclm_(const double& dhterm);
  void AllocatePointers();
  void DeallocatePointers();
  //!%-------------------------------------------------------------
  //!%-------------------------------------------------------------
  //!%-------------------------------------------------------------
  static const double debyeA_pitzer;
  static const double debyeB_pitzer;
  static const double debyeBdot_pitzer;
  static const double cwater;
  //!%-------------------------------------------------------------
  //!%-------------------------------------------------------------
  //!%-------------------------------------------------------------
  static const double bdh;
  //!%-------------------------------------------------------------
  //!% Limiting Debye-Hückel slope to 25º   0.39153  0.392
  //!%-------------------------------------------------------------
  static const double aphi25;
  //!%-------------------------------------------------------------
  //!% Limiting Debye-Hückel slope to 25º for density calculations
  //!% dAphi/dP??
  //!% cm3 kg1/2/mol3/2
  //!% taken from Monnin (1994)
  //!%-------------------------------------------------------------
  static const double aphi25vol;
  //!%-------------------------------------------------------------
  //!cprovi aphi25=0.39153d0,   & ! Limiting Debye-Hückel slope to 25º   0.39153  0.392
  //!%-------------------------------------------------------------
  static const double c0aphi;                          // Temperature depending coefficients
  static const double c1aphi;                          // Temperature depending coefficients
  static const double c2aphi;                          // Temperature depending coefficients
  static const double c3aphi;                          // Temperature depending coefficients
  static const double c4aphi;                          // Temperature depending coefficients
  static const double c5aphi;                          // Temperature depending coefficients
  static const double c6aphi;                          // Temperature depending coefficients
  static const double c7aphi;                          // Temperature depending coefficients
  static const double c8aphi;                          // Temperature depending coefficients
  static const double c9aphi;                          // Temperature depending coefficients
  //--------------------------------------------------------------
  //--------------------------------------------------------------
  //--------------------------------------------------------------
  //--------------------------------------------------------------
  std::vector<double> beta0;                           // Beta_0 virial coefficients

  std::vector<double> beta1;                           // Beta_1 virial coefficients

  std::vector<double> beta2;                           // Beta_2 virial coefficients

  std::vector<double> theta;                           // Theta virial coefficients

  std::vector<double> lamda;                           // Lambda virial coefficients

  std::vector<double> psi;                             // Psi virial coefficients

  std::vector<double> cpz;                             // C pitzer matrix (computed according equation A8)

  std::vector<double> zprod;                           // charge products [nfunj]

  std::vector<double> alpha1;                          // [nfunb]

  std::vector<double> alpha2;                          // [nfunb]

  double aphi;                                         // Debye-Hückel limiting slope

  std::vector<std::vector<int> > indnzbeta;           // Non zero indices for theta virial coefficients

  std::vector<std::vector<int> > indnzlamda;           // Non zero indices for lambda virial coefficients

  std::vector<std::vector<int> > indnzpsi;             // Non zero indices for psi virial coefficients

  std::vector<std::vector<int> > indnzcpz;             // Non zero indices for c virial coefficients

  std::vector<std::vector<int> > indnztheta;           // Non zero indices for theta virial coefficients

  int nfunb;                                           // Number of beta functions

  int nfunbvol;                                        // Number of beta functions

  int nfunj;                                           // Number of j functions

  int nnzbeta;                                         // Number of non cero Beta matrix terms

  int nnztheta;                                        // Number of non cero theta matrix terms

  int nnzcpz;                                          // Number of non cero C_ca matrix terms

  int nnzlamda;                                        // Number of non cero lambda matrix terms

  int nnzpsi;                                          // Number of non cero Psi matrix terms

  int nnzq;                                            // Number of non cero q matrix terms (q,q',q'',q-fi and q-fi')

  int ithcl;                                           // Local indice of Cl species  (usefull for macinnes convention)

  int ithw;                                            // Local indice of water species  (usefull for macinnes convention)

  int ithk;                                            // Local indice of K species   (usefull for macinnes convention)

  bool ismacinnes;                                     // ismacinnes=true, then activity coefficients will be scaled according macinnes convention

  double** g_;

  double** g_pri_;

  double** f_;

  double* j_;

  double* j_pri_;

  double* q;

  double* qphi;

  double* qpri;

  int** indnzq;

  std::vector<double> molality;

  std::vector<double> charge;

  std::vector<std::string> namesp;

  int nsp;

};
}  // namespace chemistry
}  // namespace amanzi
#endif  // AMANZI_CHEMISTRY_ACTIVITY_MODEL_PITZER_HH_
