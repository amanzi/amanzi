/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

#include "activity_model_pitzer_hwm.hh"
#include <cstdlib>
#include <math.h>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include "string_tokenizer.hh"
#include "virial_coefficient.hh"
#include "chemistry_exception.hh"
#include "exceptions.hh"


namespace amanzi {
namespace chemistry {
const double ActivityModelPitzerHWM::bdh=1.2;
const double ActivityModelPitzerHWM::cwater=55.50837e0;
//-------------------------------------------------------------
// Limiting Debye-Hückel slope to 25º   0.39153  0.392
//-------------------------------------------------------------
const double ActivityModelPitzerHWM::aphi25=0.392;
//-------------------------------------------------------------
// aphi25=0.39153d0,   & ! Limiting Debye-Hückel slope to 25º   0.39153  0.392
//-------------------------------------------------------------
const double ActivityModelPitzerHWM::c0aphi=0.13422e0;      //   Temperature depending coefficients
const double ActivityModelPitzerHWM::c1aphi=0.0368329e0;    //   Temperature depending coefficients
const double ActivityModelPitzerHWM::c2aphi=14.62718e0;     //   Temperature depending coefficients
const double ActivityModelPitzerHWM::c3aphi=1530.1474e0;    //   Temperature depending coefficients
const double ActivityModelPitzerHWM::c4aphi=80.40631e0;     //   Temperature depending coefficients
const double ActivityModelPitzerHWM::c5aphi=4.1725332e0;    //   Temperature depending coefficients
const double ActivityModelPitzerHWM::c6aphi=0.1481291e0;    //   Temperature depending coefficients
const double ActivityModelPitzerHWM::c7aphi=1.5188505e0;    //   Temperature depending coefficients
const double ActivityModelPitzerHWM::c8aphi=1.8016317e0;    //   Temperature depending coefficients
const double ActivityModelPitzerHWM::c9aphi=9.3816144e0;    //   Temperature depending coefficients
/*!
    @brief ActivityModelPitzerHWM

    @class ActivityModelPitzerHWM

    @details Create the object
*/
ActivityModelPitzerHWM::ActivityModelPitzerHWM(const std::string& namedatabase,
		                                 const std::vector<Species>& prim,
		                                 const std::vector<AqueousEquilibriumComplex>& sec)
    : ActivityModel(),
      aphi(aphi25),
      nfunb(0),
      nfunbvol(0),
      nfunj(0),
      nnzbeta(0),
      nnztheta(0),
      nnzcpz(0),
      nnzlamda(0),
      nnzpsi(0),
      nnzq(0),
      ithcl(-1),
      ithw(-1),
      ithk(-1),
      ismacinnes(false)
      {
ReadDataBase(namedatabase,prim,sec);
}  // end ActivityModelPitzer constructor
/*!
    @brief ~ActivityModelPitzerHWM

    @class ActivityModelPitzerHWM

    @details Destroy the object
*/
ActivityModelPitzerHWM::~ActivityModelPitzerHWM() {
}  // end ActivityModelPitzer destructor
/*!
    @brief Evaluate

    @class ActivityModelPitzerHWM

    @details Compute the activity coefficient (return 1)
*/
double ActivityModelPitzerHWM::Evaluate(const Species& species) {
  return 1.0;
}  // end Evaluate()
/*!
    @brief EvaluateVector

    @class ActivityModelPitzerHWM

    @details Compute the activity coefficients
*/
void ActivityModelPitzerHWM::EvaluateVector(std::vector<double>& gamma, double& actw, const std::vector<Species>& prim, const std::vector<AqueousEquilibriumComplex>& sec) {
unsigned int nsp(prim.size()+sec.size());
double gcl(1.0), gclm(1.0), osco(1.0);
ActivityModel::CalculateIonicStrength(prim,sec);
ActivityModel::CalculateSumAbsZ(prim,sec);
ActivityModel::CalculateSumC(prim,sec);
int isp(-1);
for (std::vector<Species>::const_iterator i = prim.begin(); i != prim.end(); i++) {
isp++;
molality.at(isp)=(*i).molality();
}
for (std::vector<AqueousEquilibriumComplex>::const_iterator i = sec.begin(); i != sec.end(); i++) {
isp++;
molality.at(isp)=(*i).molality();
}
ComputeQmatrices();
ComputeDH(gamma,osco,gclm);
ComputeQ(gamma,osco);
ComputeQl(osco);
ComputeQc(gamma,osco);
ComputeT(gamma,osco);
if (ismacinnes && ithcl>-1) {
   gcl=gamma.at(ithcl);
   for (int i=0; i<nsp; i++) if (i!=ithw) gamma.at(i)*= pow((gcl/gclm),charge.at(i));
}
actw=osco;
}  // end EvaluateVector()
/*!
    @brief ComputeQc

    @class ActivityModelPitzerHWM

    @details Compute the m(Qc)m product.
*/
void ActivityModelPitzerHWM::ComputeQc(std::vector<double>& gamma, double& osco) {
double qc(0.0);
double ionstrz(ActivityModel::Z_);
for (int nz=0; nz<nnzcpz; nz++) {
 int i(cpz.at(nz).GetIsp1());
 int j(cpz.at(nz).GetIsp2());
 double mi(molality.at(i));
 double mj(molality.at(j));
 double Cij(cpz.at(nz).GetVirial());
 double mjCij(mj*Cij);
 double miCij(mi*Cij);
 double mimjCij(mi*mj*Cij);
 qc+=mimjCij;
 gamma.at(i)+=mjCij*ionstrz;
 gamma.at(j)+=miCij*ionstrz;
}
osco+=ionstrz*qc;
for (int i=0; i<gamma.size(); i++) gamma.at(i)+=abs(charge.at(i))*qc;
}
/*!
    @brief ComputeQl

    @class ActivityModelPitzerHWM

    @details Compute the m(Ql)m product.
*/
void ActivityModelPitzerHWM::ComputeQl(double& osco) {
for (int nz=0; nz<nnzlamda; nz++) {
 int i(lamda.at(nz).GetIsp1());
 int j(lamda.at(nz).GetIsp2());
 double Lij(lamda.at(nz).GetVirial());
 double mi(molality.at(i));
 double mj(molality.at(j));
 osco+=mi*mj*Lij;
}
} // end ComputeQl
/*!
    @brief ComputeQ

    @class ActivityModelPitzerHWM

    @details Compute mQm product.
*/
void ActivityModelPitzerHWM::ComputeQ(std::vector<double>& gamma, double& osco) {
double q2phi(0.0),q2prim(0.0);
for (int nz=0; nz<nnzq; nz++) {
 double qij(q.at(nz));
 double qpriij(qpri.at(nz));
 double qphiij(qphi.at(nz));
 int i(indnzq.at(0).at(nz));
 int j(indnzq.at(1).at(nz));
 double mi(molality.at(i));
 double mj(molality.at(j));
 q2phi+=mi*mj*qphiij;
 q2prim+=mi*mj*qpriij;
 gamma.at(i)+=2.0*qij*mj;
 gamma.at(j)+=2.0*qij*mi;
}
osco+=q2phi;
for (int i=0; i<gamma.size(); i++) gamma.at(i)+=charge.at(i)*charge.at(i)*q2prim;
}
/*!
    @brief ComputeT

    @class ActivityModelPitzerHWM

    @details Compute mTmm product.
*/
void ActivityModelPitzerHWM::ComputeT(std::vector<double>& gamma, double& osco) {
double tril(0.0);
std::vector<double> vector;
for (int i=0; i<nsp; i++) vector.push_back(0.0);
for (int nz=0; nz<nnzpsi; nz++) {
   double Psiijk(psi.at(nz).GetVirial());
   int i(psi.at(nz).GetIsp1());
   int j(psi.at(nz).GetIsp2());
   int k(psi.at(nz).GetIsp3());
   double mi(molality.at(i));
   double mj(molality.at(j));
   double mk(molality.at(k));
   tril+=mi*mj*mk*Psiijk;
   vector.at(i)+=Psiijk*mj*mk;
   vector.at(j)+=Psiijk*mi*mk;
   vector.at(k)+=Psiijk*mj*mi;
}
osco+=tril;
osco*=2.0/(ActivityModel::M_+1.0e-20);
osco+=1.0;
osco*=-ActivityModel::M_/cwater;
osco=exp(osco);
for (int i=0; i<nsp; i++) {
 if (i==ithw) {
   gamma.at(ithw)=osco;
 } else {
   gamma.at(i)=exp(gamma.at(i)+vector.at(i));
 }
}
}
/*!
    @brief ComputeDH

    @class ActivityModelPitzerHWM

    @details Compute the Debye-Huckel term.
*/
void ActivityModelPitzerHWM::ComputeDH(std::vector<double>& gamma, double& osco, double& gclm) {
const double ionstr(ActivityModel::I_),sqri(sqrt(ActivityModel::I_));
gclm=1.0;
double den(1.0+bdh*sqri);
double dh(-aphi*sqri/den);
osco=ionstr*dh;
dh-=(2.0/bdh)*aphi*log(den);
int isp(-1);
for (std::vector<double>::iterator i=gamma.begin(); i!=gamma.end(); i++) {
 isp++;
 (*i)=charge.at(isp)*charge.at(isp)*dh;
}
if (ismacinnes) gclm=gclm_(dh);
}
/*!
    @brief gclm_

    @class ActivityModelPitzerHWM

    @details Compute the activity coefficient Cl- in a KCl-K system.
*/
double ActivityModelPitzerHWM::gclm_(const double& dhterm) {
const double ionstr(ActivityModel::I_),mtb0kcl(0.04835e0),mtb1kcl(0.2122e0),mtc0kcl(-0.00084e0);
double x(2.0*sqrt(ionstr));
double xxx(-2.0*(1.0-(1.0+x+0.5*x*x)*exp(-x))/(x*x));
xxx*=mtb1kcl/ionstr;
double yyy(2.0*(1.0-(1.0+x)*exp(-x))/(x*x));
yyy*=mtb1kcl;
yyy+=mtb0kcl;
return exp(dhterm+ionstr*ionstr*xxx+ionstr*(2.0*yyy+ionstr*mtc0kcl)+ionstr*ionstr*mtc0kcl/2.0);
}
/*!
    @brief ComputeQmatrices

    @class ActivityModelPitzerHWM

    @details Compute the Q's matrices and store
    them in sparse storage
*/
void ActivityModelPitzerHWM::ComputeQmatrices() {
for (int i=0; i<nfunb; i++) {
 g_.at(i).at(0)=0.0;
 g_.at(i).at(1)=0.0;
 g_pri_.at(i).at(0)=0.0;
 g_pri_.at(i).at(1)=0.0;
 f_.at(i).at(0)=0.0;
 f_.at(i).at(1)=0.0;
}
for (int i=0; i<nfunj; i++) {
 j_.at(i)=0.0;
 j_pri_.at(i)=0.0;
}
const double str(ActivityModel::I_);
const double stri(sqrt(str));
const double str2(str*str);
const double str3(str2*str);
const double afi_stri(2.352*stri);
ComputeFbeta();
ComputeFj();
int nz(-1);
for (int nz_loc=0; nz_loc<nnzbeta; nz_loc++) {
  nz++;
  int i(beta0.at(nz_loc).GetIsp1());
  int j(beta0.at(nz_loc).GetIsp2());
  int k(beta0.at(nz_loc).GetIfun1());
  double B0ij(beta0.at(nz_loc).GetVirial());
  double B1ij(beta1.at(nz_loc).GetVirial());
  double B2ij(beta2.at(nz_loc).GetVirial());
  double expo1(f_.at(k).at(0));
  double expo2(f_.at(k).at(1));
  double G1(g_.at(k).at(0));
  double G2(g_.at(k).at(1));
  double Gp1(g_pri_.at(k).at(0));
  double Gp2(g_pri_.at(k).at(1));
  qphi.at(nz)=B0ij+B1ij*expo1+B2ij*expo2;
  q.at(nz)=B0ij+B1ij*G1+B2ij*G2;
  qpri.at(nz)=B1ij*(Gp1/str)+B2ij*(Gp2/str);
  indnzq.at(0).at(nz)=i;
  indnzq.at(1).at(nz)=j;
}
for (int nz_loc=0; nz_loc<nnztheta; nz_loc++) {
  nz++;
  int i(theta.at(nz_loc).GetIsp1());
  int j(theta.at(nz_loc).GetIsp2());
  int funij(theta.at(nz_loc).GetIfun1());
  int funii(theta.at(nz_loc).GetIfun2());
  int funjj(theta.at(nz_loc).GetIfun3());
  double thij(theta.at(nz_loc).GetVirial());
  double jij(j_.at(funij));
  double jii(j_.at(funii));
  double jjj(j_.at(funjj));
  double jpij(j_pri_.at(funij));
  double jpii(j_pri_.at(funii));
  double jpjj(j_pri_.at(funjj));
  double zizj(zprod.at(funij));
  double zizi(zprod.at(funii));
  double zjzj(zprod.at(funjj));
  double xij(afi_stri*zizj);
  double xii(afi_stri*zizi);
  double xjj(afi_stri*zjzj);
  double eth((zizj/(4.0*str))*(jij-0.5*jii-0.5*jjj));
  q.at(nz)=thij+eth;
  double eth_i(eth/str);
  double ethpri(-eth_i+(zizj/(8.0*str2))*(xij*jpij-0.5*xii*jpii-0.5*xjj*jpjj));
  qpri.at(nz)=ethpri;
  qphi.at(nz)=thij+eth+str*ethpri;
  indnzq.at(0).at(nz)=i;
  indnzq.at(1).at(nz)=j;
}
for (int nz_loc=0; nz_loc<nnzlamda; nz_loc++) {
  nz++;
  int i(lamda.at(nz_loc).GetIsp1());
  int j(lamda.at(nz_loc).GetIsp2());
  q.at(nz)=lamda.at(nz_loc).GetVirial();
  indnzq.at(0).at(nz)=i;
  indnzq.at(1).at(nz)=j;
}
if ((nz+1)!=nnzq) {
  std::ostringstream error_stream;
  error_stream << "Warning different number non-zero terms in Q's matrices" << "\n";
  Exceptions::amanzi_throw(ChemistryInvalidInput(error_stream.str()));
}
}
/*!
    @brief ComputeFbeta

    @class ActivityModelPitzerHWM

    @details Compute the Beta functions
*/
void ActivityModelPitzerHWM::ComputeFbeta() {
const double stri(sqrt(ActivityModel::I_));
for (int j=0; j<nfunb; j++) {
	   double x1(alpha1.at(j)*stri);
	   double x1q(x1*x1);
	   double x1c(x1q*x1);
	   f_.at(j).at(0)=exp(-x1);
	   g_.at(j).at(0)=2.0*(1.0-(1.0+x1)*exp(-x1))/x1q;
	   g_pri_.at(j).at(0)=-2.0*(1.0-(1.0+x1+(x1q/2.0))*exp(-x1))/x1q;
	   if(alpha2.at(j)!=0.0) {
	     double x2(alpha2.at(j)*stri);
	     double x2q(x2*x2);
	     double x2c(x2q*x2);
	     f_.at(j).at(1)=exp(-x2);
	     g_.at(j).at(1)=2.0*(1.0-(1.0+x2)*exp(-x2))/x2q;
	     g_pri_.at(j).at(1)=-2.0*(1.0-(1.0+x2+(x2q/2.0))*exp(-x2))/x2q;
	   }
}
}
/*!
    @brief ComputeFj

    @class ActivityModelPitzerHWM

    @details Compute the J's functions.
    This constants were taken from Table III in Pitzer (1975).
*/
void ActivityModelPitzerHWM::ComputeFj() {
const double e1(4.581), e2(0.7237),
		     e3(0.012), e4(0.528), e12(7.8963);
const double afi_stri(2.352e0*sqrt(ActivityModel::I_));
double d[6];
d[1]=4.118;
d[2]=7.247;
d[3]=-4.408;
d[4]=1.837;
d[5]=-0.251;
d[6]=0.0164;
for (int i=0;i<nfunj; i++) {
 double zizj(zprod.at(i));
 double x(zizj*afi_stri);
 double x2(x*x);
 double x3(x2*x);
 double x4(x3*x);
 if (x<=0.03) {
  double s1(d[6]/x);
  double s3(6.0*s1);
  for (int k=5; k>=1; k--) {
	s1=(s1+d[k])/x;
	s3=(s3+k*d[k])/x;
  }
  s3=s3/x;
  double s1q(s1*s1);
  j_.at(i)=-(1.0/6.0)*x2*log(x)*exp(-10.0*x2)+(1.0/s1);
  j_pri_.at(i)=((10.0*x2-1.0)*log(x)-0.5)*(x/3.0)*exp(-10.0*x2)+(s3/s1q);

 } else {
  double xc4(pow(x,e4));
  double xc2(pow(x,-e2));
  double td1(e1*xc2*exp(-e3*xc4));
  double td(4.0+td1);
  j_.at(i)=x/td;
  j_pri_.at(i)=(j_[i]/x2)*(x+td1*(e2+e3*e4*xc4)*j_.at(i));
 }
}
}
/*!
    @brief Display

    @class ActivityModelPitzerHWM

    @details Write the attributes
*/
void ActivityModelPitzerHWM::Display(void) const {
  std::cout << "============================================>" << std::endl;
  std::cout << "Activity model: HWM (Harvie et al., 1984)" << std::endl;
  std::cout << "============================================>" << std::endl;
  std::cout << "Species:" << std::endl;
  std::cout << "---------------------------------------------" << std::endl;
  for (int i=0; i<nsp; i++) std::cout << namesp.at(i) << " " << charge.at(i) << " " << std::endl;
  std::cout << "============================================>" << std::endl;
  std::cout << "--------------------------------------------------------------------" << std::endl;
  std::cout << " Virial coefficients" << std::endl;
  std::cout << "--------------------------------------------------------------------" << std::endl;
  int isp1(-1),isp2(-1), isp3(-1), nvirial(0);
  if (nnzbeta>0){
  std::cout << "=================> B0 ==============>" << std::endl;
  for (int i=0; i<nnzbeta; i++){
    isp1=beta0.at(i).GetIsp1();
    isp2=beta0.at(i).GetIsp2();
    if (beta0.at(i).GetVirial()!=0.0){
  	  nvirial++;
      std::cout << namesp.at(isp1) << "  " << namesp.at(isp2) << "  " << beta0.at(i).GetVirial() << std::endl;
    }
  }
  std::cout << "=================> B1 ==============>" << std::endl;
  for (int i=0; i<nnzbeta; i++){
    isp1=beta1.at(i).GetIsp1();
    isp2=beta1.at(i).GetIsp2();
    if (beta1.at(i).GetVirial()!=0.0){
  	  nvirial++;
    	  std::cout << namesp.at(isp1) << "  " << namesp.at(isp2) << "  " << beta1.at(i).GetVirial() << std::endl;
      }
  }
  std::cout << "=================> B2 ==============>" << std::endl;
  for (int i=0; i<nnzbeta; i++){
    isp1=beta2.at(i).GetIsp1();
    isp2=beta2.at(i).GetIsp2();
    if (beta2.at(i).GetVirial()!=0.0){
  	  nvirial++;
  	  std::cout << namesp.at(isp1) << "  " << namesp.at(isp2) << "  " << beta2.at(i).GetVirial() << std::endl;
    }

  }
  }
  if (nnzcpz>0){
  std::cout << "=================> Cfi ==============>" << std::endl;
  for (int i=0; i<nnzcpz; i++){
    isp1=cpz.at(i).GetIsp1();
    isp2=cpz.at(i).GetIsp2();
    if (cpz.at(i).GetVirial()!=0.0){
  	  nvirial++;
  	  std::cout << namesp.at(isp1) << "  " << namesp.at(isp2) << "  " << cpz.at(i).GetVirial() << std::endl;
    }

  }
  }
  if (nnztheta>0){
  std::cout << "=================> Theta ==============>" << std::endl;
  for (int i=0; i<nnztheta; i++){
    isp1=theta.at(i).GetIsp1();
    isp2=theta.at(i).GetIsp2();;
    if (theta.at(i).GetVirial()!=0.0){
  	  nvirial++;
  	  std::cout << namesp[isp1] << "  " << namesp[isp2] << "  " << theta.at(i).GetVirial() << std::endl;
    }

  }
  }
  if (nnzlamda>0){
  std::cout << "=================> Lamda ==============>" << std::endl;
  for (int i=0; i<nnzlamda; i++){
    isp1=lamda.at(i).GetIsp1();
    isp2=lamda.at(i).GetIsp2();
    if (lamda.at(i).GetVirial()!=0.0){
  	  nvirial++;
  	  std::cout << namesp.at(isp1) << "  " << namesp.at(isp2) << "  " << lamda.at(i).GetVirial() << std::endl;
    }

  }
  }
  if (nnzpsi>0){
  std::cout << "=================> Psi ==============>" << std::endl;
  for (int i=0; i<nnzpsi; i++){
    isp1=psi.at(i).GetIsp1();
    isp2=psi.at(i).GetIsp2();
    isp3=psi.at(i).GetIsp3();
    if (psi.at(i).GetVirial()!=0.0){
  	  nvirial++;
  	  std::cout << namesp.at(isp1) << "  " << namesp.at(isp2) << "  " << namesp.at(isp3) << "  " << psi.at(i).GetVirial() << std::endl;
    }

  }
  }
  std::cout << "=====================================>" << std::endl;
  std::cout << "Total number of virial coefficients: " << nvirial << std::endl;
  std::cout << "=====================================>" << std::endl;
  std::cout << "Total number of Beta's functions: " << nfunb << std::endl;
  std::cout << "=====================================>" << std::endl;
  std::cout << "Total number of J's functions: " << nfunj << std::endl;
  std::cout << "--------------------------------------" << std::endl;
  for (int i=0; i<nfunj; i++) std::cout << "Zi Zj product:" << zprod[i] << std::endl;
  std::cout << "=====================================>" << std::endl;
}  // end Display()
/*!
    @brief ReadDataBase

    @class ActivityModelPitzerHWM

    @details Read the virial coefficients database
*/
void ActivityModelPitzerHWM::ReadDataBase(const std::string& namedatabase,
		                                  const std::vector<Species>& prim,
		                                  const std::vector<AqueousEquilibriumComplex>& sec) {
	const bool isdebug(true);
	const int mxlines(10000);
	const int block_b0(0), block_b1(1), block_b2(2), block_cfi(3), block_theta(4), block_lamda(5), block_psi(7), block_exit(8);
    const std::string longspace("                  ");
	// Create a vector of species

	    nsp=prim.size()+sec.size();
	    if (nsp==0) std::cout << "Error, zero number of species" << std::endl;
	    //---------------------------------------------------------
	    // Is better to use iterators
	    //---------------------------------------------------------
	    for (std::vector<Species>::const_iterator i = prim.begin();
	           i != prim.end(); i++) {
	    	 molality.push_back(0.0);
	         charge.push_back((*i).charge());
	         namesp.push_back((*i).name());
	    }
	    for (std::vector<AqueousEquilibriumComplex>::const_iterator i = sec.begin();
	    	           i != sec.end(); i++) {
	    	 molality.push_back(0.0);
	    	 charge.push_back((*i).charge());
	    	 namesp.push_back((*i).name());
	    }
//---------------------------------------------------------------------------------------------
// Open Pitzer virial coefficients database
//---------------------------------------------------------------------------------------------
if (isdebug) std::cout << "=================> Opening Pitzer Data Base ==============>" << namedatabase << std::endl;
std::ifstream database(namedatabase.c_str());
if (!database) {
	    std::ostringstream error_stream;
	    error_stream << "SimpleThermoDatabase::ReadFile(): \n";
	    error_stream << "file could not be opened.... " << namedatabase << "\n";
	    Exceptions::amanzi_throw(ChemistryInvalidInput(error_stream.str()));
} else {
		if (isdebug) std::cout << "Opening Pitzer virial coefficients was successful: " << namedatabase << std::endl;
}

std::string space(" ");
StringTokenizer no_spaces;

int count(0), iblock(-1);
while (!database.eof() && count < mxlines) {
	    count++;
	    std::string error_section("");
	    std::string line;
	    getline(database,line);
	    // Read the first character
	    char first = line[0];
	    // Check if the file was finished
	    StringTokenizer input(line,space);
	    no_spaces.tokenize(input.at(0),space);
	    std::string line1(no_spaces.at(0));
	    if (line1==space) goto exit;
	    if (first == '>') {
	    	iblock++;
	    	if (iblock==block_b0 && isdebug) std::cout << "=================> Parse B0 ==============>" << std::endl;
	    	if (iblock==block_b1 && isdebug) std::cout << "=================> Parse B1 ==============>" << std::endl;
	    	if (iblock==block_b2 && isdebug) std::cout << "=================> Parse B2 ==============>" << std::endl;
	    	if (iblock==block_cfi && isdebug) std::cout << "=================> Parse Cfi ==============>" << std::endl;
	    	if (iblock==block_theta && isdebug) std::cout << "=================> Parse Theta ==============>" << std::endl;
	    	if (iblock==block_lamda && isdebug) std::cout << "=================> Parse Lamda ==============>" << std::endl;
	    	if (iblock==block_psi && isdebug) std::cout << "=================> Parse Psi ==============>" << std::endl;
	    	getline(database,line);
	    }
	    if (iblock==block_b0) {
	    	ParseB0(line);
	    }
	    if (iblock==block_b1) {
	       	ParseB1(line);
	    }
	    if (iblock==block_b2) {
	    	ParseB2(line);
	    }
	    if (iblock==block_cfi) {
	        ParseCfi (line);
	    }
	    if (iblock==block_theta) {
	    	ParseTheta (line);
	    }
	    if (iblock==block_lamda) {
	    	ParseLamda (line);
	    }
	    if (iblock==block_psi) {
	    	ParsePsi (line);
	    }

        if (iblock==block_exit) goto exit;

}

exit:
if (isdebug) std::cout << "=================> Assign Beta's functions ==============>" << std::endl;
AssignFbeta();
if (isdebug) std::cout << "=================> Assign F's functions ==============>" << std::endl;
AssignFj();
if (isdebug) std::cout << "=================> Compute total number of non-zero terms ==============>" << std::endl;
nnzq=nnzbeta+nnztheta+nnzlamda;
PushPrivateVectors();
Update(273.15,0.0);
for (int isp=0; isp<nsp; isp++) {
	if (namesp.at(isp)=="H2O" || namesp.at(isp)=="h2o") ithw=isp;
	if (namesp.at(isp)=="Cl-" || namesp.at(isp)=="cl-") ithcl=isp;
	if (namesp.at(isp)=="K+"  || namesp.at(isp)=="k+") ithk=isp;
}
if (ithcl>-1 && ithk>-1) ismacinnes=true;
database.close();
if (isdebug) std::cout << "=================> Closing Pitzer Data Base ==============>" << std::endl;
if (isdebug) std::cout << "=================> Print virial coefficients ==============>" << std::endl;
if (isdebug) Display();
if (isdebug) std::cout << "=================> End print virial coefficients ==============>" << std::endl;
}  // end ReadDataBase()
/*!
    @brief ParseB0

    @class ActivityModelPitzerHWM

    @details Compute the J's functions
*/
void ActivityModelPitzerHWM::ParseB0 (const std::string& data) {
// Local variables and constants
const double r0(0.0e0);
bool isstored(false);
bool isspecies(false);
bool isdebug(false);
std::string semicolon(";");
std::string space(" ");
StringTokenizer b0(data,space);
StringTokenizer no_spaces;

int isp1_old(0), isp2_old(0);

// get name 1
no_spaces.tokenize(b0.at(0),space);
std::string name1(no_spaces.at(0));
// get nam2 2
no_spaces.tokenize(b0.at(1),space);
std::string name2(no_spaces.at(0));
int isp1(0), isp2(0);

for (int i=0; i<nsp; i++) {
	if (namesp[i]==name1){
		isp1=i;
		isspecies=true;
        goto exit_sp1;
	}
}

exit_sp1:

if (isspecies) {
	isspecies=false;
	for (int i=0; i<nsp; i++) {
		if (namesp[i]==name2){
			isp2=i;
			isspecies=true;
            goto exit_sp2;
		}
}}

exit_sp2:

if(isspecies) {

    no_spaces.tokenize(b0.at(2),space);
    double virial(std::atof(no_spaces.at(0).c_str()));
    std::vector<double> virial_vec;
    virial_vec.push_back(virial);
    SetVirial (virial_vec,"b0",isp1,isp2,-1);
    if (isdebug) std::cout << namesp[isp1] << namesp[isp2] << virial <<std::endl;

}
}  // ParseB0
/*!
    @brief ParseB1

    @class ActivityModelPitzerHWM

    @details Compute the J's functions
*/
void ActivityModelPitzerHWM::ParseB1 (const std::string& data) {
	// Local variables and constants
	const double r0(0.0e0);
	bool isstored(false);
	bool isspecies(false);
	bool isdebug(false);
	std::string semicolon(";");
	std::string space(" ");
	StringTokenizer b1(data,space);
	StringTokenizer no_spaces;

	int isp1_old(0), isp2_old(0);

	// get name 1
	no_spaces.tokenize(b1.at(0),space);
	std::string name1(no_spaces.at(0));
	// get nam2 2
	no_spaces.tokenize(b1.at(1),space);
	std::string name2(no_spaces.at(0));
	int isp1(0), isp2(0);

	for (int i=0; i<nsp; i++) {
		if (namesp[i]==name1){
			isp1=i;
			isspecies=true;
	        goto exit_sp1;
		}
	}

	exit_sp1:

	if (isspecies) {
		isspecies=false;
		for (int i=0; i<nsp; i++) {
			if (namesp[i]==name2){
				isp2=i;
				isspecies=true;
	            goto exit_sp2;
			}
	}}

	exit_sp2:

	if(isspecies) {

	no_spaces.tokenize(b1.at(2), space);
	double virial(std::atof(no_spaces.at(0).c_str()));
	std::vector<double> virial_vec;
	virial_vec.push_back(virial);
	SetVirial (virial_vec,"b1",isp1,isp2,-1);
	if (isdebug) std::cout << namesp[isp1] << namesp[isp2] << virial <<std::endl;
}
}
/*!
    @brief ParseB2

    @class ActivityModelPitzerHWM

    @details Compute the J's functions
*/
void ActivityModelPitzerHWM::ParseB2 (const std::string& data) {
	// Local variables and constants
	const double r0(0.0e0);
	bool isstored(false);
	bool isspecies(false);
	bool isdebug(false);
	std::string semicolon(";");
	std::string space(" ");
	StringTokenizer b2(data,space);
	StringTokenizer no_spaces;

	int isp1_old(0), isp2_old(0);

	// get name 1
	no_spaces.tokenize(b2.at(0),space);
	std::string name1(no_spaces.at(0));
	// get nam2 2
	no_spaces.tokenize(b2.at(1),space);
	std::string name2(no_spaces.at(0));
	int isp1(0), isp2(0);

	for (int i=0; i<nsp; i++) {
		if (namesp[i]==name1){
			isp1=i;
			isspecies=true;
	        goto exit_sp1;
		}
	}

	exit_sp1:

	if (isspecies) {
		isspecies=false;
		for (int i=0; i<nsp; i++) {
			if (namesp[i]==name2){
				isp2=i;
				isspecies=true;
	            goto exit_sp2;
			}
	}}

	exit_sp2:

	if(isspecies) {

	no_spaces.tokenize(b2.at(2), space);
	double virial(std::atof(no_spaces.at(0).c_str()));
	std::vector<double> virial_vec;
	virial_vec.push_back(virial);
	SetVirial (virial_vec,"b2",isp1,isp2,-1);
	if (isdebug) std::cout << namesp[isp1] << namesp[isp2] << virial <<std::endl;
}
}
/*!
    @brief CComputeFj

    @class ActivityModelPitzerHWM

    @details Compute the J's functions
*/
void ActivityModelPitzerHWM::ParseCfi (const std::string& data) {
	// Local variables and constants
	const double r0(0.0e0), r1(1.0e0), r2(2.0e0);
	bool isstored(false);
	bool isspecies(false);
	bool isdebug(false);
	std::string semicolon(";");
	std::string space(" ");
	StringTokenizer cfi(data,space);
	StringTokenizer no_spaces;


	int isp1_old(0), isp2_old(0);

	// get name 1
	no_spaces.tokenize(cfi.at(0),space);
	std::string name1(no_spaces.at(0));
	// get nam2 2
	no_spaces.tokenize(cfi.at(1),space);
	std::string name2(no_spaces.at(0));
	int isp1(0), isp2(0);

	for (int i=0; i<nsp; i++) {
		if (namesp[i]==name1){
			isp1=i;
			isspecies=true;
	        goto exit_sp1;
		}
	}

	exit_sp1:

	if (isspecies) {
		isspecies=false;
		for (int i=0; i<nsp; i++) {
			if (namesp[i]==name2){
				isp2=i;
				isspecies=true;
	            goto exit_sp2;
			}
	    }
	}

	exit_sp2:

	if(isspecies) {

	no_spaces.tokenize(cfi.at(2), space);
	double virial(std::atof(no_spaces.at(0).c_str()));
	virial/=r2*sqrt(abs(charge[isp1]*charge[isp2]));
	std::vector<double> virial_vec;
	virial_vec.push_back(virial);
	SetVirial (virial_vec,"cfi",isp1,isp2,-1);
	if (isdebug) std::cout << namesp[isp1] << namesp[isp2] << virial <<std::endl;
    }
}
/*!
    @brief CComputeFj

    @class ActivityModelPitzerHWM

    @details Compute the J's functions
*/
void ActivityModelPitzerHWM::ParseTheta (const std::string& data) {
	// Local variables and constants
	const double r0(0.0e0), r1(1.0e0), r2(2.0e0);
	const bool isdebug(false);
	bool isstored(false);
	bool isspecies(false);
	std::string semicolon(";");
	std::string space(" ");
	StringTokenizer th(data,space);
	StringTokenizer no_spaces;


	int isp1_old(0), isp2_old(0);

	// get name 1
	no_spaces.tokenize(th.at(0),space);
	std::string name1(no_spaces.at(0));
	// get name 2
	no_spaces.tokenize(th.at(1),space);
	std::string name2(no_spaces.at(0));
	int isp1(-1), isp2(-1);

	for (int i=0; i<nsp; i++) {
		if (namesp[i]==name1){
			isp1=i;
			isspecies=true;
	        goto exit_sp1;
		}
	}

	exit_sp1:

	if (isspecies) {
		isspecies=false;
		for (int i=0; i<nsp; i++) {
			if (namesp[i]==name2){
				isp2=i;
				isspecies=true;
	            goto exit_sp2;
			}
	    }
	}

	exit_sp2:

	if(isspecies) {

	no_spaces.tokenize(th.at(2), space);
	double virial(std::atof(no_spaces.at(0).c_str()));
	std::vector<double> virial_vec;
	virial_vec.push_back(virial);
	SetVirial (virial_vec,"theta",isp1,isp2,-1);
	if (isdebug){
		std::cout << namesp[isp1] << namesp[isp2] << virial <<std::endl;
		for (int i=0; i<nnztheta; i++){
		    std::cout << "zi zj" << theta[i].GetIfun1() << std::endl;
		    std::cout << "zi zi" << theta[i].GetIfun2() << std::endl;
		    std::cout << "zj zj" << theta[i].GetIfun3() << std::endl;
	     }
   }
}
}
/*!
    @brief CComputeFj

    @class ActivityModelPitzerHWM

    @details Compute the J's functions
*/
void ActivityModelPitzerHWM::ParseLamda (const std::string& data) {
	// Local variables and constants
	const double r0(0.0e0), r1(1.0e0), r2(2.0e0);
	bool isstored(false);
	bool isspecies(false);
	bool isdebug(false);
	std::string semicolon(";");
	std::string space(" ");
	StringTokenizer lam(data,space);
	StringTokenizer no_spaces;


	int isp1_old(0), isp2_old(0);

	// get name 1
	no_spaces.tokenize(lam.at(0),space);
	std::string name1(no_spaces.at(0));
	// get name 2
	no_spaces.tokenize(lam.at(1),space);
	std::string name2(no_spaces.at(0));
	int isp1(0), isp2(0);

	for (int i=0; i<nsp; i++) {
		if (namesp[i]==name1){
			isp1=i;
			isspecies=true;
	        goto exit_sp1;
		}
	}

	exit_sp1:

	if (isspecies) {
		isspecies=false;
		for (int i=0; i<nsp; i++) {
			if (namesp[i]==name2){
				isp2=i;
				isspecies=true;
	            goto exit_sp2;
			}
	    }
	}

	exit_sp2:

	if(isspecies) {

	no_spaces.tokenize(lam.at(2), space);
	double virial(std::atof(no_spaces.at(0).c_str()));
	std::vector<double> virial_vec;
	virial_vec.push_back(virial);
	SetVirial (virial_vec,"lamda",isp1,isp2,-1);
	if (isdebug) std::cout << namesp[isp1] << namesp[isp2] << virial <<std::endl;

}
}
/*!
    @brief ParsePsi

    @class ActivityModelPitzerHWM

    @details Parse Psi virial coefficients
*/
void ActivityModelPitzerHWM::ParsePsi (const std::string& data) {
	// Local variables and constants
	const double r0(0.0e0), r1(1.0e0), r2(2.0e0);
	bool isstored(false);
	bool isspecies(false);
	bool isdebug(false);
	std::string semicolon(";");
	std::string space(" ");
	StringTokenizer psi_loc(data,space);
	StringTokenizer no_spaces;


	int isp1_old(0), isp2_old(0);

	// get name 1
	no_spaces.tokenize(psi_loc.at(0),space);
	std::string name1(no_spaces.at(0));
	// get name 2
	no_spaces.tokenize(psi_loc.at(1),space);
	std::string name2(no_spaces.at(0));
	// get name 3
	no_spaces.tokenize(psi_loc.at(2),space);
	std::string name3(no_spaces.at(0));
	int isp1(0), isp2(0), isp3(0);

	for (int i=0; i<nsp; i++) {
		if (namesp[i]==name1){
			isp1=i;
			isspecies=true;
	        goto exit_sp1;
		}
	}

	exit_sp1:

	if (isspecies) {
		isspecies=false;
		for (int i=0; i<nsp; i++) {
			if (namesp[i]==name2){
				isp2=i;
				isspecies=true;
	            goto exit_sp2;
			}
	    }
	}

	exit_sp2:

	if (isspecies) {
		isspecies=false;
		for (int i=0; i<nsp; i++) {
			if (namesp[i]==name3){
				isp3=i;
				isspecies=true;
	            goto exit_sp3;
			}
	    }
	}

	exit_sp3:

	if(isspecies) {

	no_spaces.tokenize(psi_loc.at(3), space);
	double virial(std::atof(no_spaces.at(0).c_str()));
	std::vector<double> virial_vec;
	virial_vec.push_back(virial);
	SetVirial (virial_vec,"psi",isp1,isp2,isp3);
	if (isdebug) std::cout << namesp[isp1] << namesp[isp2] << namesp[isp3] << virial << std::endl;
}
}
/*!
    @brief SetVirial

    @class ActivityModelPitzerHWM

    @details Set virial coefficients
*/
void ActivityModelPitzerHWM::SetVirial (const std::vector<double>& virial, const std::string& typevirial,
		                             const int& isp1, const int& isp2, const int& isp3) {
VirialCoefficient vir;
VirialCoefficient vir0;
bool isstored(false);


if (typevirial=="b0") {

	 for (int i=0; i<nnzbeta; i++) {
	   	int isp1_old=beta0[i].GetIsp1();
	   	int isp2_old=beta0[i].GetIsp2();
	   	if ((isp1==isp1_old && isp2==isp2_old) || (isp1==isp2_old && isp2==isp1_old)) {

	   		for (std::vector<double>::const_iterator j=virial.begin();j!=virial.end();j++) beta0[i].SetPol((*j));
	   		isstored=true;
	   		break;
	   	}
	 }

	 if (!isstored) {
	    nnzbeta++;
	    for (std::vector<double>::const_iterator j=virial.begin();j!=virial.end();j++) vir.SetPol((*j));
	    vir.SetIsp1(isp1);
	    vir.SetIsp2(isp2);
	    vir0.SetIsp1(isp1);
	    vir0.SetIsp2(isp2);
	    vir.SetIfun1(0);
	    vir0.SetIfun1(0);
	    beta0.push_back(vir);
	    beta1.push_back(vir0);
	    beta2.push_back(vir0);
	  }

} else if (typevirial=="b1"){

	    for (int i=0; i<nnzbeta; i++) {
		   	int isp1_old=beta1.at(i).GetIsp1();
		   	int isp2_old=beta1.at(i).GetIsp2();;
		   	if ((isp1==isp1_old && isp2==isp2_old) || (isp1==isp2_old && isp2==isp1_old)) {
		   		for (std::vector<double>::const_iterator j=virial.begin();j!=virial.end();j++) beta1.at(i).SetPol((*j));
		   		isstored=true;
		   		break;
		   	}
		}
		if (!isstored) {
		      nnzbeta++;
		      for (std::vector<double>::const_iterator j=virial.begin();j!=virial.end();j++) vir.SetPol((*j));
		      vir.SetIsp1(isp1);
		      vir.SetIsp2(isp2);
		      vir0.SetIsp1(isp1);
		      vir0.SetIsp2(isp2);
		      vir.SetIfun1(0);
		      vir0.SetIfun1(0);
		      beta1.push_back(vir);
		      beta0.push_back(vir0);
		      beta2.push_back(vir0);
		}

} else if(typevirial=="b2"){


		for (int i=0; i<nnzbeta; i++) {
		   	int isp1_old=beta2[i].GetIsp1();
		   	int isp2_old=beta2[i].GetIsp1();
		   	if ((isp1==isp1_old && isp2==isp2_old) || (isp1==isp2_old && isp2==isp1_old)) {
		   		for (std::vector<double>::const_iterator j=virial.begin();j!=virial.end();j++) beta2.at(i).SetPol((*j));
		   		isstored=true;
		   		break;
		   	}
		}
		if (!isstored) {
		nnzbeta++;
		for (std::vector<double>::const_iterator j=virial.begin();j!=virial.end();j++) vir.SetPol((*j));
		vir.SetIsp1(isp1);
		vir.SetIsp2(isp2);
		vir0.SetIsp1(isp1);
		vir0.SetIsp2(isp2);
		vir.SetIfun1(0);
		vir0.SetIfun1(0);
		beta2.push_back(vir);
		beta0.push_back(vir0);
		beta1.push_back(vir0);
		}
} else if (typevirial=="cfi"){

   	    nnzcpz++;
   	    for (std::vector<double>::const_iterator j=virial.begin();j!=virial.end();j++) vir.SetPol((*j));
   	    vir.SetIsp1(isp1);
   	    vir.SetIsp2(isp2);
   	    cpz.push_back(vir);

} else if(typevirial=="theta"){

		nnztheta++;
		for (std::vector<double>::const_iterator j=virial.begin();j!=virial.end();j++) vir.SetPol((*j));
		vir.SetIsp1(isp1);
		vir.SetIsp2(isp2);
		vir.SetIfun1(int(charge.at(isp1)*charge.at(isp2)));
		vir.SetIfun2(int(charge.at(isp1)*charge.at(isp1)));
		vir.SetIfun3(int(charge.at(isp2)*charge.at(isp2)));
		theta.push_back(vir);


} else if(typevirial=="lamda"){

	   	nnzlamda++;
	   	for (std::vector<double>::const_iterator j=virial.begin();j!=virial.end();j++) vir.SetPol((*j));
		vir.SetIsp1(isp1);
		vir.SetIsp2(isp2);
		lamda.push_back(vir);

} else if(typevirial=="psi"){

	    nnzpsi++;
	    for (std::vector<double>::const_iterator j=virial.begin();j!=virial.end();j++) vir.SetPol((*j));
	    vir.SetIsp1(isp1);
	    vir.SetIsp2(isp2);
	    vir.SetIsp3(isp3);
	    psi.push_back(vir);

} else {
	    std::ostringstream error_stream;
	   	error_stream << "Type virial coefficient not defined" << typevirial;
	   	Exceptions::amanzi_throw(ChemistryInvalidInput(error_stream.str()));
}



}
/*!
    @brief AssignFbeta

    @class ActivityModelPitzerHWM

    @details Assign Beta's functions
*/
void ActivityModelPitzerHWM::AssignFbeta () {
// Local variables and constants
int l1(0), l2(0), l3(0), n1(0), n2(0), n3(0), nz(0);
double z1(0.0e0), z2(0.0e0);
nfunb=0;
for (nz=0; nz<nnzbeta; nz++) {

   // Take the electric charge
   z1=charge[beta0[nz].GetIsp1()];
   z2=charge[beta0[nz].GetIsp2()];
   if (abs(z1)==1.0 || abs(z2)==1.0) {
	    	       	if (l1==0) {

	    	       		   nfunb++;
						   beta0[nz].SetIfun1(nfunb-1);
						   beta1[nz].SetIfun1(nfunb-1);
						   beta2[nz].SetIfun1(nfunb-1);

						   n1 = nfunb-1;

						   alpha1.push_back(2.0);
						   alpha2.push_back(12.0);

						   l1=1;
                   	} else {

                   		beta0[nz].SetIfun1(n1);
                   		beta1[nz].SetIfun1(n1);
                   		beta2[nz].SetIfun1(n1);

				       }

	       } else {

					 if (abs(z1)!=abs(z2)) {


							 if (l2==0) {

							   nfunb++;
						       beta0[nz].SetIfun1(nfunb-1);
							   beta1[nz].SetIfun1(nfunb-1);
							   beta2[nz].SetIfun1(nfunb-1);

							   n2 = nfunb-1;

							   alpha1.push_back(2.0);
							   alpha2.push_back(50.0);

	                           l2=1;

							 } else {

								 beta0[nz].SetIfun1(n2);
								 beta1[nz].SetIfun1(n2);
								 beta2[nz].SetIfun1(n2);

	                           }

					 } else {

						     if (l3==0) {

						       nfunb++;

						       beta0[nz].SetIfun1(nfunb-1);
						       beta1[nz].SetIfun1(nfunb-1);
						       beta2[nz].SetIfun1(nfunb-1);;

						       n3=nfunb-1;

                               alpha1.push_back(1.4);
                               alpha2.push_back(12.0);

	                           l3=1;

						     } else {

						    	 beta0[nz].SetIfun1(n3);
						    	 beta1[nz].SetIfun1(n3);
						    	 beta2[nz].SetIfun1(n3);

						     }

					 }

                   	}


	       }


}  // end AssignFbeta
/*!
    @brief AssignFj

    @class ActivityModelPitzerHWM

    @details Compute the J's functions
*/
void ActivityModelPitzerHWM::AssignFj () {
int isp1(0), isp2(0);
bool isstored(false);
bool isfound(false);
VirialCoefficient vir;
//---------------------------------------------------------------
// Find other non-zero theta terms
//---------------------------------------------------------------
for (int i=0; i<nsp; i++) {
			for (int j=0; j<nsp; j++) {
			   if (i!=j && ((charge[i]>0.0 && charge[j]>0.0) || (charge[i]<0.0 && charge[j]<0.0))) {
					   for (int k=0; k<nnztheta; k++) {
					        isp1=theta[k].GetIsp1();
					        isp2=theta[k].GetIsp2();
							if ((isp1==i && isp2==j) || (isp1==j && isp2==i)) {
								isstored=true;
								goto exit_1;
							}
					   }

					   exit_1:

					   if (!isstored) {

						   nnztheta++;
						   vir.SetPol(0.0);
						   vir.SetIsp1(i);
						   vir.SetIsp2(j);
						   vir.SetIfun1(int(charge[i]*charge[j]));
						   vir.SetIfun2(int(charge[i]*charge[i]));
						   vir.SetIfun3(int(charge[j]*charge[j]));
						   theta.push_back(vir);

					   } else {
						   isstored=false;
					   }


		         }

		     }
}
//!%-----------------------------------------------------------
//!% Compute number of functions j and save
//!%-----------------------------------------------------------
nfunj=1;
zprod.push_back(double(theta[0].GetIfun1()));
double zz(0.0);
for (int nz=0; nz<nnztheta; nz++) {

	   for (int j=0; j<nfunj; j++) {
		  zz=double(theta[nz].GetIfun1());
		  if (zprod[j]==zz) {
	    	    theta[nz].SetIfun1(j);
	            isfound=true;
	            break;
	      }
	    }
	    if (isfound){
	   	   isfound=false;
	   	} else{
           zprod.push_back(zz);
    	   nfunj++;
    	   theta[nz].SetIfun1(nfunj-1);
        }
	    for (int j=0; j<nfunj; j++) {
	       zz=double(theta[nz].GetIfun2());
	       if (zprod[j]==zz) {
	      	    theta[nz].SetIfun2(j);
	            isfound=true;
	            break;
	       }
	     }
	     if (isfound){
	    	 isfound=false;
	     } else{
	   		 zprod.push_back(zz);
  	   		 nfunj++;
   	   	     theta[nz].SetIfun2(nfunj-1);
	     }
	     for (int j=0; j<nfunj; j++) {
	    	 zz=double(theta[nz].GetIfun3());
	    	 if (zprod[j]==zz) {
	    	   theta[nz].SetIfun3(j);
	    	   isfound=true;
	    	   break;
	    	 }
	     }
	     if (isfound){
	    	 isfound=false;
	     } else{
  	   		 zprod.push_back(zz);
   	   		 nfunj++;
   	   	     theta[nz].SetIfun3(nfunj-1);
         }

}
}  // end AssignFJ
/*!
    @brief PushPrivateVectors

    @class ActivityModelPitzerHWM

    @details Push back private vectors.
*/
void ActivityModelPitzerHWM::PushPrivateVectors() {
g_.resize(nfunb);
g_pri_.resize(nfunb);
g_.resize(nfunb);
f_.resize(nfunb);
for (int i=0;i<nfunb;i++) {
 g_.at(i).push_back(0.0);
 g_.at(i).push_back(0.0);
 g_pri_.at(i).push_back(0.0);
 g_pri_.at(i).push_back(0.0);
 f_.at(i).push_back(0.0);
 f_.at(i).push_back(0.0);
}
for (int i=0;i<nfunj;i++) {
 j_.push_back(0.0);
 j_pri_.push_back(0.0);
}
indnzq.resize(2);
for (int i=0; i<nnzq; i++) {
 indnzq.at(0).push_back(0);
 indnzq.at(1).push_back(0);
 q.push_back(0.0);
 qpri.push_back(0.0);
 qphi.push_back(0.0);
}
}
/*!
    @brief Update

    @class ActivityModelPitzerHWM

    @details Update virial coefficients with temperature and liquid pressure.
*/
void ActivityModelPitzerHWM::Update (const double& temp, const double& pressure) {
  for (std::vector<VirialCoefficient>::iterator i=beta0.begin(); i!=beta0.end(); i++) (*i).UpdateVirial(temp,pressure);
  for (std::vector<VirialCoefficient>::iterator i=beta1.begin(); i!=beta1.end(); i++) (*i).UpdateVirial(temp,pressure);
  for (std::vector<VirialCoefficient>::iterator i=beta2.begin(); i!=beta2.end(); i++) (*i).UpdateVirial(temp,pressure);
  for (std::vector<VirialCoefficient>::iterator i=cpz.begin(); i!=cpz.end(); i++) (*i).UpdateVirial(temp,pressure);
  for (std::vector<VirialCoefficient>::iterator i=theta.begin(); i!=theta.end(); i++) (*i).UpdateVirial(temp,pressure);
  for (std::vector<VirialCoefficient>::iterator i=psi.begin(); i!=psi.end(); i++) (*i).UpdateVirial(temp,pressure);
  for (std::vector<VirialCoefficient>::iterator i=lamda.begin(); i!=lamda.end(); i++) (*i).UpdateVirial(temp,pressure);
}
}  // namespace chemistry
}  // namespace amanzi
