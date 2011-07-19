/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

#include "activity_model_pitzer.hh"
#include <cstdlib>
#include <math.h>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include "string_tokenizer.hh"


namespace amanzi {
namespace chemistry {

// Constant parameters for Pitzer model
const double ActivityModelPitzer::debyeA_pitzer = 0.5114;     // 25C
const double ActivityModelPitzer::debyeB_pitzer = 0.3288;     // 25C
const double ActivityModelPitzer::debyeBdot_pitzer = 0.0410;  // 25C

const double ActivityModelPitzer::bdh=1.2e0;

const double ActivityModelPitzer::cwater=55.50837e0;
//!%-------------------------------------------------------------
//!% Limiting Debye-Hückel slope to 25º   0.39153  0.392
//!%-------------------------------------------------------------
const double ActivityModelPitzer::aphi25=0.392e0;
//!%-------------------------------------------------------------
//!% Limiting Debye-Hückel slope to 25º for density calculations
//!% dAphi/dP??
//!% cm3 kg1/2/mol3/2
//!% taken from Monnin (1994)
//!%-------------------------------------------------------------
const double ActivityModelPitzer::aphi25vol=1.8743e0;
//!%-------------------------------------------------------------
//!cprovi aphi25=0.39153d0,   & ! Limiting Debye-Hückel slope to 25º   0.39153  0.392
//!%-------------------------------------------------------------
const double ActivityModelPitzer::c0aphi=0.13422e0;      //   Temperature depending coefficients
const double ActivityModelPitzer::c1aphi=0.0368329e0;    //   Temperature depending coefficients
const double ActivityModelPitzer::c2aphi=14.62718e0;     //   Temperature depending coefficients
const double ActivityModelPitzer::c3aphi=1530.1474e0;    //   Temperature depending coefficients
const double ActivityModelPitzer::c4aphi=80.40631e0;     //   Temperature depending coefficients
const double ActivityModelPitzer::c5aphi=4.1725332e0;    //   Temperature depending coefficients
const double ActivityModelPitzer::c6aphi=0.1481291e0;    //   Temperature depending coefficients
const double ActivityModelPitzer::c7aphi=1.5188505e0;    //   Temperature depending coefficients
const double ActivityModelPitzer::c8aphi=1.8016317e0;    //   Temperature depending coefficients
const double ActivityModelPitzer::c9aphi=9.3816144e0;    //   Temperature depending coefficients
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
ActivityModelPitzer::ActivityModelPitzer(const std::string& namedatabase,
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
//--------------------------------------------------------------------------
// Read Pitzer virial coefficients database
//--------------------------------------------------------------------------
ReadDataBase(namedatabase,prim,sec);
}  // end ActivityModelPitzer constructor
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
ActivityModelPitzer::~ActivityModelPitzer() {
}  // end ActivityModelDebyeHuckel destructor
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
double ActivityModelPitzer::Evaluate(const Species& species) {
  const double r1(1.0e0);
  double gamma(r1);
  return gamma;
}  // end Evaluate()
//--------------------------------------------------------------------------
// Compute activity coefficients from a vector of species
// It could be more efficient computationally
//--------------------------------------------------------------------------
void ActivityModelPitzer::EvaluateVector(std::vector<double>& gamma, const std::vector<Species>* prim, const std::vector<AqueousEquilibriumComplex>* sec) {

const double r0(0.0e0), r1(1.0e0);
int nsp(prim->size()+sec->size());

double gcl(r1), gclm(r1), osco(r1);
// Compute
ActivityModel::CalculateSumAbsZ(prim,sec);
ActivityModel::CalculateSumC(prim,sec);

std::vector<double> molality;
std::vector<double> charge;
//---------------------------------------------------------
// Is better to use iterators
//---------------------------------------------------------
for (std::vector<Species>::const_iterator i = prim->begin();
          i != prim->end(); i++) {
	    charge.push_back((*i).charge());
	    molality.push_back((*i).molality());
	    }
for (std::vector<AqueousEquilibriumComplex>::const_iterator i = sec->begin();
          i != sec->end(); i++) {
	    charge.push_back((*i).charge());
	    molality.push_back((*i).molality());
  	    }
//!%-------------------------------------------------------------
//!% Resize Q's matrices
//!%-------------------------------------------------------------
std::vector<double> q, qphi, qpri;
std::vector<std::vector<int> > indnzq;
q.resize(nnzq,r0);
qphi.resize(nnzq,r0);
qpri.resize(nnzq,r0);
indnzq.resize(2);
for (int i=0; i<2; i++) indnzq[i].resize(nnzq,0);
for (int i=0; i<nnzq; i++) {
	q[i]=r0;
	qphi[i]=r0;
	qpri[i]=r0;
	indnzq[0][i]=0;
	indnzq[1][i]=0;
}
//!%-------------------------------------------------------------
//!% Compute Q's matrices
//!%-------------------------------------------------------------
ComputeQmatrices(q,qphi,qpri,indnzq);
//!%-------------------------------------------------------------
//!% Compute DH
//!%-------------------------------------------------------------
ComputeDH(gamma,osco,gclm,charge);
//!%-------------------------------------------------------------
//!% Compute q
//!%-------------------------------------------------------------
ComputeQ(gamma,osco,q,qpri,qphi,indnzq,molality,charge);
//!%-------------------------------------------------------------
//!% Compute ql
//!%-------------------------------------------------------------
ComputeQl(osco,molality);
//!%-------------------------------------------------------------
//!% Compute qc
//!%-------------------------------------------------------------
ComputeQc(gamma,osco,molality,charge);
//!%-------------------------------------------------------------
//!% Compute t
//!%-------------------------------------------------------------
ComputeT(gamma,osco,molality);
//!%-------------------------------------------------------------
//!% If ismacinnes=true, activity coefficients are scaled
//!% according to MacInnes convention (MacInnes, 1919)
//!%-------------------------------------------------------------
if (ismacinnes && ithcl>-1) {
   gcl = gamma[ithcl];
   for (int i=0; i<nsp; i++) {
	  if (i!=ithw) gamma[i] *= pow((gcl / gclm),charge[i]);
   }
}
}  // end EvaluateVector()
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
void ActivityModelPitzer::ComputeQc(std::vector<double>& gamma, double& osco, const std::vector<double>& molality,
		const std::vector<double>& charge) {
const double r0(0.0e0);
int i(0);
int j(0);
double qc(r0), mi(r0), Cij(r0), mj(r0), miCij(r0), mjCij(r0), mimjCij(r0);
double ionstrz(ActivityModel::Z_);

//!%-----------------------------------------------------------
for (int nz=0; nz<nnzcpz; nz++) {
	      i=indnzcpz[0][nz];
	      j=indnzcpz[1][nz];
	      mi=molality[i];
	      mj=molality[j];
	      Cij=cpz[nz];
	      mjCij=mj*Cij;
	      miCij=mi*Cij;
	      mimjCij=mi*mj*Cij;
	      qc+=mimjCij;
	      gamma[i]+=mjCij*ionstrz;
	      gamma[j]+=miCij*ionstrz;
}
//!%-----------------------------------------------------------
if (ithw>-1) osco+=ionstrz*qc;
//!%-----------------------------------------------------------
for (int i=0; i<gamma.size(); i++) gamma[i]+=abs(charge[i])*qc;
}
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
void ActivityModelPitzer::ComputeQl(double& osco, const std::vector<double>& molality) {
const double r0(0.0e0);
int i(0);
int j(0);
double mi(r0), mj(r0), Lij(r0);

if (ithw>-1 && nnzlamda>0) {

  for (int nz=0; nz<nnzlamda; nz++) {

      i=indnzlamda[0][nz];
      j=indnzlamda[1][nz];
      Lij=lamda[nz];
      mi=molality[i];
      mj=molality[j];
      osco+=mi*mj*Lij;

  }


}
}
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
void ActivityModelPitzer::ComputeQ(std::vector<double>& gamma, double& osco, const std::vector<double>& q,
		                           const std::vector<double>& qpri, const std::vector<double>& qphi,
		                           const std::vector<std::vector<int> >& indnzq,
		                           const std::vector<double>& molality,
		                           const std::vector<double>& charge) {
const double r0(0.0e0), r1(1.0e0), r2(2.0e0);
int i(0), j(0), nz(0);
//!%-----------------------------------------------------------

double q2phi(r0), q2prim(r0), qij(r0), qpriij(r0), qphiij(r0), mi(r0), mj(r0);
//!%-----------------------------------------------------------
for (nz=0; nz<nnzq; nz++) {

      qij=q[nz];
      qpriij=qpri[nz];
      qphiij=qphi[nz];
      i=indnzq[0][nz];
      j=indnzq[1][nz];
      mi=molality[i];
      mj=molality[j];
//!%------------------------------------------------------------
//!% q2_PHi = sUMij(mi*mj*qPHiij)
//!%------------------------------------------------------------
      q2phi+=mi*mj*qphiij;
//!%------------------------------------------------------------
//!%q'2 = sUMij(mi*mj*q'ij)
//!%------------------------------------------------------------
      q2prim+=mi*mj*qpriij;
//!%------------------------------------------------------------
//!%GAMi = GAMi + sUMj(qij*mj)
//!%------------------------------------------------------------
      gamma[i]+=r2*qij*mj;

      gamma[j]+=r2*qij*mi;

}
//!%------------------------------------------------------------
//!% Add the above terms to osco and osco = osco + q2_PHi
//!%------------------------------------------------------------
if (ithw>-1) osco+=q2phi;
//!%------------------------------------------------------------
 for (i=0; i<gamma.size(); i++) {
//!%------------------------------------------------------------
//!% GAMi = GAMi + (chargei^2)*q'2
//!%------------------------------------------------------------
     gamma[i]+=charge[i] * charge[i] * q2prim;

 }
}
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
void ActivityModelPitzer::ComputeT(std::vector<double>& gamma, double& osco, const std::vector<double>& molality) {

const double r0(0.0e0), r1(1.0e0), r2(2.0e0);
int i(0), j(0), k(0);
double mi(r0), mj(r0), mk(r0), Psiijk(r0);
double tril(r0);
const double ctw1(ActivityModel::M_/cwater), ctw2(r2/cwater);

std::vector<double> vector;
vector.resize(gamma.size(),r0);
for (int i=0; i<gamma.size(); i++) vector[i]=r0;
//!%------------------------------------------------------
//!% Compute t and vect
//!%------------------------------------------------------
for (int nz=0; nz<nnzpsi; nz++) {

   Psiijk=psi[nz];
   i=indnzpsi[0][nz];
   j=indnzpsi[1][nz];
   k=indnzpsi[2][nz];
   mi=molality[i];
   mj=molality[j];
   mk=molality[k];
   tril+=mi*mj*mk*Psiijk;
//!%----------------------------------------------------------
//!% este se utiliza para el cálculo de los coeficientes de
//!% actividad de las especies acuosas que no sean H2O
//!%----------------------------------------------------------
   vector[i]+=Psiijk*mj*mk;
   vector[j]+=Psiijk*mi*mk;
   vector[k]+=Psiijk*mj*mi;
}
//!%----------------------------------------------------
//!%
//!%----------------------------------------------------
for (int i=0; i<gamma.size(); i++) {

	if (i==ithw) {
      osco+=tril;
      osco=(ctw2*osco)+r1;
      gamma[ithw]=-osco*ctw1;
      gamma[ithw]=exp(gamma[i]);
    } else {
    	gamma[i]=exp(gamma[i]+vector[i]);
    }
}

}

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
void ActivityModelPitzer::ComputeDH(std::vector<double>& gamma, double& osco, double& gclm, const std::vector<double>& charge) {

const double r0(0.0e0), r1(1.0e0), r2(2.0e0), ionstr(ActivityModel::I_), sqri(sqrt(ActivityModel::I_));
double den(r0), dh(r0);
//!%-----------------------------------------------------------
//!% Zeroing
//!%-----------------------------------------------------------
gclm=r1;
//!%-----------------------------------------------------------
//!% Eq. A1 and A2
//!%-----------------------------------------------------------
den = r1+bdh*sqri;
dh = -aphi*sqri/den;
if (ithw>-1) osco = ionstr * dh;
dh -= (r2/bdh)*aphi*log(den);
//!%-----------------------------------------------------------
int isp(-1);
for (std::vector<double>::iterator i=gamma.begin(); i!=gamma.end(); i++) {
	isp++;
	(*i) = charge[isp] * charge[isp] * dh;
}
//!%-------------------------------------------------------------
//!% Compute activity coefficient for Cl in a KCl system
//!%-------------------------------------------------------------
if (ismacinnes) gclm = gclm_ (dh);
//!%-------------------------------------------------------------
}
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
double ActivityModelPitzer::gclm_(const double& dhterm) {

const double r0(0.0e0), rhalf(0.5e0), r1(1.0e0), r2(2.0e0), ionstr(ActivityModel::I_),
		     mtb0kcl(0.04835e0), mtb1kcl(0.2122e0), mtc0kcl(-0.00084e0);
double x(r0), xxx(r0), yyy(r0);
double gclm(r0);

x=r2*sqrt(ionstr);
xxx=-r2*(r1-(r1+x+rhalf*x*x)*exp(-x))/(x*x);
xxx*=mtb1kcl/ionstr;
yyy=r2*(r1-(r1+x)*exp(-x))/(x*x);
yyy*=mtb1kcl;
yyy+=mtb0kcl;
gclm = dhterm+ionstr*ionstr*xxx+ionstr*(r2*yyy+ionstr*mtc0kcl)+ionstr*ionstr*mtc0kcl/r2;
gclm = exp(gclm);
return gclm;
//!%------------------------------------------------------------
}
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
void ActivityModelPitzer::ComputeQmatrices(std::vector<double>& q,
										   std::vector<double>& qphi,
										   std::vector<double>& qpri,
										   std::vector<std::vector<int> >& indnzq) {

const double r0(0.0e0), rhalf(0.5e0), r2(2.0e0), r3(3.0e0), r4(4.0e0), r8(8.0e0), r16(16.0e0);
int i(0), j(0), k(0), nz_b(0), nz_lam(0), nz_th(0), funij(0), funii(0), funjj(0);
double B0ij(r0), B1ij(r0), B2ij(r0), A1_2stri(r0), A2_2stri(r0), expo1(r0), expo2(r0),
       G1(r0), G2(r0), Gp1(r0), Gp2(r0), G2p1(r0), G2p2(r0), t1(r0), t2(r0),
       thij(r0), jij(r0), jii(r0), jjj(r0), jpij(r0), jpii(r0), jpjj(r0),
       j2pij(r0), j2pii(r0), j2pjj(r0), zizj(r0), zizi(r0), zjzj(r0),
       xij(r0), eth(r0), eth_i(r0), eth_i2(r0), ethpri(r0), ethpri_i(r0),
       ttij(r0), ttii(r0), ttjj(r0), eth2pri(r0),
       xii(r0), xjj(r0);
const double ctw1(ActivityModel::M_/cwater), ctw2(r2/cwater);

std::vector<std::vector<double> > g_;
std::vector<std::vector<double> > g_pri_;
std::vector<std::vector<double> > g_2pri_;
std::vector<std::vector<double> > f_;
g_.resize(nfunb);
g_pri_.resize(nfunb);
g_2pri_.resize(nfunb);
f_.resize(nfunb);
for (int i=0; i<nfunb; i++) {
	g_[i].resize(2,r0);
	g_pri_[i].resize(2,r0);
	g_2pri_[i].resize(2,r0);
	f_[i].resize(2,r0);
	g_[i][0]=r0;
	g_[i][1]=r0;
	g_pri_[i][0]=r0;
	g_pri_[i][1]=r0;
	g_2pri_[i][0]=r0;
	g_2pri_[i][1]=r0;
	f_[i][0]=r0;
	f_[i][1]=r0;
}
// Local vectors for J's functions
std::vector<double> j_;
std::vector<double> j_pri_;
std::vector<double> j_2pri_;
j_.resize(nfunj,r0);
j_pri_.resize(nfunj,r0);
j_2pri_.resize(nfunj,r0);
for (int i=0; i<nfunj; i++) {
	j_[i]=r0;
	j_pri_[i]=r0;
	j_2pri_[i]=r0;
}
//!%--------------------------------------------------------------------
//!% Compute constant values
//!%--------------------------------------------------------------------
const double str(ActivityModel::I_);
const double stri(sqrt(str));
const double str2(str*str);
const double str3(str2*str);
const double afi_stri(2.352e0*stri);
//!%--------------------------------------------------------------------
//!% Computes g functions (eq. A6, A7 and B5)
//!%--------------------------------------------------------------------
ComputeFbeta(f_, g_, g_pri_, g_2pri_);
//!%--------------------------------------------------------------------
//!% Computes j functions (eq. A14, A15, B9, B10, B11 and B12)
//!%--------------------------------------------------------------------
ComputeFj(j_, j_pri_, j_2pri_);
//!%--------------------------------------------------------------------
int nz(-1);
//!%--------------------------------------------------------------------
//!% Computes Q, Q', Qphi, dQphi and dQ' matrices
//!%--------------------------------------------------------------------
for (nz_b=0; nz_b<nnzbeta; nz_b++) {

      nz++;

      i=indnzbeta[0][nz_b];
      j=indnzbeta[1][nz_b];
      k=indnzbeta[2][nz_b];
//!%-------------
      B0ij=beta0[nz_b];
      B1ij=beta1[nz_b];
      B2ij=beta2[nz_b];
//!%-------------
      A1_2stri=alpha[0][k]/(r2*stri);
      A2_2stri=alpha[1][k]/(r2*stri);
//!%-------------
      expo1=f_[k][0];
      expo2=f_[k][1];
//!%-------------
      G1=g_[k][0];
      G2=g_[k][1];
//!%------------
      Gp1=g_pri_[k][0];
      Gp2=g_pri_[k][1];
//!%-------------
      G2p1=g_2pri_[k][0];
      G2p2=g_2pri_[k][1];
//!%-------------------------------
//!% Eq. A5
//!%-------------------------------
      qphi[nz]=B0ij+B1ij*expo1+B2ij*expo2;
//!%-------------------------------
//!%  Eq. A3
//!%-------------------------------
      q[nz]=B0ij+B1ij*G1+B2ij*G2;
//!%-------------------------------
//!%  Eq. A4
//!%-------------------------------
      qpri[nz]=B1ij*(Gp1/str)+B2ij*(Gp2/str);

      t1=(str*G2p1*A1_2stri-Gp1)/str2;
      t2=(str*G2p2*A2_2stri-Gp2)/str2;
//!%-------------------------------
//!% Eq. B3
//!%-------------------------------
      indnzq[0][nz]=i;
      indnzq[1][nz]=j;


}
//!%--------------------------------------------------------------------
//!%--------------------------------------------------------------------
//!%--------------------------------------------------------------------
//!% Anions-Anions, Cations-Cations terms
//!%--------------------------------------------------------------------
//!%--------------------------------------------------------------------
//!%--------------------------------------------------------------------
for (nz_th=0; nz_th<nnztheta; nz_th++) {
//!%--------------------------------------------------------------------
      nz++;
//!%--------------------------------------------------------------------
      i=indnztheta[0][nz_th];
      j=indnztheta[1][nz_th];
      funij=indnztheta[2][nz_th];
      funii=indnztheta[3][nz_th];
      funjj=indnztheta[4][nz_th];
      thij=theta[nz_th];
//!%--------------------------------------------------------------------
      jij = j_[funij];
      jii=j_[funii];
      jjj=j_[funjj];
      jpij=j_pri_[funij];
      jpii=j_pri_[funii];
      jpjj=j_pri_[funjj];
      j2pij=j_2pri_[funij];
      j2pii=j_2pri_[funii];
      j2pjj=j_2pri_[funjj];
//!%--------------------------------------------------------------------
      zizj=zprod[funij];
      zizi=zprod[funii];
      zjzj=zprod[funjj];
      xij=afi_stri*zizj;
      xii=afi_stri*zizi;
      xjj=afi_stri*zjzj;
//!%-----------------------------
//!%  Eq. A12
//!%-----------------------------
      eth=(zizj/(r4*str))*(jij-rhalf*jii-rhalf*jjj);
      q[nz]=thij+eth;
      eth_i=eth/str;
      eth_i2=eth/str2;
//!%----------------------------
//!% Eq. A13
//!%----------------------------
      ethpri=-eth_i+(zizj/(r8*str2))*(xij*jpij-rhalf*xii*jpii-rhalf*xjj*jpjj);
      qpri[nz]=ethpri;
      qphi[nz]=thij+eth+str*ethpri;
      ethpri_i=ethpri/str;
//!%----------------------------
//!%  Eq. B8
//!%----------------------------
      ttij=(jpij+xij*j2pij)*xij;
      ttii=(jpii+xii*j2pii)*xii;
      ttjj=(jpjj+xjj*j2pjj)*xjj;
      eth2pri=-r3*ethpri_i-eth_i2+(zizj/(r16*str3))*(ttij-rhalf*ttii-rhalf*ttjj);
      indnzq[0][nz]=i;
      indnzq[1][nz]=j;

}
//!%--------------------------------------------------------------------
//!%--------------------------------------------------------------------
//!%--------------------------------------------------------------------
//!%--------------------------------------------------------------------
//!% Neutral ions-cations and Neutral ions-Anions terms
//!%--------------------------------------------------------------------
//!%--------------------------------------------------------------------
//!%--------------------------------------------------------------------
//!%--------------------------------------------------------------------
for (nz_lam=0; nz_lam<nnzlamda; nz_lam++) {
 nz++;
 i=indnzlamda[0][nz_lam];
 j=indnzlamda[1][nz_lam];
//!%-------------------------------
//!% Add in Q matrix (mat. 12)
//!% the Lij terms
//!%-------------------------------
 q[nz]=lamda[nz_lam];
 indnzq[0][nz]=i;
 indnzq[1][nz]=j;
}
//!%--------------------------------------------------------------------
//!% Check the number of Q'terms
//!%--------------------------------------------------------------------
if ((nz+1)!=nnzq)
std::cout << "Warning different number non-zero terms in Q's matrices" << std::endl;
//!%--------------------------------------------------------------------
}
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
void ActivityModelPitzer::ComputeFbeta(std::vector<std::vector <double> >& f_,
		                               std::vector<std::vector <double> >& g_,
		                               std::vector<std::vector <double> >& g_pri_,
		                               std::vector<std::vector <double> >& g_2pri_) {
//!%--------------------------------------------------------------------
const double r0(0.0e0), r1(1.0e0), r2(2.0e0), r3(3.0e0), r4(4.0e0);
const double stri(sqrt(ActivityModel::I_));
double alfa1(r0), alfa2(r0), x1(r0), x2(r0), x1q(r0), x2q(r0), x1c(r0), x2c(r0),
	   expo1(r0), expo2(r0);
int j(0);

for (j=0; j<nfunb; j++) {

	   alfa1=alpha[0][j];
	   alfa2=alpha[1][j];

	   x1=alfa1*stri;
	   x1q=x1*x1;
	   x1c=x1q*x1;

	   expo1=exp(-x1);

	   f_[j][0]=expo1;
//!%---------------------------------------------------------------------
//!% eq. A6
//!%---------------------------------------------------------------------
	   g_[j][0]=r2*(r1-(r1+x1)*expo1)/x1q;
//!%---------------------------------------------------------------------
//!% eq. A7
//!%---------------------------------------------------------------------
	   g_pri_[j][0]=-r2*(r1-(r1+x1+(x1q/r2))*expo1)/x1q;
//!%---------------------------------------------------------------------
//!% eq. B5
//!%---------------------------------------------------------------------
	   g_2pri_[j][0]=r4*(r1-(r1+x1+(x1q/r2)+(x1c/r4))*expo1)/x1c;
//!%---------------------------------------------------------------------
//!% when the ions aren't monovalents
//!%---------------------------------------------------------------------

	   if(alfa2!=r0) {
	     x2=alfa2*stri;
	     x2q=x2*x2;
	     x2c=x2q*x2;

	     expo2=exp(-x2);
	     f_[j][1]=expo2;
//!%---------------------------------------------------------------------
//!% eq. A6
//!%---------------------------------------------------------------------
	     g_[j][1]=r2*(r1-(r1+x2)*expo2)/x2q;
//!%---------------------------------------------------------------------
//!% eq. A7
//!%---------------------------------------------------------------------
	     g_pri_[j][1]=-r2*(r1-(r1+x2+(x2q/r2))*expo2)/x2q;
//!%---------------------------------------------------------------------
//!% eq. B5
//!%---------------------------------------------------------------------
	     g_2pri_[j][1]=r4*(r1-(r1+x2+(x2q/r2)+(x2c/r4))*expo2)/x2c;
	   }
//!%------------------------------------------------------------
}
}
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
void ActivityModelPitzer::ComputeFj(std::vector<double>& j_,
		                            std::vector<double>& j_pri_,
		                            std::vector<double>& j_2pri_) {
//!%--------------------------------------------------------------------
const double r0(0.0e0), r1(1.0e0), r2(2.0e0), r3(3.0e0), r4(4.0e0), r10(10.0e0), r20(20.0e0),
             r50(50.0e0), r200(200.0e0), rhalf(0.5e0), r1_5(1.5e0), r80(80.0e0), r6(6.0e0),
             r7(7.0e0);
const double e1(4.581e0), e2(0.7237e0),
		     e3(0.0120e0), e4(0.528e0), e12(7.8963e0),
		     e13(0.029025e0), e14(0.1957e0), afi_stri(2.352e0*sqrt(ActivityModel::I_));
double x2(r0), x3(r0), x4(r0), x(r0), expo(r0), j2p(r0), termxc4(r0), jp_j(r0),
	   jp2_j(r0), s1(r0), s2(r0), s3(r0), s1q(r0), s3q(r0), s1c(r0), zizj(r0),
	   logar(r0), xc4(r0), xc2(r0), td(r0), td1(r0);
int i(0), j(0), k(0);

double d[6];

d[1]=4.118e0;
d[2]=7.247e0;
d[3]=-4.408e0;
d[4]=1.837e0;
d[5]=-0.251e0;
d[6]=0.0164e0;
//!%------------------------------------------------------------------
//!% Allocate pointers
//!%------------------------------------------------------------------
for (i=0;i<nfunj; i++) {

	zizj = zprod[i];

	x=zizj*afi_stri;
	x2=x*x;
	x3=x2*x;
	x4=x3*x;

	//!cprovi if (x<=0.03d0) then    !For x<= 0.03
	if (x<=r80) {

	//!%------------------------------------------------------------------
	//!% Compute s1=sum1 (C(k)/x^k)
	//!% Compute s2=suma((d(k)*k*(k+1.d0)/(x**(k+2.d0))))
	//!% Compute s3=suma((d(k)*k/(x**(k+1.d0))))
	//!%------------------------------------------------------------------

	     s1=d[6]/x;
	     s3=r6*s1;
	     s2=r7*s3;

	     for (k=5; k>=1; --k) {
	       s1=(s1+d[k])/x;
	       s2=(s2+k*(k+1)*d[k])/x;
	       s3=(s3+k*d[k])/x;
	     }


	      s2=s2/x2;
	      s3=s3/x;
	      s1q=s1*s1;
	      s1c=s1q*s1;
	      s3q=s3*s3;

	      logar=log(x);
	      expo=exp(-r10*x2);
	//!%-------------------------------------
	//!% A14
	//!%-------------------------------------
	      j_[i]=-(r1/r6)*x2*logar*expo+(1/s1);
	//!%-------------------------------------
	//!% eq.B9
	//!%-------------------------------------
	      j_pri_[i]=((r10*x2-r1)*logar-rhalf)*(x/r3)*expo+(s3/s1q);
	//!%-------------------------------------
	//!% eq.B10
	//!%-------------------------------------
	      j2p=-(s2/s1q)+(r2*s3q/s1c);
	      j2p+=((r50*x2-r200*x4-r1)*logar+r20*x2-r1_5)*expo/r3;
	      j_2pri_[i]=j2p;
	//!%----------------------------------

	} else {

	      xc4=pow(x,e4);
	      xc2=pow(x,e2);
	      expo=exp(-e3*xc4);
	      td1=(e1/xc2)*expo;
	      td=r4+td1;
	//!%-------------------------------------
	//!% eq. A15
	//!%-------------------------------------
	      j_[i]=x/td;
	//!%-------------------------------------
	//!% eq. B11
	//!%-------------------------------------
	      j_pri_[i]=(j_[i]/x2)*(x+td1*(e2+e3*e4*xc4)*j_[i]);
	//!%-------------------------------------
	//!% eq. B12
	//!%-------------------------------------
	      termxc4=e3*e4*xc4;
	      jp_j=j_pri_[i]/j_[i];
	      jp2_j=jp_j*j_pri_[i];
	      j2p=e4-e2-termxc4;
	      j2p=j2p*termxc4/(e2+termxc4)+(jp_j*x-e2);
	      j2p=j2p*(jp_j*x-r1)+r1;
	      j_2pri_[i]=j2p*j_[i]/x2-r2*j_pri_[i]/x+jp2_j;
	//!%-------------------------------------

	}



}
}
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
void ActivityModelPitzer::Display(void) const {
  std::cout << "Activity model: Pitzer" << std::endl;
}  // end Display()
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
void ActivityModelPitzer::ReadDataBase(const std::string& namedatabase,
		                               const std::vector<Species>& prim,
		                               const std::vector<AqueousEquilibriumComplex>& sec) {
	const int mxlines(10000);
	const int block_b0(0), block_b1(1), block_b2(2), block_cfi(3), block_theta(4), block_lamda(5), block_psi(7), block_exit(8);
    const std::string longspace("                  ");
	// Create a vector of species

	    int npri=prim.size();
	    int naqx=sec.size();
	    std::vector<std::string> namesp;
	    std::vector<double> charge;
	    if ((npri+naqx)==0) {
	       std::cout << "Error, zero number of species" << std::endl;
	    }
	    //---------------------------------------------------------
	    // Is better to use iterators
	    //---------------------------------------------------------
	    for (std::vector<Species>::const_iterator i = prim.begin();
	           i != prim.end(); i++) {
	         charge.push_back((*i).charge());
	         namesp.push_back((*i).name());
	    }
	    for (std::vector<AqueousEquilibriumComplex>::const_iterator i = sec.begin();
	    	           i != sec.end(); i++) {
	    	         charge.push_back((*i).charge());
	    	         namesp.push_back((*i).name());
	    	    }

//---------------------------------------------------------------------------------------------
// Open Pitzer virial coefficients database
//---------------------------------------------------------------------------------------------
std::cout << "=================> Opening Pitzer Data Base ==============>" << namedatabase << std::endl;
std::ifstream database(namedatabase.c_str());
if (!database) {
	    std::ostringstream error_stream;
	    error_stream << "SimpleThermoDatabase::ReadFile(): \n";
	    error_stream << "file could not be opened.... " << namedatabase << "\n";
} else {
		std::cout << "Opening Pitzer virial coefficients database was succefully opened: " << namedatabase << std::endl;
}

std::string space(" ");
StringTokenizer no_spaces;
;

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
	    	if (iblock==block_b0) std::cout << "=================> Parse B0 ==============>" << std::endl;
	    	if (iblock==block_b1) std::cout << "=================> Parse B1 ==============>" << std::endl;
	    	if (iblock==block_b2) std::cout << "=================> Parse B2 ==============>" << std::endl;
	    	if (iblock==block_cfi) std::cout << "=================> Parse Cfi ==============>" << std::endl;
	    	if (iblock==block_theta) std::cout << "=================> Parse Theta ==============>" << std::endl;
	    	if (iblock==block_lamda) std::cout << "=================> Parse Lamda ==============>" << std::endl;
	    	if (iblock==block_psi) std::cout << "=================> Parse Psi ==============>" << std::endl;
	    	getline(database,line);
	    }
	    if (iblock==block_b0) {
	    	ParseB0(line,namesp);
	    }
	    if (iblock==block_b1) {
	       	ParseB1(line,namesp);
	    }
	    if (iblock==block_b2) {
	    	ParseB2(line,namesp);
	    }
	    if (iblock==block_cfi) {
	        ParseCfi (line,namesp,charge);
	    }
	    if (iblock==block_theta) {
	    	ParseTheta (line,namesp,charge);
	    }
	    if (iblock==block_lamda) {
	    	ParseLamda (line,namesp);
	    }
	    if (iblock==block_psi) {
	    	ParsePsi (line,namesp);
	    }

        if (iblock==block_exit) goto exit;

}

exit:

//---------------------------------------------
// Assign Beta's function indices
//---------------------------------------------
AssignFbeta(charge);
//---------------------------------------------
// Find the new non-zero theta terms and
// Assign J's function indices
//---------------------------------------------
AssignFj(charge);
//---------------------------------------------
// Compute the total number of non-zeta terms
//---------------------------------------------
nnzq=nnzbeta+nnztheta+nnzlamda;
//---------------------------------------------
// Store the indices for H2O, Cl and K
//---------------------------------------------
int isp(-1);
for (std::vector<std::string>::iterator i=namesp.begin(); i!=namesp.end(); i++) {
	isp++;
	if ((*i)=="H2O" || (*i)=="h2o") ithw=isp;
	if ((*i)=="Cl-" || (*i)=="cl-") ithcl=isp;
	if ((*i)=="K+" || (*i)=="k+") ithk=isp;
}
//---------------------------------------------
// Close Pitzer virial database
//---------------------------------------------
database.close();
std::cout << "=================> Closing Pitzer Data Base ==============>" << std::endl;
//------------------------------------------------------------
// Print the virial coeffients
//------------------------------------------------------------
std::cout << "=================> Print stored virial coefficients ==============>" << std::endl;
Print(namesp);
std::cout << "=================> End print stored virial coefficients ==============>" << std::endl;
}  // end ReadDataBase()
//---------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------
void ActivityModelPitzer::Print(const std::vector<std::string>& namesp) {
const double r0(0.0e0);
std::cout << "--------------------------------------------------------------------" << std::endl;
std::cout << " Virial coefficients" << std::endl;
std::cout << "--------------------------------------------------------------------" << std::endl;
int isp1(-1),isp2(-1), isp3(-1), nvirial(0);
if (nnzbeta>0){
std::cout << "=================> B0 ==============>" << std::endl;
for (int i=0; i<nnzbeta; i++){
  isp1=indnzbeta[0][i];
  isp2=indnzbeta[1][i];
  if (beta0[i]!=r0){
	  nvirial++;
  	  std::cout << namesp[isp1] << "  " << namesp[isp2] << "  " << beta0[i] << std::endl;
    }
}
std::cout << "=================> B1 ==============>" << std::endl;
for (int i=0; i<nnzbeta; i++){
  isp1=indnzbeta[0][i];
  isp2=indnzbeta[1][i];
  if (beta1[i]!=r0){
	  nvirial++;
  	  std::cout << namesp[isp1] << "  " << namesp[isp2] << "  " << beta1[i] << std::endl;
    }
}
std::cout << "=================> B2 ==============>" << std::endl;
for (int i=0; i<nnzbeta; i++){
  isp1=indnzbeta[0][i];
  isp2=indnzbeta[1][i];
  if (beta2[i]!=r0){
	  nvirial++;
	  std::cout << namesp[isp1] << "  " << namesp[isp2] << "  " << beta2[i] << std::endl;
  }

}
}
if (nnzcpz>0){
std::cout << "=================> Cfi ==============>" << std::endl;
for (int i=0; i<nnzcpz; i++){
  isp1=indnzcpz[0][i];
  isp2=indnzcpz[1][i];
  if (cpz[i]!=r0){
	  nvirial++;
	  std::cout << namesp[isp1] << "  " << namesp[isp2] << "  " << cpz[i] << std::endl;
  }

}
}
if (nnztheta>0){
std::cout << "=================> Theta ==============>" << std::endl;
for (int i=0; i<nnztheta; i++){
  isp1=indnztheta[0][i];
  isp2=indnztheta[1][i];
  if (theta[i]!=r0){
	  nvirial++;
	  std::cout << namesp[isp1] << "  " << namesp[isp2] << "  " << theta[i] << std::endl;
  }

}
}
if (nnzlamda>0){
std::cout << "=================> Lamda ==============>" << std::endl;
for (int i=0; i<nnzlamda; i++){
  isp1=indnzlamda[0][i];
  isp2=indnzlamda[1][i];
  if (lamda[i]!=r0){
	  nvirial++;
	  std::cout << namesp[isp1] << "  " << namesp[isp2] << "  " << lamda[i] << std::endl;
  }

}
}
if (nnzpsi>0){
std::cout << "=================> Psi ==============>" << std::endl;
for (int i=0; i<nnzpsi; i++){
  isp1=indnzpsi[0][i];
  isp2=indnzpsi[1][i];
  isp3=indnzpsi[2][i];
  if (psi[i]!=r0){
	  nvirial++;
	  std::cout << namesp[isp1] << "  " << namesp[isp2] << "  " << namesp[isp3] << "  " << psi[i] << std::endl;
  }

}
}
std::cout << "Total number of virial coefficients: " << nvirial << std::endl;
}  // Print
//---------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------
void ActivityModelPitzer::ParseB0 (const std::string& data, std::vector<std::string>& namesp) {
// Local variables and constants
const double r0(0.0e0);
bool isstored(false);
bool isspecies(false);
bool isdebug(false);
std::string semicolon(";");
std::string space(" ");
StringTokenizer b0(data,space);
StringTokenizer no_spaces;

int nsp(namesp.size());
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

if (nnzbeta>0) {
	for (int i=0; i<nnzbeta; i++) {
	   	isp1_old = indnzbeta[0][i];
	   	isp2_old = indnzbeta[1][i];
	   	if ((isp1==isp1_old && isp2==isp2_old) || (isp1==isp2_old && isp2==isp1_old)) {
	   		beta0[i]=virial;
	   		isstored=true;
	   	}
	}
    // If the virial coefficient is not stored
	    if (!isstored) {
	//!%--------------
	      nnzbeta++;
	      beta0.push_back(virial);
	      indnzbeta[0].push_back(isp1);
	      indnzbeta[1].push_back(isp2);
	//!%--------------
	      beta1.push_back(r0);
	//!%--------------
	      beta2.push_back(r0);
	    }

} else {
		nnzbeta = 1;
		beta0.push_back(virial);
		//--------------------------
		beta1.push_back(r0);
		beta2.push_back(r0);
		//--------------------------
		indnzbeta.resize(3);
		for (int i=0; i<3; i++) {
			indnzbeta[i].resize(nnzbeta,0);
			for (int j=0; j<nnzbeta; j++) indnzbeta[i][j]=0;
		}
		indnzbeta[0][nnzbeta-1]=isp1;
		indnzbeta[1][nnzbeta-1]=isp2;

}
    if (isdebug) std::cout << namesp[isp1] << namesp[isp2] << virial <<std::endl;

}
}  // ParseB0
//---------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------
void ActivityModelPitzer::ParseB1 (const std::string& data, std::vector<std::string>& namesp) {
	// Local variables and constants
	const double r0(0.0e0);
	bool isstored(false);
	bool isspecies(false);
	bool isdebug(false);
	std::string semicolon(";");
	std::string space(" ");
	StringTokenizer b1(data,space);
	StringTokenizer no_spaces;

	int nsp(namesp.size());
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

	if (nnzbeta>0) {
		for (int i=0; i<nnzbeta; i++) {
		   	isp1_old = indnzbeta[0][i];
		   	isp2_old = indnzbeta[1][i];
		   	if ((isp1==isp1_old && isp2==isp2_old) || (isp1==isp2_old && isp2==isp1_old)) {
		   		beta1[i]=virial;
		   		isstored=true;
		   	}
		}
	    // If the virial coefficient is not stored
		    if (!isstored) {
		//!%--------------
		      nnzbeta++;
		      beta1.push_back(virial);
		      indnzbeta[0].push_back(isp1);
		      indnzbeta[1].push_back(isp2);
		//!%--------------
		      beta0.push_back(r0);
		//!%--------------
		      beta2.push_back(r0);
		    }

	} else {
			nnzbeta = 1;
			beta1.push_back(virial);
			beta0.push_back(r0);
			beta2.push_back(r0);
			indnzbeta.resize(3);
			for (int i=0; i<3; i++) {
				indnzbeta[i].resize(nnzbeta,0);
				for (int j=0; j<nnzbeta; j++) indnzbeta[i][j]=0;
			}
			indnzbeta[0][nnzbeta-1]=isp1;
			indnzbeta[1][nnzbeta-1]=isp2;

}
	if (isdebug) std::cout << namesp[isp1] << namesp[isp2] << virial <<std::endl;
}
}
//---------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------
void ActivityModelPitzer::ParseB2 (const std::string& data, std::vector<std::string>& namesp) {
	// Local variables and constants
	const double r0(0.0e0);
	bool isstored(false);
	bool isspecies(false);
	bool isdebug(false);
	std::string semicolon(";");
	std::string space(" ");
	StringTokenizer b2(data,space);
	StringTokenizer no_spaces;

	int nsp(namesp.size());
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

	if (nnzbeta>0) {
		for (int i=0; i<nnzbeta; i++) {
		   	isp1_old = indnzbeta[0][i];
		   	isp2_old = indnzbeta[1][i];
		   	if ((isp1==isp1_old && isp2==isp2_old) || (isp1==isp2_old && isp2==isp1_old)) {
		   		beta2[i]=virial;
		   		isstored=true;
		   	}
		}
	    // If the virial coefficient is not stored
		    if (!isstored) {
		//!%--------------
		      nnzbeta++;
		      beta2.push_back(virial);
		      indnzbeta[0].push_back(isp1);
		      indnzbeta[1].push_back(isp2);
		//!%--------------
		      beta0.push_back(r0);
		//!%--------------
		      beta1.push_back(r0);
		    }

	} else {
			nnzbeta = 1;
			beta2.push_back(virial);
			beta0.push_back(r0);
			beta1.push_back(r0);
			indnzbeta.resize(3);
			for (int i=0; i<3; i++) {
				indnzbeta[i].resize(nnzbeta,0);
				for (int j=0; j<nnzbeta; j++) indnzbeta[i][j]=0;
			}
			indnzbeta[0][nnzbeta-1]=isp1;
			indnzbeta[1][nnzbeta-1]=isp2;
}
	if (isdebug) std::cout << namesp[isp1] << namesp[isp2] << virial <<std::endl;
}
}
//---------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------
void ActivityModelPitzer::ParseCfi (const std::string& data, std::vector<std::string>& namesp, std::vector<double>& charge) {
	// Local variables and constants
	const double r0(0.0e0), r1(1.0e0), r2(2.0e0);
	bool isstored(false);
	bool isspecies(false);
	bool isdebug(false);
	std::string semicolon(";");
	std::string space(" ");
	StringTokenizer cfi(data,space);
	StringTokenizer no_spaces;

	int nsp(namesp.size());
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



	if (nnzcpz>0) {
		//!%--------------
   	          nnzcpz++;
   	          virial/=r2*sqrt(abs(charge[isp1]*charge[isp2]));
   	          cpz.push_back(virial);
		      indnzcpz[0].push_back(isp1);
		      indnzcpz[1].push_back(isp2);


	} else {
		    nnzcpz = 1;
		    virial/=r2*sqrt(abs(charge[isp1]*charge[isp2]));
		    cpz.push_back(virial);
			indnzcpz.resize(2);
			for (int i=0; i<2; i++) {
				indnzcpz[i].resize(nnzcpz,0);
				for (int j=0; j<nnzcpz; j++) indnzcpz[i][j]=0;
			}
			indnzcpz[0][nnzcpz-1]=isp1;
			indnzcpz[1][nnzcpz-1]=isp2;
}
	if (isdebug) std::cout << namesp[isp1] << namesp[isp2] << virial <<std::endl;
}
}
//---------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------
void ActivityModelPitzer::ParseTheta (const std::string& data, std::vector<std::string>& namesp, std::vector<double>& charge) {
	// Local variables and constants
	const double r0(0.0e0), r1(1.0e0), r2(2.0e0);
	const bool isdebug(false);
	bool isstored(false);
	bool isspecies(false);
	std::string semicolon(";");
	std::string space(" ");
	StringTokenizer th(data,space);
	StringTokenizer no_spaces;

	int nsp(namesp.size());
	int isp1_old(0), isp2_old(0);

	// get name 1
	no_spaces.tokenize(th.at(0),space);
	std::string name1(no_spaces.at(0));
	// get name 2
	no_spaces.tokenize(th.at(1),space);
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

	no_spaces.tokenize(th.at(2), space);
	double virial(std::atof(no_spaces.at(0).c_str()));

	if (nnztheta>0) {

		nnztheta++;
		theta.push_back(virial);
		indnztheta[0].push_back(isp1);
		indnztheta[1].push_back(isp2);
		if (charge[isp1]==r0 || charge[isp2]==r0) std::cout << "Warning electric charge must be different than zero" << std::endl;
		indnztheta[2].push_back(charge[isp1]*charge[isp2]);
		indnztheta[3].push_back(charge[isp1]*charge[isp1]);
		indnztheta[4].push_back(charge[isp2]*charge[isp2]);



	} else {
		    nnztheta = 1;
		    theta.push_back(virial);
			indnztheta.resize(5);
			for (int i=0; i<5; i++) {
				indnztheta[i].resize(nnztheta,0);
				for (int j=0; j<nnztheta; j++) indnztheta[i][j]=0;
			}
			indnztheta[0][nnztheta-1]=isp1;
			indnztheta[1][nnztheta-1]=isp2;
			if (charge[isp1]==r0 || charge[isp2]==r0) std::cout << "Warning electric charge must be different than zero" << std::endl;
			indnztheta[2][nnztheta-1]=charge[isp1]*charge[isp2];
			indnztheta[3][nnztheta-1]=charge[isp1]*charge[isp1];
			indnztheta[4][nnztheta-1]=charge[isp2]*charge[isp2];
}
	if (isdebug){
		std::cout << namesp[isp1] << namesp[isp2] << virial <<std::endl;
		for (int i=0; i<nnztheta; i++){
		    std::cout << "zi zj" << indnztheta[2][i] << std::endl;
		    std::cout << "zi zi" << indnztheta[3][i] << std::endl;
		    std::cout << "zj zj" << indnztheta[4][i] << std::endl;
	     }
   }
}
}
//---------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------
void ActivityModelPitzer::ParseLamda (const std::string& data, std::vector<std::string>& namesp) {
	// Local variables and constants
	const double r0(0.0e0), r1(1.0e0), r2(2.0e0);
	bool isstored(false);
	bool isspecies(false);
	bool isdebug(false);
	std::string semicolon(";");
	std::string space(" ");
	StringTokenizer lam(data,space);
	StringTokenizer no_spaces;

	int nsp(namesp.size());
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



	if (nnzlamda>0) {
		//!%--------------
   	          nnzlamda++;
		      lamda.push_back(virial);
		      indnzlamda[0].push_back(isp1);
		      indnzlamda[1].push_back(isp2);


	} else {
		    nnzlamda = 1;
		    lamda.push_back(virial);
			indnzlamda.resize(2);
			for (int i=0; i<2; i++) {
				indnzlamda[i].resize(nnzlamda,0);
				for (int j=0; j<nnzlamda; j++) indnzlamda[i][j]=0;
			}
			indnzlamda[0][nnzlamda-1]=isp1;
			indnzlamda[1][nnzlamda-1]=isp2;
}
	if (isdebug) std::cout << namesp[isp1] << namesp[isp2] << virial <<std::endl;

}
}
//---------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------
void ActivityModelPitzer::ParsePsi (const std::string& data, std::vector<std::string>& namesp) {
	// Local variables and constants
	const double r0(0.0e0), r1(1.0e0), r2(2.0e0);
	bool isstored(false);
	bool isspecies(false);
	bool isdebug(false);
	std::string semicolon(";");
	std::string space(" ");
	StringTokenizer psi_loc(data,space);
	StringTokenizer no_spaces;

	int nsp(namesp.size());
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



	if (nnzpsi>0) {

   	          nnzpsi++;
   	          psi.push_back(virial);
		      indnzpsi[0].push_back(isp1);
		      indnzpsi[1].push_back(isp2);
		      indnzpsi[2].push_back(isp3);


	} else {
		    nnzpsi = 1;
		    psi.push_back(virial);
			indnzpsi.resize(3);
			for (int i=0; i<3; i++) {
				indnzpsi[i].resize(nnzpsi,0);
				for (int j=0; j<nnzpsi; j++) indnzpsi[i][j]=0;
			}
			indnzpsi[0][nnzpsi-1]=isp1;
			indnzpsi[1][nnzpsi-1]=isp2;
			indnzpsi[2][nnzpsi-1]=isp3;
}
	if (isdebug) std::cout << namesp[isp1] << namesp[isp2] << namesp[isp3] << virial << std::endl;
}
}
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
void ActivityModelPitzer::AssignFbeta (const std::vector<double>& charge) {
// Local variables and constants
int l1(0), l2(0), l3(0), n1(0), n2(0), n3(0), nz(0);
double z1(0.0e0), z2(0.0e0);
const double r0(0.0e0), r1(1.0e0), r2(2.0e0), r50(50.0e0), r1_4(1.4e0), r12(12.0e0);
const int ndim1(2), ndim2(100);
//!%-----------------------------------------------------------
//! Assign indices of Beta's functions
//!%-----------------------------------------------------------
nfunb=0;
std::vector<std::vector<double> > alphaloc;
alphaloc.resize(ndim2);
for (int i = 0; i < ndim1; ++i){
	 alphaloc[i].resize(ndim2);
	 for (int j=0; j<ndim2; j++) alphaloc[i][j]=r0;
}


for (nz=0; nz<nnzbeta; nz++) {
   // Take the electric charge
   z1=charge[indnzbeta[0][nz]];
   z2=charge[indnzbeta[1][nz]];
   if (abs(z1)==r1 || abs(z2)==r1) {
	    	       	if (l1==0) {

						   nfunb++;
						   indnzbeta[2][nz] = nfunb-1;

						   n1 = nfunb-1;

						   alphaloc[0][nfunb-1]=r2;
						   alphaloc[1][nfunb-1]=r0;

						   l1=1;
                   	} else {

						   indnzbeta[2][nz] = n1;

				       }

	       } else {

					 if (abs(z1)!=abs(z2)) {


							 if (l2==0) {
	                           nfunb++;
						       indnzbeta[2][nz] = nfunb-1;

							   n2 = nfunb-1;

							   alphaloc[0][nfunb-1]=r2;
						       alphaloc[1][nfunb-1]=r50;

	                           l2=1;

							 } else {

							         indnzbeta[2][nz]=n2;

	                           }

					 } else {

						     if (l3==0) {
						       nfunb++;
						       indnzbeta[2][nz] = nfunb-1;
                               n3=nfunb-1;

							   alphaloc[0][nfunb-1] = r1_4;
						       alphaloc[1][nfunb-1] = r12;

	                           l3=1;

						     } else {

                               indnzbeta[2][nz]=n3;

						     }

					 }

                   	}


	       }

// Resize alpha array
if (nfunb>0) {
	alpha.resize(2);
	for (int i=0; i<2; i++) alpha[i].resize(nfunb,r0);
	for (int i=0; i<nfunb; i++) {
		alpha[0][i]=alphaloc[0][i];
		alpha[1][i]=alphaloc[1][i];
	}
}

}  // end AssignFbeta
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
void ActivityModelPitzer::AssignFj (const std::vector<double>& charge) {
// Local variables and constants
const int ndim(100);
const double r0(0.0e0);
double zprodloc[ndim];
double zz;
int nn(0), isp1(0), isp2(0);
bool isstored(false);
const int nsp(charge.size());
//---------------------------------------------------------------
// Find other non-zero theta terms
//---------------------------------------------------------------
for (int i=0; i<nsp; i++) {
			for (int j=0; j<nsp; j++) {
			   if (i!=j && ((charge[i]>r0 && charge[j]>r0) || (charge[i]<r0 && charge[j]<r0))) {
					   for (int k=0; k<nnztheta; k++) {
					        isp1=indnztheta[0][k];
					        isp2=indnztheta[1][k];
							if ((isp1==i && isp2==j) || (isp1==j || isp2==i)) {
								isstored=true;
								goto exit_1;
							}
					   }

					   exit_1:

					   if (!isstored) {

						   nnztheta++;

						   if (nnztheta=1){
							   indnztheta.resize(5);
							   for (int k=0; k<5; k++) {
							     indnztheta[k].resize(nnztheta,0);
							     for (int m=0; m<nnztheta; m++) indnztheta[k][m]=0;
							   }
						   }

						   theta.push_back(r0);
						   indnztheta[0].push_back(i);
						   indnztheta[1].push_back(j);
						   if (charge[i]==r0 || charge[j]==r0) std::cout << "Warning electric charge must be different than zero" << std::endl;
						   indnztheta[2].push_back(charge[i]*charge[j]);
						   indnztheta[3].push_back(charge[i]*charge[i]);
						   indnztheta[4].push_back(charge[j]*charge[j]);

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
for (int nz=0; nz<nnztheta; nz++) {
	for (int k=2; k<=4; k++) {
	    zz=double(indnztheta[k][nz]);
	    for (int j=0; j<nfunj; j++) {
	      if (zprodloc[j]==zz) {
	            nn=j;
	            goto exit;
	      }
	    }
	    zprodloc[nfunj-1] = zz;
	    nn=nfunj-1;
	    nfunj++;
	    exit:
	    indnztheta[k][nz]=nn;
	}
}
//!%------------------------------------------------------------
//!% Resize the zprod vector
//!%------------------------------------------------------------
zprod.resize(nfunj,r0);
for (int i=0; i<nfunj; i++)
zprod[i] = zprodloc[i];
}  // end AssignFJ
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
}  // namespace chemistry
}  // namespace amanzi
