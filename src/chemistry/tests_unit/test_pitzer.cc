/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>
#include <typeinfo>

#include <UnitTest++.h>

#include "species.hh"
#include "aqueous_equilibrium_complex.hh"
#include "activity_model_factory.hh"
#include "activity_model_pitzer.hh"
#include "activity_model.hh"
#include "chemistry_exception.hh"

using std::vector;
using std::string;

SUITE(TestPitzer) {

namespace ac = amanzi::chemistry;

TEST(System_System_1) {
 ac::ActivityModelFactory amfac_;
 ac::ActivityModel* am_;
 vector<ac::Species> sp_;
 vector<ac::AqueousEquilibriumComplex> aqx_;
 ac::Species H(0, "H+", 1.0, 1.0079, 9.0);
 ac::Species OH(1, "OH-", -1.0, 17.0073, 3.5);
 ac::Species Cl(2, "Cl-", -1.0, 40.0780, 6.0);
 ac::Species Na(3, "Na+", 1.0, 96.0636, 4.0);
 ac::Species K(3, "K+", 1.0, 96.0636, 4.0);
 ac::Species Ca(3, "Ca+2", 2.0, 96.0636, 4.0);
 ac::Species Mg(3, "Mg+2", 2.0, 96.0636, 4.0);
 ac::Species CO3(3, "CO3-2", -2.0, 96.0636, 4.0);
 ac::Species CO2(3, "CO2", 0.0, 96.0636, 4.0);
 ac::Species HCO3(3, "HCO3-", -1.0, 96.0636, 4.0);
 ac::Species MgOH(3, "MgOH+", 1.0, 96.0636, 4.0);
 ac::Species MgCO3(3, "MgCO3", 0.0, 96.0636, 4.0);
 ac::Species CaCO3(3, "CaCO3", 0.0, 96.0636, 4.0);
 ac::Species H2O(3, "H2O", 0.0, 96.0636, 4.0);
 H.update(4.776e-08);
 OH.update(2.570e-07);
 Cl.update(3.0e0);
 Na.update(3.0e0);
 Ca.update(9.516e-02);
 Mg.update(9.711e-02);
 CO3.update(3.339e-03);
 HCO3.update(2.760e-01);
 CO2.update(1.299e-02);
 MgOH.update(1.131e-6);
 MgCO3.update(2.885e-03);
 CaCO3.update(4.838e-3);
 K.update(0.1e0);
 H2O.update(1.0);
 aqx_.clear();
 sp_.clear();
 sp_.push_back(Cl);
 sp_.push_back(Na);
 sp_.push_back(H);
 sp_.push_back(K);
 sp_.push_back(H2O);
 sp_.push_back(Ca);
 sp_.push_back(Mg);
 sp_.push_back(CO3);
 sp_.push_back(HCO3);
 sp_.push_back(CO2);
 sp_.push_back(MgCO3);
 sp_.push_back(CaCO3);
 sp_.push_back(MgOH);
 sp_.push_back(OH);
 vector<double> gamma;
 for (int i=0; i<sp_.size();i++) gamma.push_back(1.0);
 am_= amfac_.Create("pitzer","phreeqc_pitzer.dat",sp_, aqx_);
 am_->Display();
 am_->EvaluateVector(gamma,sp_, aqx_);
 for (int i=0; i<sp_.size();i++) {
	 gamma[i]=log10(gamma[i]);
	//std:: cout << sp_[i].name() << "  " <<gamma[i] << std::endl;
 }
 //std::cout << "Testing coeff. 1" << std::endl;
 // Results are compared with PHREEQC
 CHECK_CLOSE(-0.243, gamma[0], 1.0e-2); // Cl-
 CHECK_CLOSE(-0.027, gamma[1],1.0e-2);  // Na+
 CHECK_CLOSE(0.321, gamma[2], 1.0e-2);  // H+
 CHECK_CLOSE(-0.177, gamma[3], 1.0e-2); // K+
 CHECK_CLOSE(-0.055, gamma[4], 1.0e-2); // H2O
 CHECK_CLOSE(-0.109, gamma[5], 1.0e-2); // Ca
 CHECK_CLOSE(-0.119, gamma[6], 1.0e-2); // Mg
 CHECK_CLOSE(-1.859, gamma[7], 1.0e-2); // CO3
 CHECK_CLOSE(-0.05, gamma[12], 1.0e-2); // MgOH
 CHECK_CLOSE(0.283, gamma[9], 1.0e-2);  // CO2
 CHECK_CLOSE(-0.463, gamma[13], 1.0e-2); // OH
 CHECK_CLOSE(-0.437, gamma[8], 1.0e-2); // HCO3
}

TEST(System_System_2) {
 ac::ActivityModelFactory amfac_;
 ac::ActivityModel* am_;
 vector<ac::Species> sp_;
 vector<ac::AqueousEquilibriumComplex> aqx_;
 ac::Species H(0, "H+", 1.0, 1.0079, 9.0);
 ac::Species OH(1, "OH-", -1.0, 17.0073, 3.5);
 ac::Species Cl(2, "Cl-", -1.0, 40.0780, 6.0);
 ac::Species Na(3, "Na+", 1.0, 96.0636, 4.0);
 ac::Species K(3, "K+", 1.0, 96.0636, 4.0);
 ac::Species Ca(3, "Ca+2", 2.0, 96.0636, 4.0);
 ac::Species Mg(3, "Mg+2", 2.0, 96.0636, 4.0);
 ac::Species CO3(3, "CO3-2", -2.0, 96.0636, 4.0);
 ac::Species CO2(3, "CO2", 0.0, 96.0636, 4.0);
 ac::Species HCO3(3, "HCO3-", -1.0, 96.0636, 4.0);
 ac::Species MgOH(3, "MgOH+", 1.0, 96.0636, 4.0);
 ac::Species MgCO3(3, "MgCO3", 0.0, 96.0636, 4.0);
 ac::Species CaCO3(3, "CaCO3", 0.0, 96.0636, 4.0);
 ac::Species HSO4(3, "HSO4-", -1.0, 96.0636, 4.0);
 ac::Species SO4(3, "SO4-2", -2.0, 96.0636, 4.0);
 ac::Species H2O(3, "H2O", 0.0, 96.0636, 4.0);
 H.update(4.98e-08);
 OH.update(2.574e-07);
 Cl.update(3.0e0);
 Na.update(3.0e0);
 Ca.update(9.545e-02);
 Mg.update(9.717e-02);
 CO3.update(3.382e-03);
 HCO3.update(2.764e-01);
 CO2.update(1.283e-02);
 MgOH.update(1.224e-6);
 MgCO3.update(2.829e-03);
 CaCO3.update(4.55e-3);
 HSO4.update(3.0e-8);
 SO4.update(0.1e0);
 K.update(0.1e0);
 H2O.update(1.0);
 aqx_.clear();
 sp_.clear();
 sp_.push_back(Cl);
 sp_.push_back(Na);
 sp_.push_back(H);
 sp_.push_back(K);
 sp_.push_back(H2O);
 sp_.push_back(Ca);
 sp_.push_back(Mg);
 sp_.push_back(CO3);
 sp_.push_back(HCO3);
 sp_.push_back(CO2);
 sp_.push_back(MgCO3);
 sp_.push_back(CaCO3);
 sp_.push_back(MgOH);
 sp_.push_back(OH);
 sp_.push_back(HSO4);
 sp_.push_back(SO4);
 vector<double> gamma;
 for (int i=0; i<sp_.size();i++) gamma.push_back(1.0);
 am_= amfac_.Create("pitzer","phreeqc_pitzer.dat",sp_, aqx_);
 am_->Display();
 am_->EvaluateVector(gamma,sp_, aqx_);
 for (int i=0; i<sp_.size();i++) {
	 gamma[i]=log10(gamma[i]);
	// std:: cout << sp_[i].name() << "  " <<gamma[i] << std::endl;
 }
 //std::cout << "Testing coeff. 1" << std::endl;
 // Results are compared with PHREEQC
 CHECK_CLOSE(-0.241, gamma[0], 1.0e-2); // Cl-
 CHECK_CLOSE(-0.039, gamma[1],1.0e-2);  // Na+
 CHECK_CLOSE(0.303, gamma[2], 1.0e-2);  // H+
 CHECK_CLOSE(-0.192, gamma[3], 1.0e-2); // K+
 CHECK_CLOSE(-0.055, gamma[4], 1.0e-2); // H2O
 CHECK_CLOSE(-0.14, gamma[5], 1.0e-2); // Ca
 CHECK_CLOSE(-0.131, gamma[6], 1.0e-2); // Mg
 CHECK_CLOSE(-1.862, gamma[7], 1.0e-2); // CO3
 CHECK_CLOSE(-0.096, gamma[12], 1.0e-2); // MgOH
 CHECK_CLOSE(0.291, gamma[9], 1.0e-2);  // CO2
 CHECK_CLOSE(-0.464, gamma[13], 1.0e-2); // OH
 CHECK_CLOSE(-0.435, gamma[8], 1.0e-2); // HCO3
 CHECK_CLOSE(-0.321, gamma[14], 1.0e-2); // HSO4
 CHECK_CLOSE(-1.823, gamma[15], 1.0e-2); // SO4
}
}  // end SUITE(TestPitzer)
