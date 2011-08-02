/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include "virial_coefficient.hh"

#include <cmath>

#include <iostream>

namespace amanzi {
namespace chemistry {

VirialCoefficient::VirialCoefficient()
    :  isp1(-1),
       isp2(-1),
       isp3(-1),
       npol(0),
       ifun1(-1),
       ifun2(-1),
       ifun3(-1),
       virial(0.0e0){

}


VirialCoefficient::~VirialCoefficient(){
}


void VirialCoefficient::Display(){

}

void VirialCoefficient::UpdateVirial(const double& temp, const double& pressure)
{
	  for (int i=0;i<npol;i++){
		  if (i==0){
			  virial=pol[i];
		  } else if (i==1) {
		    virial += pol[i]*temp;
	      }	else if (i==2) {
			virial += pol[i]*temp;
		  } else if (i==3) {
			virial += pol[i]/temp;
		  }	else if (i==3) {
			virial += pol[i]*std::log10(temp);
		  }	else if (i==4) {
			virial += pol[i]/(temp*temp);
		  }

	  }

}

}  // namespace chemistry
}  // namespace amanzi
