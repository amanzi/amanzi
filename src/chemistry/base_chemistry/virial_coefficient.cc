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
		  switch(i) {
		  case 0:
			  virial=pol[i];
			  break;
	      case 1:
			  virial += pol[i]*temp;
		      break;
	      case 2:
			  virial += pol[i]/temp;
			  break;
	      case 3:
	    	  virial += pol[i]*std::log10(temp);
	    	  break;
	      case 4:
	    	  virial += pol[i]/(temp*temp);
		      break;
		  }

	  }

}

}  // namespace chemistry
}  // namespace amanzi
