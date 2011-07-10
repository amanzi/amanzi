#ifndef _PChemDriver_H_
#define _PChemDriver_H_

#include <string>
#include <map>

#include <Array.H>
#include <FArrayBox.H>

// STRINGIFY is a macro to expand its argument and surround it with strings
#define HACKSTRINGIFY(x) #x
#define STRINGIFY(x) HACKSTRINGIFY(x)

class PChemDriver 
{
public: 

  PChemDriver (int sz, Array<int> num, Array<Real> mw);	  
  
  ~PChemDriver () {}

  void solveTransient(FArrayBox&        Y,
		      FArrayBox&        P,
		      FArrayBox&        C,
		      FArrayBox&        F,
		      const Box&        box,
		      int               ncomps,
		      Real              dt,
		      FArrayBox&        kappa,
		      FArrayBox&        rock_phi);

  void solveSaturate(FArrayBox&        Y,
		     FArrayBox&        X,
		     FArrayBox&        P,
		     FArrayBox&        C,
		     const Box&        box,
		     int               ncomps,
		     FArrayBox&        kappa,
		     FArrayBox&        phi);

  void totalToInd(FArrayBox&        Y,
		  FArrayBox&        P,
		  FArrayBox&        C,
		  const Box&        box,
		  int               ncomps,
		  FArrayBox&        kappa,
		  FArrayBox&        phi); 

protected:

  void getSpeciesNames ();
  void getElementNames ();
  void getStoichCoeffs ();

private:

  static bool mInitialized;
  int         n_primary;
  int         n_secondary;
  int         n_core;
  int         n_type[4];
  Array<Real> mweight;
  Real        mEQTOL;

};

#endif /*_ChemDriver_H_*/
