#ifndef _PGslib_H_
#define _PGslib_H_

#include <string>
#include <map>

#include <Array.H>
#include <FArrayBox.H>
#include <Box.H>
#include <BoxArray.H>
#include <Geometry.H>
#include <LevelBld.H>
#include <BC_TYPES.H>
#include <BC_TYPES_PM.H>
#include <AmrLevel.H>

namespace PGslib 
{
  void seqTurn3D(Array<int> n_cell, 
		 int        twoexp, 
		 MultiFab&  mfdata);

  void parRand(Array <Real> kappaval, 
	       Real         dkappa, 
	       Array<int>   n_cell, 
	       int          twoexp, 
	       MultiFab&     mfdata);

  void seqGaussianSim(Array <Real> kappaval, 
		      Real         dkappa, 
		      Array<int>   n_cell, 
		      Array<Real>  problo,
		      Array<Real>  probhi,
		      int          twoexp,
		      MultiFab&    mfdata,
		      std::string  gsfile);

  void rdpGaussianSim(Array<Real> kappaval, 
		      Real        dkappa,  
		      Array<int>  n_cell, 
		      Array<Real> problo,
		      Array<Real> probhi,
		      int         twoexp,
		      MultiFab&   mfdata,
		      std::string gsfile);

  void cndGaussianSim(Array<Real> kappaval, 
		      Real        dkappa,  
		      Array<int>   n_cell, 
		      Array<Real>  problo,
		      Array<Real>  probhi,
		      int          twoexp,
		      MultiFab&    mfdata,
		      std::string  kcfile,
		      std::string  gsfile);
 
};

#endif /*_PGslib_H_*/
