/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#ifndef _GSLibInt_H_
#define _GSLibInt_H_

#include <Array.H>
#include <MultiFab.H>

namespace GSLibInt
{
  void parRand(Array <Real> kappaval,
	       Real         dkappa,
	       Array<int>   n_cell,
	       int          twoexp,
	       MultiFab&     mfdata);

  void rdpGaussianSim(const Array<Real>& kappaval,
                      const Box&         domain,
                      const Array<Real>& problo,
                      const Array<Real>& probhi,
                      MultiFab&          mfdata,
                      int                max_grid_size_fine_gen,
                      int                ngrow_fine_gen,
                      const std::string& gsfile);

  void seqGaussianSim(const Array <Real>& kappaval,
                      const Box&          domain,
		      const Array<Real>&  problo,
		      const Array<Real>&  probhi,
		      int                 twoexp,
		      MultiFab&           mfdata,
		      const std::string&  gsfile);

  void cndGaussianSim(const Array<Real>& kappaval,
                      const Box&         domain,
                      const Array<Real>& problo,
                      const Array<Real>& probhi,
                      MultiFab&          mfdata,
                      int                max_grid_size_fine_gen,
                      int                ngrow_fine_gen,
                      std::string&       kcfile,
                      const std::string& gsfile);
 };

#endif /*_GSLibInt_H_*/
