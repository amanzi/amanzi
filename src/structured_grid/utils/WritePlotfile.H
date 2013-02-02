#ifndef _WritePlotfile_H_
#define _WritePlotfile_H_

#include <PArray.H>
#include <MultiFab.H>

//namespace AmanziS {

  void WritePlotfile(const std::string         &pfversion,
		     const Array<Array<MultiFab*> > &data,
		     const Real                 time,
		     const Real                *probLo,
		     const Real                *probHi,
		     const Array<int>          &refRatio,
		     const Array<Box>          &probDomain,
		     const Array<Array<Real> > &dxLevel,
		     const int                  coordSys,
		     const std::string         &oFile,
		     const Array<std::string>  &names,
		     const bool                 verbose,
		     const bool                 isCartGrid,
		     const Real                *vfeps,
		     const int                 *levelSteps);
  
  void WritePlotfile(const std::string         &pfversion,
		     const Array<MultiFab*>    &data,
		     const Real                 time,
		     const Real                *probLo,
		     const Real                *probHi,
		     const Array<int>          &refRatio,
		     const Array<Box>          &probDomain,
		     const Array<Array<Real> > &dxLevel,
		     const int                  coordSys,
		     const std::string         &oFile,
		     const Array<std::string>  &names,
		     const bool                 verbose,
		     const bool                 isCartGrid,
		     const Real                *vfeps,
		     const int                 *levelSteps);
//}

#endif
