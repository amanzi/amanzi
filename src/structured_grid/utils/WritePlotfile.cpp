/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include <WritePlotfile.H>

#include <iomanip>
#include <fstream>

#include <Utility.H>
#include <VisMF.H>

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
			    const int                 *levelSteps)
{
  BL_ASSERT(data.size()>0);
  int nLevs = data[0].size();
  int nComp = data[0][0]->nComp();
  for (int i=1; i<data.size(); ++i) {
    BL_ASSERT(data[i].size()==nLevs);
    for (int lev=0; lev<nLevs; ++lev) {
      BL_ASSERT(data[i][lev]->boxArray()==data[0][lev]->boxArray());
    }
    nComp += data[i][0]->nComp();
  }
  BL_ASSERT(nComp == names.size());

  if(ParallelDescriptor::IOProcessor()) {
    if( ! BoxLib::UtilCreateDirectory(oFile,0755)) {
      BoxLib::CreateDirectoryFailed(oFile);
    }
  }
  //
  // Force other processors to wait till directory is built.
  //
  ParallelDescriptor::Barrier();

  std::string oFileHeader(oFile);
  oFileHeader += "/Header";

  VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

  std::ofstream os;
  const int finestLevel = nLevs - 1;

  //os.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

  if(ParallelDescriptor::IOProcessor()) {
    if(verbose) {
      std::cout << "Opening file = " << oFileHeader << '\n';
    }

    os.open(oFileHeader.c_str(), std::ios::out|std::ios::binary);

    if(os.fail()) {
      BoxLib::FileOpenFailed(oFileHeader);
    }

    //
    // Start writing plotfile.
    //
    os << pfversion << '\n';
    os << nComp << '\n';
    for (int n = 0; n < nComp; n++) os << names[n] << '\n';
    os << BL_SPACEDIM << '\n';
    os << std::setprecision(30) << time << '\n';
    os << finestLevel << '\n';
    for (int i = 0; i < BL_SPACEDIM; i++) os << probLo[i] << ' ';
    os << '\n';
    for (int i = 0; i < BL_SPACEDIM; i++) os << probHi[i] << ' ';
    os << '\n';
    for (int i = 0; i < finestLevel; i++) os << refRatio[i] << ' ';
    os << '\n';
    for (int i = 0; i <= finestLevel; i++) os << probDomain[i] << ' ';
    os << '\n';
    if(levelSteps != 0) {
      for (int i = 0; i <= finestLevel; i++) os << levelSteps[i] << ' ';
    } else {
      for (int i = 0; i <= finestLevel; i++) os << 0 << ' ';
    }
    os << '\n';
    for(int i = 0; i <= finestLevel; i++) {
      for(int k = 0; k < BL_SPACEDIM; k++) {
        os << dxLevel[i][k] << ' ';
      }
      os << '\n';
    }
    if(isCartGrid) {
      for(int i(0); i <= finestLevel; i++) {
        os << vfeps[i] << ' ';
      }
      os << '\n';
    }
    os << coordSys << '\n';
    os << 0 << '\n';                  // --------------- The bndry data width.
  }

  //
  // Write out level by level.
  //
  for(int iLevel(0); iLevel <= finestLevel; ++iLevel) {
    //
    // Write state data.
    //
    const BoxArray& ba = data[0][iLevel]->boxArray();
    int nGrids = ba.size();
    char buf[64];
    snprintf(buf, 64, "Level_%d", iLevel);

    if(ParallelDescriptor::IOProcessor()) {
      os << iLevel << ' ' << nGrids << ' ' << time << '\n';
      if(levelSteps != 0) {
	os << levelSteps[iLevel] << '\n';
      } else {
	os << 0 << '\n';
      }

      for(int i(0); i < nGrids; ++i) {
	const Box &b = ba[i];
	for(int n(0); n < BL_SPACEDIM; ++n) {
	  Real glo = b.smallEnd()[n] * dxLevel[iLevel][n];
	  Real ghi = (b.bigEnd()[n]+1) * dxLevel[iLevel][n];
	  os << glo << ' ' << ghi << '\n';
	}
      }
      //
      // Build the directory to hold the MultiFabs at this level.
      //
      std::string Level(oFile);
      Level += '/';
      Level += buf;

      if( ! BoxLib::UtilCreateDirectory(Level, 0755)) {
	BoxLib::CreateDirectoryFailed(Level);
      }
    }
    //
    // Force other processors to wait till directory is built.
    //
    ParallelDescriptor::Barrier();
    //
    // Now build the full relative pathname of the MultiFab.
    //
    static const std::string MultiFabBaseName("MultiFab");

    std::string PathName(oFile);
    PathName += '/';
    PathName += buf;
    PathName += '/';
    PathName += MultiFabBaseName;

    if(ParallelDescriptor::IOProcessor()) {
      //
      // The full name relative to the Header file.
      //
      std::string RelativePathName(buf);
      RelativePathName += '/';
      RelativePathName += MultiFabBaseName;
      os << RelativePathName << '\n';
    }
    MultiFab tot(data[0][iLevel]->boxArray(),nComp,0);
    int ccnt=0;
    for (int i=0; i<data.size(); ++i) {
      MultiFab::Copy(tot,*data[i][iLevel],0,ccnt,data[i][iLevel]->nComp(),0);
      ccnt += data[i][iLevel]->nComp();
    }
    VisMF::Write(tot, PathName);
  }

  os.close();
}

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
		   const Array<std::string>  &varnames,
		   const bool                 verbose,
		   const bool                 isCartGrid,
		   const Real                *vfeps,
		   const int                 *levelSteps)
{
  WritePlotfile(pfversion,Array<Array<MultiFab*> >(1,data),time,probLo,probHi,refRatio,probDomain,
		dxLevel,coordSys,oFile,varnames,verbose,isCartGrid,vfeps,levelSteps);
}

//} // namespace AmanziS
