/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include <iostream>
#include <fstream>

using std::cout;
using std::endl;

#include <MultiFabCellVector.H>


int
main (int   argc,
      char* argv[])
{
    BoxLib::Initialize(argc,argv);

    // Make a sample BoxArray
    Box domain(IntVect(D_DECL(0,0,0)),IntVect(D_DECL(3,3,3)));
    BoxArray ba(domain); ba.maxSize(2);
    if (ParallelDescriptor::IOProcessor()) {
      //cout << "BoxArray: " << ba << endl;
    }

    int nGrow = 0;
    int nComp = 3;
    MultiFab mf(ba,nComp,nGrow);
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
      for (int n=0; n<nComp; ++n) {
        mf[mfi].setVal((mfi.index()+1)*(n+1),mf[mfi].box(),n,1);
      }
    }
    MultiFabCellVector vec0(mf,nGrow,0);
    MultiFabCellVector vec1(mf,nGrow,1);
    MultiFabCellVector vec2(mf,nGrow,2);

    Real set_val = 20;
    int cnt = 0;
    for (int p=0; p<ParallelDescriptor::NProcs(); ++p) {
      ParallelDescriptor::Barrier();
      if (ParallelDescriptor::MyProc() == p) {

        for (MultiFabCellVector::iterator it=vec0.begin(), End=vec0.end(); it!=End; ++it) {
          cnt++;
          //cout << p << " " << it.index() << " " << vec0[it] << " " << vec1[it] << " " << vec2[it] << endl;
          vec0[it] = set_val;
        }
	//cout << "Visited " << cnt << " cells" << endl;
      }
      ParallelDescriptor::Barrier();
    }
    ParallelDescriptor::ReduceIntSum(cnt);

    Real tot0 = 0;
    Real tot1 = 0;
    Real tot2 = 0;
    for (int p=0; p<ParallelDescriptor::NProcs(); ++p) {
      ParallelDescriptor::Barrier();
      if (ParallelDescriptor::MyProc() == p) {
        for (MultiFabCellVector::iterator it=vec0.begin(), End=vec0.end(); it!=End; ++it) {
          tot0 += vec0[it];
          tot1 += vec1[it];
          tot2 += vec2[it];
        }
      }
      ParallelDescriptor::Barrier();
    }
    ParallelDescriptor::Barrier();
    ParallelDescriptor::ReduceRealSum(tot0);
    ParallelDescriptor::ReduceRealSum(tot1);
    ParallelDescriptor::ReduceRealSum(tot2);

#if BL_SPACEDIM==2
    if (tot0 != set_val * cnt  ||  tot1 != 80  ||  tot2 != 120) {
#else
    if (tot0 != 1280  ||  tot1 != 576  ||  tot2 != 864) {
#endif
      cout << "tot0: " << tot0 << endl;
      cout << "tot1: " << tot1 << endl;
      cout << "tot2: " << tot2 << endl;
      BoxLib::Abort();
    }

    BoxLib::Finalize();
    return 0;
}
