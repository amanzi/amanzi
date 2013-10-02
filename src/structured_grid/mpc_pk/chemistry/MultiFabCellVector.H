#ifndef MultiFabCellVector_H
#define MultiFabCellVector_H

#include <MultiFab.H>
#include <MultiFabCellIterator.H>

class MultiFabCellVector
{
public:
  MultiFabCellVector(MultiFab& mf,
		     int       nGrow,
		     int       sComp = 0,
		     int       nComp = 1)
    : mfab(mf),
      iter(mf.DistributionMap(),mf.boxArray(),nGrow),
      iterBegin(iter.begin()),
      iterEnd(iter.end()),
      src_comp(sComp),
      num_comp(nComp)
    {
      BL_ASSERT(src_comp+num_comp<=mfab.nComp());
    }

  typedef MultiFabCellIterator iterator;
  const iterator& begin() const {return iterBegin;}
  const iterator& end() const {return iterEnd;}

  Real& operator[] (iterator& it) {
    return mfab[it.boxID()](it.index(),src_comp);
  }

  const Real operator[] (iterator& it) const {
    return mfab[it.boxID()](it.index(),src_comp);
  }

protected:
  MultiFab& mfab;
  MultiFabCellIterator iter, iterBegin, iterEnd;
  int src_comp, num_comp;
};

#endif // MultiFabCellVector_H
