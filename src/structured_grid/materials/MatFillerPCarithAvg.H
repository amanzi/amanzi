#ifndef _MatFillerPCarithAvg_H_
#define _MatFillerPCarithAvg_H_

#include <MatIDFiller.H>

class MatFillerPCarithAvg
  : public MatIDFiller
{
public:

  MatFillerPCarithAvg();
  MatFillerPCarithAvg(const Array<Geometry>&  geomArray,
		      const Array<IntVect>&   refRatio,
		      const PArray<Material>& materials);

  void define(const Array<Geometry>&  geomArray,
              const Array<IntVect>&   refRatio,
              const PArray<Material>& materials);

  virtual void SetProperty(Real t, int level, MultiFab& mf, const std::string& pname, int dComp, int nGrow);
  void FillCoarseCells(Real t, int level, MultiFab& mfc, const std::string& pname, int dComp);
};

#endif
