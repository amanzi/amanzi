#ifndef _MatFillerPC_H_
#define _MatFillerPC_H_

#include <MatIDFiller.H>

class MatFillerPC
  : public MatIDFiller
{
public:

  MatFillerPC();
  MatFillerPC(const Array<Geometry>&  geomArray,
              const Array<IntVect>&   refRatio,
              const PArray<Material>& materials);
  
  void define(const Array<Geometry>&  geomArray,
              const Array<IntVect>&   refRatio,
              const PArray<Material>& materials);

  virtual void SetProperty(Real t, int level, MultiFab& mf, const std::string& pname, int dComp, int nGrow, void* ctx = 0);
  void FillCoarseCells(Real t, int level, MultiFab& mfc, const std::string& pname, int dComp, int nComp, void* ctx = 0);
  virtual void FillCellsOutsideDomain(Real t, int level, MultiFab& mf, const std::string& pname, int dComp, int nComp, int nGrow);

  virtual void Average(const FArrayBox& fineFab,
                       int              sComp,
                       FArrayBox&       crseFab,
                       const Box&       crseBox,
                       int              dComp,
                       int              nComp,
                       const IntVect&   ref) const = 0;
};

class MatFillerPCarithAvg
  : public MatFillerPC
{
public:
  MatFillerPCarithAvg();
  MatFillerPCarithAvg(const Array<Geometry>&  geomArray,
                      const Array<IntVect>&   refRatio,
                      const PArray<Material>& materials)
    : MatFillerPC(geomArray,refRatio,materials) {}

  virtual void Average(const FArrayBox& fineFab,
                       int              sComp,
                       FArrayBox&       crseFab,
                       const Box&       crseBox,
                       int              dComp,
                       int              nComp,
                       const IntVect&   ref) const;
};

#endif
