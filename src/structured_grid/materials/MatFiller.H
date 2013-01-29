#ifndef _MatFiller_H_
#define _MatFiller_H_

#include <MatIDFiller.H>

class MatFiller
  : public MatIDFiller
{
public:

  MatFiller()
    : MatIDFiller() {}
  
  MatFiller(const Array<Geometry>&  geomArray,
            const Array<IntVect>&   refRatio,
            const PArray<Material>& materials)
    : MatIDFiller(geomArray,refRatio,materials) {}
  
  void define(const Array<Geometry>&  geomArray,
              const Array<IntVect>&   refRatio,
              const PArray<Material>& materials);


  virtual void CoarsenData(const FArrayBox& fineFab,
                           int              sComp,
                           FArrayBox&       crseFab,
                           const Box&       crseBox,
                           int              dComp,
                           int              nComp,
                           const IntVect&   ref) const = 0;

  virtual void RefineData(const FArrayBox& crseFab,
                          int              sComp,
                          FArrayBox&       fineFab,
                          const Box&       crseBox,
                          int              dComp,
                          int              nComp,
                          const IntVect&   ref) const = 0;

  void SetProperty(Real               t,
                   int                level,
                   MultiFab&          mf,
                   const std::string& pname,
                   int                dComp,
                   int                nGrow,
                   void*              ctx = 0);
  
  void FillCoarseCells(Real               t,
                       int                level,
                       MultiFab&          mfc,
                       const std::string& pname,
                       int                dComp,
                       int                nComp,
                       void*              ctx = 0);

  void FillCellsOutsideDomain(Real               t,
                              int                level,
                              MultiFab&          mf,
                              const std::string& pname,
                              int                dComp,
                              int                nComp,
                              int                nGrow);
};

class MatFillerPCarithAvg
  : public MatFiller
{
public:

  MatFillerPCarithAvg()
    : MatFiller() {}

  MatFillerPCarithAvg(const Array<Geometry>&  geomArray,
                      const Array<IntVect>&   refRatio,
                      const PArray<Material>& materials)
    : MatFiller(geomArray,refRatio,materials) {}

  virtual void CoarsenData(const FArrayBox& fineFab,
                           int              sComp,
                           FArrayBox&       crseFab,
                           const Box&       crseBox,
                           int              dComp,
                           int              nComp,
                           const IntVect&   ref) const;

  virtual void RefineData(const FArrayBox& crseFab,
                          int              sComp,
                          FArrayBox&       fineFab,
                          const Box&       crseBox,
                          int              dComp,
                          int              nComp,
                          const IntVect&   ref) const;
};

#endif

