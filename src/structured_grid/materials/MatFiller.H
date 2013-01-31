#ifndef _MatFiller_H_
#define _MatFiller_H_

#include <MatIDFiller.H>

class MatFiller
  : public MatIDFiller
{
public:

  MatFiller(const Array<Geometry>&  geomArray,
            const Array<IntVect>&   refRatio,
            const PArray<Material>& materials);
  
  void SetProperty(Real               t,
                   int                level,
                   MultiFab&          mf,
                   const std::string& pname,
                   int                dComp,
                   int                nGrow,
                   void*              ctx = 0);

  Property::CoarsenRule coarsenRule (const std::string& pname) const;
  Property::RefineRule refineRule (const std::string& pname) const;
  int nComp(const std::string& property_name) const;

protected:

  virtual void CoarsenData(const FArrayBox&   fineFab,
                           int                sComp,
                           FArrayBox&         crseFab,
                           const Box&         crseBox,
                           int                dComp,
                           int                nComp,
                           const IntVect&     ref,
                           const std::string& pname) const;

  virtual void RefineData(const FArrayBox&   crseFab,
                          int                sComp,
                          FArrayBox&         fineFab,
                          const Box&         crseBox,
                          int                dComp,
                          int                nComp,
                          const IntVect&     ref,
                          const std::string& pname) const;

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

  void VerifyProperties();
  std::map<std::string,int> property_nComps;
  std::map<std::string,Property::CoarsenRule> property_cRules;
  std::map<std::string,Property::RefineRule> property_rRules;
};

#endif

