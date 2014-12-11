#ifndef _IdxRegionData_H_
#define _IdxRegionData_H_

#include <FArrayBox.H>
#include <IArrayBox.H>
#include <RegionData.H>

/*
  IdxRegionData:

  The object is constructed with a list of regions and an evaluator object.
  Its primary function is apply, which fills a fab at each cell that is within
  the given regions, based on the integer from the IArrayFab at the same IntVect.

  The intention of this class is to support the setting of initial conditions
  specified over a set of regions that may encompass more than a single material,
  but the IC itself is a constant value that is material-dependent.

  The evaluator supports lazy evaluation, so that the function is evaluated only
  for material indices that are accessed.
 */

class IdxRegionData
{
public:
  class IdxRDEval
  {
  public:
    virtual Real operator()(int i, Real t) const = 0;
  };

  IdxRegionData(const std::string&          label,
                const Array<const Region*>& regions,
                const std::string&          typeStr,
                IdxRDEval*                  eval);

  virtual void apply(FArrayBox&       fab,
                     const IArrayBox& idx,
                     const Real*      dx,
                     int              dcomp,
                     Real             time = 0) const;

protected:
  void SetRegions(const Array<const Region*>& regions);
  std::string mLabel;
  Array<const Region*> mRegions;
  std::string mType;
  IdxRDEval* mEvaluator;
};

#endif
