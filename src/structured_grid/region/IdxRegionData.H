/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

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
    IdxRDEval();
    virtual ~IdxRDEval();
    IdxRDEval(const Array<Real>&        vals,
	      const Array<Real>&        times,
	      const Array<std::string>& forms);
    IdxRDEval(Real val);

    IdxRDEval(const IdxRDEval& rhs);
    virtual IdxRDEval * clone () const;

    Array<Real> Time() const;
    virtual std::vector<Real> operator()(int i, Real t) const;
    virtual int NComp() const {return 1;}

  protected:
    TabularFunction mFunc;
    mutable std::vector<Real> mRetData;
  };

  IdxRegionData(const std::string&          label,
                const Array<const Region*>& regions,
                const std::string&          typeStr,
                const IdxRDEval&            eval);

  IdxRegionData(const std::string&          label,
                const Array<const Region*>& regions,
                const std::string&          typeStr,
                Real                        val);

  IdxRegionData(const std::string&          label,
                const Array<const Region*>& regions,
                const std::string&          typeStr,
                const Array<Real>&          vals,
		const Array<Real>&          times,
		const Array<std::string>&   forms);

  virtual ~IdxRegionData();

  virtual void apply(FArrayBox&       fab,
                     const IArrayBox& idx,
                     const Real*      dx,
                     int              dcomp,
                     Real             time = 0) const;

  const Array<const Region*>& Regions() const {return mRegions;}
  std::string Type() const {return mType;}
  std::string Label() const {return mLabel;}
  const IdxRDEval& Evaluator() const {return *mEvaluator;}
  Array<Real> Time() const;

protected:
  void SetRegions(const Array<const Region*>& regions);
  std::string mLabel;
  Array<const Region*> mRegions;
  std::string mType;
  IdxRDEval* mEvaluator;
};
#endif
