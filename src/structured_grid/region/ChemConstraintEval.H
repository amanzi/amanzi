#ifndef _ChemConstraintEval
#define _ChemConstraintEval

#include <IdxRegionData.H>
#include <RockManager.H>
#include <ChemistryHelper_Structured.H>

class ChemConstraintEval
  : public IdxRegionData::IdxRDEval
{
public:
  ChemConstraintEval(const std::string&          name,
		     int                         tracerIdx,
		     RockManager*                rockMgr,
		     ChemistryHelper_Structured* chemHelper,
                     Real&                       water_density,
                     Real&                       temperature);

  ChemConstraintEval(const std::vector<std::string>& names,
		     const std::vector<Real>&        times,
		     int                             tracerIdx,
		     RockManager*                    rockMgr,
		     ChemistryHelper_Structured*     chemHelper,
                     Real&                           water_density,
                     Real&                           temperature);

  virtual ~ChemConstraintEval() {}
  virtual ChemConstraintEval * clone () const;
  ChemConstraintEval(const ChemConstraintEval& rhs);

  virtual const std::vector<Real>& operator()(int i, Real t) const;
  virtual int NComp() const;
  bool Initialized(int i, Real t) const;

protected:
  void Initialize(int i, Real t) const;

  std::vector<std::string>    mConstraintNames;
  std::vector<Real>           mTimes;
  int                         mTracerIdx;
  RockManager*                mRockMgr;
  ChemistryHelper_Structured* mChemHelper;
  Real&                       mWaterDensity;
  Real&                       mTemperature;
  mutable std::vector<std::vector<std::vector<Real> > > mVals; // mVals[constraintIdx][materialIdx][comp]
  mutable std::vector<std::vector<bool> > mInit; // mInit[constraintIdx][materialIdx]
};

class ChemConstraint
  : public IdxRegionData
{
public:
  ChemConstraint(const std::string&          label,
                 const Array<const Region*>& regions,
                 const std::string&          typeStr,
                 const IdxRDEval&            eval)
    : IdxRegionData(label,regions,typeStr,eval) {}

  ~ChemConstraint();

  virtual void apply(FArrayBox&       fab,
                     FArrayBox&       aux,
                     const IArrayBox& idx,
                     const Real*      dx,
                     int              vcomp,
                     int              acomp,
                     const Box&       box,
                     Real             time = 0) const;
};


#endif
