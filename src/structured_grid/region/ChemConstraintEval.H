/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

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
		     ChemistryHelper_Structured* chemHelper);

  ChemConstraintEval(const std::vector<std::string>& names,
		     const std::vector<Real>&        times,
		     int                             tracerIdx,
		     RockManager*                    rockMgr,
		     ChemistryHelper_Structured*     chemHelper);

  virtual ~ChemConstraintEval() {}
  virtual ChemConstraintEval * clone () const;
  ChemConstraintEval(const ChemConstraintEval& rhs);

  virtual std::vector<Real> operator()(int i, Real t) const {
    BoxLib::Abort("ChemConstraintEval::operator() not implemented"); return std::vector<Real>(1);
  }
  std::vector<Real> operator()(const Array<Real>& primary_species_mobile,
			       const Array<Real>& auxiliary_chem_data,
			       Real               aqueous_phase_density,
			       Real               aqueous_phase_temperature,
			       Real               evaluation_time) const;

  virtual int NComp() const;

  const ChemistryHelper_Structured* ChemHelper() const {return mChemHelper;}

protected:
  std::vector<std::string>    mConstraintNames;
  std::vector<Real>           mTimes;
  int                         mTracerIdx;
  RockManager*                mRockMgr;
  ChemistryHelper_Structured* mChemHelper;
};


class ChemConstraint : public IdxRegionData
{
public:
  ChemConstraint(const std::string&          label,
                 const Array<const Region*>& regions,
                 const std::string&          typeStr,
                 const IdxRDEval&            eval)
    : IdxRegionData(label,regions,typeStr,eval) {}

  ~ChemConstraint() {};

  virtual void apply(FArrayBox&       fab, int vcomp,
		     FArrayBox&       aux, int acomp,
		     const IArrayBox& idx,
		     const Real*      dx,
		     const Box&       box,
		     const FArrayBox& density,     int dComp,
		     const FArrayBox& temperature, int TComp,
		     Real             time) const;
};


#endif
