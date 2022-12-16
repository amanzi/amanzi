/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include <ChemConstraintEval.H>

static int
which_constraint(Real time, const std::vector<Real>& times)
{
  int Ntimes = times.size();
  int i = 0;
  if (Ntimes > 1) {
    for ( ; i<Ntimes-1; ++i) {
      if (time < times[i+1]) return i;
    }
  }
  return i;
}

template< typename T>
static void
MY_RESIZE(std::vector<T>& vec, const T& def, int idx)
{
  if (idx >= vec.size()) {
    std::vector<T> TMP(idx+1);
    for (int j=0; j<vec.size(); ++j) {
      TMP[j] = vec[j];
    }
    for (int j=vec.size(); j<idx; ++j) {
      TMP[j] = def;
    }
    std::swap(vec,TMP);
  }
}

template< typename T>
static void
MY_RESIZEV(std::vector<std::vector<T> >& vec, const T& def, int idx, int size)
{
  if (idx >= vec.size()) {
    std::vector<std::vector<T> > TMP(idx+1);
    for (int j=0; j<vec.size(); ++j) {
      TMP[j] = vec[j];
    }
    for (int j=vec.size(); j<idx; ++j) {
      TMP[j].resize(size);
      for (int L=0; L<size; ++L) {
        TMP[j][L] = def;
      }
    }
    std::swap(vec,TMP);
  }
}

ChemConstraintEval::ChemConstraintEval(const std::string&          name,
				       int                         tracerIdx,
				       RockManager*                rockMgr,
				       ChemistryHelper_Structured* chemHelper)
  : mTracerIdx(tracerIdx), mRockMgr(rockMgr), mChemHelper(chemHelper)
{
  mTimes.push_back(0);
  mConstraintNames.push_back(name);
}

ChemConstraintEval::ChemConstraintEval(const std::vector<std::string>& names,
				       const std::vector<Real>&        times,
				       int                             tracerIdx,
				       RockManager*                    rockMgr,
				       ChemistryHelper_Structured*     chemHelper)
  : mConstraintNames(names), mTimes(times), mTracerIdx(tracerIdx), mRockMgr(rockMgr), mChemHelper(chemHelper)
{}

ChemConstraintEval *
ChemConstraintEval::clone () const
{
  return new ChemConstraintEval(*this);
}

ChemConstraintEval::ChemConstraintEval(const ChemConstraintEval& rhs)
{
  mConstraintNames = rhs.mConstraintNames;
  mTracerIdx = rhs.mTracerIdx;
  mRockMgr = rhs.mRockMgr;
  mChemHelper = rhs.mChemHelper;
  mTimes = rhs.mTimes;
}

std::vector<Real>
ChemConstraintEval::operator()(const Array<Real>& primary_species_mobile,
			       const Array<Real>& auxiliary_chem_data,
			       Real               aqueous_phase_density,
			       Real               aqueous_phase_temperature,
			       Real               evaluation_time) const
{
  int j = which_constraint(evaluation_time,mTimes);

  const IntVect iv(D_DECL(0,0,0));
  Box box(iv,iv);
  const std::map<std::string,int>& auxChemVariablesMap = mChemHelper->AuxChemVariablesMap();

  int Nmobile = mChemHelper->NumMobile();
  BL_ASSERT(Nmobile <= primary_species_mobile.size());
  FArrayBox prim(box,Nmobile);
  for (int i=0; i<Nmobile; ++i) {
    prim(iv,i) = primary_species_mobile[i];
  }

  int Naux = auxChemVariablesMap.size();
  BL_ASSERT(Naux <= auxiliary_chem_data.size());
  FArrayBox aux(box,Naux);
  for (int i=0; i<Naux; ++i) {
    aux(iv,i) = auxiliary_chem_data[i];
  }

  int chem_verbose = 0;
  mChemHelper->EnforceCondition(prim,0,aux,aqueous_phase_density,aqueous_phase_temperature,
				box,mConstraintNames[j],evaluation_time,chem_verbose);

  int nComp = NComp();
  Array<Real> retVals(nComp);
  BL_ASSERT(nComp >= Naux+1);
  for (int i=0; i<Naux; ++i) {
    retVals[i] = aux(iv,i);
  }
  retVals[nComp - 1] = prim.dataPtr()[mTracerIdx];
  return retVals; // Copied
}

int
ChemConstraintEval::NComp() const
{
  int Naux = mChemHelper->AuxChemVariablesMap().size();
  return Naux + 1; // Here, doing 1 mobile solute per constraint
}

void
ChemConstraint::apply(FArrayBox&       fab, int vcomp,
		      FArrayBox&       aux, int acomp,
		      const IArrayBox& idx,
		      const Real*      dx,
		      const Box&       box,
		      const FArrayBox& density,     int dComp,
		      const FArrayBox& temperature, int TComp,
		      Real             time) const
{
  const ChemConstraintEval* evaluator = dynamic_cast<const ChemConstraintEval*>(mEvaluator);
  if (evaluator == 0) {
    BoxLib::Abort("ChemConstraint not constructed with a ChemConstraintEval evaluator, unable to proceed");
  }

  if (vcomp>=fab.nComp()) BoxLib::Abort();
  FArrayBox mask(box,1); mask.setVal(-1);
  for (int j=0; j<mRegions.size(); ++j) {
    mRegions[j]->setVal(mask,1,0,dx,0);
  }

  for (IntVect iv=box.smallEnd(); iv<=box.bigEnd(); box.next(iv)) {
    if (mask(iv,0) > 0) {

      const ChemistryHelper_Structured* chemHelper = evaluator->ChemHelper();
      BL_ASSERT(chemHelper != 0);

      int Nmobile = chemHelper->NumMobile();
      Array<Real> primary_species_mobile(Nmobile,0.0);

#if 0  // FIXME: Currently, assume the constraint evaluator does not need solute concentrations, saturations and pressures on input
      for (int i=0; i<Nmobile; ++i) {
	primary_species_mobile[i] = fab(iv,vcomp+i);
      }
#endif

      int Naux = chemHelper->AuxChemVariablesMap().size();
      Array<Real> auxiliary_chem_data(Naux);
      for (int i=0; i<Naux; ++i) {
	auxiliary_chem_data[i] = aux(iv,acomp+i);
      }

      Real density_pt = density(iv,dComp);
      Real temperature_pt = temperature(iv,TComp);

      const std::vector<Real> val = (*evaluator)(primary_species_mobile, auxiliary_chem_data,
						 density_pt, temperature_pt, time);

      BL_ASSERT(val.size() == Naux + 1);
      for (int i=0; i<val.size()-1; ++i) {
        aux(iv,acomp+i) = val[i];
      }
      fab(iv,vcomp) = val[val.size()-1]; // Extract only one component

    }
  }
}
