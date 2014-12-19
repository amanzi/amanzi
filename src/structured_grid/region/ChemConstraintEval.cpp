
#include <ChemConstraintEval.H>

ChemConstraintEval::ChemConstraintEval(const std::string&          name,
				       int                         tracerIdx,
				       RockManager*                rockMgr,
				       ChemistryHelper_Structured* chemHelper,
                                       Real&                       water_density,
                                       Real&                       temperature)
  : mTracerIdx(tracerIdx), mRockMgr(rockMgr), mChemHelper(chemHelper),
    mWaterDensity(water_density), mTemperature(temperature)
{
  mTimes.push_back(0);
  mConstraintNames.push_back(name);
  int Nnames = mConstraintNames.size();
  mInit.resize(Nnames);
  mVals.resize(Nnames);
}

ChemConstraintEval::ChemConstraintEval(const std::vector<std::string>& names,
				       const std::vector<Real>&        times,
				       int                             tracerIdx,
				       RockManager*                    rockMgr,
				       ChemistryHelper_Structured*     chemHelper,
                                       Real&                           water_density,
                                       Real&                           temperature)
  : mConstraintNames(names), mTimes(times), mTracerIdx(tracerIdx), mRockMgr(rockMgr), mChemHelper(chemHelper),
    mWaterDensity(water_density), mTemperature(temperature)
{
  int Nnames = mConstraintNames.size();
  mInit.resize(Nnames);
  mVals.resize(Nnames);
}

ChemConstraintEval *
ChemConstraintEval::clone () const
{
  return new ChemConstraintEval(*this);
}

ChemConstraintEval::ChemConstraintEval(const ChemConstraintEval& rhs)
  : mWaterDensity(rhs.mWaterDensity), mTemperature(rhs.mTemperature)
{
  mConstraintNames = rhs.mConstraintNames;
  mTracerIdx = rhs.mTracerIdx;
  mRockMgr = rhs.mRockMgr;
  mChemHelper = rhs.mChemHelper;
  mVals = rhs.mVals;
  mInit = rhs.mInit;
  mTimes = rhs.mTimes;
}

static int
which_constraint(Real time, const std::vector<Real>& times)
{
  int Ntimes = times.size();
  BL_ASSERT(Ntimes>1);
  int i = 0;
  if (Ntimes > 1) {
    for ( ; i<Ntimes-1; ++i) {
      if (time < times[i+1]) return i;
    }
  }
  return i;
}

const std::vector<Real>&
ChemConstraintEval::operator()(int i, Real t) const
{
  Initialize(i,t);
  int j = which_constraint(t,mTimes);
  BL_ASSERT(j<mVals.size());
  BL_ASSERT(i<mVals[j].size());
  return mVals[j][i];
}

bool
ChemConstraintEval::Initialized(int i, Real t) const
{
  int j = which_constraint(t,mTimes);
  BL_ASSERT(j<mInit.size());
  BL_ASSERT(i<mInit[j].size());
  return mInit[j][i];
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

int
ChemConstraintEval::NComp() const
{
  int Naux = mChemHelper->AuxChemVariablesMap().size();
  //int Nmobile = mChemHelper->NumMobile();
  return Naux + 1; // Here, doing 1 mobile solute per constraint
}

void
ChemConstraintEval::Initialize(int i, Real t) const
{
  if (i<0) {
    BoxLib::Abort("VecEvaluator::Initialize: Invalid material index");
  }
  int j = which_constraint(t,mTimes);
  BL_ASSERT(j<mVals.size());
  BL_ASSERT(j<mInit.size());
  MY_RESIZE(mInit[j],false,i);
  MY_RESIZEV(mVals[j],-1.0,i,NComp());

  if (!mInit[j][i]) {
    Real cur_time = 0;
    Box boxTMP(IntVect(D_DECL(0,0,0)),IntVect(D_DECL(0,0,0)));
    const std::map<std::string,int>& auxChemVariablesMap = mChemHelper->AuxChemVariablesMap();
    int Nmobile = mChemHelper->NumMobile();
    FArrayBox primTMP(boxTMP,Nmobile); primTMP.setVal(0);
    int Naux = auxChemVariablesMap.size();
    FArrayBox auxTMP(boxTMP,Naux);

    const std::map<std::string,int>& aux_chem_variables_map = mChemHelper->AuxChemVariablesMap();
    const std::string& material_name = mRockMgr->GetMaterial(i).Name();

    mRockMgr->RockChemistryProperties(auxTMP,material_name,aux_chem_variables_map);
    mChemHelper->EnforceCondition(primTMP,0,auxTMP,mWaterDensity,mTemperature,boxTMP,mConstraintNames[j],cur_time);

    mVals[j][i].resize(NComp());
    for (int L=0; L<Naux; ++L) {
      mVals[j][i][L] = auxTMP.dataPtr(L)[0];
    }
    mVals[j][i][Naux] = primTMP.dataPtr()[mTracerIdx];
    mInit[j][i] = true;
  }
}

void
ChemConstraint::apply(FArrayBox&       fab,
                      FArrayBox&       aux,
                      const IArrayBox& idx,
                      const Real*      dx,
                      int              vcomp,
                      int              acomp,
                      const Box&       box,
                      Real             time) const
{
  if (vcomp>=fab.nComp()) BoxLib::Abort();
  FArrayBox mask(box,1); mask.setVal(-1);
  for (int j=0; j<mRegions.size(); ++j) { 
    mRegions[j]->setVal(mask,1,0,dx,0);
  }

  for (IntVect iv=box.smallEnd(); iv<=box.bigEnd(); box.next(iv)) {
    if (mask(iv,0) > 0) {
      const std::vector<Real>& val = (*mEvaluator)(idx(iv,0),time);
      if (acomp+val.size()-1>aux.nComp()) BoxLib::Abort();
      for (int i=0; i<val.size()-1; ++i) {
        aux(iv,acomp+i) = val[i];
      }
      fab(iv,vcomp) = val[val.size()-1];
    }
  }
}
