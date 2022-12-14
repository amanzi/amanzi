/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*
  A test driver for IdxRegionData

  A specialized evaluator is defined for this test, defining three materials "brown", "red" and "green",
  and evaluates to 1, 2, or 3 in each material, respectively.  The input file must define at least two
  regions named "SoilLeft" and "SoilBottom"; the union of all defined regions must cover the physical
  space from (0,0) to (31,31).

  The first test defines the two named regions above as a combination of "brown" and "red".  A fab covering
  the domain with 1x1 cells is initialized to zero, and the test sets values in the two regions.  A check
  confirms that "green" was never evaluated.  Next, the entire domain is labeled "green", and the fab
  is set accordingly.  The results are verified, and it is further ensured that "brown" and "red" were
  not evaluated in this case.
 */

#include <iostream>
#include <IdxRegionData.H>
#include <RegionManager.H>

// closing DSO objects
#include "VerboseObject_objs.hh"

class VecEvaluator : public IdxRegionData::IdxRDEval
{
public:
  VecEvaluator(const std::vector<std::string>& colors);
  VecEvaluator(const VecEvaluator& rhs);
  virtual ~VecEvaluator();
  virtual VecEvaluator * clone () const;

  virtual std::vector<Real> operator()(int i, Real time) const;
  bool Initialized(const std::string& color) const;

protected:
  void Initialize(int i) const;
  int ColorIdx(const std::string& color) const;

  std::vector<std::string> mColors;
  mutable std::vector<Real> mVals;
  mutable std::vector<bool> mInit;
  mutable std::vector<Real> mRetVal;
};

VecEvaluator::VecEvaluator(const std::vector<std::string>& colors)
  : mColors(colors)
{
  mVals.resize(mColors.size());
  mInit.resize(mColors.size());
  for (int i=0; i<mInit.size(); ++i) {
    mInit[i] = false;
  }
}

VecEvaluator::VecEvaluator(const VecEvaluator& rhs)
{
  mColors = rhs.mColors;
  mVals = rhs.mVals;
  mInit = rhs.mInit;
}

VecEvaluator *
VecEvaluator::clone () const
{
  return new VecEvaluator(*this);
}

VecEvaluator::~VecEvaluator() {}

int
VecEvaluator::ColorIdx(const std::string& color) const {
  for (int i=0; i<mColors.size(); ++i) {
    if (mColors[i] == color) {
      return i;
    }
  }
  BoxLib::Abort("VecEvaluator::Initialize: Invalid color");
  return -1;
}

void
VecEvaluator::Initialize(int i) const
{
  // A fairly trivial initializer...
  if ( i<0 || i>=mColors.size()) {
    BoxLib::Abort("VecEvaluator::Initialize: index out of bounds");
  }
  switch (i)
  {
  case 0:
  case 1:
  case 2:
    mVals[i] = i+1;
    mInit[i] = true;
    break;
  default:
    BoxLib::Abort("VecEvaluator::Initialize: Invalid color index");
  }
}

bool
VecEvaluator::Initialized(const std::string& color) const
{
  return mInit[ColorIdx(color)];
}


std::vector<Real>
VecEvaluator::operator()(int i, Real time) const
{
  if (!mInit[i]) {
    Initialize(i);
  }
  mRetVal.resize(NComp());
  for (int i=0; i<mRetVal.size(); ++i) {
    mRetVal[i] = mVals[i];
  }
  return mRetVal;
}


int
main (int   argc,
      char* argv[])
{
  BoxLib::Initialize(argc,argv);

  RegionManager rm;
  Array<Real> dx(BL_SPACEDIM,1);

  int nVals = 3;
  std::vector<std::string> v(nVals);
  v[0] = "brown";
  v[1] = "red";
  v[2] = "green";
  VecEvaluator VE(v);

  Array<std::string> reg(2);
  reg[0] = "SoilLeft";
  reg[1] = "SoilBottom";
  Array<const Region*> regarr = rm.RegionPtrArray(reg);

  /*
    Build a IdxRegionData that will fill only regions listed above, based on
    a idx array that identifies the 'material' at each cell, which is then
    unmapped by the evaluator.  So, e.g. wherever idx=0, the 'material' is
    brown, and all brown material in the regions above will get set to a value
    of 1, based on the innerds of the evaluator.
   */
  std::string label = "label"; // Remnant bit of info handy for Amanzi implementation
  std::string typeStr = "type"; // --ditto--
  IdxRegionData idxRegionData(label,regarr,typeStr,VE);

  Box domain(IntVect(D_DECL(0,0,0)),IntVect(D_DECL(31,31,0)));
  IArrayBox idx(domain,1);
  int cnt=0;
  for (IntVect iv=domain.smallEnd(), End=domain.bigEnd(); iv<=End; domain.next(iv)) {
    idx(iv,0) = (cnt++) % 2; // Since "2" is never selected, lazy eval should prevent "green" from initializing
  }

  FArrayBox fab(domain,1); fab.setVal(0);

  Real time = 0;
  idxRegionData.apply(fab,idx,dx.dataPtr(),0,time);
  if (VE.Initialized("green")) {
    BoxLib::Abort("Green should not be initialized");
  }

  VecEvaluator VE2(v);
  Array<const Region*> regarr2 = rm.RegionPtrArray(); // All of them
  IdxRegionData idxRegionData2(label,regarr2,typeStr,VE2);
  idx.setVal(2);
  fab.setVal(0);
  idxRegionData2.apply(fab,idx,dx.dataPtr(),0,time);

  // If all went well, every point in fab should have been set

  Real correct_result = VE2(2,time)[0] * fab.box().numPts();
  Real result = fab.sum(0);

  Real error = (correct_result - result)/correct_result;

  const IdxRegionData::IdxRDEval& E = idxRegionData2.Evaluator();
  const VecEvaluator* ve = static_cast<const VecEvaluator*>(&E);
  if (ve==0) {
    BoxLib::Abort("Upcast of evaluator failed");
    return 1;
  }

  if (std::abs(error) > 1.e-14 || ve->Initialized("brown") || ve->Initialized("red") ) {
    BoxLib::Abort("Something went wrong in tIdxRegionData");
    return 1;
  }

  return 0;
}
