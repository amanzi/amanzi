/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include <Utility.H>
#include <ccse-mpi.H>
#include <Region.H>
#include <MultiFab.H>

// closing DSO objects
#include "VerboseObject_objs.hh"

#ifdef _OPENMP
#include "omp.h"
#endif

static double geom_eps = 1.e-10;
static double sum_eps = 1.e-20;

int main(int argc, char* argv[])
{
  BoxLib::Initialize(argc,argv);

  Array<Real> problo(BL_SPACEDIM,0), probhi(BL_SPACEDIM,30);

  Region::geometry_eps = geom_eps;
  Region::domlo = problo;
  Region::domhi = probhi;
  bool pass = true;
  bool ioproc = ParallelDescriptor::IOProcessor();

  std::string name1 = "Fred";
  std::string purpose1 = "Test";
  Array<Real> blo1(BL_SPACEDIM); D_EXPR(blo1[0]=0,  blo1[1]=0,  blo1[2]=0);
  Array<Real> bhi1(BL_SPACEDIM); D_EXPR(bhi1[0]=10, bhi1[1]=10, bhi1[2]=10);
  BoxRegion region1(name1,purpose1,blo1,bhi1);
  Real vol1=1;
  for (int i=0; i<BL_SPACEDIM; ++i) vol1 *= bhi1[i] - blo1[i];

  std::string name2 = "Burt";
  std::string purpose2 = "Test";
  Array<Real> blo2(BL_SPACEDIM); D_EXPR(blo2[0]=10, blo2[1]=0,  blo2[2]=0);
  Array<Real> bhi2(BL_SPACEDIM); D_EXPR(bhi2[0]=20, bhi2[1]=10, bhi2[2]=10);
  BoxRegion region2(name2,purpose2,blo2,bhi2);
  Real vol2=1;
  for (int i=0; i<BL_SPACEDIM; ++i) vol2 *= bhi2[i] - blo2[i];

  std::string name3 = "Jim";
  std::string purpose3 = "Test";
  Array<Real> blo3(BL_SPACEDIM); D_EXPR(blo3[0]=20, blo3[1]=0,  blo3[2]=0);
  Array<Real> bhi3(BL_SPACEDIM); D_EXPR(bhi3[0]=30, bhi3[1]=10, bhi3[2]=10);
  BoxRegion region3(name3,purpose3,blo3,bhi3);
  Real vol3=1;
  for (int i=0; i<BL_SPACEDIM; ++i) vol3 *= bhi3[i] - blo3[i];

  std::string name4 = "Bob";
  std::string purpose4 = "Test";
  Array<Real> blo4(BL_SPACEDIM); D_EXPR(blo4[0]=7,  blo4[1]=0,  blo4[2]=0);
  Array<Real> bhi4(BL_SPACEDIM); D_EXPR(bhi4[0]=12, bhi4[1]=10, bhi4[2]=10);
  BoxRegion region4(name4,purpose4,blo4,bhi4);
  Real vol4=1;
  for (int i=0; i<BL_SPACEDIM; ++i) vol4 *= bhi4[i] - blo4[i];

  AllRegion regionA;
  Real volA=1;
  for (int i=0; i<BL_SPACEDIM; ++i) volA *= probhi[i] - problo[i];

  std::string cname = "Comp";
  std::string cpurpose = "Test1";
  Array<const Region*> set1(1);
  set1[0] = &region1;
  UnionRegion compound(cname,cpurpose,set1);

  Box box(IntVect(D_DECL(0,0,0)),
          IntVect(D_DECL(29,29,29)));
  FArrayBox fab(box,1);
  Array<Real> dx(BL_SPACEDIM,1);
  Real vol = 1;
  for (int i=0; i<BL_SPACEDIM; ++i) {
    vol *= dx[i];
  }

  fab.setVal(0);
  compound.setVal(fab,vol,0,dx.dataPtr(),0);
  Real res0 = fab.sum(0);
  pass &= std::abs(res0 - vol1) < sum_eps;
  if (ioproc && !pass) {
    std::cout << "Test 0 fail: fab.sum(): " << fab.sum(0) << " " << res0 << std::endl;
  }

  Array<const Region*> set2(2);
  set2[0] = &region1;
  set2[1] = &region2;
  UnionRegion compound1(cname,cpurpose,set2);
  fab.setVal(0);
  compound1.setVal(fab,vol,0,dx.dataPtr(),0);
  Real res1 = fab.sum(0);
  pass &= std::abs(fab.sum(0) - (vol1 + vol2)) < sum_eps;
  if (ioproc && !pass) {
    std::cout << "Test 1 fail: fab.sum(): " << fab.sum(0) << std::endl;
  }

  Array<const Region*> set3(2);
  set3[0] = &regionA;
  set3[1] = &region1;
  SubtractionRegion compound4(cname,cpurpose,set3);

  fab.setVal(0);
  compound4.setVal(fab,vol,0,dx.dataPtr(),0);
  pass &= std::abs(fab.sum(0) - (volA - vol1)) < sum_eps;
  if (ioproc && !pass) {
    std::cout << "Test 6 fail: fab.sum(): " << fab.sum(0) << std::endl;
  }

  Array<const Region*> set4(2);
  set4[0] = &region1;
  set4[1] = &region2;
  SubtractionRegion compound5(cname,cpurpose,set4);

  fab.setVal(0);
  compound5.setVal(fab,vol,0,dx.dataPtr(),0);
  pass &= std::abs(fab.sum(0) - vol1) < sum_eps;
  if (ioproc && !pass) {
    std::cout << "Test 7 fail: fab.sum(): " << fab.sum(0) << std::endl;
  }

  Real test_sum = 0;
  for (int i=0; i<10; ++i) {
    Array<Real> coord(BL_SPACEDIM,2);
    coord[0]=10*i - .001;
    PointRegion ptReg("myPt","source",coord);
    fab.setVal(0);
    ptReg.setVal(fab,1,0,dx.dataPtr(),0);
    test_sum += fab.sum(0,1);
  }
  pass &= test_sum == 3;
  if (ioproc && !pass) {
    std::cout << "Test 8 fail" << std::endl;
  }

  BoxLib::Finalize();
  return 0;
}
