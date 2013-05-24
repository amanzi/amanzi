#include <Utility.H>
#include <ccse-mpi.H>
#include <Region.H>
#include <MultiFab.H>

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
  CompoundRegion compound(cname,cpurpose,region1);

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

  CompoundRegion compound1(cname,cpurpose,region1);
  compound1.Union(region2);
  fab.setVal(0);
  compound1.setVal(fab,vol,0,dx.dataPtr(),0);
  Real res1 = fab.sum(0);
  pass &= std::abs(fab.sum(0) - (vol1 + vol2)) < sum_eps;
  if (ioproc && !pass) {
    std::cout << "Test 1 fail: fab.sum(): " << fab.sum(0) << std::endl;
  }

  CompoundRegion compound2(cname,cpurpose,region1);
  compound2.Union(region2);
  compound2.Intersect(region3);
  fab.setVal(0);
  compound2.setVal(fab,vol,0,dx.dataPtr(),0);
  Real res2 = fab.sum(0);
  pass &= std::abs(fab.sum(0)) < sum_eps;
  if (ioproc && !pass) {
    std::cout << "Test 2 fail: fab.sum(): " << fab.sum(0) << std::endl;
  }


  CompoundRegion compound3(cname,cpurpose,region1);
  compound3.Union(region2);
  compound3.Intersect(region4);
  fab.setVal(0);
  compound3.setVal(fab,vol,0,dx.dataPtr(),0);
  pass &= std::abs(fab.sum(0) - vol4) < sum_eps;
  if (ioproc && !pass) {
    std::cout << "Test 3 fail: fab.sum(): " << fab.sum(0) << std::endl;
  }

  Box box2(IntVect(D_DECL(0,0,0)),
           IntVect(D_DECL(59,59,59)));
  FArrayBox fab2(box,1);
  Array<Real> dx2(BL_SPACEDIM,0.5);
  Real volB = 1;
  for (int i=0; i<BL_SPACEDIM; ++i) {
    volB *= dx2[i];
  }

  fab2.setVal(0);
  compound3.setVal(fab2,volB,0,dx2.dataPtr(),0);
  pass &= std::abs(fab2.sum(0) - vol4) < sum_eps;
  if (ioproc && !pass) {
    std::cout << "Test 4 fail: fab2.sum(): " << fab2.sum(0) << std::endl;
  }

  BoxArray ba(box2); ba.maxSize(30);
  MultiFab mf(ba,2,0);
  mf.setVal(0);
  Real sum = 0;
  for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
    FArrayBox& fab_mf = mf[mfi];
    compound3.setVal(fab_mf,volB,1,dx2.dataPtr(),0);
    sum += fab_mf.sum(1);
  }
  ParallelDescriptor::ReduceRealSum(sum);
  pass &= std::abs(sum - vol4) < sum_eps;
  if (ioproc && !pass) {
    std::cout << "Test 5 fail: sum: " << sum << std::endl;
  }
  
  CompoundRegion compound4(cname,cpurpose,regionA);
  compound4.Subtract(region1);
  fab.setVal(0);
  compound4.setVal(fab,vol,0,dx.dataPtr(),0);
  pass &= std::abs(fab.sum(0) - (volA - vol1)) < sum_eps;
  if (ioproc && !pass) {
    std::cout << "Test 6 fail: fab.sum(): " << fab.sum(0) << std::endl;
  }

  CompoundRegion compound5(cname,cpurpose,region1);
  compound5.Subtract(region2);
  fab.setVal(0);
  compound5.setVal(fab,vol,0,dx.dataPtr(),0);
  pass &= std::abs(fab.sum(0) - vol1) < sum_eps;
  if (ioproc && !pass) {
    std::cout << "Test 7 fail: fab.sum(): " << fab.sum(0) << std::endl;
  }

  BoxLib::Finalize();  
  return 0;
}
