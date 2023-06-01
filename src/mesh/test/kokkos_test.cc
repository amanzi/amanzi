#include "UnitTest++.h"

#include "Point.hh"
#include "MeshUtils.hh"

using namespace Amanzi; 
using namespace std::chrono;

struct raggedarray{
  Kokkos::DualView<AmanziGeometry::Point*> entries;
  Kokkos::DualView<int*> rows; 
};

void test_for(int n, int np,  Amanzi::RaggedArray_DualView<AmanziGeometry::Point>& mra ){
    Kokkos::parallel_for("", n, KOKKOS_LAMBDA(int i){
      auto v = mra.getRow<MemSpace_kind::HOST>(i);
      for(int j = 0 ; j < np; ++j){
	v[j].set(i, i, i);
      }
    });
}

void test_for_2(int n, int np, raggedarray& kra){
  Kokkos::parallel_for("", n, KOKKOS_LAMBDA(int i){
      auto v = Kokkos::subview(kra.entries.h_view, Kokkos::make_pair(kra.rows.h_view[i], kra.rows.h_view[i+1]));
      for(int j = 0 ; j < np; ++j){
	v[j].set(j, j, j);
      }
    }); 
}

TEST(KOKKOS_TEST)
{
  const int n = 1000000; 
  const int np = 10;
  int val = 0;
  using timer = decltype(std::chrono::high_resolution_clock::now()); 

  {
    auto start = std::chrono::high_resolution_clock::now();
    val = 0; 
    raggedarray kra; 
    Kokkos::resize(kra.rows, n+1);
    Kokkos::resize(kra.entries, n*np);
    for(int i = 0 ; i < n+1; ++i){
      kra.rows.h_view[i] = i*np; 
    }
    auto stop = std::chrono::high_resolution_clock::now();
    std::cout<<"Kokkos Init: "<<(std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count())/1000.<<"ms"<<std::endl;
    start = std::chrono::high_resolution_clock::now(); 

    for(int i = 0; i < n; ++i){
      auto v = Kokkos::subview(kra.entries.h_view, Kokkos::make_pair(kra.rows.h_view[i], kra.rows.h_view[i+1]));
      for(int j = 0 ; j < np; ++j){
	v[j].set(val, val, val);
	++val; 
      }
    }
    stop = std::chrono::high_resolution_clock::now();
    std::cout<<"Kokkos For : "<<(std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count())/1000.<<"ms"<<std::endl; 

    start = std::chrono::high_resolution_clock::now(); 
    test_for_2(n, np, kra); 
    stop = std::chrono::high_resolution_clock::now();
    std::cout<<"Kokkos PFor: "<<(std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count())/1000.<<"ms"<<std::endl; 

  }
  std::cout<<std::endl;

  {
    auto start = std::chrono::high_resolution_clock::now();
    val = 0; 
    Amanzi::RaggedArray_DualView<AmanziGeometry::Point> mra;
    mra.resize(n, np);
    auto stop = std::chrono::high_resolution_clock::now();
    std::cout<<"Mesh   Init: "<<(std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count())/1000.<<"ms"<<std::endl;

    start = std::chrono::high_resolution_clock::now(); 
    for(int i = 0; i < n; ++i){
      auto v = mra.getRow<MemSpace_kind::HOST>(i);
      for(int j = 0 ; j < np; ++j){
	v[j].set(val, val, val);
	++val; 
      }
    }
    stop = std::chrono::high_resolution_clock::now();
    std::cout<<"Mesh   For : "<<(std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count())/1000.<<"ms"<<std::endl;

    start = std::chrono::high_resolution_clock::now();
    test_for(n, np, mra); 
    stop = std::chrono::high_resolution_clock::now();
    std::cout<<"Mesh   PFor: "<<(std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count())/1000.<<"ms"<<std::endl; 

  }

  
}
