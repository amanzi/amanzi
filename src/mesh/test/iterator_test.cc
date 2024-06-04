#include "MeshView.hh"

void
mesh_view_test()
{
  constexpr int my_size = 100;
  Kokkos::MeshView<int*, Kokkos::DefaultExecutionSpace> dv("Device View", my_size);
  Kokkos::MeshView<int*, Kokkos::HostSpace> hv("Host View", my_size);

  int val = 0;
  for (auto& v : hv) v = val++;

  Kokkos::deep_copy(dv, hv);

  for (auto v : hv) { printf("%d ", v); }

  Kokkos::parallel_for(
    my_size, KOKKOS_LAMBDA(const int i) {
      for (auto v : dv) { printf("%d ", v); }
    });
}

int
main(int argc, char* argv[])
{
  Kokkos::initialize(argc, argv);
  mesh_view_test();
  Kokkos::finalize();
}
