/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include <vector>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <float.h>
#include <TestReporterStdout.h>
#include "UnitTest++.h"
#include "hdf5.h"

#include "state_evaluators_registration.hh"
#include "VerboseObject_objs.hh"


// This computes the L2 error norm for the given component in the output file by comparing
// its values with those in the given reference file.
double
ComputeL2Error(hid_t output,
               const std::string& output_component_name,
               hid_t reference,
               const std::string& reference_component_name,
               int step,
               double time)
{
  // The name of the dataset for the Amanzi file is the time step number, and datasets
  // are stored in groups that are named after the component.
  char output_dataset_name[128];
  snprintf(output_dataset_name, 128, "%d", step);
  hid_t output_group = H5Gopen2(output, output_component_name.c_str(), H5P_DEFAULT);
  if (output_group < 0) return FLT_MAX;
  hid_t output_data = H5Dopen2(output_group, output_dataset_name, H5P_DEFAULT);
  if (output_data < 0) return FLT_MAX;
  hsize_t o_size = H5Dget_storage_size(output_data);
  std::vector<double> o_data((size_t)(o_size / sizeof(double)));
  H5Dread(output_data, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &o_data[0]);
  H5Dclose(output_data);
  H5Gclose(output_group);

  // Currently we use PFlotran as a benchmark. In PFlotran files, the Group "/" contains
  // a Group for each time, in the format "Time:  x.yzxyzE+xy y". Within a Group for the
  // Time, we search for the given component.
  hid_t slash = H5Gopen2(reference, "/", H5P_DEFAULT);
  char reference_group_name[1025];
  snprintf(reference_group_name, 1024, "Time:  %1.5E y", time);
  hid_t reference_group = H5Gopen2(slash, reference_group_name, H5P_DEFAULT);
  if (reference_group < 0) return FLT_MAX;
  hid_t reference_data = H5Dopen2(reference_group, reference_component_name.c_str(), H5P_DEFAULT);
  if (reference_data < 0) return FLT_MAX;
  hsize_t r_size = H5Dget_storage_size(reference_data);
  std::vector<double> r_data((size_t)(r_size / sizeof(double)));
  H5Dread(reference_data, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &r_data[0]);
  H5Dclose(reference_data);
  H5Gclose(reference_group);
  H5Gclose(slash);

  // Make sure the datasets are the same size (for now).
  if (o_data.size() != r_data.size()) return FLT_MAX;

  // Now that everything is in place, compute the L2 error norm.
  double L2 = 0.0;
  for (size_t i = 0; i < o_data.size(); ++i) {
    double err = o_data[i] - r_data[i];
    L2 += std::sqrt(err * err);
  }

  return L2;
}


SUITE(ChemistryBenchmarkTests)
{
  class Chemistry1DBenchmarkTest {
   public:
    Chemistry1DBenchmarkTest();
    ~Chemistry1DBenchmarkTest();

    void RunTest(const std::string& filename);

   protected:
    std::string amanzi_exe_;
    std::string benchmark_dir_;
  };

  Chemistry1DBenchmarkTest::Chemistry1DBenchmarkTest()
  {
    H5open();
    amanzi_exe_ = std::string(CMAKE_BINARY_DIR) +
                  std::string("/src/common/standalone_simulation_coordinator/amanzi");
    benchmark_dir_ =
      std::string(CMAKE_SOURCE_DIR) + std::string("/test_suites/benchmarking/chemistry");
  }

  Chemistry1DBenchmarkTest::~Chemistry1DBenchmarkTest()
  {
    H5close();
  }

  // Amanzi U Calcite benchmarks
  TEST_FIXTURE(Chemistry1DBenchmarkTest, AmanziUCalciteA)
  {
    // Construct the Calcite benchmark directory.
    char test_dir[1025];
    snprintf(test_dir, 1024, "%s/calcite_1d", benchmark_dir_.c_str());

    // Copy the contents of the directory to cwd.
    char command[1025];
    snprintf(command, 1024, "cp -R %s/* .", test_dir);
    int status = std::system(command);

    // Run Amanzi.
    snprintf(
      command, 1024, "%s --xml_file=test/chemistry_benchmarks_1d_a.xml", amanzi_exe_.c_str());
    std::cout << command << "\n";
    status = std::system(command);

    // Fetch the newly-created output and the reference data.
    hid_t output = H5Fopen("calcite_data.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
    char reference_file[1025];
    snprintf(reference_file, 1024, "%s/pflotran/1d-calcite.h5", test_dir);
    hid_t reference = H5Fopen(reference_file, H5F_ACC_RDONLY, H5P_DEFAULT);

    // Compute the L2 error norm for the Calcite concentration by reading data from
    // the HDF5 files.
    double conc_L2 = ComputeL2Error(
      output, "total_component_concentration.Ca++", reference, "Total_Ca++ [M]", 27, 9.0);
    std::cout << "Ca++ concentration L2 norm: " << conc_L2 << std::endl;
    CHECK(conc_L2 < 0.025);

    // Compute the L2 error norm for the Calcite volume fraction.
    double VF_L2 =
      ComputeL2Error(output, "mineral_volume_fractions.Calcite", reference, "Calcite_VF", 27, 9.0);
    std::cout << "Ca++ volume fraction L2 norm: " << VF_L2 << std::endl;
    CHECK(VF_L2 < 0.0002);

    // Close the files.
    H5Fclose(output);
    H5Fclose(reference);
  }

  TEST_FIXTURE(Chemistry1DBenchmarkTest, AmanziUCalciteB)
  {
    char test_dir[1025];
    snprintf(test_dir, 1024, "%s/calcite_1d", benchmark_dir_.c_str());

    char command[1025];
    snprintf(
      command, 1024, "%s --xml_file=test/chemistry_benchmarks_1d_b.xml", amanzi_exe_.c_str());
    std::cout << command << "\n";
    int status = std::system(command);

    // Fetch the newly-created output and the reference data.
    hid_t output = H5Fopen("calcite_data.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
    char reference_file[1025];
    snprintf(reference_file, 1024, "%s/pflotran/1d-calcite.h5", test_dir);
    hid_t reference = H5Fopen(reference_file, H5F_ACC_RDONLY, H5P_DEFAULT);

    // Compute the L2 error norm for the Calcite concentration by reading data from
    // the HDF5 files.
    double conc_L2 = ComputeL2Error(
      output, "total_component_concentration.Ca++", reference, "Total_Ca++ [M]", 27, 9.0);
    std::cout << "Ca++ concentration L2 norm: " << conc_L2 << std::endl;
    CHECK(conc_L2 < 0.025);

    // Compute the L2 error norm for the Calcite volume fraction.
    double VF_L2 =
      ComputeL2Error(output, "mineral_volume_fractions.Calcite", reference, "Calcite_VF", 27, 9.0);
    std::cout << "Ca++ volume fraction L2 norm: " << VF_L2 << std::endl;
    CHECK(VF_L2 < 0.0002);

    // Close the files.
    H5Fclose(output);
    H5Fclose(reference);
  }
}

int
main(int argc, char* argv[])
{
  return UnitTest::RunAllTests();
}
