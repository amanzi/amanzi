/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include <UnitTest++.h>
#include <TestReporterStdout.h>
#include "Teuchos_GlobalMPISession.hpp"
#include "hdf5.h"

// This computes the L2 error norm for the given component in the output file by comparing 
// its values with those in the given reference file.
double ComputeL2Error(hid_t output, hid_t reference, const std::string& component)
{

  return 0.0;
}

SUITE(ChemistryBenchmarkTests) {

  class Chemistry1DBenchmarkTest {
   public:
    Chemistry1DBenchmarkTest();
    ~Chemistry1DBenchmarkTest();
 
    void RunTest(const std::string name, double* gamma);

   private:
    std::string amanzi_exe_;
  };  

  Chemistry1DBenchmarkTest::Chemistry1DBenchmarkTest() {

    // Initialize the HDF5 library.
    H5open();

    // Figure out where Amanzi lives.

  }

  Chemistry1DBenchmarkTest::~Chemistry1DBenchmarkTest() {

    // Close the HDF5 library.
    H5close();
  }

  void Chemistry1DBenchmarkTest::RunTest(const std::string name, double * gamma) {
  }  // end Chemistry1DBenchmarkTest::RunTest()

  // Amanzi U Calcite benchmark.
  TEST_FIXTURE(Chemistry1DBenchmarkTest, AmanziUCalcite) {

    // Construct the Calcite benchmark directory.
    char test_dir[1024];
    // FIXME

    // Run Amanzi.
    char command[1024];
    snprintf(command, 1024, "%s --xml_file=%s/amanzi-u-1d-calcite.xml", amanzi_exe, test_dir);
    system(command);

    // Fetch the newly-created output and the reference data.
    hid_t output = H5Fopen("calcite_data.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
    char reference_file[1024];
    snprintf(reference_file, 1024, "%s/pflotran/1d-calcite.h5");
    hid_t reference = H5Fopen(reference_file, H5F_ACC_RDONLY, H5P_DEFAULT);

    // Compute the L2 error norm for the Calcite concentration by reading data from 
    // the HDF5 files.
    double conc_L2 = ComputeL2Error(output, reference, "Total Ca++ [M]");

    // Compute the L2 error norm for the Calcite volume fraction.
    double VF_L2 = ComputeL2Error(output, reference, "Calcite_VF");

    // Close the files.
    H5Fclose(output);
    H5Fclose(reference);

  }  // end TEST_FIXTURE()

}  // end SUITE(ChemistryBenchmarkTests)

int main(int argc, char* argv[]) {
  Teuchos::GlobalMPISession mpiSession(&argc,&argv);
  return UnitTest::RunAllTests();
}
