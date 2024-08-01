/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#ifndef AMANZI_MPC_UTILS_HH
#define AMANZI_MPC_UTILS_HH

#include <iostream>
#include "stdlib.h"
#include "math.h"

/* ******************************************************************
* TBW
****************************************************************** */
std::vector<std::pair<double, double>>
getObservations(const std::string& filename, int skip)
{
  std::vector<std::pair<double, double>> data;

  std::ifstream ifs(filename.c_str(), std::ios::in);

  // skip two lines
  char line[80];
  ifs.getline(line, 80);
  ifs.getline(line, 80);

  do {
    std::string word;
    ifs >> word;
    if (ifs.eof()) break;

    for (int i = 0; i < skip; ++i) ifs >> word;
    double val1 = std::atof(word.c_str());
    ifs >> word;
    double val2 = std::atof(word.c_str());
    data.push_back(std::make_pair(val1, val2));
    if (ifs.eof()) break;
  } while (true);

  return data;
}


/* ******************************************************************
* Write data to an HDF5 file.
****************************************************************** */
void CreateApertureFile(int ncells, double time)
{
  hid_t hout = H5Fcreate("test/aperture_dynamic_test.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  hid_t gout = H5Gcreate(hout, "fracture-aperture.cell.0", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // times for aperture are hardcoded
  std::vector<double> aperture0(ncells, 4.6416e-4);
  std::vector<double> aperture1(ncells, 5.6416e-4);

  std::vector<double> times = { 0.0, time };

  hsize_t dims[2] = { (hsize_t)ncells, 1 };
  hid_t dataspace = H5Screate_simple(2, dims, NULL);
  hid_t dataset0 = H5Dcreate(gout, "0", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset0, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, aperture0.data());
  H5Dclose(dataset0);

  hid_t dataset1 = H5Dcreate(gout, "1", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, aperture1.data());
  H5Dclose(dataset1);

  hsize_t time_dims[1] = {2};
  hid_t time_dataspace = H5Screate_simple(1, time_dims, NULL);
  hid_t time_dataset = H5Dcreate(hout, "time", H5T_NATIVE_DOUBLE, time_dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(time_dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, times.data());
  H5Dclose(time_dataset);

  H5Sclose(dataspace);
  H5Sclose(time_dataspace);
  H5Gclose(gout);
  H5Fclose(hout);
}
#endif
