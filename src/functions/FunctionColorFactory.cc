#include <iostream>
#include <sstream>
#include <iomanip>

#include "AmanziComm.hh"
#include "FunctionColorFactory.hh"
#include "FunctionColor.hh"
#include "FunctionGridColor.hh"
#include "errors.hh"

namespace Amanzi {

std::unique_ptr<FunctionColor>
FunctionColorFactory::Create(std::string &filename, const Comm_type &comm) const
{
  int error;
  std::unique_ptr<FunctionColor> f;

  // Open the input file.
  std::fstream infile;
  if (comm.MyPID() == 0) {
    infile.open(filename.c_str(), std::ios::in);
    error = !infile.good();
  }
  comm.Broadcast(&error, 1, 0);
  if (error) {
    Errors::Message m;
    m << "unable to open file " << filename.c_str();
    Exceptions::amanzi_throw(m);
  }

  // Read the DATATYPE record and broadcast.
  int datatype;
  if (comm.MyPID() == 0) {
    infile >> datatype;
    error = !infile.good();
  }
  comm.Broadcast(&error, 1, 0);
  if (error) {
    Errors::Message m;
    m << "error reading DATATYPE record";
    Exceptions::amanzi_throw(m);
  }
  comm.Broadcast(&datatype, 1, 0);

  // Verify data is of int type.
  if (datatype != 0) {
    Errors::Message m;
    m << "require DATATYPE == 0";
    Exceptions::amanzi_throw(m);
  }

  // Read the GRIDTYPE record.
  std::string gridtype;
  if (comm.MyPID() == 0) {
    infile >> gridtype;
    error = !infile.good();
  }
  comm.Broadcast(&error, 1, 0);
  if (error) {
    Errors::Message m;
    m << "error reading GRIDTYPE record";
    Exceptions::amanzi_throw(m);
  }

  // Broadcast gridtype; this is painful.
  int n = gridtype.size();
  comm.Broadcast(&n, 1, 0);
  char *data = new char[n];
  if (comm.MyPID() == 0) gridtype.copy(data, n);
  comm.Broadcast(data, (int) n, 0);
  if (comm.MyPID() != 0) gridtype.assign(data, n);
  delete [] data;

  // Proceed with the rest of the file based on the value of gridtype.
  if (gridtype == "1DCoRectMesh" || gridtype == "2DCoRectMesh" || gridtype == "3DCoRectMesh") {
    std::stringstream ss(gridtype.substr(0,1));
    int dim;
    ss >> dim;
    f = create_grid_color_function(dim, infile, comm);
  } else {
    Errors::Message m;
    m << "unknown GRIDTYPE: " << gridtype.c_str();
    Exceptions::amanzi_throw(m);
  }

  infile.close();

  return f;
}

std::unique_ptr<FunctionColor>
FunctionColorFactory::create_grid_color_function(int dim, std::fstream &infile, const Comm_type& comm) const
{
  int error;

  // Read and broadcast the NXNYNZ record.
  std::vector<int> count(dim);
  if (comm.MyPID() == 0) {
    for (int k = 0; k < dim; ++k) {
      infile >> count[k];
      if ((error = !infile.good())) break;
    }
  }
  comm.Broadcast(&error, 1, 0);
  if (error) {
    Errors::Message m;
    m << "error reading NXNYNZ record";
    Exceptions::amanzi_throw(m);
  }
  comm.Broadcast(&count[0], dim, 0);

  // Check NXNYNZ values.
  for (int k = 0; k < dim; ++k) {
    if (count[k] < 1) {
      Errors::Message m;
      m << "invalid NXNYNZ value";
      Exceptions::amanzi_throw(m);
    }
  }

  // Read and broadcast the CORNERLO record.
  std::vector<double> x0(dim);
  if (comm.MyPID() == 0) {
    for (int k = 0; k < dim; ++k) {
      infile >> x0[k];
      if ((error = !infile.good())) break;
    }
  }
  comm.Broadcast(&error, 1, 0);
  if (error) {
    Errors::Message m;
    m << "error reading CORNERLO record";
    Exceptions::amanzi_throw(m);
  }
  comm.Broadcast(&x0[0], dim, 0);

  // Read the CORNERHI record; generate and broadcast the grid spacing data.
  std::vector<double> dx(dim);
  if (comm.MyPID() == 0) {
    for (int k = 0; k < dim; ++k) {
      infile >> dx[k];
      if ((error = !infile.good())) break;
      dx[k] = (dx[k] - x0[k]) / count[k];
    }
  }
  comm.Broadcast(&error, 1, 0);
  if (error) {
    Errors::Message m;
    m << "error reading CORNERHI record";
    Exceptions::amanzi_throw(m);
  }
  comm.Broadcast(&dx[0], dim, 0);

  // Check DX values.
  for (int k = 0; k < dim; ++k) {
    if (dx[k] == 0.0) {
      Errors::Message m;
      m << "invalid CORNERHI value";
      Exceptions::amanzi_throw(m);
    }
  }

  // Read the DATALOC record and broadcast.
  int dataloc;
  if (comm.MyPID() == 0) {
    infile >> dataloc;
    error = !infile.good();
  }
  comm.Broadcast(&error, 1, 0);
  if (error) {
    Errors::Message m;
    m << "error reading DATALOC record";
    Exceptions::amanzi_throw(m);
  }
  comm.Broadcast(&dataloc, 1, 0);

  // Check DATALOC value.
  if (dataloc != 0 && dataloc != 1) {
    Errors::Message m;
    m << "invalid DATALOC value";
    Exceptions::amanzi_throw(m);
  }

  // If data is point-based, tweak grid to make it cell-based.
  if (dataloc == 1) {
    for (int k = 0; k < dim; ++k) {
      x0[k] -= 0.5 * dx[k];
      count[k] += 1;
    }
  }

  // Read the DATACOL record and broadcast.
  int ncol;
  if (comm.MyPID() == 0) {
    infile >> ncol;
    error = !infile.good();
  }
  comm.Broadcast(&error, 1, 0);
  if (error) {
    Errors::Message m;
    m << "error reading DATACOL record";
    Exceptions::amanzi_throw(m);
  }
  comm.Broadcast(&ncol, 1, 0);

  // Check DATACOL value.
  if (ncol < 1) {
    Errors::Message m;
    m << "invalid DATACOL value";
    Exceptions::amanzi_throw(m);
  }
  if (ncol > 1) {
    Errors::Message m;
    m << "DATACOL > 1 not yet implemented";
    Exceptions::amanzi_throw(m);
  }

  // Read the DATAVAL record and broadcast.
  int n = 1; // product of the counts
  for (int k = 0; k < dim; ++k) n *= count[k];
  std::vector<int> array(n);
  if (comm.MyPID() == 0) {
    for (int i = 0; i < n; ++i) {
      infile >> array[i];
      if (i == n-1) { // okay to see an EOF on the last value
        if ((error = infile.fail())) break;
      } else {
        if ((error = !infile.good())) break;
      }
    }
  }
  comm.Broadcast(&error, 1, 0);
  if (error) {
    Errors::Message m;
    m << "error reading DATAVAL record";
    Exceptions::amanzi_throw(m);
  }
  comm.Broadcast(&array[0], n, 0);

  return std::make_unique<FunctionGridColor>(dim, count, x0, dx, array);
}

} // namespace Amanzi
