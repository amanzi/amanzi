/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Tabulated equations of state: rows/columns are p/T
*/

#include <fstream>
#include <string>

#include "errors.hh"
#include "LookupTable.hh"

namespace Amanzi {
namespace AmanziEOS {

/* ******************************************************************
* Populate internal structures
****************************************************************** */
LookupTable::LookupTable(Teuchos::ParameterList& plist)
{
  int nT, nP, ndata;
  std::string filename = plist.get<std::string>("table name");
  std::string field = plist.get<std::string>("field name");

  std::ifstream ifs;
  ifs.open(filename, std::ifstream::in);

  Errors::Message msg;
  msg << "\nFailed to open/read data from file: " << filename;
  if (ifs.fail()) Exceptions::amanzi_throw(msg);

  // primary variables are T and p
  ReadMetaData_(ifs, "temperature", &nT, &scaleT_, &shiftT_);
  ReadMetaData_(ifs, "pressure", &nP, &scaleP_, &shiftP_);

  // density [kg/m3] internal_energy [kJ/mol] enthalpy [kJ/mol]
  // cv [J/mol*K] cp [J/mol*K] viscosity [Pa*s]	thermal_conductivity [W/m*K]
  // phase (l or g)
  double scaleD, scaleIE, scaleE, scaleCv, scaleCp, scaleV, scaleTC;
  double shiftD, shiftIE, shiftE, shiftCv, shiftCp, shiftV, shiftTC;

  ReadMetaData_(ifs, "density", &ndata, &scaleD, &shiftD);
  ReadMetaData_(ifs, "internal_energy", &ndata, &scaleIE, &shiftIE);
  ReadMetaData_(ifs, "enthalpy", &ndata, &scaleE, &shiftE);
  ReadMetaData_(ifs, "cv", &ndata, &scaleCv, &shiftCv);
  ReadMetaData_(ifs, "cp", &ndata, &scaleCp, &shiftCp);
  ReadMetaData_(ifs, "viscosity", &ndata, &scaleV, &shiftV);
  ReadMetaData_(ifs, "thermal_conductivity", &ndata, &scaleTC, &shiftTC);

  if (field == "density") {
    scaleF_ = scaleD;
    shiftF_ = shiftD;
  } else if (field == "internal_energy") {
    scaleF_ = scaleIE;
    shiftF_ = shiftIE;
  } else if (field == "viscosity") {
    scaleF_ = scaleV;
    shiftF_ = shiftV;
  }

  // data blocks
  ReadBlock_(ifs, field, nP, nT, scaleF_, shiftF_);
}

 
/* ******************************************************************
* Main function evaluations
****************************************************************** */
double LookupTable::Function(double T, double p)
{
  int nT = axisT_.size();
  int nP = axisP_.size();

  if (T < axisT_[0] || T > axisT_[nT - 1] ||
      p < axisP_[0] || p > axisP_[nP - 1]) {
    Errors::Message msg;
    msg << "\nOut of bounds input values of T=" << T << " or P=" << p;
    amanzi_throw(msg);
  }

  int ip(0), jp(0);  // right-top corner of the table block
  for (int i = 1; i < nP; ++i) {
    if (axisP_[i] >= p) {
      ip = i;
      break;
    } 
  } 
  for (int i = 1; i < nT; ++i) {
    if (axisT_[i] >= T) {
      jp = i;
      break;
    } 
  } 

  // bilinear interpolation
  double a = (axisP_[ip] - p) / (axisP_[ip] - axisP_[ip - 1]);
  double b = (axisT_[jp] - T) / (axisT_[jp] - axisT_[jp - 1]);

  double val = (1.0 - a) * ((1.0 - b) * F_[ip][jp] + b * F_[ip][jp - 1]) 
             + a * ((1.0 - b) * F_[ip - 1][jp] + b * F_[ip - 1][jp - 1]);
  return val;
}

double LookupTable::DFunctionDT(double T, double p) {
  return 0.0;
}

double LookupTable::DFunctionDp(double T, double p) {
  return 0.0;
}


/* ******************************************************************
* Read a meta data block 
****************************************************************** */
void LookupTable::ReadMetaData_(std::ifstream& ifs, const std::string& label,
                                int* n, double* scale, double* shift)
{
  Errors::Message msg;
  msg << "\nFailed to read meta data for label: " << label;

  ifs >> *n;
  if (ifs.fail()) Exceptions::amanzi_throw(msg);

  ifs >> *scale;
  if (ifs.fail()) Exceptions::amanzi_throw(msg);
  ifs >> *shift;
  if (ifs.fail()) Exceptions::amanzi_throw(msg);

  std::string tmp;
  ifs >> tmp;
  if (ifs.fail()) Exceptions::amanzi_throw(msg);
  if (tmp != label) Exceptions::amanzi_throw(msg);
}


/* ******************************************************************
* Read a simple block of data: table + vector of doubles
****************************************************************** */
void LookupTable::ReadBlock_(std::ifstream& ifs, const std::string& field,
                             int nP, int nT, double scale, double shift)
{
  Errors::Message msg;
  msg << "\nFailed to read a data block from input stream";

  axisP_.resize(nP);
  axisT_.resize(nT);

  F_.resize(nP);
  for (int i = 0; i < nP; ++i) F_[i].resize(nT);

  double values[9]; 
  std::string tmp;

  for (int i = 0; i < nP; ++i) {
    for (int j = 0; j < nT; ++j) {
      for (int k = 0; k < 9; ++k) {
        ifs >> values[k];
        if (ifs.fail()) Exceptions::amanzi_throw(msg);
      }

      ifs >> tmp;
      if (ifs.fail()) Exceptions::amanzi_throw(msg);

      axisP_[i] = shiftP_ + scaleP_ * values[1];
      axisT_[j] = shiftT_ + scaleT_ * values[0];

      if (field == "density")
        F_[i][j] = shift + scale * values[2];
      else if (field == "internal_energy")
        F_[i][j] = shift + scale * values[3];
      else if (field == "viscosity")
        F_[i][j] = shift + scale * values[7];
    }
  }
} 


}  // namespace AmanziEOS
}  // namespace Amanzi
