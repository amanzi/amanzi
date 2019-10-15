/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//! <MISSING_ONELINE_DOCSTRING>

#ifndef AMANZI_OBSERVATION_DATA_HH
#define AMANZI_OBSERVATION_DATA_HH

#include <list>
#include <map>
#include <ostream>
#include <string>
#include <vector>

namespace Amanzi {

class ObservationData {
 public:
  struct DataQuadruple {
    DataQuadruple() : time(-1), value(-1), is_valid(false), unit(""){};
    void print(std::ostream& os) const
    {
      os << "is_valid=" << is_valid << ", time=" << time << ", data=" << value
         << ", unit=" << unit << std::endl;
    }
    double time, value;
    bool is_valid;
    std::string unit;
  };

  ObservationData(){};

  std::vector<ObservationData::DataQuadruple>
  operator[](const std::string& label) const
  {
    // If label not found, returns zero-length vector
    std::map<std::string, std::vector<DataQuadruple>>::const_iterator it =
      data.find(label);
    if (it == data.end()) {
      return std::vector<DataQuadruple>();
    } else {
      return it->second;
    }
  }

  std::vector<std::string> observationLabels() const
  {
    std::vector<std::string> result(data.size());
    int cnt = 0;
    for (std::map<std::string, std::vector<DataQuadruple>>::const_iterator it =
           data.begin();
         it != data.end();
         ++it) {
      result[cnt++] = it->first;
    }
    return result;
  }

  std::vector<ObservationData::DataQuadruple>&
  operator[](const std::string& label)
  {
    std::map<std::string, std::vector<DataQuadruple>>::const_iterator it =
      data.find(label);
    if (it == data.end()) {
      std::vector<DataQuadruple> dt;
      data[label] = dt;
      return data[label];
    } else {
      return data[label];
    }
  }

  std::vector<std::string> observationLabels()
  {
    std::vector<std::string> result(data.size());
    int cnt = 0;
    for (std::map<std::string, std::vector<DataQuadruple>>::const_iterator it =
           data.begin();
         it != data.end();
         ++it) {
      result[cnt++] = it->first;
    }
    return result;
  }

  void print(std::ostream& os) const
  {
    for (std::map<std::string, std::vector<DataQuadruple>>::const_iterator it =
           data.begin();
         it != data.end();
         it++) {
      os << it->first << std::endl;

      for (std::vector<DataQuadruple>::const_iterator jt = it->second.begin();
           jt != it->second.end();
           jt++)
        jt->print(os);
    }
  }

 private:
  std::map<std::string, std::vector<DataQuadruple>> data;
};

} // namespace Amanzi

#endif
