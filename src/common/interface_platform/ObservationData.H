#ifndef OBSERVATIONDATA_H
#define OBSERVATIONDATA_H

#include "mpi.h"
#include "list"
#include "vector"

namespace Amanzi
{
  class ObservationData
  {
  public:
    struct DataTriple
    {
      DataTriple() : time(-1), value(-1), is_valid(false) {}
      double time, value;
      bool is_valid;
    };
    
    ObservationData() {}
    
    std::vector<ObservationData::DataTriple> operator[](const std::string& label) const
    {
      // If label not found, returns zero-length vector
      std::map<std::string, std::vector<DataTriple> >::const_iterator it = data.find(label);
      if (it == data.end()) 
        {
          return std::vector<DataTriple>();
        }
      else
        {
          return it->second;
        }
    }
    
  private:
    
    std::map<std::string, std::vector<DataTriple> > data;
  };

} // end namespace Amanzi

#endif
