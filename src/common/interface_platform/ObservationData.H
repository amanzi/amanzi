#ifndef OBSERVATIONDATA_H
#define OBSERVATIONDATA_H

#include "mpi.h"
#include "list"
#include "vector"

#include "ostream"

namespace Amanzi
{
  class ObservationData
  {
  public:
    struct DataTriple
    {
      DataTriple() : time(-1), value(-1), is_valid(false) {}
      void print(std::ostream &os) const 
      { 
	os << "is_valid = " << is_valid << ", time = " << time << ", data = " << value << std::endl; }
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
    
    std::vector<std::string> observationLabels() const
    {
      std::vector<std::string> result(data.size());
      int cnt=0;
      for (std::map<std::string, std::vector<DataTriple> >::const_iterator it = data.begin(); it!=data.end(); ++it)
        {
          result[cnt++] = it->first;
        }
      return result;
    }
    
    std::vector<ObservationData::DataTriple>& operator[](const std::string& label)
    {
      std::map<std::string, std::vector<DataTriple> >::const_iterator it = data.find(label);
      if (it == data.end()) 
        {
    	  std::vector<DataTriple> dt;
    	  data[label] = dt;
    	  return data[label];
        }
      else
        {
          return data[label];
        }
    }

    
    void print (std::ostream& os) const
    {
      for (std::map<std::string, std::vector<DataTriple> >::const_iterator it = data.begin();
	   it != data.end(); it++) 
	{
	  os << it->first << std::endl;
	  
	  for (std::vector<DataTriple>::const_iterator jt = it->second.begin(); jt != it->second.end(); jt++)
	    jt->print(os);
	}
	  
    }



  private:
    
    std::map<std::string, std::vector<DataTriple> > data;
  };

} // end namespace Amanzi

#endif
