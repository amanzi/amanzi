#ifndef __TransportBCs_hpp__
#define __TransportBCs_hpp__

#include <vector>

class Transport_BCs {
 public:
  Transport_BCs() {};
  Transport_BCs(int ssid, int ntcc) { 
    side_set_id = ssid;  
    values.resize(ntcc); 
  }
  ~Transport_BCs() {};

 public:
  int side_set_id;
  int type;
  std::vector<double> values;
  std::vector<unsigned int> faces;

  std::vector<double>   influx; // accumulated influx and outflux
  std::vector<double>  outflux;
};

#endif
