#ifndef __TransportBCs_hpp__
#define __TransportBCs_hpp__


using namespace std;


class Transport_BCs {

public:
  Transport_BCs() {};
  Transport_BCs( int ssid, int ntcc ) { side_set_id = ssid; 
                                        values.resize( ntcc ); }
  ~Transport_BCs() {};


public:
  int  side_set_id;
  vector<double>  values;
  vector<unsigned int>  faces;
};

#endif
