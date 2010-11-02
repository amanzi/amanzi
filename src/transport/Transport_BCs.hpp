#ifndef __TransportBCs_hpp__
#define __TransportBCs_hpp__


using namespace std;


class Transport_BCs {

public:
  Transport_BCs() {};
  Transport_BCs( int ssid, int ntcc ) { side_set_id = ssid; 
                                        bc_values.resize( ntcc ); }
  ~Transport_BCs() {};

  inline
  int get_side_set_id() { return side_set_id; }
  
  vector<int>& get_BCs() { return bc_values; }

private:
  int  side_set_id;
  vector<int>  bc_values;
};

#endif
