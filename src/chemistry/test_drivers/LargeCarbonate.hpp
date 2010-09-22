#ifndef __LargeCarbonate_hpp__
#define __LargeCarbonate_hpp__

#include "Geochemistry.hpp"

// Driver class for evalating geochemical related processes at a
// single computational node

class LargeCarbonate : public Geochemistry {

public:
  LargeCarbonate();
  ~LargeCarbonate();

  void setup(std::vector<double> *total);


private:

};

#endif
