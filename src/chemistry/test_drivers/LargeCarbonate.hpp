#ifndef __LargeCarbonate_hpp__
#define __LargeCarbonate_hpp__

#include "Beaker.hpp"

// Driver class for evalating geochemical related processes at a
// single computational node

class LargeCarbonate : public Beaker {

public:
  LargeCarbonate();
  ~LargeCarbonate();

  void setup(std::vector<double> *total);


private:

};

#endif
