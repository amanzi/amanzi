#ifndef __Block_hpp__
#define __Block_hpp__

#include <math.h>
#include <iostream>
#include <iomanip>
using namespace std;

// Boost may provide us with a more optimal matrix implementation - Glenn

class Block {
  
public:
  Block();
  Block(int);
  virtual ~Block();

  int getSize(void) const { return this->size; };
  double **getValues(void) { return this->A; };

  double getRowAbsMax(int);

  void setValue(int, int, double);
  void setValues(double **);
  void setValues(int, int, Block *);

  void addValue(int, int, double);
  void addValues(double **);
  void addValues(int, int, Block *);

  void scaleRow(int, double);
  void scaleColumn(int, double);
  void scale(double);

  void zero(void);
  void setDiagonal(double d);

  void print(void);

  
private:

  int size;
  double **A;

};

#endif /*__Block_hpp__*/

