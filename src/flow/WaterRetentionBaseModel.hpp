#ifndef _WATERRETENTIONBASEMODEL_HH_
#define _WATERRETENTIONBASEMODEL_HH_

#include <string>

class WaterRetentionBaseModel 
{
public:
  
  virtual double k_relative(double p) = 0;
  virtual double saturation(double p) = 0;
  virtual double d_saturation(double p) = 0;
  
  const std::string region() { return reg; };

  virtual double pressure(double sl) = 0;

protected:
  void set_region( const std::string region_ ) { reg = region_; };
  std::string reg;

};
  
#endif
  
