#ifndef _WATERRETENTIONBASEMODEL_HH_
#define _WATERRETENTIONBASEMODEL_HH_


class WaterRetentionBaseModel 
{
public:
  
  virtual double k_relative(double p) = 0;
  virtual double saturation(double p) = 0;
  virtual double d_saturation(double p) = 0;
  
  const int mesh_block() { return meshblock; };

protected:
  void set_mesh_block(const int meshblock_) { meshblock = meshblock_; };
  int meshblock; 

};
  
#endif
  
