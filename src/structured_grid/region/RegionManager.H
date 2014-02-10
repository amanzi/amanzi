#ifndef _RegionManager_H_
#define _RegionManager_H_

#include <Array.H>

#include <Region.H>

#include <map>
#include <string>

class RegionManager
{
public:
  RegionManager();
  ~RegionManager();

  Array<const Region*> RegionPtrArray() const;
  Array<const Region*> RegionPtrArray(const Array<std::string>& region_names) const;

protected:
  void Clear();
  Array<Region*> regions;
  std::map<std::string,int> name_to_region_idx;
};

#endif
