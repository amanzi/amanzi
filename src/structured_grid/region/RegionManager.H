#ifndef _RegionManager_H_
#define _RegionManager_H_

#include <Array.H>

#include <Region.H>

#include <map>
#include <string>

class RegionManager
{
public:
  RegionManager(); // Will search for geometry plo and phi data in input
  RegionManager(const Array<Real>& plo,
                const Array<Real>& phi);
  ~RegionManager();

  Array<const Region*> RegionPtrArray() const;
  Array<const Region*> RegionPtrArray(const Array<std::string>& region_names) const;

protected:
  void Init(const Array<Real>& plo,
            const Array<Real>& phi);
  void Clear();
  Array<Region*> regions;
  std::map<std::string,int> name_to_region_idx;
};

std::ostream& operator<<(std::ostream& os, const RegionManager& rm);

#endif
