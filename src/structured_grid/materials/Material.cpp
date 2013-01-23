#include <Material.H>

Property*
ConstantProperty::clone() const
{
  const ConstantProperty* t = dynamic_cast<const ConstantProperty*>(this);
  BL_ASSERT(t!=0);
  ConstantProperty* ret = new ConstantProperty(t->value);
  return ret;
}

void 
ConstantProperty::eval (Real t, int level, const Box& box, FArrayBox& fab, int dComp) const
{
  fab.setVal(value,box,dComp,1);
}

Property*
TabularInTimeProperty::clone() const
{
  const TabularInTimeProperty* t = dynamic_cast<const TabularInTimeProperty*>(this);
  BL_ASSERT(t!=0);
  TabularInTimeProperty* ret = new TabularInTimeProperty(t->func);
  return ret;
}

void 
TabularInTimeProperty::eval (Real t, int level, const Box& box, FArrayBox& fab, int dComp) const
{
  fab.setVal(func(t),box,dComp,1);
}


Material::Material(const std::string&    _name, 
		   const PArray<Region>& _regions,
		   const std::map<std::string,Property*>& _property_map)
  : name(_name)
{
  regions.resize(_regions.size(),PArrayNoManage);
  for (int i=0; i<regions.size(); ++i) {
    regions.set(i,(Region*)&(_regions[i]));
  }
  for (std::map<std::string,Property*>::const_iterator it=_property_map.begin(), 
         End=_property_map.end(); it!=End; ++it) {
    property_map[it->first] = (it->second)->clone();
  }
}

Material::Material(const Material& rhs) 
{
  name = rhs.name;
  if (rhs.regions.size()>0) {
    regions.resize(rhs.regions.size(),PArrayNoManage);
    for (int i=0; i<regions.size(); ++i) {
      regions.set(i,(Region*)&(rhs.regions[i]));
    }
  }
  for (std::map<std::string,Property*>::const_iterator it=rhs.property_map.begin(), 
         End=rhs.property_map.end(); it!=End; ++it) {
    property_map[it->first] = (it->second)->clone();
  }
}

Material::~Material()
{
  for (std::map<std::string,Property*>::iterator it=property_map.begin(), End=property_map.end(); it!=End; ++it) {
    delete it->second;
  }
  property_map.clear();
}

void 
Material::setVal(FArrayBox& fab,Real val,int comp, const Real* dx) const
{

  int ng = 0;
  for (int k=0; k<regions.size(); ++k) {
    regions[k].setVal(fab,val,comp,dx,ng);
  }
}

