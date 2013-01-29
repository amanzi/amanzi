#include <Material.H>

Property*
ConstantProperty::clone() const
{
  const ConstantProperty* t = dynamic_cast<const ConstantProperty*>(this);
  BL_ASSERT(t!=0);
  ConstantProperty* ret = new ConstantProperty(t->Name(), t->values);
  return ret;
}

void 
ConstantProperty::eval (Real t, int level, const Box& box, FArrayBox& fab, int dComp, void* ctx) const
{
  for (int i=0, N=values.size(); i<N; ++i) {
    fab.setVal(values[i],box,dComp+i,1);
  }
}


Property*
TabularInTimeProperty::clone() const
{
  const TabularInTimeProperty* t = dynamic_cast<const TabularInTimeProperty*>(this);
  BL_ASSERT(t!=0);
  TabularInTimeProperty* ret = new TabularInTimeProperty(t->Name(), t->Functions());
  return ret;
}

void 
TabularInTimeProperty::eval (Real t, int level, const Box& box, FArrayBox& fab, int dComp, void* ctx) const
{
  for (int i=0, N=funcs.size(); i<N; ++i) {
    Real val = funcs[i](t);
    fab.setVal(val,box,dComp+i,1);
  }
}

Material::Material(const std::string&    _name, 
		   const PArray<Region>& _regions,
		   const std::vector<Property*>& _properties)
  : name(_name)
{
  regions.resize(_regions.size(),PArrayNoManage);
  for (int i=0; i<regions.size(); ++i) {
    regions.set(i,(Region*)&(_regions[i]));
  }
  ClearProperties();
  int nprop = _properties.size();
  properties.resize(nprop);
  for (int i=0; i<nprop; ++i) {
    properties[i] = _properties[i]->clone();
  }
}

const Property* 
Material::Prop(const std::string& pname) const 
{
  int ip=-1;
  for (int i=0; i<properties.size()&&ip<0; ++i) {
    if (properties[i]->Name() == pname) {
      return properties[i];
    }
  }
  std::string str = "No property registered with name: " + pname;
  BoxLib::Abort(str.c_str());
  return 0;
}

Array<std::string> 
Material::PropertyNames() const
{
  int nprop = properties.size();
  Array<std::string> names(nprop);
  for (int i=0; i<nprop; ++i) {
    names[i] = properties[i]->Name();
  }
  return names;
}

void
Material::ClearProperties()
{
  for (int i=0, nprop = properties.size(); i<nprop; ++i) {
    delete properties[i];
  }
  properties.clear();
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
  ClearProperties();
  int nprop = rhs.properties.size();
  properties.resize(nprop);
  for (int i=0; i<nprop; ++i) {
    properties[i] = rhs.properties[i]->clone();
  }
}

Material::~Material()
{
  ClearProperties();
}

void 
Material::setVal(FArrayBox& fab,Real val,int comp, const Real* dx) const
{

  int ng = 0;
  for (int k=0; k<regions.size(); ++k) {
    regions[k].setVal(fab,val,comp,dx,ng);
  }
}

