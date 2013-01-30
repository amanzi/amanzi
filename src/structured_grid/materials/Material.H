#ifndef _Material_H_
#define _Material_H_

#include <PArray.H>
#include <Region.H>
#include <TabularFunction.H>

class Property
{
public:
  Property(const std::string& _name) {name=_name;}
  virtual ~Property() {}
  virtual Property* clone() const = 0;
  virtual int nComp() const = 0;
  virtual void eval (Real t, int level, const Box& box, FArrayBox& fab, int dComp, void* ctx = 0) const = 0;
  const std::string& Name() const {return name;}
protected:
  std::string name;
};


class ConstantProperty
  : public Property
{
public:
  ConstantProperty(const std::string& _name, Real value)
    : Property(_name) {values = Array<Real>(1,value);}
  ConstantProperty(const std::string& _name, const Array<Real>& _values)
    : Property(_name) {values = _values;}
  virtual Property* clone() const;
  virtual int nComp() const {return values.size();}
  virtual void eval (Real t, int level, const Box& box, FArrayBox& fab, int dComp, void* ctx = 0) const;
  const Array<Real> Values() const {return values;}
protected:
  Array<Real> values;
};



class TabularInTimeProperty
  : public Property
{
public:
  TabularInTimeProperty(const std::string& _name, const TabularFunction& _func)
    : Property(_name), funcs(Array<TabularFunction>(1,_func)) {}
  TabularInTimeProperty(const std::string& _name, const Array<TabularFunction>& _funcs)
    : Property(_name), funcs(_funcs) {}
  virtual Property* clone() const;
  virtual int nComp() const {return funcs.size();}
  virtual void eval (Real t, int level, const Box& box, FArrayBox& fab, int dComp, void* ctx = 0) const;
  const Array<TabularFunction>& Functions() const {return funcs;}
  
protected:
  Array<TabularFunction> funcs;
};

class Material
{
public:
  Material(const std::string&    name, 
	   const PArray<Region>& regions,
	   const std::vector<Property*>& properties = std::vector<Property*>());
  
  Material(const Material& rhs);
  ~Material();
  
  void setVal(FArrayBox& fab,Real val,int comp, const Real* dx) const;
  
  const PArray<Region>& Regions() const {return regions;}
  const std::string& Name() const {return name;}
  const Property* Prop(const std::string& pname) const;
  Array<std::string> PropertyNames() const;

protected:
  void ClearProperties();

  std::string name;
  PArray<Region> regions;
  std::vector<Property*> properties;
};
#endif

