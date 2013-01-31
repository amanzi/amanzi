#ifndef _Material_H_
#define _Material_H_

#include <PArray.H>
#include <Region.H>
#include <TabularFunction.H>

class Property
{
public:

  enum CoarsenRule { Arithmetic = 0, ComponentHarmonic = 1, INVALID_CR = 2 };
  enum RefineRule  { PiecewiseConstant = 0, INVALID_RR = 3};

  Property(const std::string&    _name,
           Property::CoarsenRule _coarsen_rule = Property::Arithmetic,
           Property::RefineRule  _refine_rule = Property::PiecewiseConstant)
    : name(_name), coarsen_rule(_coarsen_rule), refine_rule(_refine_rule) {}

  virtual void eval (Real t, int level, const Box& box, FArrayBox& fab, int dComp, void* ctx = 0) const = 0;

  virtual ~Property() {}
  virtual Property* clone() const = 0;
  virtual int nComp() const = 0;
  const std::string& Name() const {return name;}
  CoarsenRule coarsenRule() const {return coarsen_rule;}
  RefineRule refineRule() const {return refine_rule;}
protected:
  std::string name;
  CoarsenRule coarsen_rule;
  RefineRule refine_rule;
};


class ConstantProperty
  : public Property
{
public:
  ConstantProperty(const std::string&    _name,
                   Real                  _value,
                   Property::CoarsenRule _coarsen_rule = Property::Arithmetic,
                   Property::RefineRule  _refine_rule = Property::PiecewiseConstant)
    : Property(_name,_coarsen_rule,_refine_rule) {values = Array<Real>(1,_value);}

  ConstantProperty(const std::string&    _name,
                   const Array<Real>&    _values,
                   Property::CoarsenRule _coarsen_rule = Property::Arithmetic,
                   Property::RefineRule  _refine_rule = Property::PiecewiseConstant)
    : Property(_name,_coarsen_rule,_refine_rule) {values = _values;}

  virtual void eval (Real t, int level, const Box& box, FArrayBox& fab, int dComp, void* ctx = 0) const;

  virtual ~ConstantProperty() {}
  virtual Property* clone() const;
  virtual int nComp() const {return values.size();}
  const Array<Real> Values() const {return values;}
protected:
  Array<Real> values;
};



class TabularInTimeProperty
  : public Property
{
public:
  TabularInTimeProperty(const std::string&     _name,
                        const TabularFunction& _func,
                        Property::CoarsenRule  _coarsen_rule = Property::Arithmetic,
                        Property::RefineRule   _refine_rule = Property::PiecewiseConstant)
    : Property(_name,_coarsen_rule,_refine_rule), funcs(Array<TabularFunction>(1,_func)) {}

  TabularInTimeProperty(const std::string&            _name, 
                        const Array<TabularFunction>& _funcs,
                        Property::CoarsenRule         _coarsen_rule = Property::Arithmetic,
                        Property::RefineRule          _refine_rule = Property::PiecewiseConstant)
    : Property(_name,_coarsen_rule,_refine_rule), funcs(_funcs) {}

  virtual void eval (Real t, int level, const Box& box, FArrayBox& fab, int dComp, void* ctx = 0) const;

  virtual ~TabularInTimeProperty() {}
  virtual Property* clone() const;
  virtual int nComp() const {return funcs.size();}
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

