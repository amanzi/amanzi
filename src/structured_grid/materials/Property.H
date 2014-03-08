#ifndef _Property_H_
#define _Property_H_

#include <PArray.H>
#include <Region.H>
#include <TabularFunction.H>

#include <AmrData.H>

class Property
{
public:

  enum CoarsenRule { Arithmetic = 0, ComponentHarmonic = 1, INVALID_CR = 2 };
  enum RefineRule  { PiecewiseConstant = 0, INVALID_RR = 3};

  Property(const std::string&    _name,
           Property::CoarsenRule _coarsen_rule = Property::Arithmetic,
           Property::RefineRule  _refine_rule = Property::PiecewiseConstant)
    : name(_name), coarsen_rule(_coarsen_rule), refine_rule(_refine_rule) {}

  virtual bool isTimeDependent() const = 0;
  virtual bool Evaluate(Real t, Array<Real>& result) const = 0;

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

  virtual bool isTimeDependent() const {return false;}
  virtual bool Evaluate(Real t, Array<Real>& result) const;

  virtual ~ConstantProperty() {}
  virtual Property* clone() const;
  virtual int nComp() const {return values.size();}
protected:
  Array<Real> values;
};


class GSLibProperty
  : public ConstantProperty
{
public:
  GSLibProperty(const std::string&    _name,
                Real                  _value,
                AmrData*              _amrd,
                Property::CoarsenRule _coarsen_rule = Property::Arithmetic,
                Property::RefineRule  _refine_rule = Property::PiecewiseConstant)
    : ConstantProperty(_name,_value,_coarsen_rule,_refine_rule) {amrData = _amrd;}

  GSLibProperty(const std::string&    _name,
                const Array<Real>&    _values,
                AmrData*              _amrd,
                Property::CoarsenRule _coarsen_rule = Property::Arithmetic,
                Property::RefineRule  _refine_rule = Property::PiecewiseConstant)
    : ConstantProperty(_name,_values,_coarsen_rule,_refine_rule) {amrData = _amrd;}

  virtual bool isTimeDependent() const {return false;}
  virtual bool Evaluate(Real t, Array<Real>& result) const;

  virtual ~GSLibProperty() {}
  virtual Property* clone() const;
  virtual int nComp() const {return values.size();}

  const AmrData* GetAmrData() const {return amrData;}
protected:
  AmrData* amrData;
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

  virtual bool isTimeDependent() const {return true;}
  virtual bool Evaluate(Real t, Array<Real>& result) const;

  virtual ~TabularInTimeProperty() {}
  virtual Property* clone() const;
  virtual int nComp() const {return funcs.size();}
  const Array<TabularFunction>& Functions() const {return funcs;}
  
protected:
  Array<TabularFunction> funcs;
};


#endif

