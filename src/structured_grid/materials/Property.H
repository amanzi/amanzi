/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#ifndef _Property_H_
#define _Property_H_

#include <PArray.H>
#include <Geometry.H>
#include <DataServices.H>

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
                Real                  _avg,
                const std::string&    _param_file,
                const std::string&    _data_file,
                const Array<Real>&    _shift,
                Property::CoarsenRule _coarsen_rule = Property::Arithmetic,
                Property::RefineRule  _refine_rule = Property::PiecewiseConstant)
    : ConstantProperty(_name,_avg,_coarsen_rule,_refine_rule),
      dataServices(0),
      avg(_avg),
      param_file(_param_file), data_file(_data_file) {}

  virtual bool isTimeDependent() const {return false;}
  virtual bool Evaluate(Real t, Array<Real>& result) const;

  virtual ~GSLibProperty() {delete dataServices;}
  virtual Property* clone() const;
  virtual int nComp() const {return num_comps;}

  const AmrData* GetAmrData() const;
  Real Avg() const {return avg;}
  const DataServices* GetDataServices() const {return dataServices;}

  void BuildDataFile(const Array<Geometry>& geom_array,
                     const Array<IntVect>&  ref_ratio,
                     int                    num_grow,
                     int                    max_grid_size_fine_gen,
                     Property::CoarsenRule  crule,
		     const std::string&     varname,
		     bool                   restart);

  const Array<std::string>& PlotfileVars() const {return varnames;}

protected:
  void BuildGSLibFile(Real                   avg,
		      const std::string&     gslib_param_file,
		      const std::string&     gslib_data_file,
		      const Array<Geometry>& geom_array,
		      const Array<IntVect>&  ref_ratio,
		      int                    num_grow,
		      int                    max_grid_size_fine_gen,
		      Property::CoarsenRule  crule);

  DataServices* dataServices;
  Real avg;
  std::string param_file, data_file;
  int num_comps;
  Array<std::string> varnames;
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
