#ifndef _Material_H_
#define _Material_H_

#include <PArray.H>
#include <Region.H>
#include <TabularFunction.H>

class Property
{
public:
  virtual ~Property() {}
  virtual Property* clone() const = 0;
  virtual void eval (Real t, int level, const Box& box, FArrayBox& fab, int dComp) const = 0;
};


class ConstantProperty
  : public Property
{
public:
  ConstantProperty(Real _value) {value = _value;}
  virtual Property* clone() const;
  virtual void eval (Real t, int level, const Box& box, FArrayBox& fab, int dComp) const;
protected:
  Real value;
};

class TabularInTimeProperty
  : public Property
{
public:
  TabularInTimeProperty(const TabularFunction& _func) {func=_func;}
  virtual Property* clone() const;
  virtual void eval (Real t, int level, const Box& box, FArrayBox& fab, int dComp) const;
protected:
  TabularFunction func;
};

// Extend the map class to include helpers to deal with keys
template<class map_type>
class key_iterator : public map_type::iterator
{
public:
    typedef typename map_type::iterator map_iterator;
    typedef typename map_iterator::value_type::first_type key_type;

    key_iterator(const map_iterator& other) : map_type::iterator(other) {} ;

    key_type& operator *()
    {
        return map_type::iterator::operator*().first;
    }
};

class Material
{
public:
  Material(const std::string&    name, 
	   const PArray<Region>& regions,
	   const std::map<std::string,Property*>& property_map = std::map<std::string,Property*>());
  
  Material(const Material& rhs);
  ~Material();
  
  void setVal(FArrayBox& fab,Real val,int comp, const Real* dx) const;
  
  const PArray<Region>& Regions() const {return regions;}
  const std::string& Name() const {return name;}
  const Property* Prop(const std::string& pname) const {
    std::map<std::string,Property*>::const_iterator it = property_map.find(pname);
    if (it==property_map.end()) {
      std::string str = "No property registered with name: " + pname;
      BoxLib::Abort(str.c_str());
    }
    return it->second;
  }

protected:
  std::string name;
  PArray<Region> regions;
  std::map<std::string,Property*> property_map;
};
#endif

