#include <Property.H>

Property*
ConstantProperty::clone() const
{
  const ConstantProperty* t = dynamic_cast<const ConstantProperty*>(this);
  BL_ASSERT(t!=0);
  ConstantProperty* ret = new ConstantProperty(t->Name(), t->values, t->coarsen_rule, t->refine_rule);
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
  TabularInTimeProperty* ret = new TabularInTimeProperty(t->Name(), t->Functions(), t->coarsen_rule, t->refine_rule);
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

