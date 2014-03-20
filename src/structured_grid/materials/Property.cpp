#include <Property.H>

Property*
ConstantProperty::clone() const
{
  const ConstantProperty* t = dynamic_cast<const ConstantProperty*>(this);
  BL_ASSERT(t!=0);
  ConstantProperty* ret = new ConstantProperty(t->Name(), t->values, t->coarsen_rule, t->refine_rule);
  return ret;
}

bool
ConstantProperty::Evaluate(Real t, Array<Real>& result) const
{
  int N = values.size();
  result.resize(N);
  for (int i=0; i<N; ++i) {
    result[i] = values[i];
  }
  return true;
}

Property*
GSLibProperty::clone() const
{
  const GSLibProperty* t = dynamic_cast<const GSLibProperty*>(this);
  BL_ASSERT(t!=0);
  GSLibProperty* ret = new GSLibProperty(t->Name(), t->values, t->amrData, t->coarsen_rule, t->refine_rule);
  return ret;
}

bool
GSLibProperty::Evaluate(Real t, Array<Real>& result) const
{
  int N = values.size();
  result.resize(N);
  for (int i=0; i<N; ++i) {
    result[i] = values[i];
  }
  return false;
}


Property*
TabularInTimeProperty::clone() const
{
  const TabularInTimeProperty* t = dynamic_cast<const TabularInTimeProperty*>(this);
  BL_ASSERT(t!=0);
  TabularInTimeProperty* ret = new TabularInTimeProperty(t->Name(), t->Functions(), t->coarsen_rule, t->refine_rule);
  return ret;
}

bool
TabularInTimeProperty::Evaluate(Real t, Array<Real>& result) const
{
  int N = funcs.size();
  result.resize(N);
  for (int i=0; i<N; ++i) {
    result[i] = funcs[i](t);
  }
  return true;
}
