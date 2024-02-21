/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov, Ethan Coon
*/

namespace Amanzi {

PK_BDF_Default::PK_BDF_Default(const Comm_ptr_type& comm,
        Teuchos::ParameterList& pk_tree,
        const Teuchos::RCP<Teuchos::ParameterList>& glist,
        const Teuchos::RCP<State>& S)
  : PK_Default(comm, pk_tree, glist, S) {}

Teuchos::RCP<Operators::Operator>
getOperator(const Operators::Operator_kind& type)
{
  return Teuchos::null;
}

double getDt() override;
  virtual void setDt(double dt) override;



} // namespace Amanzi
