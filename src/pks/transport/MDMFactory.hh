/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
*/

/*!

The mechanical dispersion coefficient may be described by a number of different
models.  Each is defined on its own region.

.. _mechanical-dispersion-typed-spec:
.. admonition:: mechanical-dispersion-typed-spec

   * `"region`" ``[string]`` Region on which the model is valid.
   * `"mechanical dispersion type`" ``[string]`` Name of the model, see below for options.
   * `"_mechanical_dispersion_type_ parameters`"
     ``[_mechanical_dispersion_type_-spec]`` See below for the required
     parameter spec for each type.

*/

#ifndef PK_TRANSPORT_MDM_FACTORY_HH_
#define PK_TRANSPORT_MDM_FACTORY_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "MDM.hh"
#include "Factory.hh"

namespace Amanzi {
namespace Transport {

class MDMFactory : public Utils::Factory<MDM> {
 public:
  Teuchos::RCP<MDM> Create(Teuchos::ParameterList& plist);
};

} // namespace Transport
} // namespace Amanzi

#endif
