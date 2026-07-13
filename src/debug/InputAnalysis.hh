/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*!

Diagnostics that print entity counts and geometric measures (area/volume) for
mesh regions.  Placed under ``mesh`` -> ``MESHNAME`` -> ``region analysis``.

The following parameters specify explicit region lists.  Use ``{*}`` as the
sole entry to print all regions that have already been resolved on that entity
kind (i.e. all regions actually used by PKs after Setup()).

.. admonition:: region-analysis-spec

  * `"cell regions`" ``[Array(string)]`` Regions to analyse on cells (count + volume).
    Use ``{*}`` for all resolved cell regions.

  * `"face regions`" ``[Array(string)]`` Regions to analyse on faces (count + area).
    Use ``{*}`` for all resolved face regions.

  * `"node regions`" ``[Array(string)]`` Regions to analyse on nodes (count only).
    Use ``{*}`` for all resolved node regions.

  The following parameters are also accepted (used primarily by Amanzi input translators):

  * `"used source regions`" ``[Array(string)]`` Analysed on cells.
  * `"used boundary condition regions`" ``[Array(string)]`` Analysed on faces.
  * `"used observation regions`" ``[Array(string)]`` Analysed on cells, falling back to faces.

.. code-block:: xml

  <ParameterList name="mesh">
    <ParameterList name="domain">
      <ParameterList name="region analysis">
        <Parameter name="cell regions" type="Array(string)" value="{*}"/>
        <Parameter name="face regions" type="Array(string)" value="{_bottom_face,_top_face}"/>
      </ParameterList>
    </ParameterList>
  </ParameterList>

*/

#ifndef AMANZI_INPUT_ANALYSIS_HH_
#define AMANZI_INPUT_ANALYSIS_HH_

#include <string>

#include "Teuchos_Array.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "InputAnalysis.hh"
#include "Mesh.hh"
#include "VerboseObject.hh"

namespace Amanzi {

class InputAnalysis {
 public:
  InputAnalysis(Teuchos::RCP<const AmanziMesh::Mesh> mesh, const std::string& domain)
    : mesh_(mesh), domain_(domain), vo_(NULL) {};
  ~InputAnalysis()
  {
    if (vo_ != NULL) delete vo_;
  };

  // main members
  void Init(Teuchos::ParameterList& plist);
  void RegionAnalysis();
  void OutputBCs();

  // supporting members
  template<class Iterator>
  Iterator SelectUniqueEntries(Iterator first, Iterator last);

 private:
  Teuchos::ParameterList* plist_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  std::string domain_;
  VerboseObject* vo_;
};

} // namespace Amanzi

#endif
