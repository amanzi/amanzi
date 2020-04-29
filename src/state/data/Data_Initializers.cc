/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//! Data_Initializers are functions that know how to initialize data.

/*
  Initialize from parameter list.
*/

#include "Data_Initializers.hh"

namespace Amanzi {
namespace Data_Initializers {

template <>
bool
Initialize<double>(Teuchos::ParameterList& plist,
                   const Teuchos::ParameterList& attrs, double& t)
{
  return InitializePrimitiveByValue(plist, attrs, t);
}

template <>
bool
Initialize<int>(Teuchos::ParameterList& plist,
                const Teuchos::ParameterList& attrs, int& t)
{
  return InitializePrimitiveByValue(plist, attrs, t);
}

template <>
bool
Initialize<std::string>(Teuchos::ParameterList& plist,
                        const Teuchos::ParameterList& attrs, std::string& t)
{
  return InitializePrimitiveByValue(plist, attrs, t);
}


template <>
bool
Initialize<CompositeVector_<int> >(Teuchos::ParameterList& plist, const Teuchos::ParameterList& attrs,
           CompositeVector_<int>& t)
{
  if (InitializeVectorByValue<CompositeVector_, int>(plist, attrs, t))
    return true;
  return false;
}


template <>
bool
Initialize<CompositeVector>(Teuchos::ParameterList &plist,
                            const Teuchos::ParameterList& attrs,
                            CompositeVector &t)
{
  bool partially_initialized = false;

  // set a default value
  if (InitializeVectorByValue<CompositeVector_, double>(plist, attrs, t))
    return true;
  
  // First try all initialization method which set the entire data structure.
  // ------ Try to set values from a restart file -----
  if (plist.isParameter("restart file")) {
    auto chkp_plist = Teuchos::rcp(new Teuchos::ParameterList());
    chkp_plist->set("file name", plist.get<std::string>("restart file"));
    Checkpoint chkp(chkp_plist, t.Comm());
    chkp.Read(attrs, t);
    return true;
  }

  // // ------ Try to set values from a file -----
  // if (plist.isSublist("exodus file initialization")) {
  //   // data must be pre-initialized to zero in case Exodus file does not
  //   // provide all values.
  //   t.putScalar(0.0);

  //   Teuchos::ParameterList& file_list = plist.sublist("exodus file
  //   initialization"); Functions::ReadExodusIIMeshFunction(file_list, t);
  //   return true;
  // }

  // Next try all partial initialization methods -- typically cells.
  // ------ Try to set cell values from a restart file -----
  // FIXME EPETRA TO TPETRA: read file
  // if (plist.isParameter("cells from file")) {
  //   auto filename = plist.get<std::string>("cells from file");
  //   Checkpoint chkp(filename, t.Comm());

  //   // read just the cells
  //   auto vec_c = t.ViewComponent("cell", false);
  //   for (int i = 0; i != t.getNumVectors("cell"); ++i) {
  //     std::stringstream name;
  //     name << fieldname << ".cell." << i;
  //     // FIXME EPETRA TO TPETRA
  //     // chkp.Read(name.str(), *vec_c(i));
  //   }
  //   chkp.Finalize();
  // }

  // ------ Set values from 1D solution -----
  // if (plist.isSublist("initialize from 1D column")) {
  //   Teuchos::ParameterList &init_plist =
  //       plist.sublist("initialize from 1D column");
  //   if (!init_plist.isParameter("f header"))
  //     init_plist.set("f header", std::string("/") + fieldname);
  //   Functions::ReadColumnMeshFunction(init_plist, t);
  //   partially_initialized = true;
  // }

  // ------ Set values using a function -----
  if (plist.isSublist("function")) {
    Teuchos::ParameterList func_plist = plist.sublist("function");
    auto func = Functions::createCompositeVectorFunction(func_plist, t.getMap()->Mesh());
    func->Compute(0.0, t);
    partially_initialized = true;
  }

  //   if (partially_initialized) {
  //     if ((t.HasComponent("face") || t.HasComponent("boundary_face")) &&
  //         t.HasComponent("cell") &&
  //         plist.get<bool>("initialize faces from cells", false)) {
  //       // FIXME EPETRA TO TPETRA
  //       // DeriveFaceValuesFromCellValues(t);
  //       return true;
  //     }
  //   }
  return partially_initialized;
}


} // namespace Data_Initializers
} // namespace Amanzi
