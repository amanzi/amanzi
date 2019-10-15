/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

#include "UnitTest++.h"
#include "VerboseObject.hh"

#include "Teuchos_ParameterList.hpp"


TEST(VERBOSE_OBJECT_DEFAULTS)
{
  Teuchos::ParameterList plist;
  Amanzi::VerboseObject vo("my verb object", plist);
  CHECK_EQUAL(true, vo.os_OK(Teuchos::VERB_MEDIUM));
  CHECK_EQUAL(false, vo.os_OK(Teuchos::VERB_HIGH));

  Teuchos::OSTab tab = vo.getOSTab();
  if (vo.os_OK(Teuchos::VERB_MEDIUM)) {
    *vo.os() << "object defaults are okay!" << std::endl;
  }
}


TEST(VERBOSE_OBJECT_GLOBAL)
{
  Amanzi::VerboseObject::global_default_level = Teuchos::VERB_HIGH;

  Teuchos::ParameterList plist;
  Amanzi::VerboseObject vo("my verb object", plist);
  CHECK_EQUAL(true, vo.os_OK(Teuchos::VERB_HIGH));
  CHECK_EQUAL(false, vo.os_OK(Teuchos::VERB_EXTREME));

  Teuchos::OSTab tab = vo.getOSTab();
  if (vo.os_OK(Teuchos::VERB_HIGH)) {
    *vo.os() << "object globals are okay!" << std::endl;
  }
}


TEST(VERBOSE_OBJECT_LOCAL)
{
  Teuchos::ParameterList plist;
  Teuchos::ParameterList& vo_plist = plist.sublist("verbose object");
  vo_plist.set("verbosity level", "high");

  Amanzi::VerboseObject vo("my verb object", plist);
  CHECK_EQUAL(true, vo.os_OK(Teuchos::VERB_HIGH));
  CHECK_EQUAL(false, vo.os_OK(Teuchos::VERB_EXTREME));

  Teuchos::OSTab tab = vo.getOSTab();
  if (vo.os_OK(Teuchos::VERB_HIGH)) {
    *vo.os() << "local object vebosity is okay!" << std::endl;
  }
}


TEST(VERBOSE_OBJECT_GET_VERBOSITY)
{
  Teuchos::ParameterList plist;
  Teuchos::ParameterList& vo_plist = plist.sublist("verbose object");
  vo_plist.set("verbosity level", "extreme");

  Amanzi::VerboseObject vo("my verb object", plist);
  CHECK_EQUAL(Teuchos::VERB_EXTREME, vo.getVerbLevel());
  CHECK_EQUAL(false, (Teuchos::VERB_LOW == vo.getVerbLevel()));

  Teuchos::OSTab tab = vo.getOSTab();
  if (vo.os_OK(Teuchos::VERB_HIGH)) {
    *vo.os() << "local object has vebosity EXTREME" << std::endl;
  }
}


TEST(VERBOSE_TWO_OBJECTS)
{
  Teuchos::ParameterList plist;
  Amanzi::VerboseObject vo1("my verb object 1", plist);
  Amanzi::VerboseObject vo2("my verb object 2", plist);

  {
    Teuchos::OSTab tab1 = vo1.getOSTab();
    if (vo1.os_OK(Teuchos::VERB_MEDIUM)) {
      *vo1.os() << "object 1 at level 0" << std::endl;
    }
  }

  {
    Teuchos::OSTab tab1 = vo2.getOSTab();
    {
      Teuchos::OSTab tab2 = vo2.getOSTab();
      if (vo2.os_OK(Teuchos::VERB_HIGH)) {
        *vo2.os() << "object 2 at level 1" << std::endl;
      }
    }

    {
      Teuchos::OSTab tab1 = vo1.getOSTab();
      if (vo1.os_OK(Teuchos::VERB_MEDIUM)) {
        *vo1.os() << "object 1 at level 1" << std::endl;
      }
    }
  }
}


TEST(VERBOSE_TWO_OBJECTS_NOTAB)
{
  Amanzi::VerboseObject::global_hide_line_prefix = true;

  Teuchos::ParameterList plist;
  Amanzi::VerboseObject vo1("my verb object 1", plist);
  Amanzi::VerboseObject vo2("my verb object 2", plist);

  if (vo1.os_OK(Teuchos::VERB_MEDIUM)) {
    *vo1.os() << "test object 1" << std::endl;
  }

  if (vo2.os_OK(Teuchos::VERB_MEDIUM)) {
    *vo2.os() << "test object 2" << std::endl;
  }

  if (vo1.os_OK(Teuchos::VERB_MEDIUM)) {
    *vo1.os() << "test object 1" << std::endl;
  }
}
