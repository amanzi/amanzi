/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@ornl.gov)
*/

#include "Teuchos_ParameterList.hpp"
#include "UnitTest++.h"

#include "errors.hh"
#include "Key.hh"

using namespace Amanzi;
using namespace Amanzi::Keys;

SUITE(UTILS_KEYS) {
TEST(BASICS) {

  CHECK(starts_with("hello", "hel"));
  CHECK(starts_with("*!.,hello", "*!.,"));
  CHECK(starts_with(" ", " "));
  CHECK(starts_with("123", "12"));
  CHECK(!starts_with("hello", "bye"));
  CHECK(!starts_with("hel", "hello"));
  CHECK(!starts_with(" hel", "hel"));

  CHECK(ends_with("hello", "o"));
  CHECK(ends_with("hello", "lo"));
  CHECK(ends_with(" ", " "));
  CHECK(!ends_with("hello", "hel"));
  CHECK(!ends_with("hel", "hello"));
  CHECK(!ends_with("hello ", "hello"));

  CHECK(in("hello", 'h'));
  CHECK(in("hello", 'l'));
  CHECK(in("h*llo", '*'));
  CHECK(!in("hello", '*'));
  CHECK(!in("hello", '2'));

  CHECK_EQUAL("hi-bye", merge("hi", "bye", '-'));
  CHECK_EQUAL("hi2bye", merge("hi", "bye", '2'));
  CHECK_EQUAL("hi:bye", merge("hi", "bye", ':'));
  CHECK_EQUAL("hi|bye", merge("hi", "bye", '|'));
}

TEST(KEY_OPERATIONS) {

  CHECK_EQUAL("surface-temp", getKey("surface", "temp"));
  CHECK_EQUAL("surface:4-temp", getKey("surface:4", "temp"));
  CHECK_EQUAL("temp", getKey("", "temp"));
  CHECK_EQUAL("temp", getKey("domain", "temp"));
  CHECK_THROW(getKey("domain", ""), Errors::Message);
  CHECK_THROW(getKey("bad-name",""), Errors::Message);
  // this can't throw because of things like getVarName("dsurface-energy|dsurface-temp") which must return "surface"
  //  CHECK_THROW(getKey("", "bad-name"), Errors::Message);

  auto result1 = KeyPair{"surface", "temp"};
  CHECK(result1 == splitKey("surface-temp"));

  auto result2 = KeyPair{"surface:4", "temp"};
  CHECK(result2 == splitKey("surface:4-temp"));

  auto result3 = KeyPair{"","temp"};
  CHECK(result3 == splitKey("temp"));

  CHECK_EQUAL("surface", getDomain("surface-temp"));
  CHECK_EQUAL("surface:4", getDomain("surface:4-temp"));
  CHECK_EQUAL("domain", getDomain("temp"));
  CHECK_EQUAL("domain", getDomain("domain-temp"));
  CHECK_EQUAL("surface", getDomain("dsurface-energy|dsurface-temp"));
  CHECK_THROW(getDomain("bad-name-is-bad"), Errors::Message);

  CHECK_EQUAL("temp", getVarName("temp"));
  CHECK_EQUAL("temp", getVarName("surface-temp"));
  CHECK_EQUAL("temp", getVarName("domain-temp"));
  CHECK_EQUAL("temp", getVarName("ddomain-temp|dsurface-temp"));
  CHECK_THROW(getVarName("bad-name-is-bad"), Errors::Message);

  CHECK_EQUAL("surface:4", getDomainInSet("surface", "4"));
  CHECK_EQUAL("surface:4", getDomainInSet("surface", 4));
  CHECK_EQUAL("surface:*", getDomainInSet("surface", "*"));
  CHECK_EQUAL("surface:upstream", getDomainInSet("surface", "upstream"));

  CHECK_EQUAL("surface", getDomainSetName("surface:4"));
  CHECK_EQUAL("surface", getDomainSetName("surface:upstream"));
  CHECK_THROW(getDomainSetName("surface"), Errors::Message);

  CHECK_EQUAL("upstream", getDomainSetIndex("surface:upstream"));
  CHECK_EQUAL("4", getDomainSetIndex("surface:4"));
  CHECK_EQUAL(4, getDomainSetIndex<int>("surface:4"));

  CHECK(isDomainSet("surface:4"));
  CHECK(isDomainSet("surface:upstream"));
  CHECK(isDomainSet("surface:*"));
  CHECK(!isDomainSet("surface"));
  CHECK(!isDomainSet(""));
  CHECK(!isDomainSet("*"));
  CHECK(!isDomainSet("surface-temp"));
  CHECK_THROW(isDomainSet("surface:-temp"), Errors::Message);

  {
    KeyTriple ds;
    bool result = splitDomainSet("surface:4-temp", ds);
    CHECK(result);
    CHECK_EQUAL("surface", std::get<0>(ds));
    CHECK_EQUAL("4", std::get<1>(ds));
    CHECK_EQUAL("temp", std::get<2>(ds));
  }
  {
    KeyTriple ds;
    bool result = splitDomainSet("surface:*-temp", ds);
    CHECK(result);
    CHECK_EQUAL("surface", std::get<0>(ds));
    CHECK_EQUAL("*", std::get<1>(ds));
    CHECK_EQUAL("temp", std::get<2>(ds));
  }
  {
    KeyTriple ds;
    bool result = splitDomainSet("surface:*", ds);
    CHECK(result);
    CHECK_EQUAL("surface", std::get<0>(ds));
    CHECK_EQUAL("*", std::get<1>(ds));
    CHECK_EQUAL("", std::get<2>(ds));
  }

  CHECK_EQUAL("surface:4-temp", getKey("surface","4","temp"));
  CHECK_EQUAL("surface:4-temp", getKey("surface",4,"temp"));
  CHECK_EQUAL("surface:*-temp", getKey("surface","*","temp"));

  CHECK(matchesDomainSet("surface", "surface:4-temp"));

  CHECK_EQUAL("dwater_content|dpressure", getDerivKey("water_content", "pressure"));
  CHECK_EQUAL("dsurface-water_content|dpressure", getDerivKey("surface-water_content", "pressure"));
}


TEST(READS_DOMAINS) {
  {
    // cleans names
    Teuchos::ParameterList plist;
    CHECK_EQUAL("isgood", cleanPListName(plist.sublist("my_sub").sublist("isgood")));
  }
  {
    // cleans names
    Teuchos::ParameterList plist;
    CHECK_EQUAL("domain-isgood", cleanPListName(plist.sublist("my_sub").sublist("domain-isgood")));
  }

  {
    // domain names
    Teuchos::ParameterList plist;
    plist.set("domain name", "surface");
    CHECK_EQUAL("surface", readDomain(plist));
    CHECK_THROW(readDomain(plist, "snow"), std::exception);
  }
  {
    Teuchos::ParameterList plist;
    plist.set("surface domain name", "my_surface");
    CHECK_EQUAL("my_surface", readDomain(plist, "surface"));
    CHECK_THROW(readDomain(plist), std::exception);
  }

  {
    Teuchos::ParameterList plist;
    CHECK_EQUAL("my_surface", readDomain(plist, "surface", "my_surface"));
    CHECK_EQUAL("my_surface", plist.get<std::string>("surface domain name")); // default set into list
  }


  // identity
  {
    Teuchos::ParameterList plist;
    CHECK_EQUAL("surface", readDomainHint(plist, "surface", "surface", "surface"));
  }
  {
    Teuchos::ParameterList plist;
    CHECK_EQUAL("surface:4", readDomainHint(plist, "surface:4", "surface", "surface"));
  }


  {
    Teuchos::ParameterList plist;
    CHECK_EQUAL("snow", readDomainHint(plist, "surface", "surface", "snow"));
  }
  {
    Teuchos::ParameterList plist;
    CHECK_EQUAL("domain", readDomainHint(plist, "surface", "surface", "subsurface"));
  }
  {
    Teuchos::ParameterList plist;
    CHECK_EQUAL("domain", readDomainHint(plist, "surface", "surface", "subsurface"));
  }
  {
    Teuchos::ParameterList plist;
    CHECK_EQUAL("domain", readDomainHint(plist, "surface", "surface", "domain"));
  }
  {
    Teuchos::ParameterList plist;
    CHECK_EQUAL("domain", readDomainHint(plist, "surface", "surface", "domain"));
  }
  {
    Teuchos::ParameterList plist;
    CHECK_EQUAL("surface", readDomainHint(plist, "", "domain", "surface"));
  }
  {
    Teuchos::ParameterList plist;
    CHECK_EQUAL("surface:4", readDomainHint(plist, "domain:4", "domain", "surface"));
  }
  {
    Teuchos::ParameterList plist;
    CHECK_EQUAL("surface_column:4", readDomainHint(plist, "column:4", "domain", "surface"));
  }
  {
    Teuchos::ParameterList plist;
    CHECK_EQUAL("column:4", readDomainHint(plist, "surface_column:4", "surface", "subsurface"));
  }
  {
    Teuchos::ParameterList plist;
    CHECK_EQUAL("column:4", readDomainHint(plist, "surface_column:4", "surface", "subsurface"));
  }


  {
    Teuchos::ParameterList plist;
    CHECK_EQUAL("snow:4", readDomainHint(plist, "surface:4", "surface", "snow"));
  }
  {
    Teuchos::ParameterList plist;
    CHECK_EQUAL("domain:4", readDomainHint(plist, "surface:4", "surface", "subsurface"));
  }
  {
    Teuchos::ParameterList plist;
    CHECK_EQUAL("domain:4", readDomainHint(plist, "surface:4", "surface", "subsurface"));
  }
  {
    Teuchos::ParameterList plist;
    CHECK_EQUAL("domain:4", readDomainHint(plist, "surface:4", "surface", "domain"));
  }
  {
    Teuchos::ParameterList plist;
    CHECK_EQUAL("domain:4", readDomainHint(plist, "surface:4", "surface", "domain"));
  }

  {
    Teuchos::ParameterList plist;
    plist.set<std::string>("snow domain name", "user_provided");
    CHECK_EQUAL("user_provided", readDomainHint(plist, "surface:4", "surface", "snow"));
  }
  {
    Teuchos::ParameterList plist;
    plist.set<std::string>("subsurface domain name", "user_provided");
    CHECK_EQUAL("user_provided", readDomainHint(plist, "surface:4", "surface", "subsurface"));
  }
}

TEST(READS_KEYS) {
  {
    // reads keys
    Teuchos::ParameterList plist;
    CHECK_EQUAL("surface-pressure", readKey(plist, "surface", "pressure", "pressure"));
    CHECK_EQUAL("temperature", readKey(plist, "domain", "temperature", "temperature"));
    CHECK_EQUAL("surface-air_temperature", readKey(plist, "surface", "air temperature", "air_temperature"));
  }

  {
    // reads keys
    Teuchos::ParameterList plist;
    plist.set<std::string>("pressure key", "notsurface-notpressure");
    plist.set<std::string>("temperature key suffix", "nottemperature");

    CHECK_EQUAL("notsurface-notpressure", readKey(plist, "surface", "pressure", "pressure"));
    CHECK_EQUAL("surface-nottemperature", readKey(plist, "surface", "temperature", "temperature"));
  }

  {
    // key takes precendence
    Teuchos::ParameterList plist;
    plist.set<std::string>("pressure key", "notsurface-notpressure");
    plist.set<std::string>("pressure key suffix", "spew");
    CHECK_EQUAL("notsurface-notpressure", readKey(plist, "surface", "pressure", "pressure"));
  }
}

}
