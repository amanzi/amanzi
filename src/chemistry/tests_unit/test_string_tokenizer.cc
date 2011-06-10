/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>
#include <string>

#include <UnitTest++.h>

#include "string_tokenizer.hh"

SUITE(GeochemistryTests_StringTokenizer) {
  TEST(TestStringTokenizer_tokenize_default_delimiters) {
    // default delimiters are space characters, " \t\n"
    std::string test_string("This is a string ; Say \"Hi\". \nHello, World!\tWe'll use it for testing.");
    StringTokenizer st;
    st.tokenize(test_string);
    CHECK_EQUAL(st.at(0), "This");
    CHECK_EQUAL(st.at(4), ";");
    CHECK_EQUAL(st.at(6), "\"Hi\".");
    CHECK_EQUAL(st.at(7), "Hello,");
    CHECK_EQUAL(st.at(9), "We'll");
  }  // end TEST()

  TEST(TestStringTokenizer_tokenize_delimiters) {
    std::string test_string("This is a string ; We'll use<it for testing.");
    StringTokenizer st;
    st.tokenize(test_string, ";<");
    CHECK_EQUAL(st.at(0), "This is a string ");
    CHECK_EQUAL(st.at(1), " We'll use");
    CHECK_EQUAL(st.at(2), "it for testing.");
  }  // end TEST()

  TEST(TestStringTokenizer_constructor_tokenize) {
    std::string test_string("This is a string ; We'll use<it for testing.");
    StringTokenizer st(test_string, ";<");
    CHECK_EQUAL(st.at(0), "This is a string ");
    CHECK_EQUAL(st.at(1), " We'll use");
    CHECK_EQUAL(st.at(2), "it for testing.");
  }  // end TEST()

  TEST(TestStringTokenizer_no_delimiter_in_string) {
    // no delimiter, don't break anything up.
    StringTokenizer st("This is a string.", ";<");
    CHECK_EQUAL(st.size(), 1);
    CHECK_EQUAL(st.at(0), "This is a string.");
  }  // end TEST()
}  // end SUITE()
