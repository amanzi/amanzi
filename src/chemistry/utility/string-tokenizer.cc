/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/*******************************************************************************
**
**  File Name: StringTokenizer.cpp
**
**  Description: This class is a simple string tokenizer.
**
**  Notes:
**    - Found July 2003, unencumbered, on the web at:
**
**      http://www.thecodezone.com/diary/archives/000057.html
**
**    - no longer available from origional source, but archived at:
**
**      http://web.archive.org/web/20030810163805/http://www.thecodezone.com/diary/archives/000057.html
**
*******************************************************************************/
#include "string-tokenizer.hh"

#include <string>

StringTokenizer::StringTokenizer(void) {
} /* end StringTokenizer() */

StringTokenizer::StringTokenizer(const std::string &source,
                                 const std::string &delimiters) {
  tokenize(source, delimiters);
} /* StringTokenizer(source, delimiters) */

void StringTokenizer::tokenize(const std::string &source,
                               const std::string &delimiters) {
  clear();
  std::string::size_type spos(source.find_first_not_of(delimiters, 0));
  std::string::size_type epos(source.find_first_of(delimiters, spos));

  while (std::string::npos != epos || std::string::npos != spos) {
    push_back(source.substr(spos, epos - spos));
    spos = source.find_first_not_of(delimiters, epos);
    epos = source.find_first_of(delimiters, spos);
  }
} /* end tokenize(source, delimitiers) */
