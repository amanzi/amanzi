/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef STRING_TOKENIZER_HH_

#define STRING_TOKENIZER_HH_

/*******************************************************************************
**
**  File Name: StringTokenizer.h
**
**  Source Control: $Id$
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

#include <string>
#include <vector>


class StringTokenizer : public std::vector<std::string> {
 public:

  StringTokenizer(void);
  StringTokenizer(const std::string &source,
                  const std::string &delimiters = " \t\n");
  void tokenize(const std::string &source,
                const std::string &delimiters = " \t\n");
};

#endif     /* STRING_TOKENIZER_HH_ */

