/* -*-  mode: c++; c-default-style: "google-c-style"; indent-tabs-mode: nil -*- */
#ifndef __STRING_TOKENIZER_HPP__

#define __STRING_TOKENIZER_HPP__

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


/*
** Include Files:
*/

#include <string>
#include <vector>


class StringTokenizer : public std::vector<std::string>
{

public:

    StringTokenizer(void);
    StringTokenizer(const std::string &source,
                    const std::string &delimiters = " \t\n");
    void tokenize(const std::string &source,
                  const std::string &delimiters = " \t\n");

};


#endif     /* __STRING_TOKENIZER_HPP__ */

