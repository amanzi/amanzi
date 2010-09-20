#ifndef _DBC_H_
#define _DBC_H_

#include <stdio.h>
#include <stdexcept>
#include <iostream>

// Default values for _DEBUG and _THROW

#ifndef _DEBUG
#define _DEBUG 1
#endif

#ifndef _THROW
#define _THROW 0
#endif


// Assertion classes.

class assertion : public std::logic_error
{
public:
    explicit assertion (std::string const & message) : std::logic_error (message) {  }
};


class DBC_assertion : public assertion
{
public:
    DBC_assertion (std::string const & condition, std::string const & file, int line);
};


// Method to abort running code.
void abort_program (std::string const & cond, std::string const & file,  int const line);


#ifdef _DEBUG

#define ASSERT(bool_expression) if (!(bool_expression)) abort_program( #bool_expression, __FILE__, __LINE__ )
#define REQUIRE(bool_expression) ASSERT (bool_expression)
#define ENSURE(bool_expression) ASSERT (bool_expression)

#else

#define ASSERT(a)
#define REQUIRE(a)
#define ENSURE(a)

#endif /* _DEBUG */

#endif /* _DBC_H_ */

