#ifndef _DBC_H_
#define _DBC_H_

#include "exceptions.hh"

namespace DBC
{


/* DBC_assertion
 *
 * An exception class for DBC asserion volations.
 * 
 */
class Assertion : public Exceptions::Amanzi_exception
{
    const char* assertion_;
    const char* filename_;
    unsigned int line_number_;

public:
    Assertion (const char* condition, const char* file, unsigned int line);
    const char* what () const throw ();
};

void amanzi_assert (const char * cond, const char * file,  unsigned int line);

} // namespace DBC


// ---------------------------------
// Macros for DBC checks in the code
// ---------------------------------

#ifdef ENABLE_DBC
#define ASSERT(bool_expression) if (!(bool_expression)) DBC::amanzi_assert( #bool_expression, __FILE__, __LINE__ )
#else
#define ASSERT(a)

#endif /* ENABLE_DBC */


#endif /* _DBC_H_ */

