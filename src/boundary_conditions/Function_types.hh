#ifndef _FUNCTION_TYPES_H_
#define _FUNCTION_TYPES_H_


namespace Boundary
{

/* Define some parent types for boundary condition functions of
   various arguments. This will allow us to write boundary
   condition evaluators which specialize on the kinds of functions
   they are give.
*/

typedef double* Space;
typedef double Time;
typedef double* Output;
typedef unsigned int Face;

// These abstract base classes are required only for implementations
// which require run-time dispatch. Any type which conforms the the
// same interface can subsistute for these base classes as a template
// parameter.


struct Space_function      { virtual void operator () (Space x,         Output c) = 0; };
struct Time_function       { virtual void operator () (         Time t, Output c) = 0; };
struct Space_time_function { virtual void operator () (Space x, Time t, Output c) = 0; };

}


#endif /* _FUNCTION_TYPES_H_ */
