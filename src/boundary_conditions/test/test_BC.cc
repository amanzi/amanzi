#include <UnitTest++.h>

#include "../Evaluator_types.hh"

struct simple_space_time
{
    void operator () (Space x, Time t, Output c) { *c = *x + t; }
};


/* A very simple face evaluator, which converts a face index f into
   f+1/2 (like a cell center) */

struct simple_face_evaluator
{
    template <typename F>
    void operator () (F& bc, Face f, Time t, Output c) const {  double x = f+0.5; bc ( &x, t, c); } 
};


TEST (Space_time_test)
{

    simple_space_time sbc;
    simple_face_evaluator sfe;

    BC_evaluator *bc = build_space_time_ev (sbc, sfe);

    double result;
    (*bc) (1, 1.0, &result);

    CHECK_EQUAL (result, 2.5);

}
