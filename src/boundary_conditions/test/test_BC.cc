#include <UnitTest++.h>

#include "../Evaluator_types.hh"

struct simple_time
{
    void operator () (Boundary::Time t, Boundary::Output c) { *c = 10*t; }
};

struct simple_space
{
    void operator () (Boundary::Space x, Boundary::Output c) { *c = *x + 1; } 
};


struct simple_space_time
{
    void operator () (Boundary::Space x, Boundary::Time t, Boundary::Output c) { *c = *x + t; }
};




/* A very simple face evaluator, which converts a face index f into
   f+1/2 (like a cell center) */

struct simple_face_evaluator
{
    template <typename F>
    void operator () (F& bc, Boundary::Face f, Boundary::Time t, Boundary::Output c) const 
    {  double x = f+0.5; bc ( &x, t, c); } 
};

struct timeless_face_evaluator
{
    template <typename F>
    void operator () (F& bc, Boundary::Face f, Boundary::Output c) const
    {
        double x = f+0.5; bc (&x, c);
    }
};


SUITE (Structs)
{
    
    TEST (Time_test)
    {
        simple_time time_bc;
        
        Boundary::BC_evaluator *bc = Boundary::build_time_BC (time_bc);

        double result;

        (*bc) (NULL, 0.0, &result);
        CHECK_EQUAL (result, 0.0);

        (*bc) (NULL, 1.0, &result);
        CHECK_EQUAL (result, 10.0);
    }


    TEST (Space_test)
    {
        simple_space space_bc;
        timeless_face_evaluator face_ev;

        Boundary::BC_evaluator *bc = Boundary::build_space_BC (space_bc, face_ev);

        double result1, result2;

        (*bc) (1, 0.0, &result1);
        CHECK_EQUAL (result1, 2.5);

        (*bc) (1, 10.0, &result2);
        CHECK_EQUAL (result2, 2.5);

        CHECK_EQUAL (result1, result2);
    }



    TEST (Space_time_test)
    {
        simple_space_time space_time_bc;
        simple_face_evaluator sfe;
        
        Boundary::BC_evaluator *bc = Boundary::build_space_time_BC (space_time_bc, sfe);
        
        double result;

        (*bc) (1, 1.0, &result);
        CHECK_EQUAL (result, 2.5);

        (*bc) (2, 2.0, &result);
        CHECK_EQUAL (result, 4.5);
    }

}
