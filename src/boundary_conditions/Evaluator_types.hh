#ifndef _EVALUATOR_TYPES_H_
#define _EVALUATOR_TYPES_H_

#include "Function_types.hh"

/* Define an abstract parent evaluator class and some specialized
   implementations.
*/

namespace Boundary
{

struct BC_evaluator
{
    virtual void operator () (Face f, Time t, Output c) = 0;
};


template <typename F>
struct BC_impl : public BC_evaluator
{

    BC_impl (F& bc) : bc_ (bc) {  }

protected:

    F &bc_;

};


// -----------------------------------------------------------------------
template <typename E>
class Face_evaluator
{
    E &face_evaluator_;

public:

    Face_evaluator (E &face_evaluator) : face_evaluator_ (face_evaluator) {  }

    template <typename F>
    void eval (F &bc, Face f, Output c) { face_evaluator_ (bc, f, c); }
    
    template <typename F>
    void eval (F &bc, Face f, Time t, Output c) { face_evaluator_  (bc, f, t, c);  }

};




// -----------------------------------------------------------------------
template <typename F, typename E>
struct Space_time_evaluator : public BC_impl<F>, Face_evaluator<E>
{

    Space_time_evaluator (F& bc, E& e) : BC_impl<F>(bc) , Face_evaluator<E> (e) {  }
    void operator () (Face f, Time t, Output c) { Face_evaluator<E>::eval(BC_impl<F>::bc_, f, t, c); }

};

// -----------------------------------------------------------------------
template <typename F, typename E>
struct Space_evaluator: public BC_impl<F>, Face_evaluator<E>
{

    Space_evaluator (F &bc, E& e) : BC_impl<F>(bc), Face_evaluator<E> (e) {  }
    void operator () (Face f, Time t, Output c) { Face_evaluator<E>::eval(BC_impl<F>::bc_, f, c); }

};
    

// -----------------------------------------------------------------------
template <typename F>
struct Time_evaluator : public BC_impl<F>
{

    Time_evaluator (F &bc) : BC_impl<F> (bc) {  }
    void operator () (Face f, Time t, Output c) { BC_impl<F>::bc_ (t, c); }

};





template <typename F, typename E>
BC_evaluator* build_space_time_BC(F &bc, E &ev) { return new Space_time_evaluator<F,E>(bc, ev); }

template <typename F, typename E>
BC_evaluator* build_space_BC (F &bc, E &ev) { return new Space_evaluator<F,E>(bc, ev); }

template <typename F>
BC_evaluator* build_time_BC (F &bc) { return new Time_evaluator<F>(bc); }



}



#endif /* _EVALUATOR_TYPES_H_ */
