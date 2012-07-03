# ifndef _MY_MACRO_HH
# define _MY_MACRO_HH

#include <cassert>

// to get compatibiliy with MAC implementations
//#define LONG_MAX 2147483647

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// INPUT/OUTPUT Macros for debugging codes
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

// to get namespace std
# include<iostream> 
# define OUTPUT std::cout
# define INPUT  std::cin
# define ERRPUT std::cerr

# define PRT(A)      ( OUTPUT << (#A) << " = " << (A) << std::endl )
# define VAL(A)      ( OUTPUT << (#A) << " = " << (A) << " " << std::flush )
# define VALV(I,A)   ( OUTPUT << "  " << (#A) << "(" << I << ") =" << A((I)) << " " <<std::flush )
# define VALA(I,A)   ( OUTPUT << "  " << (#A) << "[" << I << "] =" << A[(I)] << " " <<std::flush )
# define VALM(I,J,A) ( OUTPUT << "  " << (#A) << "(" << I << "," << J << ") =" << A((I),(J)) << std::flush )
# define PRTA(I,A)   ( OUTPUT << "  " << (#A) << "[" << I << "] =" << ( abs(A[(I)])<1.e-14? 0. : A[(I)] )<< std::endl )
# define PRTV(I,A)   ( OUTPUT << "  " << (#A) << "(" << I << ") =" << ( abs(A((I)))<1.e-14? 0. : A((I)) )<< std::endl )
# define PRTM(I,J,A) ( OUTPUT << "  " << (#A) << "(" << I << "," << J << ") =" << ( abs(A((I),(J)))<1.e-14? 0. : A((I),(J)) )<< std::endl )

# define DBGF(A)     do { OUTPUT << #A << std::endl << std::flush ; } while(0)
# define MSG(A)      do { OUTPUT << A ; } while(0)
# define MSGF(A)     do { OUTPUT << A << std::flush << std::endl ; } while(0)
# define INSERT(A)   do { OUTPUT << "insert " << (#A) << " = " ; INPUT >> A ; } while(0)

# define PAUSE								\
  do {									\
    char dummy=' ' ;							\
    ERRPUT << "pause (0 to continue)" ;					\
    while ( dummy != '0' ) { INPUT >> dummy ; }				\
  } while(0)

# define indexType int
# define PRT_VEC(V)     for( indexType i=0; i<V.size(); ++i ) { PRTV(i,V) ; }
# define PRT_ARR(V)     for( indexType i=0; i<V.size(); ++i ) { PRTA(i,V) ; }
# define PRT_ARRAY(N,V) for( indexType i=0; i<N; ++i )        { PRTA(i,V) ; }
# define PRT_MATRIX(MAT)				     \
  do {							     \
    for ( indexType j=0 ; j<MAT.size_j() ; ++j ) {	     \
      MSG("print col --> ") ; PRT(j) ;			     \
      for ( indexType i=0 ; i<MAT.size_i() ; ++i ) {	     \
	PRTM(i,j,MAT) ;					     \
      }							     \
    }							     \
  } while(0)

# define PRT_MATROW(MAT)				     \
  do {							     \
    for ( indexType i=0 ; i<MAT.size_i() ; ++i ) {	     \
      MSG("print row --> ") ; PRT(i) ;			     \
      for ( indexType j=0 ; j<MAT.size_j() ; ++j ) {	     \
	PRTM(i,j,MAT) ;					     \
      }							     \
    }							     \
  } while(0)

# define LINE(A)							\
  do {									\
    for( indexType i=0; i<20; ++i ) { OUTPUT << (#A) ; }		\
    OUTPUT << std::endl ;						\
  } while(0)

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
# endif // end of _MY_MACRO_HH
