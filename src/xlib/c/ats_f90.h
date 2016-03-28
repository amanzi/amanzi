/*----------------------------------------------------------------------------*
 * ATS
 *
 * License: see $ATS_DIR/COPYRIGHT
 *
 * ATS C Language interface implementation.
 *----------------------------------------------------------------------------*/

#if !defined(_ats_source) && !defined(_include_ats_h)
#error "Error: do not include this file directly, use #include <ats.h>"
#endif

#ifndef ats_f90_h
#define ats_f90_h

int32_t ats_init_f90(int32_t comm);
int32_t ats_finalize_f90();

#endif // ats_f90_h
