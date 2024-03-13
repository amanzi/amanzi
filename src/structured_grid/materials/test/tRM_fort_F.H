/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#ifndef _tRM_fort_F_H_
#define _tRM_fort_F_H_

#ifdef BL_LANG_FORT
#    define FORT_TEST_PCAP    testpcap
#    define FORT_TEST_INVPCAP testinvpcap
#    define FORT_TEST_RELPERM testrelperm
#    define FORT_TEST_DSDPCAP testdsdpcap
#else
#  if defined(BL_FORT_USE_UPPERCASE)
#    define FORT_TEST_PCAP    TESTPCAP
#    define FORT_TEST_INVPCAP TESTINVPCAP
#    define FORT_TEST_RELPERM TESTRELPERM
#    define FORT_TEST_DSDPCAP TESTDSDPCAP
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_TEST_PCAP    testpcap
#    define FORT_TEST_INVPCAP testinvpcap
#    define FORT_TEST_RELPERM testrelperm
#    define FORT_TEST_DSDPCAP testdsdpcap
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_TEST_PCAP    testpcap_
#    define FORT_TEST_INVPCAP testinvpcap_
#    define FORT_TEST_RELPERM testrelperm_
#    define FORT_TEST_DSDPCAP testdsdpcap_
#  endif

#include <ArrayLim.H>

extern "C" {
  void FORT_TEST_PCAP (const int* rmid, const Real* time, const int* lo, const int* hi,
		       const Real *  s, ARLIM_P(s_lo), ARLIM_P(s_hi),
		       const int * mat, ARLIM_P(m_lo), ARLIM_P(m_hi),
		       Real *        p, ARLIM_P(p_lo), ARLIM_P(p_hi));

  void FORT_TEST_INVPCAP (const int* rmid, const Real* time, const int* lo, const int* hi,
			  const Real *  p, ARLIM_P(p_lo), ARLIM_P(p_hi),
			  const int * mat, ARLIM_P(m_lo), ARLIM_P(m_hi),
			  Real *        s, ARLIM_P(s_lo), ARLIM_P(s_hi));

  void FORT_TEST_RELPERM (const int* rmid, const Real* time, const int* lo, const int* hi,
			  const Real *  s, ARLIM_P(s_lo),  ARLIM_P(s_hi),
			  const int * mat, ARLIM_P(m_lo),  ARLIM_P(m_hi),
			  Real *       kr, ARLIM_P(kr_lo), ARLIM_P(kr_hi));

  void FORT_TEST_DSDPCAP (const int* rmid, const Real* time, const int* lo, const int* hi,
			  const Real *  s, ARLIM_P(s_lo),  ARLIM_P(s_hi),
			  const int * mat, ARLIM_P(m_lo),  ARLIM_P(m_hi),
			  Real *     dsdp, ARLIM_P(ds_lo), ARLIM_P(ds_hi));
}

#endif
#endif /*_tRM_fort_H_*/
