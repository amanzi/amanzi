/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#ifndef _RockManager_F_H_
#define _RockManager_F_H_

#ifdef BL_LANG_FORT
#    define ROCK_MANAGER_PCAP     rmgrpcap
#    define ROCK_MANAGER_INVPCAP  rmgrinvpcap
#    define ROCK_MANAGER_RELPERM  rmgrrelperm
#    define ROCK_MANAGER_DSDPCAP  rmgrdsdpcap
#    define ROCK_MANAGER_RESIDSAT rmgrsr
#    define FORT_FILLPMAT         rmgrfillpmat
#else
#  if defined(BL_FORT_USE_UPPERCASE)
#    define ROCK_MANAGER_PCAP     RMGRPCAP
#    define ROCK_MANAGER_INVPCAP  RMGRINVPCAP
#    define ROCK_MANAGER_RELPERM  RMGRRELPERM
#    define ROCK_MANAGER_DSDPCAP  RMGRDSDPCAP
#    define ROCK_MANAGER_RESIDSAT RMGRSR
#    define FORT_FILLPMAT         RMGRFILLPMAT
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define ROCK_MANAGER_PCAP     rmgrpcap
#    define ROCK_MANAGER_INVPCAP  rmgrinvpcap
#    define ROCK_MANAGER_RELPERM  rmgrrelperm
#    define ROCK_MANAGER_DSDPCAP  rmgrdsdpcap
#    define ROCK_MANAGER_RESIDSAT rmgrsr
#    define FORT_FILLPMAT         rmgrfillpmat
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define ROCK_MANAGER_PCAP     rmgrpcap_
#    define ROCK_MANAGER_INVPCAP  rmgrinvpcap_
#    define ROCK_MANAGER_RELPERM  rmgrrelperm_
#    define ROCK_MANAGER_DSDPCAP  rmgrdsdpcap_
#    define ROCK_MANAGER_RESIDSAT rmgrsr_
#    define FORT_FILLPMAT         rmgrfillpmat_
#  endif

#include <ArrayLim.H>

extern "C" {

  void FORT_FILLPMAT (Real       * q, ARLIM_P(q_lo),  ARLIM_P(q_hi),
                      const int * id, ARLIM_P(id_lo), ARLIM_P(id_hi),
                      const Real * v, ARLIM_P(v_lo),  ARLIM_P(v_hi),
                      const int* myid, const int* lo, const int* hi,
                      const int* nc);
}

#endif
#endif /*_RockManager_F_H_*/
