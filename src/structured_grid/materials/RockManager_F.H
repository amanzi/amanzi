#ifndef _RockManager_F_H_
#define _RockManager_F_H_

#ifdef BL_LANG_FORT
#    define ROCK_MANAGER_PCAP     rmgrpcap
#    define ROCK_MANAGER_INVPCAP  rmgrinvpcap
#    define ROCK_MANAGER_RELPERM  rmgrrelperm
#    define ROCK_MANAGER_DSDPCAP  rmgrdsdpcap
#else
#  if defined(BL_FORT_USE_UPPERCASE)
#    define ROCK_MANAGER_PCAP     RMGRPCAP
#    define ROCK_MANAGER_INVPCAP  RMGRINVPCAP
#    define ROCK_MANAGER_RELPERM  RMGRRELPERM
#    define ROCK_MANAGER_DSDPCAP  RMGRDSDPCAP
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define ROCK_MANAGER_PCAP     rmgrpcap
#    define ROCK_MANAGER_INVPCAP  rmgrinvpcap
#    define ROCK_MANAGER_RELPERM  rmgrrelperm
#    define ROCK_MANAGER_DSDPCAP  rmgrdsdpcap
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define ROCK_MANAGER_PCAP     rmgrpcap_
#    define ROCK_MANAGER_INVPCAP  rmgrinvpcap_
#    define ROCK_MANAGER_RELPERM  rmgrrelperm_
#    define ROCK_MANAGER_DSDPCAP  rmgrdsdpcap_
#  endif

#endif
#endif /*_RockManager_F_H_*/
