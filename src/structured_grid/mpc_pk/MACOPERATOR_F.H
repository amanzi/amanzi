#ifndef _MACOPERATOR_F_H_
#define _MACOPERATOR_F_H_

#ifdef BL_LANG_FORT
#    define FORT_MACCOEF       maccoef
#    define FORT_MACRHS        macrhs
#    define FORT_MACUPDATE     macupdate
#    define FORT_MACSYNCRHS    macsyncrhs
#    define FORT_RHOGBC        rhogbc
#    define FORT_MACNORML1     macnorml1
#    define FORT_MAC_SYNC      mac_sync
#else
#  if  defined(BL_FORT_USE_UPPERCASE)
#    define FORT_MACCOEF       MACCOEF
#    define FORT_MACRHS        MACRHS
#    define FORT_MACUPDATE     MACUPDATE
#    define FORT_MACSYNCRHS    MACSYNCRHS
#    define FORT_RHOGBC        RHOGBC
#    define FORT_MACNORML1     MACNORML1
#    define FORT_MAC_SYNC      MAC_SYNC
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_MACCOEF       maccoef
#    define FORT_MACRHS        macrhs
#    define FORT_MACUPDATE     macupdate
#    define FORT_MACSYNCRHS    macsyncrhs
#    define FORT_RHOGBC        rhogbc
#    define FORT_MACNORML1     macnorml1
#    define FORT_MAC_SYNC      mac_sync
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_MACCOEF       maccoef_
#    define FORT_MACRHS        macrhs_
#    define FORT_MACUPDATE     macupdate_
#    define FORT_MACSYNCRHS    macsyncrhs_
#    define FORT_RHOGBC        rhogbc_
#    define FORT_MACNORML1     macnorml1_
#    define FORT_MAC_SYNC      mac_sync_
#endif

#include <ArrayLim.H>

extern "C"
{
#if (BL_SPACEDIM == 2)
    void FORT_MACCOEF (Real* cx, ARLIM_P(cxlo), ARLIM_P(cxhi),
                       Real* cy, ARLIM_P(cylo), ARLIM_P(cyhi),
                       const Real* lx, ARLIM_P(lxlo), ARLIM_P(lxhi),
                       const Real* ly, ARLIM_P(lylo), ARLIM_P(lyhi),
                       const Real* ax, ARLIM_P(axlo), ARLIM_P(axhi),
                       const Real* ay, ARLIM_P(aylo), ARLIM_P(ayhi),
                       const int* lo, const int* hi,
                       const int* domlo, const int* domhi,
                       const Real* dx, const int* lo_bc, const int* hi_bc);

    void FORT_MACRHS  (const Real* ux, ARLIM_P(uxlo), ARLIM_P(uxhi),
                       const Real* uy, ARLIM_P(uylo), ARLIM_P(uyhi),
                       const Real* rx, ARLIM_P(rxlo), ARLIM_P(rxhi),
                       const Real* ry, ARLIM_P(rylo), ARLIM_P(ryhi),
                       const Real* ax, ARLIM_P(axlo), ARLIM_P(axhi),
                       const Real* ay, ARLIM_P(aylo), ARLIM_P(ayhi),
                       const Real* vol, ARLIM_P(vlo), ARLIM_P(vhi),
                       Real* rhs, ARLIM_P(rlo), ARLIM_P(rhi),
                       const int* lo, const int* hi,
                       const int* domlo, const int* domhi,
                       const int* lo_bc, const int* hi_bc, const Real* scale);

    void FORT_MACUPDATE (const int *init,
                         Real* ux, ARLIM_P(uxlo), ARLIM_P(uxhi),
                         Real* uy, ARLIM_P(uylo), ARLIM_P(uyhi),
			 const Real* lx, ARLIM_P(lxlo), ARLIM_P(lxhi),
			 const Real* ly, ARLIM_P(lylo), ARLIM_P(lyhi),
                         const Real* phi, ARLIM_P(p_lo), ARLIM_P(p_hi),
                         const int* lo, const int* hi,
                         const int* domlo, const int* domhi,
                         const Real* dx, const int* lo_bc, const int* hi_bc,
                         const Real* scale);

    void FORT_MACNORML1(Real* normval, Real* ux, ARLIM_P(uxlo), ARLIM_P(uxhi),
                        Real* uy, ARLIM_P(uylo), ARLIM_P(uyhi),
			const int* lo, const int* hi,const Real* dx);

    void FORT_MAC_SYNC  (Real* ux, ARLIM_P(uxlo), ARLIM_P(uxhi),
                         Real* uy, ARLIM_P(uylo), ARLIM_P(uyhi),
			 const Real* rx, ARLIM_P(rxlo), ARLIM_P(rxhi),
                         const Real* ry, ARLIM_P(rylo), ARLIM_P(ryhi),
			 const Real* ax, ARLIM_P(axlo), ARLIM_P(axhi),
                         const Real* ay, ARLIM_P(aylo), ARLIM_P(ayhi),
			 const Real* lx, ARLIM_P(lxlo), ARLIM_P(lxhi),
			 const Real* ly, ARLIM_P(lylo), ARLIM_P(lyhi),
                         const Real* phi, ARLIM_P(p_lo), ARLIM_P(p_hi),
                         const int* lo, const int* hi,
                         const Real* dx, const Real* scale);

    void FORT_MACSYNCRHS (Real* rhs, ARLIM_P(rlo), ARLIM_P(rhi),
			  const Real* ux, ARLIM_P(uxlo), ARLIM_P(uxhi),
			  const Real* uy, ARLIM_P(uylo), ARLIM_P(uyhi),
			  const Real* rx, ARLIM_P(rxlo), ARLIM_P(rxhi),
			  const Real* ry, ARLIM_P(rylo), ARLIM_P(ryhi),
			  const Real* ax, ARLIM_P(axlo), ARLIM_P(axhi),
			  const Real* ay, ARLIM_P(aylo), ARLIM_P(ayhi),
			  const Real* vol, ARLIM_P(vlo), ARLIM_P(vhi),
                          const int* lo, const int* hi,
                          const int* domlo, const int* domhi,
                          const int* lo_bc, const int* hi_bc,
                          const Real* rhs_scale);

#endif

#if (BL_SPACEDIM == 3)
    void FORT_MACCOEF (Real* cx, ARLIM_P(cxlo), ARLIM_P(cxhi),
                       Real* cy, ARLIM_P(cylo), ARLIM_P(cyhi),
                       Real* cz, ARLIM_P(czlo), ARLIM_P(czhi),
                       const Real* lx, ARLIM_P(lxlo), ARLIM_P(lxhi),
                       const Real* ly, ARLIM_P(lylo), ARLIM_P(lyhi),
                       const Real* lz, ARLIM_P(lzlo), ARLIM_P(lzhi),
                       const Real* ax, ARLIM_P(axlo), ARLIM_P(axhi),
                       const Real* ay, ARLIM_P(aylo), ARLIM_P(ayhi),
                       const Real* az, ARLIM_P(azlo), ARLIM_P(azhi),
                       const int* lo, const int* hi,
                       const int* domlo, const int* domhi,
                       const Real* dx, const int* lo_bc, const int* hi_bc);

    void FORT_MACRHS  (const Real* ux, ARLIM_P(uxlo), ARLIM_P(uxhi),
                       const Real* uy, ARLIM_P(uylo), ARLIM_P(uyhi),
                       const Real* uz, ARLIM_P(uzlo), ARLIM_P(uzhi),
		       const Real* rx, ARLIM_P(rxlo), ARLIM_P(rxhi),
                       const Real* ry, ARLIM_P(rylo), ARLIM_P(ryhi),
                       const Real* rz, ARLIM_P(rzlo), ARLIM_P(rzhi),
                       const Real* ax, ARLIM_P(axlo), ARLIM_P(axhi),
                       const Real* ay, ARLIM_P(aylo), ARLIM_P(ayhi),
                       const Real* az, ARLIM_P(azlo), ARLIM_P(azhi),
                       const Real* vol, ARLIM_P(vlo), ARLIM_P(vhi),
                       Real* rhs, ARLIM_P(rlo), ARLIM_P(rhi),
                       const int* lo, const int* hi,
                       const int* domlo, const int* domhi,
                       const int* lo_bc, const int* hi_bc, const Real* scale);

    void FORT_MACUPDATE (const int *init,
                         Real* ux, ARLIM_P(uxlo), ARLIM_P(uxhi),
                         Real* uy, ARLIM_P(uylo), ARLIM_P(uyhi),
                         Real* uz, ARLIM_P(uzlo), ARLIM_P(uzhi),
			 const Real* lx, ARLIM_P(lxlo), ARLIM_P(lxhi),
			 const Real* ly, ARLIM_P(lylo), ARLIM_P(lyhi),
			 const Real* lz, ARLIM_P(lzlo), ARLIM_P(lzhi),
                         const Real* phi, ARLIM_P(p_lo), ARLIM_P(p_hi),
                         const int* lo, const int* hi,
                         const int* domlo, const int* domhi,
                         const Real* dx, const int* lo_bc, const int* hi_bc, 
                         const Real* scale); 

    void FORT_MACNORML1(Real* normval, Real* ux, ARLIM_P(uxlo), ARLIM_P(uxhi),
                        Real* uy, ARLIM_P(uylo), ARLIM_P(uyhi),
                        Real* uz, ARLIM_P(uzlo), ARLIM_P(uzhi),
			const int* lo, const int* hi,const Real* dx);

    void FORT_MAC_SYNC  (Real* ux, ARLIM_P(uxlo), ARLIM_P(uxhi),
                         Real* uy, ARLIM_P(uylo), ARLIM_P(uyhi),
                         Real* uz, ARLIM_P(uzlo), ARLIM_P(uzhi),
			 const Real* rx, ARLIM_P(rxlo), ARLIM_P(rxhi),
			 const Real* ry, ARLIM_P(rylo), ARLIM_P(ryhi),
			 const Real* rz, ARLIM_P(rzlo), ARLIM_P(rzhi),
			 const Real* ax, ARLIM_P(axlo), ARLIM_P(axhi),
			 const Real* ay, ARLIM_P(aylo), ARLIM_P(ayhi),
			 const Real* az, ARLIM_P(azlo), ARLIM_P(azhi),
			 const Real* lx, ARLIM_P(lxlo), ARLIM_P(lxhi),
			 const Real* ly, ARLIM_P(lylo), ARLIM_P(lyhi),
			 const Real* lz, ARLIM_P(lzlo), ARLIM_P(lzhi),
                         const Real* phi, ARLIM_P(p_lo), ARLIM_P(p_hi),
                         const int* lo, const int* hi,
                         const Real* dx, const Real* scale);

    void FORT_MACSYNCRHS (Real* rhs, ARLIM_P(rlo), ARLIM_P(rhi),
			  const Real* ux, ARLIM_P(uxlo), ARLIM_P(uxhi),
			  const Real* uy, ARLIM_P(uylo), ARLIM_P(uyhi),
			  const Real* uz, ARLIM_P(uzlo), ARLIM_P(uzhi),
			  const Real* rx, ARLIM_P(rxlo), ARLIM_P(rxhi),
			  const Real* ry, ARLIM_P(rylo), ARLIM_P(ryhi),
			  const Real* rz, ARLIM_P(rzlo), ARLIM_P(rzhi),
			  const Real* ax, ARLIM_P(axlo), ARLIM_P(axhi),
			  const Real* ay, ARLIM_P(aylo), ARLIM_P(ayhi),
			  const Real* az, ARLIM_P(azlo), ARLIM_P(azhi),
			  const Real* vol, ARLIM_P(vlo), ARLIM_P(vhi),
                          const int* lo, const int* hi,
                          const int* domlo, const int* domhi,
                          const int* lo_bc, const int* hi_bc,
                          const Real* rhs_scale);

#endif

    void FORT_RHOGBC(const Real* rho, ARLIM_P(rho_lo), ARLIM_P(rho_hi),
                     const Real* phi, ARLIM_P(phi_lo), ARLIM_P(phi_hi),
                     int * face, const Real* gravity, const Real* dx,
                     const int* domlo, const int* domhi,
                     const int* lo_bc, const int* hi_bc,
		     const Real* press_lo, const Real* press_hi);

}
#endif

#endif /*_MACOPERATOR_F_H_*/
