import numpy as np
from scipy import optimize

import wrm_vangenuchten
import capillary_pressure

class PermafrostModel(object):
    def __init__(self, wrm=None):
        if wrm is None:
            wrm = wrm_vangenuchten.VanGenuchten(0.205, 5.5e-4)
        self.wrm = wrm

    def _si_functor(self, pcliq, pcice):
        def functor(si):
            tmp = (1.0 - si) * self.wrm.saturation(pcliq);
            return tmp - self.wrm.saturation(pcice + self.wrm.capillaryPressure(tmp + si))
        return functor

    def _si_frozen_unsaturated(self, pcliq, pcice):
        func = self._si_functor(pcliq, pcice)
        s_i,r = optimize.brentq(func, 0., 1., full_output=True)
        assert r.converged
        return s_i

    def _dsi_dpcliq_frozen_unsaturated(self, pcliq, pcice, si):
        sstar =  self.wrm.saturation(pcliq);
        sstarprime = self.wrm.d_saturation(pcliq);
        tmp = (1.0 - si) * sstar;
        tmpprime = (1.0 - si) * sstarprime

        G = -self.wrm.d_saturation( pcice + self.wrm.capillaryPressure( tmp + si)) \
          * self.wrm.d_capillaryPressure( tmp + si )

        numer = tmpprime * (1 + G)
        denom = - sstar + G * ( 1.0 - sstar)
        return - numer / denom

    def _dsi_dpcice_frozen_unsaturated(self, pcliq, pcice, si):
        sstar =  self.wrm.saturation(pcliq)
        tmp = (1.0 - si) * sstar

        G1 = self.wrm.d_saturation( pcice + self.wrm.capillaryPressure( tmp + si))
        G2 = self.wrm.d_capillaryPressure( tmp + si )

        return -G1 / (sstar + G1*G2*(1-sstar))

    #-----------

    def _sats_frozen_unsaturated(self, pcliq, pcice):
        si = self._si_frozen_unsaturated(pcliq, pcice)
        sl = (1.0 - si)*self.wrm.saturation(pcliq);
        sg = 1. - si - sl
        return sg, sl, si

    def _dsats_dpcliq_frozen_unsaturated(self, pcliq, pcice):
        si = self._si_frozen_unsaturated(pcliq, pcice)

        dsi = self._dsi_dpcliq_frozen_unsaturated(pcliq, pcice, si)
        dsl = (1.0 - si) * self.wrm.d_saturation(pcliq) - dsi * self.wrm.saturation(pcliq)
        dsg = - dsi - dsl
        return dsg, dsl, dsi

    def _dsats_dpcice_frozen_unsaturated(self, pcliq, pcice):
        si = self._si_frozen_unsaturated(pcliq, pcice)

        dsi = self._dsi_dpcice_frozen_unsaturated(pcliq, pcice, si)
        dsl = - dsi * self.wrm.saturation(pcliq)
        dsg = - dsi - dsl
        return dsg, dsl, dsi

    #-----------

    def saturations_pc(self, pcliq, pcice):
        if pcice <= 0.:
            # saturations if unfrozen
            si = 0.
            sl = self.wrm.saturation(pcliq)
            sg = 1. - sl
            assert 0. <= sg <= 1.
            assert 0. <= sl <= 1.
            assert 0. <= si <= 1.

            return sg, sl, si

        elif pcliq <= 0.:
            # saturations if saturated
            sl = self.wrm.saturation(pcice)
            si = 1. - sl
            sg = 0.
            return sg, sl, si

        else:
            # saturations if frozen and unsaturated
            return self._sats_frozen_unsaturated(pcliq, pcice)

    def saturations_Tp(self, T, p):
        pcliq = capillary_pressure.pc_liq(T,p)
        pcice = capillary_pressure.pc_ice(T,p)
        return self.saturations_pc(pcliq,pcice)

    def dsaturations_dpcliq_pc(self, pcliq, pcice):
        if pcice <= 0.:
            # derivatives if unfrozen
            dsi = 0.
            dsl = self.wrm.d_saturation(pcliq)
            dsg = - dsl
            return dsg, dsl, dsi

        elif pcliq <= 0.:
            # derivatives if saturated
            dsl = 0.
            dsi = 0.
            dsg = 0.
            return dsg, dsl, dsi

        else:
            return self._dsats_dpcliq_frozen_unsaturated(pcliq, pcice)

    def dsaturations_dpcliq_Tp(self, T, p):
        pcliq = capillary_pressure.pc_liq(T,p)
        pcice = capillary_pressure.pc_ice(T,p)
        return self.dsaturations_dpcliq_pc(pcliq, pcice)

    def dsaturations_dpcice_pc(self, pcliq, pcice):
        if pcice <= 0.:
            # derivatives if unfrozen
            dsi = 0.
            dsl = 0.
            dsg = 0.
            return dsg, dsl, dsi

        elif pcliq <= 0.:
            # derivatives if saturated
            dsl = self.wrm.d_saturation(pcice)
            dsi = - dsl
            dsg = 0.
            return dsg, dsl, dsi

        else:
            return self._dsats_dpcice_frozen_unsaturated(pcliq, pcice)

    def dsaturations_dpcice_Tp(self, T, p):
        pcliq = capillary_pressure.pc_liq(T,p)
        pcice = capillary_pressure.pc_ice(T,p)
        return self.dsaturations_dpcice_pc(pcliq, pcice)
