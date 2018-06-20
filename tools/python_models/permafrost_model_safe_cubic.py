import numpy as np

import wrm_vangenuchten
import capillary_pressure
import permafrost_model

class SafePermafrostModel(permafrost_model.PermafrostModel):
    def _determine_cutoff(self, pcliq, pcice):
        # this model fails when (1-si)*sstar + si == 1
        cutoff = np.exp((np.floor(np.log(pcliq))))
        done = False
        while not done:
            si = super(SafePermafrostModel,self)._si_frozen_unsaturated(cutoff, pcice)
            sstar = self.wrm.saturation(cutoff)
            if (1 - si)*sstar + si < (1-1.e-16):
                done = True
            else:
                cutoff = np.exp(np.log(cutoff) + 1.)

        return cutoff, si

    def _fit_spline(self, pcice, cutoff, si_cutoff):
        dsi_cutoff = super(SafePermafrostModel,self)._dsi_dpcliq_frozen_unsaturated(cutoff, pcice, si_cutoff)
        si_saturated = self.saturations_pc(-1., pcice)[2]
        dsi_saturated = (si_cutoff - si_saturated) / cutoff

        # fit a cubic spline to value, deriv at the cutoff and value at the
        # saturated point, of the form:  si = a*pcliq^3 + b*pcliq^2 + c
        # note this ensures dsi/dpc_liq(0) = 0

        # -- at the saturated point, si(0) = si_saturated --> d = si_saturated
        d = si_saturated

        # -- at the saturated point, dsi_dpcliq(0) = dsi_saturated, c = dsi_saturated
        c = dsi_saturated

        # -- at the cutoff point:
        # a * cutoff^3 + b * cutoff^2 + c * cutoff + d = si_cutoff
        # d(a * cutoff^3 + b * cutoff^2 + c * cutoff + d)/dcutoff = dsi_cutoff
        M = np.zeros((2,2),'d')
        rhs = np.zeros((2,),'d')

        M[0,0] = 3*cutoff*cutoff
        M[0,1] = 2*cutoff
        rhs[0] = dsi_cutoff - c

        M[1,0] = cutoff**3
        M[1,1] = cutoff**2
        rhs[1] = si_cutoff - (c*cutoff + d)

        detM = M[0,0]*M[1,1] - M[0,1]*M[1,0]
        a = M[1,1]/detM * rhs[0] + -M[0,1]/detM * rhs[1]
        b = M[0,0]/detM * rhs[1] + -M[1,0]/detM * rhs[0]

        return a,b,c,d

    def _si_frozen_unsaturated(self, pcliq, pcice):
        # check if we are in the spline region
        cutoff, si_cutoff = self._determine_cutoff(pcliq, pcice)
        if pcliq > cutoff:
            si = super(SafePermafrostModel,self)._si_frozen_unsaturated(pcliq, pcice)
            return si
        else:
            # in the spline region, form the spline
            a,b,c,d = self._fit_spline(pcice, cutoff, si_cutoff)
            return ((a * pcliq + b)*pcliq + c) * pcliq + d

    def _dsi_dpcliq_frozen_unsaturated(self, pcliq, pcice, si):
        # check if we are in the spline region
        cutoff, si_cutoff = self._determine_cutoff(pcliq, pcice)
        if pcliq > cutoff:
            dsi = super(SafePermafrostModel,self)._dsi_dpcliq_frozen_unsaturated(pcliq, pcice, si)
            return dsi
        else:
            # in the spline region, form the spline
            a,b,c = self._fit_spline(pcice, cutoff, si_cutoff)
            return (3 * a * pcliq + 2 * b) * pcliq + c

    def _dsi_dpcice_frozen_unsaturated(self, pcliq, pcice, si):
        # check if we are in the spline region
        cutoff, si_cutoff = self._determine_cutoff(pcliq, pcice)
        if pcliq > cutoff:
            dsi = super(SafePermafrostModel,self)._dsi_dpcice_frozen_unsaturated(pcliq, pcice, si)
            return dsi
        else:
            # in the spline region, form the spline
            a1,b1,c1,d1 = self._fit_spline(pcice, cutoff, si_cutoff)

            pcice2 = pcice + 10.
            si_cutoff2 = super(SafePermafrostModel,self)._si_frozen_unsaturated(cutoff, pcice2)
            a2,b2,c2,d2 = self._fit_spline(pcice2, cutoff, si_cutoff2)

            da = (a2 - a1)/10.
            db = (b2 - b1)/10.
            dc = (c2 - c1)/10.
            dd = (d2 - d1)/10.
            return ((da * pcliq + db)*pcliq + dc)*pcliq + dd

