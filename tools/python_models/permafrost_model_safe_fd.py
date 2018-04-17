import numpy as np

import wrm_vangenuchten
import capillary_pressure
import permafrost_model

class SafePermafrostModel(permafrost_model.PermafrostModel):
    def _dsi_dpcliq_frozen_unsaturated(self, pcliq, pcice, si):
        # derivatives fail when (1-si)*sstar + si == 1
        sstar = self.wrm.saturation(pcliq)
        if (1 - si)*sstar + si > (1-1.e-16):
            delta = pcliq / 100.
            si2 = self._si_frozen_unsaturated(pcliq + delta, pcice)
            return (si2 - si) / delta
        else:
            return super(SafePermafrostModel,self)._dsi_dpcliq_frozen_unsaturated(pcliq, pcice, si)

    def _dsi_dpcice_frozen_unsaturated(self, pcliq, pcice, si):
        # derivatives fail when (1-si)*sstar + si == 1
        sstar = self.wrm.saturation(pcliq)
        if (1 - si)*sstar + si > (1-1.e-16):
            delta = pcice / 100.
            si2 = self._si_frozen_unsaturated(pcliq, pcice + delta)
            return (si2 - si) / delta
        else:
            return super(SafePermafrostModel,self)._dsi_dpcliq_frozen_unsaturated(pcliq, pcice, si)
