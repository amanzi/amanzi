import numpy as np
#from scipy import optimize

import wrm_vangenuchten
import capillary_pressure

class PermafrostModel(object):
    def __init__(self, wrm=None, dT=1., pc_ice=capillary_pressure.pc_ice, pc_liq=capillary_pressure.pc_liq):
        if wrm is None:
            wrm = wrm_vangenuchten.VanGenuchten(0.8, 1.5e-4)
        self.wrm = wrm
        self.pc_ice = pc_ice
        self.pc_liq = pc_liq
        self.dp = self.pc_ice(273.15-dT, 101325.)

    #-----------

    def saturations_pc(self, pcliq, pcice):
        if pcice <= max(pcliq, 0.):
            # unfrozen
            si = 0.
            sl = self.wrm.saturation(pcliq)
            sg = 1. - sl
        else:
            # frozen
            sstar_ice = self.wrm.saturation(pcice)
            sstar_liq = self.wrm.saturation(pcliq)
            sr = self.wrm._sr
            
            if pcliq <= 0.:
                # saturated            
                sl_sm = (sstar_liq - sr)*np.exp( -pcice/self.dp) + sr
                sl = max(sstar_ice, sl_sm)
                si = 1. - sl
                sg = 0.
                
            else:
                sl_sm = (sstar_liq - sr)*np.exp( (pcliq - pcice) / self.dp) + sr
                sl = max(sstar_ice, sl_sm)
                si = 1. - sl / sstar_liq
                sg = 1. - sl - si
            
        return sg, sl, si


    def saturations_Tp(self, T, p):
        pcliq = self.pc_liq(T,p)
        pcice = self.pc_ice(T,p)
        return self.saturations_pc(pcliq,pcice)

    def dsaturations_dpc_liq(self, pcliq, pcice):
        dsats = np.zeros((3,),'d')
        if (pcice <= max(pcliq,0.)):
            dsats[2] = 0.
            dsats[1] = self.wrm.d_saturation(pcliq)
            dsats[0] = - dsats[1]

        else:
            sstar_ice = self.wrm.saturation(pcice)
            sstar_liq = self.wrm.saturation(pcliq)
            sr = self.wrm._sr

            if (pcliq <= 0.):
                sl_sm = (sstar_liq - sr)*np.exp(-pcice/self.dp) + sr
                if (sstar_ice > sl_sm):
                    dsats[1] = 0.
                    dsats[2] = 0.
                    dsats[0] = 0.
                else:
                    sstarprime_liq = self.wrm.d_saturation(pcliq)
                    dsats[1] =  sstarprime_liq * np.exp(-pcice/self.dp)
                    dsats[2] = -dsats[1]
                    dsats[0] = 0.

            else:
                sl_sm = (sstar_liq - sr)*np.exp( (pcliq - pcice) / self.dp) + sr
                if (sstar_ice > sl_sm):
                    dsats[1] = 0.
                    dsats[2] = sstar_ice / pow(sstar_liq,2) * self.wrm.d_saturation(pcliq)
                    dsats[0] = - dsats[2]

                else:
                    sstarprime_liq = self.wrm.d_saturation(pcliq)
                    dsats[1] =  sstarprime_liq * np.exp( (pcliq - pcice)/self.dp ) + (sstar_liq - sr) * np.exp((pcliq - pcice)/self.dp) / self.dp
                    dsats[2] = -dsats[1] / sstar_liq + sl_sm / pow(sstar_liq, 2) * sstarprime_liq
                    dsats[0] = -dsats[1] - dsats[2]
        return dsats[0], dsats[1], dsats[2]

    def dsaturations_dpc_ice(self, pc_liq, pc_ice):
        dsats = np.zeros((3,),'d')
        if (pc_ice <= max(pc_liq,0.)):
            dsats[0] = 0.
            dsats[1] = 0.
            dsats[2] = 0.
            
        else:
            sstar_ice = self.wrm.saturation(pc_ice)
            sstar_liq = self.wrm.saturation(pc_liq)
            sr = self.wrm._sr

            if (pc_liq <= 0.):
                sl_sm = (sstar_liq - sr)*np.exp(-pc_ice/self.dp) + sr
                if (sstar_ice > sl_sm):
                    dsats[1] = self.wrm.d_saturation(pc_ice)
                else:
                    dsats[1] = -(sstar_liq - sr) * np.exp(-pc_ice/self.dp) / self.dp

                dsats[2] = -dsats[1]
                dsats[0] = 0.

            else:
                sl_sm = (sstar_liq - sr)*np.exp( (pc_liq - pc_ice)/self.dp) + sr
                if (sstar_ice > sl_sm):
                    dsats[1] = self.wrm.d_saturation(pc_ice)
                else:
                    dsats[1] = -(sstar_liq - sr) * np.exp( (pc_liq - pc_ice)/self.dp) / self.dp

                dsats[2] = -dsats[1] / sstar_liq
                dsats[0] = -dsats[1] - dsats[2]
        return dsats[0], dsats[1], dsats[2]
