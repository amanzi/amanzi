import numpy as np
import eos

def pc_liq(T,p=None):
    if p is None:
        return pc_liq(p,T)
    return eos.p_atm - p

def pc_ice(T,p):
    if T >= 273.15:
        return 0
    else:
        Mv = 0.0180153
        gamma = 2.7/33.1 * 3.34e5 * Mv
        return gamma * (273.15 - T)/273.15 * eos.n_water(T,max(p,eos.p_atm))

def pc_ice_regularized(halfwidth):
    def pc_ice(T,p):
        Mv = 0.0180153
        gamma = 2.7/33.1 * 3.34e5 * Mv        
        alpha = gamma * eos.n_water(T,max(p,eos.p_atm)) / 273.15
        a = alpha * halfwidth / 4.
        b = -alpha / 2.
        c = alpha / (4*halfwidth)
        x = T - 273.15
        
        if x < -halfwidth:
            pc = -alpha * x
        elif x > halfwidth:
            pc = 0.
        else:
            pc = a + b*x + c*x*x
        return pc
    return pc_ice
