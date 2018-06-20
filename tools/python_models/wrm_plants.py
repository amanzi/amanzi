"""Implementation of plant hydraulics WRMs.

Unsmoothed version"""

import numpy as np
import scipy.optimize
import os,sys
#import ridders as rd

p_atm = 101325.


class WRMPlantHydraulics(object):
    """Class that implements WRM for plant hydraulics

    Parameters:
      psi_0     | water potential at saturation [MPa]
                |   default = -0.08
      s_r       | residual saturation [-]
      s_tlp     | saturation at turgor loss [-]
      s_ft      | saturation at full turgor [-]
      pi_0      | osmotic potential at full turgor [MPa]
                |   default = -1.93116, corresponds to
                |       sapwood with wood density of 0.65
      eps       | bulk elastic modulus [MPa]
                |   default = 15.90459885634534, corresponds
                |       sapwood with wood density of 0.65
    """

    def __init__(self, s_r, s_tlp,
                 eps=15.90459885634534,
                 psi_0=-0.08,
                 pi_0=-1.93116,
                 psi_cap=-0.39
        ):
        assert 0.0 <= s_r <= s_tlp <= 1.0
        # traits, from wood density
        self.pi_0 = pi_0
        self.eps = eps
        
        # inflection points in s
        # -- left endpoint
        self.s_r = s_r  # (psi = -inf)

        # -- right endpoint
        self.psi_0 = psi_0 # (s = 1)

        # -- loss of turgor
        self.s_tlp = s_tlp

        # -- full-turgor is not a true inflection point -- instead we
        # fit a line to the capillary data and take the min of that
        # line and the turgor curve (which goes from turgor loss point
        # to full turgor).
        star = 1. - abs(pi_0) / eps
        self.s_ft = (s_tlp - (1-star)*s_r) / star
        
        # -- psi(s = 1) = psi_0 and psi(s = s_cap) = psi_cap are the
        # fits to data for the capillary line.  Note that the cap
        # point is NOT an inflection point as this line crosses the
        # turgor curve and we take the min of those two
        self.s_cap = 1.0 - 0.61*(1-self.s_tlp) 
        self.psi_cap = psi_cap

        # -- and the psi at tlp, ft
        s_tlp_e = (s_tlp - s_r) / (self.s_ft - s_r)
        self.psi_tlp = - abs(pi_0) / s_tlp_e

        # -- now, the slope of the capillary curve line fit to data
        self.m_cap = (self.psi_0 - self.psi_cap)/(1 - self.s_cap)
        self.m_cap_star = self.m_cap * (self.s_ft - self.s_r)

        # -- calculate s_cap_ft_trans, which is the true inflection
        # point, at which the turgor curve equals the capillary line.
        assert self.s_ft > self.s_cap
        def func(s):
            return self.potentialLinear(s) - (self.potentialP(s) + self.potentialSol(s))
        self.s_cap_ft_trans = scipy.optimize.brentq(func, self.s_tlp, self.s_ft)
        self.psi_cap_ft_trans = self.potentialLinear(self.s_cap_ft_trans)

        assert 0.0 <= self.s_r <= self.s_tlp <= self.s_cap <= self.s_ft <= 1.0
        assert 0.0 <= self.s_r <= self.s_tlp <= self.s_cap_ft_trans <= self.s_ft <= 1.0
        

    def capillaryPressure(self, s):
        return - self.potential(s)*1e6        
        
    def potential(self, s):
        assert 1.0 >= s >= self.s_r
        #print "Ss:", s, self.s_ft, self.s_tlp, self.s_r
        if s > self.s_ft:
            psi = self.potentialLinear(s)
            #print "Cap:", s, psi
        elif s > self.s_tlp:
            psi1 = self.potentialSol(s) + self.potentialP(s)
            psi2 = self.potentialLinear(s)
            psi = min(psi1,psi2)
            #print "Turgor:", s, psi, psi1, psi2
        else:
            psi = self.potentialSol(s)
            #print "Embolysm:", s, psi, self.potentialSol(s)
        return psi

    def potentialLinear(self, s):
        return self.psi_0 - self.m_cap*(1-s)
    
    def potentialSol(self, s):
        sstar = (s - self.s_r) / (self.s_ft - self.s_r)
        return -abs(self.pi_0) / sstar

    def potentialP(self, s):
        sstar = (s - self.s_r) / (self.s_ft - self.s_r)
        return abs(self.pi_0) - self.eps*(1-sstar)

    def saturation(self, pc):
        psi = - pc*1e-6

        if psi > self.psi_0:
            # saturated
            s1e = (1.0 - self.s_r) / (self.s_ft - self.s_r)
            se = s1e

        elif psi > self.psi_cap_ft_trans:
            # linear branch
            s1e = (1.0 - self.s_r) / (self.s_ft - self.s_r)
            se = s1e - (self.psi_0 - psi) / self.m_cap_star

        elif psi > self.psi_tlp:
            # turgor branch
            b = abs(self.pi_0) - psi - self.eps
            se = (-b + np.sqrt(b**2 + 4*self.eps*abs(self.pi_0))) / (2*self.eps)
            
        else:
            # embolism, turgor lost
            se = - abs(self.pi_0)/psi

        s = se * (self.s_ft - self.s_r) + self.s_r
        return s
        
def plot(wrm, fmt, ax):
    from matplotlib import pyplot as plt
    sats = np.linspace(wrm.s_r + 1.e-6, 1.0, 1000)[1:]
    pcs = np.array([wrm.capillaryPressure(s) for s in sats])

    # pcs_linear = p_atm - 1.e6*np.array([wrm.potentialLinear(s) for s in sats])
    # pcs_p_plus_sol = p_atm - 1.e6*(np.array([wrm.potentialP(s) for s in sats])+np.array([wrm.potentialSol(s) for s in sats]))
    # pcs_sol = p_atm - 1.e6*np.array([wrm.potentialSol(s) for s in sats])

    sats_inv = np.array([wrm.saturation(pc) for pc in pcs])

    ax.semilogy(sats, pcs, fmt)
    ax.semilogy(sats_inv, pcs, fmt)
    # ax.plot(sats, pcs_linear, 'b')
    # ax.plot(sats, pcs_p_plus_sol, 'r')
    # ax.plot(sats, pcs_sol, 'g')
    # ax.plot(sats, pcs, 'k--')

if __name__ == "__main__":
    from matplotlib import pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(111)

    # typical values: s_r = 0.4, s_tlp = 0.85 (both stem and leaf)
    wrm = WRMPlantHydraulics(0.4, 0.85)
    plot(wrm,'b',ax)
    
    import wrm_vangenuchten
    wrmvg = wrm_vangenuchten.VanGenuchten(m=.3,alpha=1.e-3, sr=0.1)
    plot(wrmvg,'g',ax)
    wrmvg = wrm_vangenuchten.VanGenuchten(n=1.27,alpha=5.e-3, sr=0.13)
    plot(wrmvg,'r',ax)
    plt.show()
