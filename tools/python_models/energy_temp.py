import numpy as np
import wrm_vangenuchten
from matplotlib import pyplot as plt

from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset


_p_atm = 101325.0

def u_water(T):
    return 76. * (T - 273.15)

def u_ice(T):
    dT = T - 273.15
    return -6007.86 + 37.7841*dT + 0.0659661*dT*dT

def u_gas(T,omega):
    dT = T - 273.15
    return (1.0 + 0.622*omega)*13.0*dT + omega*40650.0

def mol_frac_gas(T):
    ka0 = 16.635764
    ka = -6096.9385
    kb = -2.7111933e-2
    kc = 1.673952e-5
    kd = 2.433502
    return 100.0*np.exp(ka0 + ka/T + (kb + kc*T)*T + kd*np.log(T)) / _p_atm

def n_water(T,p=_p_atm):
    Mv = 0.0180153
    # ka = 999.915
    # kb = 0.0416516
    # kc = 0.0100836
    # kd = 0.000206355
    # kT0 = 273.15
    # kalpha = 5.0e-10
    # kp0 = 1.0e5

    # dT = T - kT0
    # rho1bar = ka + (kb + (kc + kd*dT)*dT)*dT
    # rho = rho1bar * (1.0 + kalpha*(p - kp0))
    rho = 1000.
    return rho / Mv

def n_gas(T,p=_p_atm):
    return p / 8.3144621 / T

def n_ice(T,p=_p_atm):
    Mv = 0.0180153
    ka = 916.724
    kb = -0.147143
    kc = -0.000238095
    kT0 = 273.15
    kalpha = 1.0e-10
    kp0 = 1.0e5

    dT = T - kT0
    rho1bar = ka + (kb + kc*dT)*dT
    rho = rho1bar * (1.0 + kalpha*(p - kp0))
    return rho / Mv

def saturations_A(T):
    wrm = wrm_vangenuchten.VanGenuchten(0.8, 1.5e-4)
    Mv = 0.0180153
    gamma = 72.7/33.1 * 3.34e5 * Mv
    A = 1.0 / wrm.saturation(gamma * (273.15 - T)/273.15 * n_water(T))
    return A

def saturations_A2(T):
    wrm = wrm_vangenuchten.VanGenuchten(0.8, 1.5e-4)
    Mv = 0.0180153
    gamma = 72.7/33.1 * 3.34e3 * Mv
    A = 1.0 / wrm.saturation(gamma * (273.15 - T)/273.15 * n_water(T))
    return A

def saturations_A3(T):
    wrm = wrm_vangenuchten.VanGenuchten(0.8, 1.5e-4)
    Mv = 0.0180153
    gamma = 72.7/33.1 * 3.34e4 * Mv
    A = 1.0 / wrm.saturation(gamma * (273.15 - T)/273.15 * n_water(T))
    return A

def saturations_AB(T,p,afunc=saturations_A):
    wrm = wrm_vangenuchten.VanGenuchten(0.8, 1.5e-4)
    B = 1.0 / wrm.saturation(_p_atm - p)
    A = afunc(T)
    return A,B

def saturations(T,p,afunc=saturations_A):
    A,B = saturations_AB(T,p,afunc)
    print "A,B = ", A,B
    s_l = 1.0/(A+B-1.0)
    s_g = s_l*(B-1.0)
    s_i = s_l*(A-1.0)
    return s_l, s_i, s_g

def saturation_with_wc(T,wc, afunc=saturations_A):
    nl = n_water(T)
    ni = n_ice(T)
    ng = n_gas(T)
    omega = mol_frac_gas(T)

    A = afunc(T)
    wrm = wrm_vangenuchten.VanGenuchten(0.8, 1.5e-4)
    B = (nl + ni*(A-1) - omega*ng - (A-1)*wc) / (wc - omega*ng)
    pc = wrm.capillaryPressure(1.0/B)
    p = _p_atm - pc

    s_l = 1.0/(A+B-1.0)
    s_g = s_l*(B-1.0)
    s_i = s_l*(A-1.0)
    return (s_l, s_i, s_g), p

def water_of_T_sg(T,sg,afunc=saturations_A):
    nl = n_water(T)
    ni = n_ice(T)
    ng = n_gas(T)
    omega = mol_frac_gas(T)

    A = afunc(T)
    B = (1. - sg + sg*A) / (1. - sg)

    wrm = wrm_vangenuchten.VanGenuchten(0.8, 1.5e-4)
    pc = wrm.capillaryPressure(1.0/B)
    p = _p_atm - pc

    s_l = 1.0/(A+B-1.0)
    s_i = s_l*(A-1.0)
    return (s_l, s_i, s_g), p


def water_content(T,p,afunc=saturations_A):
    nl = n_water(T)
    ni = n_ice(T)
    ng = n_gas(T)
    omega = mol_frac_gas(T)
    sl,si,sg = saturations(T,p,afunc)

    wc = (sl*nl + si*ni + sg*ng*omega)
    return wc


def internal_energy(T,p,afunc=saturations_A):
    nl = n_water(T)
    ni = n_ice(T)
    ng = n_gas(T)
    omega = mol_frac_gas(T)
    sl,si,sg = saturations(T,p,afunc)
    ul = u_water(T)
    ui = u_ice(T)
    ug = u_gas(T,omega)

    return (ul*sl*nl + ui*si*ni + ug*sg*omega*ng)

def internal_energy_per_mol(T,p,afunc=saturations_A):
    return internal_energy(T,p,afunc) / water_content(T,p,afunc)

def plot_A_of_T():
    Ts1 = np.arange(273.15-10, 273.15-5, 1.)
    Ts2 = np.arange(273.15-5, 273.15-1, .1)
    Ts3 = np.arange(273.15-1, 273.15+1, .0001)
    Ts4 = np.arange(273.15+1, 273.15+5, .1)
    Ts5 = np.arange(273.15+5, 273.15+50, 1.)

    Ts = np.concatenate((Ts1, Ts2, Ts3, Ts4, Ts5))

    # water content at T=273.65, p=85000.
    A = []
    A2 = []
    A3 = []
    for T in Ts:
        A.append(saturations_A(T))
        A2.append(saturations_A2(T))
        A3.append(saturations_A3(T))

    plt.plot(Ts,np.array(A), 'b-x')
    plt.plot(Ts,np.array(A2), 'r-x')
    plt.plot(Ts,np.array(A3), 'g-x')
    plt.xlabel("temperature [K]")
    plt.ylabel("A (saturation parameter)")
    plt.gca().set_ylim(0,100)
    plt.show()


def get_reasonable_Ts():
    Ts1 = np.arange(273.15-40, 273.15-5, 1.)
    Ts2 = np.arange(273.15-5, 273.15-1, .1)
    Ts3 = np.arange(273.15-1, 273.15+1, .0001)
    Ts4 = np.arange(273.15+1, 273.15+5, .1)
    Ts5 = np.arange(273.15+5, 273.15+30, 1.)
    Ts = np.concatenate((Ts1, Ts2, Ts3, Ts4, Ts5))
    return Ts

def plot_e_wc():
    """Plots e and wc at a range of saturations and temperatures"""
    Ts = get_reasonable_Ts()

    sgs = [0., 0.5, 0.9, 0.99, 0.999]
    colors = ['red', 'orange', 'green', 'blue', 'violet']

    fig = plt.figure(figsize=(7,7))
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)
    ax1.ylabel("water content [mol/m^3]")
    ax2.ylabel("internal energy [J/m^3]")

    for sg,color in zip(sgs,colors):
        ps = []
        wcs = []
        ies = []

        for T in Ts:
            (sl,si,sg2),p = water_of_T_sg(T,sg)
            ps.append(p)
            wcs.append(water_content(T,p))
            ies.append(internal_energy(T,p))




def plot_e_of_T():
    Ts = get_reasonable_Ts()

    # water content at T=273.65, p=85000.
    wc = water_content(273.65, 85000., saturations_A)
    ie = []
    for T in Ts:
        # for each T, determine the p required to maintain same mass
        (sl,si,sg),p = saturation_with_wc(T,wc, saturations_A)
        ie.append(internal_energy_per_mol(T,p, saturations_A))

    plt.plot(Ts - 273.15,np.array(ie), 'b')

    # water content at T=273.65, p=100000.
    wc = water_content(273.65, 97000.)
    ie = []
    for T in Ts:
        # for each T, determine the p required to maintain same mass
        (sl,si,sg),p = saturation_with_wc(T,wc)
        ie.append(internal_energy_per_mol(T,p))
    plt.plot(Ts,np.array(ie), 'r-x')

    # # water content at T=273.65, p=65000.
    # wc = water_content(273.65, 65000.)
    # ie = []
    # for T in Ts:
    #     # for each T, determine the p required to maintain same mass
    #     (sl,si,sg),p = saturation_with_wc(T,wc)
    #     ie.append(internal_energy_per_mol(T,p))
    # plt.plot(Ts,np.array(ie), 'g-x')

    # # water content at T=273.65, p=95000.
    # wc = water_content(273.65, 95000.)
    # ie = []
    # for T in Ts:
    #     # for each T, determine the p required to maintain same mass
    #     (sl,si,sg),p = saturation_with_wc(T,wc)
    #     ie.append(internal_energy_per_mol(T,p))
    # plt.plot(Ts,np.array(ie), 'k-x')

    plt.xlabel("temperature [C]")
    plt.ylabel("energy density (water) [J/mol]")
    plt.show()

def nice_plot_e_of_T():
    Ts1 = np.arange(273.15-10, 273.15-5, 1.)
    Ts2 = np.arange(273.15-3, 273.15-1, .1)
    Ts3 = np.arange(273.15-1, 273.15+1, .0001)
    Ts4 = np.arange(273.15+1, 273.15+5, .1)
    Ts5 = np.arange(273.15+5, 273.15+50, 1.)

    Ts = np.concatenate((Ts1, Ts2, Ts3, Ts4, Ts5))

    # water content at T=273.65, p=85000.
    wc = water_content(273.65, 85000., saturations_A)
    ie = []

    for T in Ts:
        # for each T, determine the p required to maintain same mass
        (sl,si,sg),p = saturation_with_wc(T,wc, saturations_A)
        ie.append(internal_energy_per_mol(T,p, saturations_A))

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(Ts - 273.15,np.array(ie), 'b')
    ax.set_xlim(-10,20)
    ax.set_ylim(-8000,2000)
    plt.xlabel("temperature [C]")
    plt.ylabel("energy density (water) [J/mol]")

    axins = inset_axes(ax, 3, 3, loc=4)
    axins.plot(Ts - 273.15,np.array(ie), 'b')
    axins.set_xlim(-.01,.01)
    axins.set_ylim(-7000,500)
    axins.set_xticks([-0.01, 0, 0.01])
    axins.tick_params(labeltop=1, labelbottom=0)
    axins.set_yticks([-6000,-4000,-2000,0])

    mark_inset(ax, axins, loc1=1, loc2=3, fc="none", ec="0.5")
    plt.show()

def nice_plot_e_of_T_simple():
    Ts1 = np.arange(273.15-10, 273.15-5, 1.)
    Ts2 = np.arange(273.15-3, 273.15-1, .1)
    Ts3 = np.arange(273.15-1, 273.15+1, .0001)
    Ts4 = np.arange(273.15+1, 273.15+5, .1)
    Ts5 = np.arange(273.15+5, 273.15+50, 1.)

    Ts = np.concatenate((Ts1, Ts2, Ts3, Ts4, Ts5))

    # water content at T=273.65, p=85000.
    wc = water_content(273.65, 85000., saturations_A)
    ie = []

    for T in Ts:
        # for each T, determine the p required to maintain same mass
        (sl,si,sg),p = saturation_with_wc(T,wc, saturations_A)
        ie.append(internal_energy_per_mol(T,p, saturations_A))

    fig = plt.figure(figsize=(3,6))
    ax = fig.add_subplot(111)
    ax.plot(Ts - 273.15,np.array(ie), 'b')
    ax.set_xlim(-5,5)
    ax.set_ylim(-6200,1000)
    plt.xlabel("temperature [C]")
    plt.ylabel("energy density (water) [J/mol]")

    return ax


def calc_T_of_e_at_p():
    def func(e,p):
        def e_of_T(T):
            return internal_energy_per_mol(T,p) - e
        return e_of_T

    # bisection?
    from scipy.optimize import bisect
    def do_bisect(p,e,Tmin=220,Tmax=320,verbose=True):
        T, result = bisect(func(e,p), Tmin, Tmax, full_output=True,disp=False)
        if verbose:
            print "Bisection:  Energy = %g J/mol, pressure = %g Pa"%(e,p)
            print "  T = %g"%T
            if result.converged:
                print "  converged in %d function calls"%result.function_calls
            else:
                print "DID NOT CONVERGE"
                print "------------------------------------------------------------"
                print ""
        return result

    # secant method?
    from secant import secant
    def do_secant(p,e,T0=220.,T1=320., verbose=True):
        T, result = secant(func(e,p), T0, T1, maxiter=1000, full_output=True,disp=False)
        if verbose:
            print "Secant method:  Energy = %g J/mol, pressure = %g Pa"%(e,p)
            print "  T = %g"%T
            if result.converged:
                print "  converged in %d function calls"%result.function_calls
            else:
                print "DID NOT CONVERGE"
                print "------------------------------------------------------------"
                print ""
        return result

    # regula falsi method?
    from secant import regula_falsi
    def do_rf(p,e,T0=220.,T1=320., verbose=True):
        T, result = regula_falsi(func(e,p), T0, T1, maxiter=1000, full_output=True,disp=False)
        if verbose:
            print "Regula falsi:  Energy = %g J/mol, pressure = %g Pa"%(e,p)
            print "  T = %g"%T
            if result.converged:
                print "  converged in %d function calls"%result.function_calls
            else:
                print "DID NOT CONVERGE"
                print "------------------------------------------------------------"
                print ""
        return result

    # brent method?
    from scipy.optimize import brentq
    def do_brent(p,e,T0=220.,T1=320., verbose=True):
        T, result = brentq(func(e,p), T0, T1, maxiter=1000, full_output=True,disp=False)
        if verbose:
            print "Brent:  Energy = %g J/mol, pressure = %g Pa"%(e,p)
            print "  T = %g"%T
            if result.converged:
                print "  converged in %d function calls"%result.function_calls
            else:
                print "DID NOT CONVERGE"
                print "------------------------------------------------------------"
                print ""
        return result

    # hybrid approaches
    def do_hybrid_secant(p,e, verbose=True):
        if verbose: print "Hybrid approach:  Energy = %g J/mol, pressure = %g Pa"%(e,p)
        e_above = internal_energy_per_mol(273.15,p)
        e_below = internal_energy_per_mol(273.0,p)
        if e > e_above:
            result = do_secant(p, e, 273.15, 320, verbose)
        elif e < e_below:
            result = do_secant(p, e, 220, 273.0, verbose)
        else:
            result = do_secant(p,e, 273.0, 273.15, verbose)
        return result

    def do_hybrid_bisect(p,e, verbose=True):
        if verbose: print "Hybrid approach:  Energy = %g J/mol, pressure = %g Pa"%(e,p)
        e_above = internal_energy_per_mol(273.15,p)
        e_below = internal_energy_per_mol(273.0,p)
        if e > e_above:
            result = do_secant(p, e, 273.15, 320, verbose)
        elif e < e_below:
            result = do_secant(p, e, 220, 273.0, verbose)
        else:
            result = do_bisect(p,e, 273.0, 273.15, verbose)
        return result

    def do_hybrid_rf(p,e, verbose=True):
        if verbose: print "Hybrid approach:  Energy = %g J/mol, pressure = %g Pa"%(e,p)
        e_above = internal_energy_per_mol(273.15,p)
        e_below = internal_energy_per_mol(273.0,p)
        if e > e_above:
            result = do_secant(p, e, 273.15, 320, verbose)
        elif e < e_below:
            result = do_secant(p, e, 220, 273.0, verbose)
        else:
            result = do_rf(p,e, 273.0, 273.15, verbose)
        return result

    def do_hybrid_brent(p,e, verbose=True):
        if verbose: print "Hybrid approach:  Energy = %g J/mol, pressure = %g Pa"%(e,p)
        e_above = internal_energy_per_mol(273.15,p)
        e_below = internal_energy_per_mol(273.0,p)
        if e > e_above:
            result = do_secant(p, e, 273.15, 320, verbose)
        elif e < e_below:
            result = do_secant(p, e, 220, 273.0, verbose)
        else:
            result = do_brent(p,e, 273.0, 273.15, verbose)
        return result


    def test_all(p,e):
        print "Comparison of Algorithms"
        print ""
        res = do_bisect(p,e,verbose=False)
        print "  Bisection:", res.converged, res.root, res.function_calls

        res = do_secant(p,e,verbose=False)
        print "  Secant:", res.converged, res.root, res.function_calls

        res = do_rf(p,e,verbose=False)
        print "  Regula Falsi:", res.converged, res.root, res.function_calls

        res = do_brent(p,e,verbose=False)
        print "  Brent:", res.converged, res.root, res.function_calls

        res = do_hybrid_bisect(p,e,verbose=False)
        print "  Hybrid Bisect:", res.converged, res.root, res.function_calls

        res = do_hybrid_secant(p,e,verbose=False)
        print "  Hybrid Secant:", res.converged, res.root, res.function_calls

        res = do_hybrid_rf(p,e,verbose=False)
        print "  Hybrid RF:", res.converged, res.root, res.function_calls

        res = do_hybrid_brent(p,e,verbose=False)
        print "  Hybrid Brent:", res.converged, res.root, res.function_calls

        print "============================================================"
        print ""
        print ""

    test_all(85000,2000)
    test_all(85000,1.344)
    test_all(85000,-2)
    test_all(85000,-3000)
    test_all(85000,-5800)
    test_all(85000,-6018)
    test_all(85000,-6300)
