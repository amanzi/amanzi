import numpy as np

p_atm = 101325.

def u_water(T):
    return 76. * (T - 273.15)

def u_ice(T):
    dT = T - 273.15
    return -6007.86 + 37.7841*dT + 0.0659661*dT*dT

def u_gas(T,omega):
    dT = T - 273.15
    return (1.0 + 0.622*omega)*13.0*dT + omega*40650.0

def u_rock(T):
    return 620. * (T - 273.15)

def mol_frac_gas(T):
    ka0 = 16.635764
    ka = -6096.9385
    kb = -2.7111933e-2
    kc = 1.673952e-5
    kd = 2.433502
    return 100.0*np.exp(ka0 + ka/T + (kb + kc*T)*T + kd*np.log(T)) / p_atm

def n_water(T,p):
    p = max(p,p_atm)
    Mv = 0.0180153
    ka = 999.915
    kb = 0.0416516
    kc = -0.0100836
    kd = 0.000206355
    kT0 = 273.15
    kalpha = 5.0e-10
    kp0 = 1.0e5

    dT = T - kT0
    rho1bar = ka + (kb + (kc + kd*dT)*dT)*dT
    rho = rho1bar * (1.0 + kalpha*(p - kp0))
    return rho / Mv

def n_gas(T,p):
    p = max(p,p_atm)
    return p / 8.3144621 / T

def n_ice(T,p):
    p = max(p,p_atm)
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

def visc_water(T):
    kav1_ = 998.333
    kbv1_ = -8.1855
    kcv1_ = 0.00585
    kbv2_ = 1.3272
    kcv2_ = -0.001053
    kT1_ = 293.15

    dT = kT1_ - T

    if (T < kT1_):
        A = kav1_ + (kbv1_ + kcv1_*dT)*dT
        xi = 1301.0 * (1.0/A - 1.0/kav1_)
    else:
        A = (kbv2_ + kcv2_*dT)*dT
        xi = A/(T - 168.15)

    return 0.001 * pow(10.0, xi)





