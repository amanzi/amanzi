import eos
import capillary_pressure
import permafrost_model_safe_cubic

def compressible_poro(p,poro):
    if p > eos.p_atm:
        comp_poro = 1.e-9 * (p - eos.p_atm) + poro
    else:
        comp_poro = poro
    return comp_poro


def water_content(T,p,poro,pm=None):
    nl = eos.n_water(T,p)
    ni = eos.n_ice(T,p)
    ng = eos.n_gas(T,p)
    omega = eos.mol_frac_gas(T)

    if pm is None:
        pm = permafrost_model_safe_cubic.SafePermafrostModel()
    sg,sl,si = pm.saturations_Tp(T,p)

    cporo = compressible_poro(p,poro)
    return (ng*omega*sg + nl*sl + ni*si) * cporo

def energy(T,p,poro,pm=None):
    nl = eos.n_water(T,p)
    ni = eos.n_ice(T,p)
    ng = eos.n_gas(T,p)
    omega = eos.mol_frac_gas(T)

    if pm is None:
        pm = permafrost_model_safe_cubic.SafePermafrostModel()
    sg,sl,si = pm.saturations_Tp(T,p)

    ul = eos.u_water(T)
    ui = eos.u_ice(T)
    ug = eos.u_gas(T,omega)
    u_rock = eos.u_rock(T)
    rho_rock = 2170.

    cporo = compressible_poro(p,poro)
    return (ng*omega*sg*ug + nl*sl*ul + ni*si*ui) * cporo + u_rock*rho_rock*(1-cporo)


def ewc(T,p,poro,pm=None):
    if pm is None:
        pm = permafrost_model_safe_cubic.PermafrostModel()
    return energy(T,p,poro,pm), water_content(T,p,poro,pm)
