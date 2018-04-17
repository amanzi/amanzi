# scipy's secant not an all that useful secant method...
from scipy.optimize.zeros import RootResults

def secant(func, x0, x1, args=(), tol=1.48e-8, maxiter=50, full_output=False, disp=True):
    # Secant method
    p0 = x0
    p1 = x1

    calls = 2
    q0 = func(*((p0,) + args))
    q1 = func(*((p1,) + args))

    converged = False
    for iteration in range(maxiter):
        if q1 == q0:
            if p1 != p0:
                msg = "Tolerance of %s reached" % (p1 - p0)
                warnings.warn(msg, RuntimeWarning)
            root = (p1 + p0)/2.0
            converged = True
            break
        else:
            p = p1 - q1*(p1 - p0)/(q1 - q0)
        if abs(p - p1) < tol:
            root = p
            converged = True
            break
        p0 = p1
        q0 = q1
        p1 = p
        q1 = func(*((p1,) + args))
        calls = calls + 1
        root = p

    if not converged and disp:
        msg = "Failed to converge after %d iterations, value is %s" % (maxiter, p)
        raise RuntimeError(msg)

    if full_output:
        if converged:
            res = RootResults(root, iteration, calls, 0)
        else:
            res = RootResults(root, iteration, calls, -1)
        return root, res
    return root


def regula_falsi(func, x0, x1, args=(), tol=1.48e-8, maxiter=50, full_output=False, disp=True):
    p0 = x0
    p1 = x1

    calls = 2
    q0 = func(*((p0,) + args))
    q1 = func(*((p1,) + args))

    converged = False
    for iteration in range(maxiter):
        if q1 == q0:
            if p1 != p0:
                msg = "Tolerance of %s reached" % (p1 - p0)
                warnings.warn(msg, RuntimeWarning)
            root = (p1 + p0)/2.0
            converged = True
            break
        else:
            p = (p0*q1 - p1*q0) / (q1 - q0)
        if abs(p - p1) < tol:
            root = p
            converged = True
            break

        q = func(*((p,) + args))
        calls = calls + 1

        if q*q0 < 0.:
            p1 = p
            q1 = q
        else:
            p0 = p
            q0 = q

        root = p

    if not converged and disp:
        msg = "Failed to converge after %d iterations, value is %s" % (maxiter, p)
        raise RuntimeError(msg)

    if full_output:
        if converged:
            res = RootResults(root, iteration, calls, 0)
        else:
            res = RootResults(root, iteration, calls, -1)
        return root, res
    return root
