import numpy

class SteadyLinear(object):
    """Solves the system:

        \div \frac{\rho}{\mu} K \grad (p + \rho g z) = 0  on the domain [x_0, x_1] \cross [z_0,z_1]

    Boundary conditions are given by:
       \left. q \cdot \hat{n} \right|_{z=z_0,z_1} = 0
       \right. q \cdot \hat{n} \right|_{x=x_0} = q_left  [kg / m^2 s]
       \right. h + (z_1 - z) \left|_{x=x_1} = head_right  [ m ]

    Parameters are in units of:
       \rho : density, [kg/m^3]
       \mu  : viscosity, [kg / m s^2]
       K    : absolute permeability, [ m^2 ]
       g    : gravity, used in converting head to pressure, [ m / s^2 ]
    """

    def __init__(self, params=None):
        if params is None:
            params = dict()

        params.setdefault("x0",0)
        params.setdefault("x1",1)
        params.setdefault("z0",0)
        params.setdefault("z1",1)
        params.setdefault("g",9.81)
        params.setdefault("K",1.e-12)
        params.setdefault("rho",1000)
        params.setdefault("mu",1.e-3)
        params.setdefault("q_left",1.95e-2)
        params.setdefault("h_right",120)
        self.__dict__.update(params)

    def run(self, coords):
        """Evaluates the solution at (x,z)-coordinates in the [n_points, 2] array."""
        pres = numpy.zeros((len(coords),), 'd')

        p_right = self.rho*self.g*self.h_right + 101325.0

        grad_p = self.q_left / self.K / self.rho * self.mu

        for i,(x,z) in enumerate(coords):
            xi = (x - self.x1) / (self.x1 - self.x0)
            p_upper_boundary = p_right + xi*grad_p
            pres[i] = p_upper_boundary + self.rho*self.g*(self.z1 - z)

        return pres

def createFromXML(filename):
    # grab params from input file
    params = dict()

    import amanzi_xml.utils.io
    xml = amanzi_xml.utils.io.fromFile(filename)
    import amanzi_xml.utils.search as search
    # this line does not work... yet
    #    params["mu"] = search.find(xml, "..../viscosity").get("value")

    # instantiate the class
    return SteadyLinear(params)



if __name__ == "__main__":
    sl = SteadyLinear()
    x = ...
    c1 = numpy.array()
    # run with generic values

    # plot
    plt.show()
