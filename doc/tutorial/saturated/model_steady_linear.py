import numpy
import matplotlib.pylab as plt

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
        params.setdefault("g",9.80665)
        params.setdefault("K",1.e-12)
        params.setdefault("rho",1000)
        params.setdefault("mu",1.e-3)
        params.setdefault("q_Upstream",1.95e-2)
        params.setdefault("h_Downstream",120)
        self.__dict__.update(params)

    def run(self, coords):
        """Evaluates the solution at (x,z)-coordinates in the [n_points, 2] array."""
        
        pres = numpy.zeros((len(coords),), 'd')

        p_right = self.rho*self.g*self.h_Downstream + 101325.0

        grad_p = - self.q_Upstream / self.K / self.rho * self.mu

        for i,(x,z) in enumerate(coords):
            xi = self.x1 - x
            p_lower_boundary = p_right + xi*grad_p
            pres[i] = p_lower_boundary - self.rho*self.g*z

        return pres

def createFromXML(filename):

    # grab params from input file
    params = dict()

    import amanzi_xml.utils.io
    xml = amanzi_xml.utils.io.fromFile(filename)
    import amanzi_xml.utils.search as search
   
    params["x0"] = search.getElementByPath(xml, "/Main/Mesh/Unstructured/Generate Mesh/Uniform Structured/Domain Low Corner").value[0]
    params["z0"] = search.getElementByPath(xml, "/Main/Mesh/Unstructured/Generate Mesh/Uniform Structured/Domain Low Corner").value[2]
    params["x1"] = search.getElementByPath(xml, "/Main/Mesh/Unstructured/Generate Mesh/Uniform Structured/Domain High Corner").value[0]
    params["z1"] = search.getElementByPath(xml, "/Main/Mesh/Unstructured/Generate Mesh/Uniform Structured/Domain High Corner").value[2] 
    params["K"] = search.getElementByPath(xml, "/Main/Material Properties/Soil/Intrinsic Permeability: Uniform/Value").value
    params["mu"]= search.getElementByPath(xml, "/Main/Phase Definitions/Aqueous/Phase Properties/Viscosity: Uniform/Viscosity").value
    params["rho"]= search.getElementByPath(xml, "/Main/Phase Definitions/Aqueous/Phase Properties/Density: Uniform/Density").value
    params["q_Upstream"] = -search.getElementByPath(xml, "/Main/Boundary Conditions/Upstream BC/BC: Flux/Inward Mass Flux").value[0]
    params["h_Downstream"] = search.getElementByPath(xml, "/Main/Boundary Conditions/Downstream BC/BC: Hydrostatic/Water Table Height").value[0] 
    params.setdefault("g",9.80665)
   

    # instantiate the class
    return SteadyLinear(params)



if __name__ == "__main__":

    sl = SteadyLinear()

    x = numpy.linspace(0,1,11)

    coords = numpy.zeros((11,2))
    coords[:,0]=x

    coords[:,1]=0.3
    p1 = sl.run(coords)

    coords[:,1]=0.7
    p2 = sl.run(coords)

    plt.plot(x,p1)
    plt.plot(x,p2)

    # plot
    plt.show()
