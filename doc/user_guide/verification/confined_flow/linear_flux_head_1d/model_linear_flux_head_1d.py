import numpy
import matplotlib.pylab as plt

class LinearFluxHead(object):

    """
    Solves the system:

        \div \frac{\rho}{\mu} k \grad (p + \rho g z) = 0  on the domain [x_0, x_1] \cross [z_0,z_1]

    Boundary conditions are given by:

       U(x_0,z,t) = U_0 [m/s]  =>  p(x_0,z,t)=(U_0 L / K + h_L - z) \rho g
       h(x_1,z,t) = h_L [m]    =>  p(x_1,z,t)=(h_L-z) \rho g

    Parameters are in units of:

       \rho : density, [kg/m^3]
       \mu  : viscosity, [kg / m s^2]
       K    : absolute permeability, [ m^2 ]
       g    : gravity, used in converting head to pressure, [ m / s^2 ]

    """

    def __init__(self, params=None):
        if params is None:
            params = dict()

        params.setdefault("x_0",0)
        params.setdefault("x_1",100)
        params.setdefault("z_0",0)
        params.setdefault("z_1",10)

        params.setdefault("k",1.0)
        params.setdefault("rho",998.2)
        params.setdefault("mu",1.002e-3)

        params.setdefault("U_0",0.01)
        params.setdefault("h_1",19.0)

        params.setdefault("g",9.80665)
        params.setdefault("p_atm",101325.0)

        self.__dict__.update(params)

        self.K = self.k*self.g*self.rho / self.mu

    def head(self, coords):

        """
        Compute the head at the x-values given by coords[:]

          h(x) = \frac{U_0 L}{K}(1 - \frac{x}{L}) + h_1

        """

        head = numpy.zeros(len(coords))
        head[:] = self.U_0 * self.x_1 / self.K * (1.0 - coords[:,0]/self.x_1) + self.h_1

        return head

    def pressure(self, coords):

        """
        Compute the pressure at (x,z)-coordinates in the coords[:,:] array.
        Note: coords has dimension len(coords) x 2.
        """
        
        pressure = numpy.zeros((len(coords),),'d')
        head = self.head(coords)

        pressure[:]=self.p_atm+(head[:]-coords[:,1])*self.rho*self.g

        return pressure

def createFromXML(filename):

    # grab params from input file
    params = dict()

    import amanzi_xml.utils.io
    xml = amanzi_xml.utils.io.fromFile(filename)
    import amanzi_xml.utils.search as search
   
    #
    #  Domain Size
    #
    xyz = search.find_tag_path(xml, ["amanzi_input","mesh","generate","box",]).get("low_coordinates")
    params["x_0"] = float(xyz.split(',')[0])
    params["z_0"] = float(xyz.split(',')[2])
    xyz = search.find_tag_path(xml, ["amanzi_input","mesh","generate","box",]).get("high_coordinates")
    params["x_1"] = float(xyz.split(',')[0])
    params["z_1"] = float(xyz.split(',')[2])

    #
    #  Material Properties
    #
    strK = search.find_tag_path(xml, ["amanzi_input","materials","material","permeability",]).get('x')
    params["k"] = float(strK)
    strMu = search.find_tag_path(xml, ["amanzi_input","phases","liquid_phase","viscosity",]).text
    params["mu"] = float(strMu)
    strRho = search.find_tag_path(xml, ["amanzi_input","phases","liquid_phase","density",]).text
    params["rho"] = float(strRho)

    #
    #  Boundary Conditions
    #
    strU0 = search.find_tag_path(xml, ["amanzi_input","boundary_conditions","boundary_condition","liquid_phase","liquid_component","inward_mass_flux",]).get("value")
    params["U_0"] = float(strU0) / params["rho"]
    strhL = search.find_tag_path(xml, ["amanzi_input","boundary_conditions","boundary_condition,RightBC","liquid_phase","liquid_component","hydrostatic",]).get("value")
    params["h_L"] = float(strhL)

    #
    #  Standard Gravity
    #
    params.setdefault("g",9.80665)
   
    # instantiate the class
    return LinearFluxHead(params)

if __name__ == "__main__":

    # Instantiate the class 
    lfh = LinearFluxHead()

    # Get 11 equally spaced points: dx=(x_1-x_0)/10
    x = numpy.linspace(lfh.x_0,lfh.x_1,11)

    # Create space for a set of (x,z) points
    coords = numpy.zeros((11,2))
    # set x
    coords[:,0]=x
    # set z
    coords[:,1]=3
    
    # compute heads and pressures 
    h1 = lfh.head(coords)
    p1 = lfh.pressure(coords)

    # reset z
    coords[:,1]=7

    # compute heads and pressures
    h2 = lfh.head(coords)
    p2 = lfh.pressure(coords)

    # plot 
    plt.plot(x,h1)
    plt.plot(x,h2)
    plt.xlabel('x-coordinate [m]')
#    plt.ylabel('Pressure [Pa]')

    # show the plot
    plt.show()
