import numpy
import matplotlib.pylab as plt

class LinearMaterialsSerial(object):

    """
    Solves the system:

        \div \frac{\rho}{\mu} k \grad (p + \rho g z) = 0  on the domain [x_0, x_1] \cross [z_0,z_1]

    Boundary conditions are given by:

       h(x_0,z,t) = h_0 [m]
       h(x_1,z,t) = h_L [m]

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
        params.setdefault("L_i",50)
        params.setdefault("z_0",0)
        params.setdefault("z_1",10)

        params.setdefault("k1",1.1847e-12)
        params.setdefault("k2",1.1847e-11)
        params.setdefault("rho",998.2)
        params.setdefault("mu",1.002e-3)

        params.setdefault("h_0",20.0)
        params.setdefault("h_L",19.0)

        params.setdefault("g",9.80665)
        params.setdefault("p_atm",101325.0)

        self.__dict__.update(params)

    def head(self, coords):

        """
        Compute the head at the x-values given by coords[:]

          h(x) = (h_i - h_0)*(x/L_i) + h0, x <= L_i
               = (h_L - h_i)*((x-L_i)/(L-L_i)) + h_i, x > L_i

        """

        head = numpy.zeros(len(coords))
        # Compute the head at the interface.
        h0 = self.h_0
        hL = self.h_L
        L = self.x_1 - self.x_0
        Li = self.L_i
        K1 = self.rho * self.g * self.k1 / self.mu
        K2 = self.rho * self.g * self.k2 / self.mu
        hi = (K1*(L - Li)*h0 + K2*Li*hL) / (K1*(L-Li) + K2*Li)
        for i in range(len(coords[:,0])):
            x = coords[i,0]
            if x <= Li:
                head[i] = (hi - h0)*(x/Li) + h0
            else:
                head[i] = (hL - hi)*((x-Li)/(L-Li)) + hi

        return head

def createFromXML(filename):

    # grab params from input file
    params = dict()

    import amanzi_xml.utils.io
    xml = amanzi_xml.utils.io.fromFile(filename)
    import amanzi_xml.utils.search as search
   
    #
    #  Domain Size
    #
    xyz = search.getElementByTags(xml, "/amanzi_input/mesh/generate/box").get("low_coordinates")
    params["x_0"] = float(xyz.split(',')[0])
    params["z_0"] = float(xyz.split(',')[2])
    xyz = search.getElementByTags(xml, "/amanzi_input/mesh/generate/box").get("high_coordinates")
    params["x_1"] = float(xyz.split(',')[0])
    params["z_1"] = float(xyz.split(',')[2])

    #
    #  Material Properties
    #
    strK = search.getElementByTags(xml, "/amanzi_input/materials/material,Left half/permeability").get('x')
    params["k1"] = float(strK)
    strK = search.getElementByTags(xml, "/amanzi_input/materials/material,Right half/permeability").get('x')
    params["k2"] = float(strK)
    strMu = search.getElementByTags(xml, "/amanzi_input/phases/liquid_phase/viscosity").text
    params["mu"] = float(strMu)
    strRho = search.getElementByTags(xml, "/amanzi_input/phases/liquid_phase/density").text
    params["rho"] = float(strRho)

    #
    #  Boundary Conditions
    #
    strh0 = search.getElementByTags(xml, "/amanzi_input/boundary_conditions/boundary_condition,LeftBC/" +
                                         "liquid_phase/liquid_component/hydrostatic").get("value")
    params["h_0"] = float(strh0)
    strhL = search.getElementByTags(xml, "/amanzi_input/boundary_conditions/boundary_condition,RightBC/" +
                                         "liquid_phase/liquid_component/hydrostatic").get("value")
    params["h_L"] = float(strhL)

    #
    #  Standard Gravity
    #
    params.setdefault("g",9.80665)
   
    # instantiate the class
    return LinearMaterialsSerial(params)

if __name__ == "__main__":

    # Instantiate the class 
    lhh = LinearMaterialsSerial()

    # Get 11 equally spaced points: dx=(x_1-x_0)/10
    x = numpy.linspace(lhh.x_0,lhh.x_1,11)

    # Create space for a set of (x,z) points
    coords = numpy.zeros((11,2))
    # set x
    coords[:,0]=x
    # set z
    coords[:,1]=3
    
    # compute head.
    h1 = lhh.head(coords)

    # reset z
    coords[:,1]=7

    # compute head.
    h2 = lhh.head(coords)

    # plot 
    plt.plot(x,h1)
    plt.plot(x,h2)
    plt.xlabel('x-coordinate [m]')
    plt.ylabel('Head [m]')

    # show the plot
    plt.show()
