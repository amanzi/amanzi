import numpy
import matplotlib.pylab as plt

class UnconfinedSeepageHead(object):

    ft2m = 0.3048        # feet to meters 
    m2ft = 3.28084       # meters to feet
    d2s  = 86400         # days to seconds
    y2s  = 31536000      # years to seconds

    """
    Solves the system:  MUST HAVE RICHARDS EQUATION HERE?

        \div \frac{\rho}{\mu} k \grad (p + \rho g z) = 0  on the domain [x_0, x_1] \cross [z_0,z_1]

    Boundary conditions are given by:

       h(x_0,z,t) = h_0 [m]  =>  p(x_0,z,t)=(h_0-z) \rho g
       h(x_1,z,t) = h_L [m]  =>  p(x_1,z,t)=(h_L-z) \rho g

    Parameters are in units of:

       \rho : density, [kg/m^3]
       \mu  : viscosity, [kg / m s^2]
       K    : absolute permeability, [ m^2 ]
       g    : gravity, used in converting head to pressure, [ m / s^2 ]

    """

    def __init__(self, params=None):
        if params is None:
            params = dict()

        """ all numbers in metric system """

        params.setdefault("x_0",0)
        params.setdefault("x_1",self.ft2m*1000)  
        params.setdefault("z_0",self.ft2m*100)  
        params.setdefault("z_1",self.ft2m*50)  

        params.setdefault("k",self.ft2m*self.ft2m*1)   
        params.setdefault("rho",998.2)  
        params.setdefault("mu",1.002e-3)

        params.setdefault("h_0",self.ft2m*80)   
        params.setdefault("h_1",self.ft2m*50)    
        params.setdefault("Q_src",self.ft2m*1/self.y2s)  # 1 ft/year to m/s

        params.setdefault("g",9.80665)
        params.setdefault("p_atm",101325.0)

        self.__dict__.update(params)

        self.K = self.k*self.g*self.rho / self.mu

        """
        We should really write the expression for this 
        """

        self.L_s = self.ft2m*829.0
        self.h_s = (self.ft2m*50.0)*(2-self.L_s/self.x_1)

    def head(self, coords):

        """
        Compute the head at the x-values given by coords[:]

          h(x)^2 = h_0^2 + (h_1^2 - h_0^2)(x/L) + 
                   ((Q_src)(L^2)/K)(x/L)(1-x/L)

        where L=x_1-x_0

        """

        head = []
        Ls = self.L_s
        h0 = self.h_0
        h1 = self.h_1
        hs = self.h_s

        for x in coords[:,0]:
            if x < Ls:
                a = x / Ls
                tmp = numpy.sqrt(h0*h0 + (hs*hs - h0*h0) * a + (self.Q_src*Ls*Ls/self.K) * a * (1 - a))
            else:
                a = x / self.x_1
                tmp = h1 * (2 - a)
            head.append(tmp)

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
    strH0 = search.find_tag_path(xml, ["amanzi_input","boundary_conditions","boundary_condition,LeftBC","liquid_phase","liquid_component","hydrostatic",]).get("value")
    params["h_0"] = float(strH0)

    strH1 = search.find_tag_path(xml, ["amanzi_input","boundary_conditions","boundary_condition,RightBC","liquid_phase","liquid_component","hydrostatic",]).get("value")
    params["h_1"] = float(strH1)

    #
    #  Standard Gravity
    #
    params.setdefault("g",9.80665)
   
    # instantiate the class
    return UnconfinedSeepageHead(params)

if __name__ == "__main__":

    # Instantiate the class 
    ush = UnconfinedSeepageHead()

    # Get 11 equally spaced points: dx=(x_1-x_0)/10
    x = numpy.linspace(ush.x_0,ush.x_1,11)

    # Create space for a set of (x,z) points
    coords = numpy.zeros((11,2))
    # set x
    coords[:,0]=x
    # set z
    coords[:,1]=0.524
    
    # compute heads (converted to feet) and pressures 
    h1 = ush.head(coords)
    p1 = ush.pressure(coords)
    
    # plot 
    plt.plot(x,h1)

    plt.xlabel('x-coordinate [m]')
    plt.ylabel('Hydraulic head [m]')

    # show the plot
    plt.show()
