import math 
import numpy
import matplotlib.pyplot as plt
drawdown_x = []
drawdown_y = []
drawdown_xy = []
class TransientAnisotropic(object):
    """Governing Equation:
    \frac{\partial{T_x \frac{\partial{h}}{\partial{x}}}}}{\partial{x}}  +  \frac{\partial{T_y \frac{\partial{h}}{\partial{y}}}}}{\partial{y}}  =  S \frac{\partial{h}}{\partial{t}} + Q \delta(x) \delta(y)

    -\infty < x < \infty and  -\infty < y < \infty

    Initial Condition:
    
    h(x,y,0)=h_0

    Solution to Governing Equation (Hantush and Thomas):
    
    s = h_0 - h(r,t) = \frac{Q}{4*pi*(T_x*T_y)**(1/2)} \int_u'^infty \frac{exp[-/tao]}{\tao} d\tao

    where,

    u' = \frac{(x**2*T_y + y**2*T_x)S}{4*T_x*T_y*t}

    """

    def __init__(self,params=None):
       #### Default permeability values are set to satisfy Theis ####
       ####                                                      #### 
        if params is None:
            params = dict()
        params.setdefault("z",10)
        params.setdefault("x",[55,0,55])
        params.setdefault("y",[0,55,55])
        params.setdefault("S_s",7.5e-5)
        params.setdefault("pi",math.pi)
        params.setdefault("g",9.80665)
        params.setdefault("K_x",2.35e-11)
        params.setdefault("K_y",2.35e-10)
        params.setdefault("rho",1000)
        params.setdefault("mu",1.002e-3)
        params.setdefault("Q",-4.0)
        params.setdefault("h_0",15)
       
        self.__dict__.update(params)

        self.x
        print self.x

        self.y
        print self.y

        self.h_0
  
        self.Q_vol = -self.Q / self.rho
        
        self.K_hx = self.K_x*self.g*self.rho / self.mu

        self.K_hy = self.K_y*self.g*self.rho / self.mu

        self.T_x =self.K_hx*self.z
        
        self.T_y=self.K_hy*self.z

        self.var = self.Q_vol / 4 / self.pi / self.T_x**0.5 / self.T_y**0.5
        
        self.S = self.S_s*self.z

    def runAnisotropic(self,times):
        #### This method evaluates Hantush and Thomas solution for coordinates (55,0), (0,55), and (55,55) (z coordinate will be equal for each point) at multiple progressions in time.  If different coordinates are input spec the if statements need to match the graphs sought to be made. ####
        for xa,ya in zip(self.x,self.y):
            if xa == 55 and ya == 0:
                print xa,ya
                for t in times:
                    u = (xa ** 2 * self.T_y + ya**2 * self.T_x) * self.S / 4 / (self.T_x * self.T_y) / t
                    W = getWellFunction(u)
                    s = self.var * W
                    drawdown_x.append(s)
                
            if xa == 0 and ya == 55:
                print xa,ya
                for t in times:
                    u = (xa ** 2 * self.T_y + ya**2 * self.T_x) * self.S / 4 / (self.T_x * self.T_y) / t
                    W = getWellFunction(u)
                    s = self.var * W
                    drawdown_y.append(s)
                        
            if xa == ya:
                print xa,ya
                for t in times:
                    u = (xa ** 2 * self.T_y + ya**2 * self.T_x) * self.S / 4 / (self.T_x * self.T_y) / t
                    W = getWellFunction(u)
                    s = self.var * W
                    drawdown_xy.append(s)
        

def getWellFunction(U):
    u = float(U)
    a_0 =-0.57722
    a_1 = 0.99999
    a_2 =-0.24991
    a_3 = 0.05520
    a_4 =-0.00976
    a_5 = 0.00108
    b_1 = 8.57333
    b_2 = 18.05902
    b_3 = 8.63476
    b_4 = 0.26777
    c_1 = 9.57332
    c_2 = 25.63296
    c_3 = 21.09965
    c_4 = 3.95850

    try:
        u < 0
    except KeyError:
        raise RunTimeError("u cannot be negative check input times and radius!")
    if u < 1:
        W = -math.log(u) + a_0 + a_1*u + a_2*u**2 + a_3*u**3 + a_4*u**4 + a_5*u**5
        return W
    if u >= 1:
        W = (math.exp(-u) / u)*((u**4 +b_1*(u**3) + b_2*(u**2) + b_3*u + b_4)/(u**4 + c_1*(u**3) + c_2*(u**2) + c_3*u + c_4))
        return W
    
    
def createFromXML(filename):
        
    #-- grab parameters from input XML file.
    import amanzi_xml.observations.ObservationXML as Obs_xml
    observations = Obs_xml.ObservationXML(filename)
    xml = observations.xml
    coords = observations.getAllCoordinates()
    import amanzi_xml.utils.search as search
    params = dict()
   
    params["x"] = []
    params["y"] = []
    for coord in coords.itervalues():
        params["x"].append(coord[0])
    print x

    for coord in coords.itervalues():
        params["y"].append(coord[1])
    print y

    params.setdefault("g",9.80665)
    params.setdefault("pi",math.pi)  
    params["z"] = search.getElementByPath(xml, "/Main/Regions/Well/Region: Box/High Coordinate").value[2]
    params["K_x"] = search.getElementByPath(xml, "/Main/Material Properties/Soil/Intrinsic Permeability: Anisotropic Uniform/Horizontal").value
    params["K_y"] = search.getElementByPath(xml, "/Main/Material Properties/Soil/Intrinsic Permeability: Anisotropic Uniform/Vertical").value
    params["mu"] = search.getElementByPath(xml, "/Main/Phase Definitions/Aqueous/Phase Properties/Viscosity: Uniform/Viscosity").value
    params["rho"] = search.getElementByPath(xml, "/Main/Phase Definitions/Aqueous/Phase Properties/Density: Uniform/Density").value
    params["Q"] = search.getElementByPath(xml, "/Main/Sources/Pumping Well/Source: Volume Weighted/Values").value[0]
    params["S_s"] = search.getElementByPath(xml, "/Main/Material Properties/Soil/Specific Storage: Uniform/Value").value

if __name__ == "__main__":
    #instaniate the class
    TA = TransientAnisotropic()
    tindex=numpy.arange(.1,1.0e5)
    times = []

    for i in tindex:
        times.append(i)

    TA.runAnisotropic(times)
    fig1 = plt.figure()

    ax1 = fig1.add_axes([.1,.1,.8,.8])
    ax1.set_xlabel('t [s]')
    ax1.set_ylabel('Drawdown [m]')
    ax1.plot(times,drawdown_x, label = '$x=55$ $y=0$')

    ax1.plot(times,drawdown_y,label = '$x=0$ $y=55$')
    
    ax1.plot(times,drawdown_xy,label = '$x=y=55$')

    ax1.legend(title='Hantush and Thomas',loc = 'lower right', fancybox=True, shadow = True)

    plt.show()
