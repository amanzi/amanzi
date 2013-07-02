import numpy
import matplotlib.pylab as plt
import math

class TransientTheis(object):
    """ Solves the system
         \frac{\partial^2{h}}{\partial^2{r}} + \frac{1}{r} \frac{\partial{h}}{\partial{r}} = \frac{S}{T} \frac{\partial{h}}{partial{r}}
         
         Initial Condition:
         h(r,0) = h_0 where h_0 is the initial heigth of the water table at t = 0.

         Boundary Conditions:
         h(\infin,t) = h_0 and constant pumping rate Q [volume/time].

         Theis Solution:
         s = \frac{Q}{4 \pi T} W(u)
         where,
         u = \frac{r^2 S}{4 T t} this quantity is a measure of aquifer response time. 
         W(u) = \int_u^\infty \frac{exp[-u]}{u} du = -0.5772 - ln(u) + u - \frac{u^2}{2*2!} + \frac{u^3}{3*3!} - ...

         Parameters are in units of:
         Q   : Pumping Rate [m^3/s]
         s   : Drawdown [m]
         h_0 : Initial height of water table [m]
         T   : Transmissivity [m^2/s]
         r   : radial distance measuresd outward from well [m]
         S   : Storage coefficient (unitless)
         t   : duration of pumping [s]
    """

    def __init__(self, params=None):
        if params is None:
            params = dict()
        
        params.setdefault("x_low",-1)
        params.setdefault("x_high",1)
        params.setdefault("z_low",0)
        params.setdefault("z_high",1)
        params.setdefault("y_low",-1)
        params.setdefault("y_high",1)
        params.setdefault("S_s",0.00001)
        params.setdefault("pi",math.pi)
        params.setdefault("g",9.80665)
        params.setdefault("K",1.e-12)
        params.setdefault("rho",1000)
        params.setdefault("mu",1.e-3)
        params.setdefault("h_0", 10)
        params.setdefault("Q",-1.e-3)
        self.__dict__.update(params)
        
        self.K_h = self.K*self.g*self.rho / self.mu
        
        self.T =self. K_h*self.h_0
        
        self.var = -self.Q / 4 / self.pi /self.T
        
        self.S = self.S_s*self.h_0

    def runForFixedTime(self, radi, time):
        #This method evaluates Theis solution for a given time at multiple radial distances from the well 
        drawdown_r = []
        for r in radi:
            u = r ** 2 * self.S / 4 / self.T / time 
            W = -0.577216 - math.log(u) + u - (u**2 /( 2 * math.factorial(2))) + (u**3 /( 3 * math.factorial(3))) - (u**4/(4*math.factorial(4))) + (u**5/(5*math.factorial(5))) - (u**6/(6*math.factorial(6))) + (u**7/(7*math.factorial(7))) 
            s = self.var*W
            drawdown_r.append(s)
        return drawdown_r

    def runForFixedRadius(self, times,radius):
        #This method evaluates Theis solution for a given radius at multiple progressions in time
        drawdown_t = []
        for t in times:
            u = radius ** 2 * self.S / 4 / self.T / t
            W = -0.577216 - math.log(u) + u - (u**2 /(2*math.factorial(2))) + (u**3/(3*math.factorial(3))) - (u**4/(4*math.factorial(4))) + (u**5/(5*math.factorial(5))) - (u**6/(6*math.factorial(6))) + (u**7/(7*math.factorial(7))) 
            s = self.var*W
            drawdown_t.append(s)
        return drawdown_t
 
def paramsFromXML(filename):
        
        #grab parameters from input XML file.
    import amanzi_xml.utils.io
    xml = amanzi_xml.utils.io.fromFile(filename)
    import amanzi_xml.utils.search as search

    params = dict()
    obs_lists = []
    obs_list = search.getElementByPath(xml, "/Main/Regions/Region: Point")
    for el in obs_list.getchildren():
        if el.tag == "Parameter":
            obs_lists.append(el)
    print obs_list

    params.setdefault("g",9.80665)
    params.setdefault("pi",math.pi)
    params["K"] = search.getElementByPath(xml, "/Main/Material Properties/Soil/Intrinsic Permeability: Uniform/Value").value
    params["mu"]= search.getElementByPath(xml, "/Main/Phase Definitions/Aqueous/Phase Properties/Viscosity: Uniform/Viscosity").value
    params["rho"]= search.getElementByPath(xml, "/Main/Phase Definitions/Aqueous/Phase Properties/Density: Uniform/Density").value
    params["h_0"]= search.getElementByPath(xml, "/Main/Boundary Conditions/Far Field Head/BC: Hydrostatic/Water Table Height").value[0]
    params["Q"]= search.getElementByPath(xml, "/Main/Sources/Pumping Well/Source: Volume Weighted/Values").value[0]
    params["S_s"]= search.getElementByPath(xml, "/Main/Material Properties/Soil/Specific Storage: Uniform/Value").value
    params["n"]= search.getElementByPath(xml, "/Main/Material Properties/Soil/Porosity: Uniform/Value").value
    params["x_low"]=search.getElementByPath(xml, "/Main/Regions/Entire Domain/Region: Box/Low Coordinate").value[0] 
    params["y_low"]=search.getElementByPath(xml, "/Main/Regions/Entire Domain/Region: Box/Low Coordinate").value[1]    
    params["z_low"]=search.getElementByPath(xml, "/Main/Regions/Entire Domain/Region: Box/Low Coordinate").value[2]
    params["x_high"]=search.getElementByPath(xml, "/Main/Regions/Entire Domain/Region: Box/High Coordinate").value[0]  
    params["y_high"]=search.getElementByPath(xml, "/Main/Regions/Entire Domain/Region: Box/High Coordinate").value[1] 
    params["z_high"]=search.getElementByPath(xml, "/Main/Regions/Entire Domain/Region: Box/High Coordinate").value[2]

    return TransientTheis(params)

if __name__ == "__main__":
    #instaniate the class
    Tt = TransientTheis()
    
    radi = []
    rindex = numpy.arange(.1,50,.3)
    time = [1500,3600,86400,360000] #100 hours
    
    tindex=numpy.arange(100)
    times=[]

    for i in rindex:
        radi.append(i)

    for i in tindex:
        times.append(50+math.exp(float(tindex[i])*(i+1)/(10*len(tindex))))

    radius = [1,5,10,20]

#s_fixed_time = Tt.runForFixedTime(radi, time)
# s_fixed_radius = Tt.runForFixedRadius(times, radius)

    fig1 = plt.figure()
    fig2 = plt.figure()
    ax1 = fig1.add_axes([.1,.1,.8,.8]) 
    ax2 = fig2.add_axes([.1,.1,.8,.8])
    ax1.set_xlabel('Radial Distance from Well [m]')
    ax1.set_ylabel('Drawdown [m]')
    ax2.set_xlabel('Time after pumping [s]')
    ax2.set_ylabel('Drawddown [m]')
    ax1.set_title('Drawdown vs Radial Distance')
    ax2.set_title('Drawdown vs Time after Pumping')

    for r in radius:
        s_fixed_radius = Tt.runForFixedRadius(times, r)
        ax2.plot(times , s_fixed_radius, label ='$r=%0.1f m$'%r )
    

    for t in time: 
        s_fixed_time = Tt.runForFixedTime(radi, t)
        ax1.plot(radi , s_fixed_time, label = '$t=%0.1f s$'%t )
    

    ax1.legend(title='Theis Solution',loc = 'upper right', fancybox=True, shadow = True)
    ax2.legend(title = 'Theis Solution',loc = 'lower right', fancybox=True, shadow = True)

    plt.show()
    

    
    
