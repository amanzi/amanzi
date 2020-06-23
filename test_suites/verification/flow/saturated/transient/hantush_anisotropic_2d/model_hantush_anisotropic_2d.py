import numpy
import matplotlib.pylab as plt
import math
import matplotlib.gridspec as gridspec

class TransientTheis(object):
    """ Solves the system
         \frac{\partial^2{h}}{\partial^2{r}} + \frac{1}{r} \frac{\partial{h}}{\partial{r}} = \frac{S}{T} \frac{\partial{h}}{partial{r}}
         
         Initial Condition:
         h(r,0) = h_0 where h_0 is the initial heigth of the water table at t = 0.

         Boundary Conditions:
         h(\infin,t) = h_0 and constant pumping rate Q [volume/time].

         Theis Solution:
         s = \frac{Q}{4 \pi \sqrt{T_x T_y}} W(u)
         where,
         u = \frac{x^2 T_y + y^2 T_x}{4 T_x T_y t} this quantity is a measure of aquifer response time. 
         W(u) = \int_u^\infty \frac{exp[-u]}{u} du = -0.5772 - ln(u) + u - \frac{u^2}{2*2!} + \frac{u^3}{3*3!} - ...

         Parameters are in units of:
         Q_vol : Pumping Rate [m^3/s]
         s     : Drawdown [m]
         z     : Thickness of Confined Aquifer [m]
         T_x   : Transmissivity in x direction [m^2/s]
         T_y   : Transmissivity in y direction [m^2/s]
         S     : Storage coefficient (unitless)
         t     : duration of pumping [s]
    """

    def __init__(self, params=None):
        if params is None:
            params = dict()
        params.setdefault("z",5)
        params.setdefault("S_s",7.5e-5)
        params.setdefault("pi",math.pi)
        params.setdefault("g",9.80665)
        params.setdefault("Kx",2.35e-11)
        params.setdefault("Ky",2.35e-12)
        params.setdefault("rho",998.2)
        params.setdefault("mu",1.002e-3)
        params.setdefault("Q",-2.0)
       
        self.__dict__.update(params)

        self.Q_vol = -self.Q / self.rho
       
        self.Kx_h = self.Kx * (self.g*self.rho/self.mu) * self.z
        self.Ky_h = self.Ky * (self.g*self.rho/self.mu) * self.z
        
        self.T = math.sqrt(self.Kx_h * self.Ky_h)
        
        self.var = self.Q_vol / 4 / self.pi / self.T
        
        self.S = self.S_s * self.z
       
    def runForFixedPoint(self, times, x, y):
        drawdown_t = []
        for t in times:
            if t == 0.0000:
                drawdown_t.append(0.0)
            else:
                u = (x*x * self.Ky_h + y*y * self.Kx_h) * self.S / 4 / self.T / self.T / t
                W = getWellFunction(u)
                s = self.var * W
                drawdown_t.append(s)
        return drawdown_t
 

def createFromXML(filename):
    #-- grab parameters from input XML file.
    import amanzi_xml.observations.ObservationXMLv2 as Obs_xml
    observations = Obs_xml.ObservationXMLv2(filename)
    xml = observations.xml
    coords = observations.getAllCoordinates()
    import amanzi_xml.utils.search as search
    params = dict()

    params["times"] = []
    for i in search.getElementByTags(xml, "/amanzi_input/definitions/macros/time_macro"):
        params["times"].append(float(i.text))

    params.setdefault("g", 9.80665)
    params.setdefault("pi", math.pi)

    for i in search.getElementByTags(xml, "/amanzi_input/regions"):
        if (i.get("name") == "Well"):
            params["z"] = float(i.find("box").get("high_coordinates").split(',')[2])

    strKx = search.getElementByTags(xml, "/amanzi_input/materials/material/permeability").get('x')
    params["Kx"] = float(strKx)

    strKy = search.getElementByTags(xml, "/amanzi_input/materials/material/permeability").get('y')
    params["Ky"] = float(strKy)

    strMu = search.getElementByTags(xml, "/amanzi_input/phases/liquid_phase/viscosity").text
    params["mu"] = float(strMu)

    strRho = search.getElementByTags(xml, "/amanzi_input/phases/liquid_phase/density").text
    params["rho"] = float(strRho)

    strQ = search.getElementByTags(xml, "/amanzi_input/sources/source/liquid_phase" + 
                                        "/liquid_component/volume_weighted").get("value")
    params["Q"] = float(strQ)

    strSs = search.getElementByTags(xml, "/amanzi_input/materials/material" + 
                                         "/mechanical_properties/specific_storage").get("value")
    params["S_s"] = float(strSs)
    
    return TransientTheis(params)

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
        raise RuntimeError("u cannot be negative check input times and radius!")
    if u < 1:
        W = -math.log(u) + a_0 + a_1*u + a_2*u**2 + a_3*u**3 + a_4*u**4 + a_5*u**5
        return W
    if u >= 1:
        W = (math.exp(-u) / u)*((u**4 +b_1*(u**3) + b_2*(u**2) + b_3*u + b_4)/(u**4 + c_1*(u**3) + c_2*(u**2) + c_3*u + c_4))
        return W

    
