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
         s = \frac{Q}{4 \pi T} W(u)
         where,
         u = \frac{r^2 S}{4 T t} this quantity is a measure of aquifer response time. 
         W(u) = \int_u^\infty \frac{exp[-u]}{u} du = -0.5772 - ln(u) + u - \frac{u^2}{2*2!} + \frac{u^3}{3*3!} - ...

         Parameters are in units of:
         Q_vol : Pumping Rate [m^3/s]
         s     : Drawdown [m]
         z     : Thickness of Confined Aquifer [m]
         T     : Transmissivity [m^2/s]
         r     : radial distance measuresd outward from well [m]
         S     : Storage coefficient (unitless)
         t     : duration of pumping [s]
    """

    def __init__(self, params=None):
        if params is None:
            params = dict()
        params.setdefault("z",10)
        params.setdefault("r", [20,30,55])
        params.setdefault("S_s",7.5e-5)
        params.setdefault("pi",math.pi)
        params.setdefault("g",9.80665)
        params.setdefault("K",2.35e-11)
        params.setdefault("rho",1000)
        params.setdefault("mu",1.002e-3)
        params.setdefault("Q",-4.0)
        params.setdefault("times",[.1,.1,100,200,500,1000,2000,3600])
       
        self.__dict__.update(params)

        self.r
  
        self.Q_vol = -self.Q / self.rho
       
        self.K_h = self.K*self.g*self.rho / self.mu
        
        self.T =self.K_h*self.z
        
        self.var = self.Q_vol / 4 / self.pi / self.T
        
        self.S = self.S_s*self.z
       
    def runForFixedTime(self, radi, time):
        #This method evaluates Theis solution for a given time at multiple radial distances from the well 
        drawdown_r = []
        for r in radi:
            u = r ** 2 * self.S / 4 / self.T / time 
            W = getWellFunction(u)
            s = self.var*W
            drawdown_r.append(s)
        return drawdown_r

    def runForFixedRadius(self, times,radius):
        #This method evaluates Theis solution for a given radius at multiple progressions in time
        drawdown_t = []
        for t in times:
            if t == 0.0000:
                drawdown_t.append(0.0)
            else:
                u = radius ** 2 * self.S / 4 / self.T / t
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

    params["r"] = []
    for (coord) in coords.itervalues():
        params["r"].append(coord[0]) 
    
    params["times"] = []
    for i in search.getElementByTags(xml, "/amanzi_input/definitions/macros/time_macro"):
        params["times"].append(float(i.text))

    params.setdefault("g", 9.80665)
    params.setdefault("pi", math.pi)

    for i in search.getElementByTags(xml, "/amanzi_input/regions"):
        if (i.get("name") == "Well"):
            params["z"] = float(i.find("box").get("high_coordinates").split(',')[2])

    strK = search.getElementByTags(xml, "/amanzi_input/materials/material/permeability").get('x')
    params["K"] = float(strK)

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

if __name__ == "__main__":
    #instaniate the class
    Tt = TransientTheis()
    
    radi = []
    rindex = numpy.arange(.1,50,.3)
    time = [1500,3600,86400,360000, 860000] #100 hours
    
    tindex=numpy.arange(1.2e2,1.0e5,120)
    times=[]
    table_values=[]
    error=[]
    PORFLOW =[6.2764e-3,2.9013E-02,5.3206E-02,7.5200E-02,9.4700E-02,1.1203E-01,1.2757E-01,1.4161E-01,1.5440E-01,1.6614E-01,1.7699E-01,1.8706E-01,1.9645E-01,2.0526E-01,2.1355E-01,2.2137E-01,2.2877E-01,2.3581E-01,2.4250E-01,2.4889E-01,2.5500E-01,2.6085E-01,2.6647E-01,2.7186E-01,2.7706E-01,2.8207E-01,2.8690E-01,2.9157E-01,2.9609E-01,3.0046E-01]

    for i in rindex:
        radi.append(i)

    for i in tindex:
        times.append(i)
#1+math.exp(float(i)*(i+1)/(8.5*len(tindex))))

    radius = [55]

    col_labels = ['time [s]','Theis [m]','PORFLOW [m]','Difference']
    fig1 = plt.figure()
    fig2 = plt.figure()
    fig3 = plt.figure()
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
        if r == 55:
            for s,t,p in zip(s_fixed_radius, times,PORFLOW):
                
                table_values.append([format(t,"0.1f") ,format(s,"0.7f") , p,format(s-p,"0.7f")])
        
    
    for t in time: 
        s_fixed_time = Tt.runForFixedTime(radi, t)
        ax1.plot(radi , s_fixed_time, label = '$t=%0.1f s$'%t )
    
    ax1.legend(title='Theis Solution',loc = 'upper right', fancybox=True, shadow = True)
    ax2.legend(title = 'Theis Solution',loc = 'lower right', fancybox=True, shadow = True)

    cellText = table_values
    ax = plt.subplot(111,frame_on = False)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    the_table = plt.table(cellText = table_values ,  colWidths = [.1 , .3 , .3,.3 ], colLabels = col_labels,  cellLoc = 'center', colLoc = 'center', loc = 'best')
    properties = the_table.properties()
    cells = properties['child_artists']
    for cell in cells: 
        cell.set_height(.045)
    
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(12)
   
    plt.show()
    

    
    
