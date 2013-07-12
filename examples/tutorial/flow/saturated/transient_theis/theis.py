import matplotlib.pyplot as plt
import numpy
import math
from amanzi_xml.observations.ObservationXML import ObservationXML as ObsXML
from amanzi_xml.observations.ObservationData import ObservationData as ObsDATA
import amanzi_xml.utils.search as search

# load input xml file
#  -- create an ObservationXML object
def loadInputXML(filename):
    Obs_xml = ObsXML(filename)
    return Obs_xml

# load the data file
#  -- use above xml object to get observation filename
#  -- create an ObservationData object
#  -- load the data using this object

def loadDataFile(Obs_xml):
    output_file =  Obs_xml.getObservationFilename()
    Obs_data = ObsDATA(output_file)
    Obs_data.getObservationData()
    coords = Obs_xml.getAllCoordinates()

    for obs in Obs_data.observations.itervalues():
        region = obs.region
        obs.coordinate = coords[region]
    
    return Obs_data

def plotTheisObservations(Obs_xml, Obs_data, axes1):
    #=== SPECIAL CODE === Theis EXAMPLE
    r_vals =[]
    for obs in Obs_data.observations.itervalues():
        r_vals.append(obs.coordinate[0])
    r_vals = set(r_vals)
    r_vals=list(r_vals)
    r_vals.sort(reverse=True)
    colors = ['r','g','b']
    cmap = dict((rval,color) for (rval,color) in zip(r_vals, colors))

    for obs in Obs_data.observations.itervalues():
        color = cmap[obs.coordinate[0]]
        pres0 = 101325 - 9793.5192 * (obs.coordinate[2] - 200)
        pres_drop = (pres0 - numpy.array([obs.data]))
        drawdown = pres_drop / 9793.5192
        axes1.scatter(obs.times, drawdown, marker='s', s=25, c=color)
        print  pres_drop

    return cmap

def plotTheisAnalytic(filename, cmap, axes1, Obs_xml ,Obs_data):
    import model.theis
    mymodel = model.theis.createFromXML(filename)
    tindex = numpy.arange(100)
    times = []
    table_values = []
    axes1.set_ylabel('Drawdown [m]')
    axes1.set_xlabel('Time after pumping [s]')
    axes1.set_title('Drawdown vs Time after Pumping')
    
    for i in tindex:
        times.append(1+math.exp(float(i)*(i+1)/(8.5*len(tindex))))
    
    for r in mymodel.r:
       
        drawdown=mymodel.runForFixedRadius(times, r)
        axes1.plot(times,drawdown, label ='$r=%0.1f m$'%r )
        
 # if r == 55:
        #     for t in times:
        #         table_values.append([format(t,"0.2f")
      
        
    axes1.legend(title = 'Theis Solution',loc = 'lower right', fancybox=True, shadow = True)

if __name__ == "__main__":

    import os
    #import run_amanzi

    input_filename = "theis.xml"
    CWD = os.getcwd()

    try:
        obs_xml=loadInputXML("theis.xml")
        obs_data=loadDataFile(obs_xml)

        fig1= plt.figure()
        axes1=fig1.add_axes([.1,.1,.8,.8])
       
        cmap = plotTheisObservations(obs_xml,obs_data, axes1)
        plotTheisAnalytic(input_filename, cmap, axes1, obs_xml, obs_data)
        
        plt.show()

    finally:
        os.chdir(CWD)


