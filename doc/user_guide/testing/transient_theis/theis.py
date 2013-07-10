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
    Obs_xml.getAllCoordinates()
    return Obs_xml
            
# load the data file
#  -- use above xml object to get observation filename
#  -- create an ObservationData object
#  -- load the data using this object

def loadDataFile(Obs_xml):
    output_file =  Obs_xml.getObservationFilename()
    Obs_data = ObsDATA(output_file)
    Obs_data.getObservationData()
    return Obs_data

def plotTheisObservations(Obs_xml, Obs_data, axes1):
    #=== SPECIAL CODE === Theis EXAMPLE
    # Collect the used z-values
    r_vals = [coord[0] for coord in Obs_xml.coordinates]
    r_vals = set(r_vals)
    colors = ['r','b','g']
    cmap = dict((rval,color) for (rval,color) in zip(r_vals, colors))
    h_0 = 100
    
    #-- grab the correct coordinate and corresponding time 
    coord = Obs_xml.coordinates 
    for pressure in Obs_data.pressure_value:
        for time in Obs_data.times:
            head = pressure / 9793.5192
            drawdown = h_0 - head
            axes1.scatter(time, drawdown, marker='s', s=25)
    return cmap

def plotTheisAnalytic(filename, cmap, axes1, Obs_xml ,Obs_data):
    import model.theis
    mymodel = model.theis.createFromXML(filename)
    tindex = numpy.arange(100)
    times =[]
    axes1.set_ylabel('Drawdown [m]')
    axes1.set_xlabel('Time after pumping [s]')
    axes1.set_title('Drawdown vs Time after Pumping')
    
    for time in tindex:
        times.append(1+math.exp(float(time)*(time+1)/(10.5*len(tindex))))
    
    for r in mymodel.r:
        drawdown=mymodel.runForFixedRadius(times, r)
        axes1.plot(times,drawdown, label ='$r=%0.1f m$'%r )
        
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


