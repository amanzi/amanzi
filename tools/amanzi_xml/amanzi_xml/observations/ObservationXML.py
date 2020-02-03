import amanzi_xml.utils.io as io
import amanzi_xml.utils.search as search


class ObservationXML(object):
    def __init__(self, filename):
        self.filename = filename
        self.xml = io.fromFile(filename)
        self.obs_list = self.getObservationList()
        self.obs_lists = []
        for el in self.obs_list.getchildren():
            if el.tag == "ParameterList":
                self.obs_lists.append(el)

        self.obs_file = self.getObservationFilename()
        self.coordinates = []
        self.names = []

    def getObservationList(self):
        return search.getElementByTagPath(self.xml, "/Main/Output/Observation Data")

    def getRegionList(self):
        return search.getElementByTagPath(self.xml, "/Main/Regions")

    def getAllNames(self):
        self.names = []
        for obs in self.obs_lists:
            self.names.append(obs.get("name"))
        return self.names

    def getAllCoordinates(self):
        self.coordinates = {}
        for item in self.obs_lists:
            well_name = item.getElement("Region").value
            region = search.getElementByTagPath(self.xml, "/Main/Regions/"+ well_name)
            local = region.sublist("Region: Point")
            location = search.getElementByTagPath(local, "/Region: Point/Coordinate")
            coordinate = location.value
            self.coordinates[well_name] = coordinate
        return self.coordinates

    def getCoordinateFromList(self, one_list):
        well_name = one_list.getElement("Region").value
        region = search.getElementByTagPath(self.xml, "/Main/Regions/"+ well_name)
        local = region.sublist("Region: Point")
        location = search.getElementByTagPath(local, "/Region: Point/Coordinate")
        coordinate = location.value
        return coordinate

    def getObservationFilename(self):
        obs = search.getElementByTagPath(self.xml, "/Main/Output/Observation Data/Observation Output Filename")
        obs_file = obs.value.strip()
        return  obs_file

    def printSummary(self):
        print("Read input file:", self.filename)
        print("  Found", len(self.obs_lists), "observations")


if __name__ == "__main__":
    import optparse
    p = optparse.OptionParser()
    p.add_option("--input_file",action="store", type='string', dest="input_file")
    p.set_defaults(obs_file = "f.xml") #need to get observation.out from input 
    (options, args) = p.parse_args()

    print("Found input filename", options.input_file)
    print(" creating ObservationXML object")
    obs = ObservationXML(options.input_file)
    obs.printSummary()
