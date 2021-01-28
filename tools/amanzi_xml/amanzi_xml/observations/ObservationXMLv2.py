import amanzi_xml.utils.io as io
import amanzi_xml.utils.search as search


class ObservationXMLv2(object):
    def __init__(self, filename):
        self.filename = filename
        self.xml = io.fromFile(filename)
        self.obs_list = self.getObservationList()
        self.obs_lists = []
        for el in self.obs_list.getchildren():
            self.obs_lists.append(el)

        self.obs_file = self.getObservationFilename()
        self.coordinates = []
        self.names = []

    def getObservationList(self):
        return search.find_tag_path(self.xml, ["amanzi_input","output","observations","liquid_phase",])
        #return search.find_path(self.xml, ["amanzi_input","output","observations","liquid_phase",], False)

    def getRegionList(self):
        return search.find_path(self.xml, ["amanzi_input","regions",], False)

    def getAllNames(self):
        self.names = []
        for obs in self.obs_lists:
            self.names.append(obs.get("name"))
        return self.names

    def getAllCoordinates(self):
        self.coordinates = {}
        for item in self.obs_lists:
            children = [el for el in search.children_by_tag(item, "assigned_regions")]
            well_name = children[0].text
            regions = search.find_tag_path(self.xml, ["amanzi_input","regions",])
            children = [el for el in search.children_by_name(regions, well_name)]
            coords = children[0].get("coordinate")
            coordinate = []
            for el in coords.strip(",").split(","):
                coordinate.append(float(el))
            self.coordinates[well_name] = coordinate
        return self.coordinates

    def getCoordinateFromList(self, one_list):
        well_name = one_list.getElement("Region").value
        region = search.child_by_name(search.find_path(self.xml, ["amanzi_input","regions",], False), well_name)

        try:
            coordinate = region.get("coordinate")
        except KeyError:
            raise RuntimeError("Region is not of type point")
        return coordinate

    def getObservationFilename(self):
        obs = search.find_tag_path(self.xml, ["amanzi_input","output","observations","filename",])
        obs_file = obs.text.strip()
        return obs_file

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
    print(" creating ObservationXMLv2 object")
    obs = ObservationXMLv2(options.input_file)
    obs.printSummary()
