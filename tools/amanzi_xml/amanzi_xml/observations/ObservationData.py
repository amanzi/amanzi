import optparse


class ObservationData(object):
    class Data(object):
        def __init__(self, var_ref, region, obs_type, var_name):
            self.var_ref = var_ref
            self.region = region
            self.obs_type = obs_type
            self.var_name = var_name
            self.times = []
            self.data = []
            self.coordinate = None
            self.head = None

    def __init__(self, obs_file):
        self.obs_file = obs_file
        self.observations={}

    def getObservationData(self):
        obs_fid = open(self.obs_file)
        obs_fid.next() # pop header line "name, region,..."
        obs_fid.next() # pop header line "===="

        for line in obs_fid:
            [var_ref, region, obs_type, var_name, time, value] = line.rstrip().split(",")
            var_ref = var_ref.strip()
            region = region.strip()
            obs_type = obs_type.strip()
            var_name = var_name.strip()

            key = (var_ref,region)
            if key not in self.observations:
                self.observations[key] = self.Data(var_ref, region, obs_type, var_name)

            self.observations[key].times.append(float(time.strip()))
            self.observations[key].data.append(float(value.strip()))

    def printSummary(self):
        print("Read observation data file:", self.obs_file)
        for key in self.observations.keys():
            print("  obs:", key[0], "on region", key[1])



if __name__ == "__main__":
    p = optparse.OptionParser()
    p.add_option("--obs_file",action="store", type='string', dest="obs_file")
    p.set_defaults(obs_file = "observation.out") #need to get observation.out from input 
    (options, args) = p.parse_args()

    print("Found filename", options.obs_file)
    print(" creating ObservationData object")
    obs = ObservationData(options.obs_file)
    obs.getObservationData()
    obs.printSummary()

