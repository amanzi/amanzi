import os
import sys
import h5py
import numpy as np

def get_observation_data(obs_file, observation):

    # open amanzi observations file
    datafile = open(obs_file,'r')

    # parse file to extra location and value
    # throw away 2 header lines
    datafile.readline()
    datafile.readline()
    # read the rest
    obs = []
    for line in datafile.readlines():
        vals = line.split(",")
        obs.append(vals)

    # get observation regions and names
    obs_regions = []
    obs_names = []
    for ob in obs:
        if (ob[1] not in obs_regions):
            obs_regions.append(ob[1].strip())
        if (ob[3] not in obs_names):
            obs_names.append(ob[3].strip())
    if observation not in obs_names:
        print('%s: observation %s was not found in file %s.'%(sys.argv[0], observation, obs_file))
        exit(0)
    obs_coords = {}
    obs_values = {}
    for region in obs_regions:
        tmp_coords = []
        tmp_values = []
        for ob in obs:
            if (ob[1].strip() == region):
                tmp_coords.append(float(ob[len(ob)-2]))
                tmp_values.append(float(ob[len(ob)-1]))
        obs_coords[region] = np.array(tmp_coords)
        obs_values[region] = np.array(tmp_values)

    datafile.close()

    return (obs_coords, obs_values)

def find_output_observation_path(input_xml):
    import xml.etree.ElementTree as ET
    tree = ET.parse(input_xml)
    root = tree.getroot()
    if root.tag != 'amanzi_input':
        raise RuntimeError('The given XML file is not an Amanzi input.')
    output = root.find('output')
    observations = output.find('observations')
    obs_file = observations.find('filename').text
    return obs_file

def compute_error_norms(x_output, c_output, x_ref, c_ref, relative_tol, absolute_tol):
    errors_tols = {}
    for region in x_ref.keys():
        diff = c_ref[region] - c_output[region]
        # Save a tuple of errors and tolerance
        errors_tols[region] = (np.linalg.norm(diff),
                               relative_tol*np.linalg.norm(c_ref[region]) + absolute_tol)

    return errors_tols

if __name__ == "__main__":

    # Check argument number.
    if len(sys.argv) != 6:
        print("%s: usage:"%sys.argv[0])
        print("%s observation input.xml reference.obs relative_tolerance absolute_tolerance"%sys.argv[0])
        exit(0)

    # Get arguments.
    observation = sys.argv[1]
    input_xml = sys.argv[2]
    ref_file = sys.argv[3]
    relative_tol = float(sys.argv[4])
    absolute_tol = float(sys.argv[5])

    # Check arguments.
    if not os.path.exists(input_xml):
        print('%s: input file %s does not exist.'%(sys.argv[0], input_xml))
        exit(0)
    if not os.path.exists(ref_file):
        print('%s: reference file %s does not exist.'%(sys.argv[0], ref_file))
        exit(0)
    if (relative_tol < 0.0) or (absolute_tol < 0.0):
        print('%s: tolerance must be non-negative.'%sys.argv[0])
        exit(0)

    # Root through the input file and find out where the output lives.
    obs_file = find_output_observation_path(input_xml)
    if not os.path.exists(obs_file):
        print('%s: output file %s does not exist.'%(sys.argv[0], obs_file))

    # Now fetch the data from the output and reference files.
    x_output, c_output = get_observation_data(obs_file, observation)
    x_ref, c_ref = get_observation_data(ref_file, observation)

    # Compute the error norms for each of the observations in their respective regions.
    errors_tols = compute_error_norms(x_output, c_output, x_ref, c_ref, relative_tol, absolute_tol)
    msg = ""

    # Report.
    total_errors = 0
    error_lines = []
    for region, value in errors_tols.items():
        error, tolerance = value[0], value[1]
        if (error > tolerance):
            total_errors += 1
            error_lines.append('  error %g > tolerance %g in region %s'%(error, tolerance, region))

    if total_errors == 0:
        sys.exit('Comparison Passed')
    else:
        sys.exit('\n'.join(['Comparison Failed for observation %s:'%observation] + error_lines))
