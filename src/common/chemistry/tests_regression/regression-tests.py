#!/bin/env python
#
# cross platform regression testing script...?
#
# read a configuration file with a list of program options, test
# suites and individual tests.
#
# Run the tests/suites specified on the command line
#
# Generate a SHA1 hash of the test output and compare it to a known value.
#
# to generate the hash from an existing output file:
# $ shasum --algorithm 1 --portable file/name
#

import argparse
import ConfigParser
import subprocess
import hashlib

def parse_options():
    parser = argparse.ArgumentParser(description='Run regression tests, generate a SHA1 hash of the output and compare it to a known value.')
    parser.add_argument('-c', '--config-file', nargs=1, default='regression-tests.cfg',
                        help='test configuration file to use')
    parser.add_argument('--debug', action='store_true',
                        help='extra debugging output')
    parser.add_argument('-d', '--dry-run', dest='do_tests',
                        default='store_true', action='store_false',
                        help='perform a dry run, setup the test commands but don\'t run them')
    parser.add_argument('-e', '--executable', nargs=1, default=['/usr/bin/false'],
                        help='executable to use for testing')
    parser.add_argument('--list-suites', action='store_true',
                        help='print the list of test suites from the config file and exit')
    parser.add_argument('--list-tests', action='store_true',
                        help='print the list of tests from the config file and exit')
    parser.add_argument('-s', '--suites', nargs="+", default=False,
                        help='space seperated list of test suite names')
    parser.add_argument('-t', '--tests', nargs="+", default=False,
                        help='space seperated list of test names')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='verbose output')

    options = parser.parse_args()
    return options

def get_configuration(options):
    config = ConfigParser.SafeConfigParser()
    config.read(options.config_file)

    # all sections except 'setup' and 'suites' are tests
    if not config.has_section('setup'):
        raise Exception("""
Config file must contain a \'setup\' section with the following fields:
[setup]
input arg :
input dir :
input suffix :
results dir : 
""")
    setup = config.items('setup')
    setup = convert_list_to_dictionary(setup)

    if not config.has_section('suites'):
        config.add_section('suites')
        config.set('suites', 'none', "")
    available_suites = config.items('suites')
    available_suites = convert_list_to_dictionary(available_suites)

    # extract the test sections
    test_names = config.sections()
    test_names.remove('setup')
    if config.has_section('suites'):
        test_names.remove('suites')

    available_tests = {}
    for t in test_names:
        available_tests[t] = config.items(t)
        available_tests[t] = convert_list_to_dictionary(available_tests[t])

    return setup, available_suites, available_tests

def convert_list_to_dictionary(input_list):
    output_dict = {}
    for item in input_list:
        output_dict[item[0]] = item[1]
    return output_dict

def display_available_tests(available_tests):
    print("Available tests: ")
    for t in sorted(available_tests.keys()):
        print("    {0}".format(t))

def display_available_suites(suites):
    print("Available test suites: ")
    for s in suites.keys():
        print("    {0} :".format(s))
        for t in suites[s].split():
            print("        {0}".format(t))

def identify_tests(options, available_suites, available_tests):
    tests_to_run = {}
    # If no tests or suites were specified on the command line, then run
    # all the tests
    if not options.tests and not options.suites:
        for t in available_tests.keys():
            tests_to_run[t] = available_tests[t]

    # if a test or list of tests in given, run them.
    if options.tests:
        for t in options.tests:
            if t in available_tests.keys():
                tests_to_run[t] = available_tests[t]
            else:
                print("Unknown test specified on command line: {0}".format(t))

    # if suite or list
    # of suites is given, convert to a list of tests and run them.
    if options.suites:
        for s in options.suites:
            if s in available_suites.keys():
                for t in available_suites[s].split():
                    if t in available_tests.keys():
                        tests_to_run[t] = available_tests[t]
                    else:
                        print("Unknown test \'{0}\' specified in suite \'{1}\'".format(t, s))
            else:
                print("Unknown test suite specified on command line: {0}".format(s))

    if options.debug:
        print("Found tests:")
        for t in tests_to_run.keys():
            print("    {0} : ".format(t))
            for k in tests_to_run[t].keys():
                print("        {0} : {1}".format(k, tests_to_run[t][k]))

    return tests_to_run

def run_tests(options, tests_to_run):
    num_failed = 0
    for r in tests_to_run.keys():
        test_info = tests_to_run[r]
        input_file = setup['input dir'] + r + "." + setup['input suffix']
        cmdline = "{0} ".format(options.executable[0])
        cmdline += "{0} {1} ".format(setup['input arg'], input_file)
        if options.verbose or not options.do_tests:
            print(80*'-')
            print("Running test \'{0}\' with the command:\n\t{1}".format(r, cmdline))
            verification_name = None
            if 'verification' in test_info:
                verification_name = test_info['verification']
            print("    verification problem : {0}".format(verification_name))
        else:
            print("{0} ....".format(r))
        if options.do_tests:
            results_name = r + '.stdout'
            subprocess.Popen(['mv', '-f', results_name, results_name+'.old'], 
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
            results = open(results_name, 'w')
            subprocess.Popen(cmdline.split(), 
                             stdout=results, stderr=subprocess.STDOUT).wait()
            results.close()
            results_hash = hashlib.sha1()
            with open(results_name, 'r') as results:
                for line in results:
                    results_hash.update(line)
            #print results_hash.hexdigest()
            if test_info['hash'] == results_hash.hexdigest():
                print("Passed")
            else:
                num_failed += 1
                print("Failed:")
                if options.verbose:
                    print("Expected hash: {0}".format(test_info['hash']))
                    print("Received hash: {0}".format(results_hash.hexdigest()))
                    verified_results_name = setup['results dir'] + r + '.test'
                    diff_command = 'diff {0} {1}'.format(results_name, 
                                                         verified_results_name)
                    print(diff_command)

    print(80*'-')
    if not options.do_tests:
        print("End of dry run.")
    elif num_failed:
        print("FAILURE: {0} of {1} tests failed.".format(num_failed, len(tests_to_run)))
    else:
        print("SUCCESS: {0} tests passed.".format(len(tests_to_run)))


if __name__ == "__main__":
    options = parse_options()
    setup, available_suites, available_tests = get_configuration(options)

    if options.list_tests:
        display_available_tests(available_tests)

    if options.list_suites:
        display_available_suites(available_suites)

    if options.list_suites or options.list_tests:
        exit()

    if options.executable[0] == '/usr/bin/false':
        options.do_tests = False

    tests_to_run = identify_tests(options, available_suites, available_tests)

    run_tests(options, tests_to_run)
