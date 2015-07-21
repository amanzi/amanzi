import sys,os
try:
    sys.path.append(os.path.join(os.environ['ATS_SRC_DIR'], 'tools', 'utils'))
except KeyError:
    pass

import parse_ats
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract a subset of h5 files")
    parser.add_argument("infile", help="input filename")
    parser.add_argument("-d", "--directory", dest="directory", default=".",
                        help="directory containing input")
    parser.add_argument("-o", "--outfile", default=None, help="output filename")
    parser.add_argument("--start", dest="start", type=float, default=0.0, help="start time [yrs]")
    parser.add_argument("--end", dest="end", type=float, default=-1.0, help="end time [yrs], -1 for end of simulation")
    parser.add_argument("--interval", dest="interval", type=int, default=1, help="subsampling interval")
    parser.add_argument("--names", dest="names", type=str, nargs="+", help="subsampling interval")
    args = parser.parse_args()

    assert args.end < 0 or args.end > args.start
    assert args.infile.endswith(".h5")
    
    if args.outfile == None:
        args.outfile = args.infile[:-3]+"_%d_%d_%i.h5"%(args.start,args.end,args.interval)

    if args.end < 0:
        args.end = 1.e99
        
    parse_ats.subsetFile(args.directory, args.infile, args.outfile, interval=args.interval,
                         time_range=(args.start,args.end), names=args.names)
    sys.exit(0)

        
    
