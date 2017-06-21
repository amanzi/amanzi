"""Renames a set of variables in an output file for use in cases where
variable names change but we still need to run regression tests.
2
"""

import h5py

aliases = dict(ponded_depth="surface-ponded_depth")

def rename(datin, datout):
    for var in datin.keys():
        try:
            varout = aliases[var]
        except KeyError:
            varout = var

        grp = datout.create_group(varout)
        for k in datin[var].keys():
            grp.create_dataset(k, data=datin[var][k][:])

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("FILENAME_IN", type=str)
    parser.add_argument("FILENAME_OUT", type=str)

    args = parser.parse_args()

    with h5py.File(args.FILENAME_IN, 'r') as fin, h5py.File(args.FILENAME_OUT, 'w') as fout:
        rename(fin, fout)
    
        
    
