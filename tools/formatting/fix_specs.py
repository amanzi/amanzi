import sys, os
import re



def fixString(string):
    """Returns line with dashes as underscores"""
    for pattern in re.finditer('[a-zA-Z][a-zA-Z\-_]*-spec', string):
        start, end = pattern.span()
        end = end - len('-spec')

        # a few special cases
        if string[start:end].endswith('-typedinline'):
            end = end - len('-typedinline')
        if string[start:end].endswith('-typed'):
            end = end - len('-typed')

        pattern_match = string[start:end]
        string = string[0:start] + pattern_match.replace('-','_') + string[end:]
    return string

def fixFile(fname):
    fname_out = fname+'_tmp'
    with open(fname, 'r') as fin, open(fname_out, 'w') as fout:
        fout.write(fixString(fin.read()))
    os.rename(fname_out, fname)

if __name__ == '__main__':
    arg = sys.argv[-1]
    if os.path.isdir(arg):
        for d, _, files in os.walk(arg):
            for f in files:
                if f.endswith('.hh'):
                    fixFile(os.path.join(d, f))
    elif os.path.isfile(arg):
        fixFile(arg)
    
