#!/usr/bin/env python

def fix_file(fin,fout):
    for line in fin:
        sline = line.split()
        for i,word in enumerate(sline):
            if "Array" in word:
                nextword = sline[i+1].strip('"')
                line = line.replace("Array "+nextword, "Array("+nextword+")")
        fout.write(line)

def fout_name(fin_name, suffix):
    return fin_name[:-4] + suffix + ".xml"


if __name__ == "__main__":
    import sys,os
    import optparse
    opts = optparse.OptionParser()
    opts.add_option("-o","--out",
                    dest="out",
                    help="output xml file")
    opts.add_option("-s","--suffix",
                    dest="suffix",default="-mod",
                    help="Creates files of name FNAME_IN-mod.xml")
    opts.add_option("-i","--inplace",
                    action="store_true",
                    dest="inplace", default=False,
                    help="Modifies files in-place")

    argv = sys.argv[1:]
    if len(argv) < 1 or not os.path.exists(argv[-1]):
        print opts.get_usage()
        sys.exit()

    infile_or_dir = argv.pop()
    user_opt, other_args = opts.parse_args(argv)
    if len(other_args) != 0:
        print opts.get_usage()
        sys.exit()

    if os.path.isfile(infile_or_dir):
        if not infile_or_dir.endswith(".xml"):
            print opts.get_usage()
            raise RunTimeError("File does not end with .xml")
        if user_opt.outfile is None:
            user_opt.outfile = fout_name(infile_or_dir, user_opt.suffix)

        with open(infile_or_dir,'rU') as fin, open(user_opt.outfile,'wU') as fout:
            fix_file(fin,fout)

        if user_opt.inplace:
            os.rename(user_opt.outfile,infile_or_dir)

    elif os.path.isdir(infile_or_dir):
        for fin_name in os.listdir(infile_or_dir):
            if fin_name.endswith('.xml'):
                fin_full = os.path.join(infile_or_dir, fin_name)
                fout_full = fout_name(fin_full, user_opt.suffix)
                with open(fin_full,'rU') as fin, open(fout_full,'w') as fout:
                    fix_file(fin,fout)

                if user_opt.inplace:
                    os.rename(fout_full, fin_full)

    else:
        print opts.get_usage()
        sys.exit()
