"""Cleans code according to my best guess at Amanzi standards."""

import os
import argparse
import subprocess

def process_exclude_patterns(patterns):
    """Formats exclude patterns into a form easier used with walk."""
    new_patterns = []
    for pattern in patterns:
        split = pattern.split('/')
        new_patterns.append(('/'.join(split[:-1]), split[-1]))
    return new_patterns


def generate_all(dirname, exclude_patterns=None, extensions=None):
    """Generator for files to clean."""
    if exclude_patterns is None:
        exclude_patterns = []
    else:
        exclude_patterns = process_exclude_patterns(exclude_patterns)
        
    if extensions is None:
        extensions = ['.h', '.hh', '.H', '.hpp', '.c', '.cc', '.C', 'cpp']

    for dirname, subdirs, files in os.walk('.'):
        # strip excludes
        for exclude in exclude_patterns:
            if len(exclude[0]) == 0 or dirname == exclude[0]:
                if exclude[1] in subdirs:
                    subdirs.remove(exclude[1])
                elif exclude[1] in files:
                    files.remove(exclude[1])

        # generate
        for f in files:
            if any(f.endswith(ext) for ext in extensions) and \
                (not f.startswith('.')) and (not f.startswith('#')):
                yield os.path.join(dirname, f)
                    
                
def clean_modeline(lines):
    """Removes emacs mode-lines"""
    if lines[0].startswith("/* -*-"):
        print('  removing modeline')
        lines.pop(0)
    return lines

def clean_opening_whitespace(lines):
    """Removes blank lines to start the document."""
    if len(lines) == 0:
        return lines
    while lines[0].strip() == '':
        lines.pop(0)
        print('  cleaning opening whitespace line')
        if len(lines) == 0:
            break
    return lines


def strip_words(string, patterns):
    done = False
    while not done:
        stripped_any = False
        for pattern in patterns:
            if string.startswith(pattern):
                stripped_any = True
                string = string[len(pattern):]
            if string.endswith(pattern):
                stripped_any = True
                string = string[:-len(pattern)]
        if not stripped_any:
            done = True
    return string

copyright = """/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
{}  
*/
"""
def clean_copyright(lines):
    i = 0

    author_lines = []
    writing = False
    while i < len(lines):
        line = lines[i]
        if 'Author:' in line or 'author:' in line or 'authors:' in line or 'Authors:' in line:
            writing = True
            author_lines.append(line)
        elif writing:
            if '*/' == line.strip():
                break
            elif '*/' in line.strip():
                author_lines.append(line.replace('*/',''))
                break
            elif line.strip() == '':
                break
            else:
                author_lines.append(line)
        i += 1
    authors = [strip_words(l.strip(),['author:', 'authors:', 'Author:', 'Authors:']).strip() for l in author_lines]
    authors = [l.strip() for l in authors if len(l) > 0]
    new_copyright = copyright.format('\n'.join([6*' '+a for a in authors]))
    
    as_string = ''.join(lines)
    if 'Copyright' in as_string:
        print('  fixing Copyright')
        split = as_string.split('Copyright')
        if len(split) != 2:
            raise ValueError('More than one "Copyright" in file')
    elif 'COPYRIGHT' in as_string:
        print('  fixing COPYRIGHT')
        split = as_string.split('COPYRIGHT')
        if len(split) != 2:
            raise ValueError('More than one "COPYRIGHT" in file')
    else:
        print('  adding Copyright')
        new_copyright = copyright.format('')
        return new_copyright, lines

    before = '/*'.join(split[0].split('/*')[:-1])
    after = '*/'.join(split[1].split('*/')[1:])
    as_string = before+after
    
    lines = as_string.split('\n')
    for i,l in enumerate(lines):
        lines[i] = l+'\n'
    return new_copyright, lines


def find_oneline_doc(lines):
    try:
        print('  found one-line docstring')
        onelines = [(i,l) for (i,l) in enumerate(lines) if l.strip().startswith('//!')]
    except StopIteration:
        print('  adding one-line docstring')
        oneline = '<MISSING_ONELINE_DOCSTRING>\n'
    else:
        if (len(onelines) == 1):
            lines.pop(onelines[0][0])
            oneline = strip_words(onelines[0][1].strip(), ['//!',]).strip()
        else:
            # pop in reverse to avoid modification issues
            for oneline in reversed(onelines):
                lines.pop(oneline[0])
            oneline = ' '.join([strip_words(l[1].strip(), ['//!',]).strip() for l in onelines])
                
    return '//! '+oneline, lines


def clang_format(filename):
    print('  running clang-format')
    return subprocess.run(['clang-format', '-i', '-style=file', filename])

def clean(filename):
    with open(filename, 'r') as fid:
        lines = fid.readlines()

    print("      ", type(lines), len(lines))
    oneline_doc, lines = find_oneline_doc(lines)
        
    print("      ", type(lines), len(lines))
    lines = clean_opening_whitespace(lines)

    print("      ", type(lines), len(lines))
    lines = clean_modeline(lines)

    print("      ", type(lines), len(lines))
    lines = clean_opening_whitespace(lines)

    print("      ", type(lines), len(lines))
    try:
        new_copyright, lines = clean_copyright(lines)
    except ValueError as err:
        print("  FAILED on Copyright:", err)
        new_copyright = ''
 
    print("      ", type(lines), len(lines))
    lines = clean_opening_whitespace(lines)

    print("      ", type(lines), len(lines))

    new_lines = [l+'\n' for l in new_copyright.split('\n')] + \
        [oneline_doc,'\n', '\n'] + \
        lines

    with open(filename, 'w') as fid:
        fid.write(''.join(new_lines))

    clang_format(filename)


def get_args():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--exclude', type=str, action='append',
                      help='Pattern to exclude from cleaning.')
    parser.add_argument('--extension', type=str, nargs='+', action='append',
                      help='File extensions to clean.')
    parser.add_argument('DIRECTORY', type=str,
                      help='Start searching for files in DIRECTORY')
    return parser.parse_args()

        
if __name__ == '__main__':
    args = get_args()
    if args.exclude is not None:
        print('Excluding:')
        print('-----------')
        for exc in args.exclude:
            print(exc)
    
    for fname in generate_all(args.DIRECTORY, args.exclude, args.extension):
        print('')
        print('Cleaning:')
        print('--------')
        print(fname)
        clean(fname)
