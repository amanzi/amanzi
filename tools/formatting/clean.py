import os

def generate_all(dirname='.'):
    for dirname, _, files in os.walk('.'):
        for f in files:
            if (f.endswith('.hh') or f.endswith('.cc')) and not f.startswith('.'):
                yield os.path.join(dirname, f)

                
def clean_modeline(lines):
    """Removes emacs mode-lines"""
    if lines[0].startswith("/* -*-"):
        lines.pop(0)
    return lines

def clean_opening_whitespace(lines):
    """Removes emacs mode-lines"""
    while lines[0].strip() == '':
        lines.pop(0)
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
    new_copyright = copyright.format('\n'.join([6*' '+a for a in authors]))
    
    as_string = ''.join(lines)
    split = as_string.split('Copyright')
    if len(split) is not 2:
        raise ValueError('More than one "Copyright" in file')

    before = '/*'.join(split[0].split('/*')[:-1])
    after = '*/'.join(split[1].split('*/')[1:])
    as_string = before+after
    
    lines = as_string.split('\n')
    for i,l in enumerate(lines):
        lines[i] = l+'\n'
    return new_copyright, lines


def find_oneline_doc(lines):
    try:
        i, oneline = next((i,l) for (i,l) in enumerate(lines) if l.strip().startswith('//!'))
    except StopIteration:
        oneline = '//! <MISSING_ONELINE_DOCSTRING>\n'
    else:
        lines.pop(i)
    return oneline, lines
        

def do_it(lines):
    print(type(lines), len(lines))
    oneline_doc, lines = find_oneline_doc(lines)
        
    print(type(lines), len(lines))
    lines = clean_opening_whitespace(lines)

    print(type(lines), len(lines))
    lines = clean_modeline(lines)

    print(type(lines), len(lines))
    lines = clean_opening_whitespace(lines)

    print(type(lines), len(lines))
    new_copyright, lines = clean_copyright(lines)
 
    print(type(lines), len(lines))
    lines = clean_opening_whitespace(lines)

    print(type(lines), len(lines))

    new_lines = [l+'\n' for l in new_copyright.split('\n')] + \
        [oneline_doc,'\n',] + \
        lines
    return new_lines
    
        
if __name__ == '__main__':
    for fname in generate_all():
        print('Running: ', fname)
        with open(fname, 'r') as fid:
            lines = fid.readlines()

        try:
            lines = do_it(lines)
        except ValueError:
            pass
        else:
            with open(fname, 'w') as fid:
                fid.write(''.join(lines))
