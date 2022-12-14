"""Cleans code according to my best guess at Amanzi standards."""

import os
import argparse
import subprocess
import warnings


def process_exclude_patterns(patterns):
    """Formats exclude patterns into a form easier used with walk."""
    new_patterns = []
    for pattern in patterns:
        split = pattern.split('/')
        new_patterns.append(('/'.join(split[:-1]), split[-1]))
    return new_patterns


def generate_all(top, exclude_patterns=None, extensions=None):
    """Generator for files to clean."""
    if os.path.isfile(top):
        yield top

    else:
        if exclude_patterns is None:
            exclude_patterns = []
        else:
            exclude_patterns = process_exclude_patterns(exclude_patterns)

        if extensions is None:
            extensions = ['.h', '.hh', '.H', '.hpp', '.c', '.cc', '.C', 'cpp']

        for dirname, subdirs, files in os.walk(top):
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


def clean_trailing_whitespace(lines):
    """Removes blank lines to end the document."""
    if len(lines) == 0:
        return lines
    while lines[-1].strip() == '':
        lines.pop(-1)
        print('  cleaning opening whitespace line')
        if len(lines) == 0:
            break
    return lines


def strip_words(string, patterns):
    """Similar to string.strip, but with words instead of characters"""
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


def find_block_comment_bounds(lines, i):
    """Finds the begin and end of a block comment in lines that includes line i.

    Note this assumes that the comment is a BLOCK comment.

    Returns first_line, last_line+1 e.g. standard python indexing conventions.
    """
    begin_line = i
    end_line = i

    if lines[i].strip().startswith('//'):
        # block comment of lines starting with //
        while lines[begin_line].strip().startswith('//'):
            begin_line = begin_line - 1
        while lines[end_line].strip().startswith('//'):
            end_line = end_line + 1

    else:
        # comment that starts with /*
        while '/*' not in lines[begin_line]:
            begin_line = begin_line - 1
        begin_line = begin_line - 1
        while '*/' not in lines[end_line]:
            end_line = end_line + 1
        end_line = end_line + 1
    return begin_line+1, end_line


def strip_comments_line(l):
    """Removes comment strings, along with ---- or ====="""
    l = l.strip()
    l = strip_words(l, ['/*', '*/', '//', 'ATS', '*'])
    l = l.strip()

    # also remove lines /* --------  and -------*/
    if l == '-'*len(l):
        l = ''
    if l == '='*len(l):
        l = ''
    return l


def clean_block_comment(lines):
    lines = [strip_comments_line(l) for l in lines]
    lines = clean_opening_whitespace(lines)
    lines = clean_trailing_whitespace(lines)
    return ['/*',] + lines + ['','*/']


def is_empty_block(lines):
    for l in lines:
        sc = strip_comments_line(l)
        if sc != '' and sc != 'ATS':
            return False
    return True


def find_remove_authors(lines):
    """Searches for block comments containing Author info, returning the
    author list and the modified lines with the author list removed."""

    author_lines = [i for (i,l) in enumerate(lines) if any((a in l.lower()) for a in ['author:', 'authors:'])]
    ret_lines = []

    if len(author_lines) > 0:
        author_block_extent = find_block_comment_bounds(lines, author_lines[0])
        end_authors = author_lines[0] + 1
        while end_authors < author_block_extent[1] and strip_comments_line(lines[end_authors]) != '':
            end_authors = end_authors + 1

        ret_lines = [strip_comments_line(lines[al]) for al in range(author_lines[0], end_authors)]
        ret_lines[0] = strip_words(ret_lines[0], ['Author', 'AUTHOR', 'author'])
        if ret_lines[0][0] in ['s', 'S']:
            ret_lines[0] = ret_lines[0][1:]
        if ret_lines[0].startswith(':'):
            ret_lines[0] = ret_lines[0][1:]
        ret_lines[0] = ret_lines[0].strip()
        if ret_lines[0] == '':
            ret_lines.pop(0)

        remaining_comment_block = lines[author_block_extent[0]:author_lines[0]] + lines[end_authors:author_block_extent[1]]

        if is_empty_block(remaining_comment_block):
            lines = lines[0:author_block_extent[0]] + lines[author_block_extent[1]:]
        else:
            if lines[end_authors].strip() == '' or lines[end_authors].strip() == "//":
                end_authors = end_authors+1
            lines = lines[0:author_lines[0]] + lines[end_authors:]
        
    return ret_lines, lines


def find_paragraph_bounds(lines, paragraph_line, block_extent):
    begin_paragraph = paragraph_line
    while begin_paragraph >= block_extent[0] and \
          strip_comments_line(lines[begin_paragraph-1]) != '':
        begin_paragraph = begin_paragraph - 1

    end_paragraph = paragraph_line + 1
    while end_paragraph < block_extent[1] and \
          strip_comments_line(lines[end_paragraph]) != '':
        end_paragraph = end_paragraph + 1
    return begin_paragraph, end_paragraph


    
copyright = """/*
  Copyright 2010-202x held jointly by participating institutions.
  {} is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:{}
*/
"""

TPL_COPYRIGHTS = ['Copyright (c) 2009',]

def find_remove_copyright(lines, ats=False):
    """Searches for block comments containing COPYRIGHT info, returning
    the modified lines with the copyright info removed."""

    # find the authors
    authors, lines = find_remove_authors(lines)

    # construct a new copyright comment block with authors
    if ats: code = 'ATS'
    else: code = 'Amanzi'
    author_block = '\n           '.join(authors)
    if len(author_block) > 0:
        author_block = ' '+author_block
    new_copyright = copyright.format(code, author_block).split('\n')

    # find all old copyright comment block(s)
    copyright_lines = [j for j,l in enumerate(lines) if \
                       any(c in l.lower() for c in ['copyright', 'license: bsd', 'bsd license', 'licencse: bsd'])]

    # Note that these may overlap.
    copyright_blocks = []
    for copyright_line in copyright_lines:
        # check if this is a THIRD PARTY COPYRIGHT!
        if any(tpl_language in lines[copyright_line] for tpl_language in TPL_COPYRIGHTS):
            # stop now!  more will be found, and these should stay!
            break
        
        # search for the block of comments containing this copyright
        # line.  note that this ASSUMES the copyright is an a comment
        # block
        comment_block_extent = find_block_comment_bounds(lines, copyright_line)

        # find the bounds of the paragraph in the comment block
        copyright_block_extent = find_paragraph_bounds(lines, copyright_line, comment_block_extent)
                    
        copyright_blocks.append((comment_block_extent, copyright_block_extent))

    # make unique, and sort by block
    blocks = dict()
    for cb in copyright_blocks:
        if cb[0] not in blocks:
            blocks[cb[0]] = set()
        blocks[cb[0]].add(cb[1])

    lines_to_remove = []
        
    for comment_block in blocks.keys():
        print('FOUND COMMENT BLOCK:')
        comment_lines = list(range(comment_block[0], comment_block[1]))
        
        if len(blocks[comment_block]) > 1:
            print('DOUBLE:', comment_block, blocks[comment_block])
        for paragraph in blocks[comment_block]:
            assert(paragraph[0] >= comment_block[0])
            assert(paragraph[1] <= comment_block[1])

            print('PARAGRAPH:')
            print('\n'.join(lines[i] for i in range(paragraph[0], paragraph[1])))
            
            for i in range(paragraph[0], paragraph[1]):
                comment_lines.remove(i)
            if paragraph[1] in comment_lines and (lines[paragraph[1]].strip() == '' or lines[paragraph[1]].strip() == '//'):
                comment_lines.remove(paragraph[1])

        remaining_comment_block = [lines[l] for l in comment_lines]
        print('REMAINING:')
        print('\n'.join(remaining_comment_block))
        if is_empty_block(remaining_comment_block):
            print(' IS EMPTY')
            # kill the whole block
            comment_lines = []

        # need to get the ones removed
        for line in range(comment_block[0], comment_block[1]):
            if line not in comment_lines:
                lines_to_remove.append(line)
                
    # now remove the lines, in reverse order so the numbers stay right
    new_lines = lines[:]
    for line in reversed(sorted(lines_to_remove)):
        new_lines.pop(line)
    return new_copyright, new_lines


def find_oneline_doc(lines, token='//!'):
    """Searches for one-line docstrings that start with //!"""
    print('  found one-line docstring')
    onelines = [(i,l) for (i,l) in enumerate(lines) if l.strip().startswith(token)]

    if len(onelines) == 0:
        return None, lines
    elif (len(onelines) == 1):
        lines.pop(onelines[0][0])
        oneline = strip_words(onelines[0][1].strip(), [token,]).strip()
        return '//! '+oneline, lines
    else:
        # pop in reverse to avoid modification issues
        for oneline in reversed(onelines):
            lines.pop(oneline[0])
        oneline = ' '.join([strip_words(l[1].strip(), [token,]).strip() for l in onelines])
        return '//! '+oneline, lines


def clang_format(filename):
    """Runs clang-format, using a style-file"""
    print('  running clang-format')
    return subprocess.run(['clang-format', '-i', '-style=file', filename])


def clean(filename, ats=False):
    with open(filename, 'r') as fid:
        lines = fid.readlines()
    lines = [l.rstrip() for l in lines]

    print("      ", type(lines), len(lines))
    oneline_doc, lines = find_oneline_doc(lines)
        
    print("      ", type(lines), len(lines))
    lines = clean_opening_whitespace(lines)

    print("      ", type(lines), len(lines))
    lines = clean_modeline(lines)

    print("      ", type(lines), len(lines))
    lines = clean_opening_whitespace(lines)

    print("      ", type(lines), len(lines))
    new_copyright, lines = find_remove_copyright(lines, ats)
 
    print("      ", type(lines), len(lines))
    lines = clean_opening_whitespace(lines)

    print("      ", type(lines), len(lines))
    lines = clean_trailing_whitespace(lines)

    print("      ", type(lines), len(lines))

    new_lines = new_copyright
    if oneline_doc is not None:
        new_lines.append(oneline_doc)
    new_lines.extend(lines)

    # one trailing empty line to force a newline on the last line of code
    new_lines.append('')
    
    with open(filename, 'w') as fid:
        fid.write('\n'.join(new_lines))

    #clang_format(filename)


def get_args():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--exclude', type=str, action='append',
                      help='Pattern to exclude from cleaning.')
    parser.add_argument('--extension', type=str, nargs='+', action='append',
                      help='File extensions to clean.')
    parser.add_argument('--ats', action='store_true',
                        help='On ATS source code')
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
        clean(fname, args.ats)
