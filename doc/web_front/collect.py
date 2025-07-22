#!/usr/bin/env python

import os
import shutil
import optparse
import fnmatch

# Support routine to find files matching a pattern
def find(pattern, path):
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return result



# Guides dictionary defining where the developer_guide.pdf should be copied
guides = {}
guides['index'] = {
    'index_title': 'Guides',
    'index_file': 'guides/index.rst',
    'index_list': ['developer_guide', 'theory_guide', 'input_spec_s', 'input_spec_u']
}

guides['developer_guide'] = {
    'from_file': '../developer_guide/developer_guide.pdf',
    'dest_file': 'guides/developer_guide.pdf',
    'index_entry': 'developer_guide.pdf'
}

guides['theory_guide'] = {
    'from_file': '../theory_guide/theory.pdf',
    'dest_file': 'guides/theory.pdf',
    'index_entry': 'theory.pdf'
}

guides['input_spec'] = {
    'from_dir': '../input_spec/build/html',
    'dest_dir': 'guides/input_spec',
    'index_entry': 'AmanziInputSpec'
}

# =======================================================================================

# Create parser and options
p = optparse.OptionParser()
p.add_option('--full-guide', help='Build the full User Guide', default=False, dest='full_guide', action='store_true')

(opts, args) = p.parse_args()

# Function to copy the developer guide PDF
def copy_developer_guide():
    src = guides['developer_guide']['from_file']
    dest = guides['developer_guide']['dest_file']
    
    # Check if source file exists before copying
    if os.path.exists(src):
        # Ensure the destination directory exists
        os.makedirs(os.path.dirname(dest), exist_ok=True)
        
        # Copy the file from src to dest
        shutil.copyfile(src, dest)
        print(f"Copied {src} to {dest}")
    else:
        print(f"\nSource file {src} does not exist! Not copied to {dest}\n")

# Function to copy the theory guide PDF
def copy_theory_guide():
    src = guides['theory_guide']['from_file']
    dest = guides['theory_guide']['dest_file']
    
    # Check if source file exists before copying
    if os.path.exists(src):
        # Ensure the destination directory exists
        os.makedirs(os.path.dirname(dest), exist_ok=True)
        
        # Copy the file from src to dest
        shutil.copyfile(src, dest)
        print(f"Copied {src} to {dest}")
    else:
        print(f"\nSource file {src} does not exist! Not copied to {dest}\n")

# Function to copy the input guide directory
def copy_input_guide():
    src = guides['input_spec']['from_dir']
    dest = guides['input_spec']['dest_dir']
        # Check if source direc exists before copying
    if os.path.exists(src):
        # Ensure the destination directory exists
        os.makedirs(os.path.dirname(dest), exist_ok=True)
        
        # Copy the file from src to dest
        shutil.copytree(src, dest)
        print(f"Copied {src} to {dest}")
    else:
        print(f"\nSource directory {src} does not exist! Not copied to {dest}\n")

# Main execution starts here
if opts.full_guide:
    # Call function to copy the developer guide if the --full-guide option is passed
    copy_developer_guide()
    copy_theory_guide()
    copy_input_guide()
