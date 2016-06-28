#
#  IndexWriter:
#
#
#  section  -  a directory is copied and entry in the top-level index is added
#  subsection - a directory is copied and an entry is added in the corresponding section-level index
#  pages    -  a single page and accompanying files get copied and no index entries are added
#
#  Create Files:
#  -------------
#
#  - section or subsection directories that do not exist
#  - create directories
#  - create entry in an existing index.rst 
#  - add its own index.rst
#
#  Copy Files:
#  -------------
#
#  - from directory
#  - dest directory
#  - exclude files or directories 
#  - only (option to grab only certain files and/or subdirectories)
#
#  Index pages
#  -------------
#
#  - location of index
#  - entry and position?  
#

### Missing/incomplete:
###
###  Insert required lines in top-level index automatically
###  Correct paths in the .. plot:: directives

import shutil, distutils.dir_util, os, sys, glob

def IndexCreate(amanzi_home,data,logfile):

    index_data=data['index']
    local_path=os.path.dirname(index_data['index_file'])
    full_path=amanzi_home+os.sep+local_path

    logfile.write('  %s\n' % (amanzi_home+os.sep+index_data['index_file']) )
    
    #
    #  Create path if it doesn't exist
    #
    try: 
        os.makedirs(full_path)
    except OSError:
        if not os.path.isdir(full_path):
            raise

    #
    #  Overwrite and we don't care?
    #
    f=open(amanzi_home+os.sep+index_data['index_file'],"w")
    
    #
    # Reference
    #
    if ( 'index_label' in index_data ):
        f.write('.. _%s:\n\n' % (index_data['index_label']))

    #
    # Title
    #
    f.write('%s\n' % (index_data['index_title']))
    f.write('%s\n\n' % ('='*len(index_data['index_title'])))

    #
    # Table of Contents
    #
    f.write('%s\n' % ('.. toctree::') )
    f.write('   %s\n\n' % (':maxdepth: 2'))

    #
    # Add entries from index_list
    #
    for le in index_data['index_list']:
        f.write('   %s\n' % data[le]['index_entry'] )

    f.close()
            
def IndexInsert(index_file,new_lines):

    # make a backup
    shutil.copyfile(index_file,index_file+"~")
    # add new lines after the toctree directive
    infile=open(index_file+"~")
    outfile=open(index_file,"w")

    skip_toctree_options=False
    old_entries_skipped=False
    new_lines_inserted=False
    for line in infile:
        if "toctree::" in line:
            skip_toctree_options=True
            outfile.write(line)
        elif ( skip_toctree_options and line.strip() ):
            outfile.write(line)
        elif ( skip_toctree_options and not line.strip() ):
            outfile.write(line)
            skip_toctree_options=False
            for nline in new_lines:
                outfile.write("   "+nline+"\n")
            new_lines_inserted=True
        elif ( new_lines_inserted and not old_entries_skipped ):
            if ( not line.strip() ):
                outfile.write("\n")
                old_entries_skipped=True
        elif ( old_entries_skipped ):
            outfile.write(line)
        else:
            outfile.write(line)
       
               
def RecurseIndex(amanzi_home,content,level,logfile):

    if not isinstance(content,dict):

        level=level-1

        return

    else:
        
        if ( level == 1 ):
            logfile.write('\n%s\n\n' % ("Creating index files with Table of Contents") )

        if ('index' in content.keys()):
            # Create the file (overwite existing)
            IndexCreate(amanzi_home,content,logfile)
        
        level=level+1

        for c in content:
            RecurseIndex(amanzi_home,content[c],level,logfile)

        level=level-1

        return


def RecurseCopy(amanzi_home,content,level,logfile):

    if not isinstance(content,dict):

        return

    elif ('from_dir' in content.keys()):
        
        content_from=amanzi_home+os.sep+content['from_dir']
        content_dest=amanzi_home+os.sep+content['dest_dir']
        logfile.write('  %s\n  %s\n' % (content_from, content_dest) )
        distutils.dir_util.copy_tree(content_from, content_dest)

        level=level-1

    elif ('from_file' in content.keys()):
        
        content_from=amanzi_home+os.sep+content['from_file']
        content_dest=amanzi_home+os.sep+content['dest_file']
        if not os.path.exists(os.path.dirname(content_dest)):
            os.makedirs(os.path.dirname(content_dest))

        logfile.write('  %s\n  %s\n' % (content_from, content_dest) )
        shutil.copyfile(content_from, content_dest)

        level=level-1
        
    else:

        if ( level == 1):
            logfile.write('\n%s\n\n' % ("Copying input/script/rst files") )

        level=level+1
        
        for c in content:

            if ( level == 2 ):
                logfile.write('%s  %s\n' % ("Section:", c) )

            RecurseCopy(amanzi_home,content[c],level,logfile)

        level=level-1

    return


def AmanziHome(logfile):

    CWD = os.getcwd()
    CWD_list=CWD.split(os.sep)
    amanzi_home=os.sep.join(CWD_list[0:len(CWD_list)-2])

    logfile.write("\n%s  %s\n" % ("Amanzi Home Directory:", amanzi_home) )

    return amanzi_home

def UpdatePlotPath(amanzi_home,rst_dir,rst_base,logfile):

    rst_file=os.path.join(rst_dir,rst_base)
    
    # Ensure this will work with a full path as well as local path
    k = rst_dir.rfind("doc"+os.sep+"user_guide"+os.sep)
    if ( k == -1 ):
        plt_path=rst_dir
    else:
        plt_path = rst_dir[k+15:len(rst_dir)]

    logfile.write("  %s\n" % rst_file )

    # make a backup
    shutil.copyfile(rst_file,rst_file+"~")
    # add new lines after the toctree directive
    infile=open(rst_file+"~")
    outfile=open(rst_file,"w")

    for line in infile:
        if (".. plot::" in line):
            plt_py=line.rstrip().split(":: ")
            newline=plt_py[0] + ":: " + plt_path + os.sep + os.path.basename(plt_py[1]) + "\n"
            outfile.write(newline)
        else:
            outfile.write(line)

    # remove the backup
    os.remove(rst_file+"~")

    return

def WalkRstFiles(amanzi_home,content,logfile):

    logfile.write("\nFix plot directive script paths \n\n")

    for c in content:
        
        logfile.write("Section : %s\n" % c )

        if os.path.isdir(c):
            for dirpath, dirnames, filenames in os.walk(c):
                for filename in filenames:
                    # Fix directives in rst files (none are in the index pages).
                    if ( filename.endswith(".rst") and filename != "index.rst" ):
                        UpdatePlotPath(amanzi_home,dirpath,filename,logfile)                        
                        


    return

if __name__ == "__main__":

    # Data to test User Guide utility functions

    # Create a dictionary for the verification directory 
    verification={}
    verification['index']={'index_title':'Verification',
                           'index_file':'doc/user_guide/verification/index.rst',
                           'index_list':['confined_flow','unconfined_flow','transport'],
                       }
    
    verification['confined_flow']={'index_entry' : 'confined_flow/index.rst',
                                   'index' : 
                                   {'index_title' : 'Confined Flow Tests',
                                    'index_file' : 'doc/user_guide/verification/confined_flow/index.rst',
                                    'index_list' : ['linear_head_head', 'linear_flux_head'],
                                    },
                                   'linear_head_head' :
                                   {'from_dir' : 'testing/verification/flow/saturated/steady-state/linear_head_head_1d',
                                    'dest_dir' : 'doc/user_guide/verification/confined_flow/linear_head_head',
                                    'index_entry' : 'linear_head_head/amanzi_linear_head_head_1d.rst'
                                    },
                                   'linear_flux_head' :
                                   {'from_dir' : 'testing/verification/flow/saturated/steady-state/linear_flux_head_1d',
                                    'dest_dir' : 'doc/user_guide/verification/confined_flow/linear_flux_head',
                                    'index_entry' : 'linear_head_head/amanzi_linear_flux_head_1d.rst'
                                    },
                               }


    verification['unconfined_flow']={'index_entry': 'unconfined_flow/index.rst'}
    verification['transport']={'index_entry': 'transport/index.rst'}


    # Tutorials
    tutorial={}
    tutorial['index']={'index_title' : 'Tutorial',
                        'index_file' : 'doc/user_guide/tutorial/index.rst',
                        'index_list' : ['flow_1','flow_2'],
                   }
    
    tutorial['flow_1']={'from_dir' : 'doc/tutorial/saturated',
                        'dest_dir' : 'doc/user_guide/tutorial/steady_confined',
                        'index_entry' : 'steady_confined/Tutorial_steady_linear.rst',
                    }
    
    tutorial['flow_2']={'from_dir' : 'doc/tutorial/richards',
                        'dest_dir' : 'doc/user_guide/tutorial/transient_infiltration',
                        'index_entry' : 'transient_infiltration/Tutorial_dvz.rst',
                    }



    level=1
    sections={'tutorial' : tutorial,
              'verification' : verification,
          }


    logfile=sys.stdout

    # Set Amanzi source directory
    amanzi_home=AmanziHome(logfile)

    RecurseIndex(amanzi_home,sections,level,logfile)
    RecurseCopy(amanzi_home,sections,level,logfile)

    # Test Plot Path correction
    rsttest=amanzi_home+os.sep+"doc/user_guide/verification/confined_flow/"+verification['confined_flow']['linear_head_head']['index_entry']
    rst_dir = os.path.dirname(rsttest)
    rst_base = os.path.basename(rsttest)
    UpdatePlotPath(amanzi_home,rst_dir,rst_base,logfile)

    # Apply it to all the files
    WalkRstFiles(amanzi_home,sections,logfile)
    
