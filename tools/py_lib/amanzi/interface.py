# ############################################################################ #

import os, sys


from amanzi import input_tree

# ############################################################################ #

class AmanziDataOutput:

  CYCLE_REGEX_PATTERN=r'.h5.(\d+).xmf'
  XMF_TIME_TAG='Time'
  XMF_TIME_VALUE_ATTR='Value'
  XMF_DATASET_TAG='Attribute'
  XMF_DATASET_NAME_ATTR='Name'

  def __init__(self,xmf_filename,directory=None):

    from xml.etree import ElementTree as ETree
    import re
   
    
    # Set to current directory if not set
    try:
      self.output_dir=os.path.abspath(directory)
    except:
      self.output_dir=os.getcwd()

    # Set filename to the full path
    self.filename=os.path.abspath(self.output_dir+os.sep+xmf_filename)

    # Define the cycle
    base_file=os.path.basename(self.filename)
    self.cycle=int(re.search(self.CYCLE_REGEX_PATTERN,base_file).group(1))

    # Search for the XMF file need to open it with ElementTree
    try:
      tree=ETree.parse(self.filename)
    except IOError:
      print('Failed to open:' + self.filename)
      raise
    except:
      raise

    # Find the time
    time_elem=tree.getiterator(tag=self.XMF_TIME_TAG)[0]
    self.time=float(time_elem.get(self.XMF_TIME_VALUE_ATTR))

    # Now find all the datasets
    self.datasets = []
    dataset_iter=tree.getiterator(tag=self.XMF_DATASET_TAG)
    for elem in tree.getiterator(tag=self.XMF_DATASET_TAG):
      self.datasets.append(elem.get(self.XMF_DATASET_NAME_ATTR))

    

class AmanziInterface:

  AMANZI_CHECKPOINT_REGEX=r'checkpoint\d+.h5'
  AMANZI_DATAFILE_REGEX=r'_data.h5'

  def __init__(self,binary=None,input=None,output=None,error=None,mpiexec='mpiexec'):

    self.binary=binary
    self.mpiexec=mpiexec
    self.input=input
    self.output=output
    self.error=error

    self.nprocs=0
    self.exit_code=None
    self.args=[]

    self.data_files=[] 


  def run(self):
    import time
    import signal
    try:
      import subprocess
    except ImportError:
      print('Failed to import module subprocess')
      print('Module is avail in Python 2.4 and higher')
      print('Current Python version ' + str(sys.version))
      raise
    except:
      print("Unknown error:", sys.exc_info())
      raise

    # Check the files
    if not self._check_binary():
      raise ValueError('Can not run binary')

    if not self._check_input():
      raise ValueError('Invalid input')

    if self.nprocs > 0:
      if not self._check_mpiexec:
        raise ValueError('Can not launch a parallel executable')

    # Open files to redirect STDOUT and STDERR
    try:
      stdout_fh = open(self.output,'w')
    except (ValueError,TypeError):
      stdout_fh=None
    except IOError:
      raise IOError('Failed to open ' + self.output)
    except:
      print("Unknown exception:", sys.exc_info()[0])
      raise

    try:
      stderr_fh = open(self.error,'w')
    except (ValueError,TypeError):
      stderr_fh=None
    except IOError:
      raise IOError('Failed to open ' + self.error)
    except:
      print('Unknown exception')
      raise

    # Define the arguments and execuables
    if self.nprocs == 0:
      executable=self.binary
      args=[executable, self.input_option()]
    else:
      executable=self.mpiexec
      args=[executable, '-n', str(self.nprocs), self.binary, self.input_option()]
  
    # Open a pipe, catch a CTRL-C shut down the pipe
    try: 
      pipe = subprocess.Popen(args,executable=executable,bufsize=-1,stdout=stdout_fh,stderr=stderr_fh)
    except ValueError:
      raise ValueError('Popen called with incorrect arguments')
    except OSError:
      print('Failed to run ' + str(args))
      raise OSError('Failed to run Amanzi binary ' + self.binary)
    else:
      try:
          pipe.wait()
      except KeyboardInterrupt:
          try:
              kill_signal=signal.SIGKILL
          except AttributeError:
              kill_signal=signal.SIGTERM
              print('Sending SIGTERM to Amanzi run PID='+str(pipe.pid))
              os.kill(pipe.pid,signal.SIGTERM)
              time.sleep(1)
              while pipe.poll() == None:
                  print('PID='+str(pipe,pid)+' still alive. Sending signal (SIGKILL) '+str(kill_signal))
                  os.kill(pipe.pid,kill_signal)
                  time.sleep(5)

    # Set the return code 
    self.exit_code=pipe.returncode

    # Flush and close the files
    try:
      stdout_fh.flush()
    except:
      pass
    else:
      stdout_fh.close()
    
    try:
      stderr_fh.flush()
    except:
      pass
    else:
      stderr_fh.close()

      

    return pipe.returncode

  def input_option(self):
    return '--xml_file='+self.input

  def _check_mpiexec(self):
    from os import path,access, R_OK, X_OK
    
    try:
      exists=os.path.exists(self.mpiexec)
    except TypeError:
      print('MPI launch command is not defined')
      raise
   
    ok_flag=False
    if exists and path.isfile(self.mpiexec) and access(self.mpiexec,R_OK|X_OK):
      ok_flag=True
    else:
      print(amanzi.mpiexec + ' does not exist or is not readable or is not executable')

    return ok_flag

  def _check_binary(self):
    from os import path,access, R_OK, X_OK
    
    try:
      exists=os.path.exists(self.binary)
    except TypeError:
      print('Amanzi binary is not defined')
      raise
   
    ok_flag=False
    if exists and path.isfile(self.binary) and access(self.binary,R_OK|X_OK):
      ok_flag=True
    else:
      print(self.binary + ' does not exist or is not readable')

    return ok_flag  

  def _check_input(self):
    from os import path,access, R_OK
    
    try:
      exists=os.path.exists(self.input)
    except TypeError:
      print('Input file is not defined')
      raise
   
    ok_flag=False
    if exists and path.isfile(self.input) and access(self.input,R_OK):
      ok_flag=True
    else:
      print(amanzi.input+ ' does not exist or is not readable')

    return ok_flag

  def find_data_files(self,basename,directory=None):
    import glob
    
    try:
      output_xmf_regex=basename+'_data.h5'+'.[0-9]*'+'.xmf'
    except: 
      print('Failed to build XMF output regular expression pattern')
      raise

    try:
      output_xmf_regex=directory+os.sep+output_xmf_regex
    except:
      pass

    print('output_xmf_regex='+output_xmf_regex)

    for output in glob.glob(output_xmf_regex):
     print('File='+output)
     self.data_files.append(AmanziDataOutput(output))

    return self.data_files

  def sort_data_files(self,files=[]):
    if len(files) == 0:
      if len(self.data_files) == 0:
        basename=self.plot_output_basename()
        self.find_data_files(basename)
      files=self.data_files

    search_dict = {}
    for data_file in files:
      search_dict[data_file.cycle]=data_file
  
    sorted_cycles=search_dict.keys()
    sorted_cycles.sort()
    sorted_files=[]
    for n in sorted_cycles:
      sorted_files.append(search_dict.get(n))

    return sorted_files

  def find_mesh_file(self,input=None):
    import glob
    basename = self.plot_output_basename(input)

    mesh_regex=basename+'_mesh.h5'
    try:
      mesh_file=glob.glob(mesh_regex)[0]
    except:
      print("Could not find mesh file matching \'" +mesh_regex+"\' pattern")
      mesh_file=None

    return mesh_file


  def plot_output_basename(self,input=None):

    from amanzi.trilinos import InputList as AmanziInput

    if input == None:
      input_tree = AmanziInput(self.input)
    else:
      input_tree = AmanziInput(input)

    output_ctrl = input_tree.find_sublist('Output')
    viz_ctrl    = output_ctrl.find_sublist('Visualization Data')

    return viz_ctrl.find_parameter('File Name Base').get_value()

if __name__ == '__main__':

  try:
    amanzi=AmanziInterface('duh.xml')
  except ValueError:
    print('Caught the invalid input error')
