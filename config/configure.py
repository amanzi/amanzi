#!/usr/bin/env python

# This script uses the BuildSystem classes (originally developed by Matt Knepley for configuring/building PETSc) to configure 
# Amanzi and its third party libraries.
# The ConfigurationManager class defined here is largely borrowed from Matt's "simpleConfigure" example.  Much thanks to Matt.

import os
import sys

class ConfigurationManager(object):
  def __init__(self, project = 'Amanzi'):
    self.Project = project
    self.PROJECT = project.upper()
    self.checkPythonVersion()
    self.normalizeLanguage()
    self.setupLog()
    self.setupBuildSystem()
    return

  def checkPythonVersion(self):
    if not hasattr(sys, 'version_info') or not sys.version_info[0] == 2 or not sys.version_info[1] >= 3:
      print('*** You must have Python2 version 2.3 or higher to run ./configure        *****')
      print('*          Python is easy to install for end users or sys-admin.              *')
      print('*                  http://www.python.org/download/                            *')
      print('*******************************************************************************')
      sys.exit(4)
    return

  def normalizeLanguage(self):
    '''Use en_US as language so that BuildSystem parses compiler messages in english'''
    if 'LC_LOCAL' in os.environ and os.environ['LC_LOCAL'] != '' and os.environ['LC_LOCAL'] != 'en_US' and os.environ['LC_LOCAL']!= 'en_US.UTF-8':
      os.environ['LC_LOCAL'] = 'en_US.UTF-8'
    if 'LANG' in os.environ and os.environ['LANG'] != '' and os.environ['LANG'] != 'en_US' and os.environ['LANG'] != 'en_US.UTF-8':
      os.environ['LANG'] = 'en_US.UTF-8'
    return

  def untar(self, tar, path = '.', leading = ''):
    if leading:
      entries = [t.name for t in tar.getmembers()]
      prefix = os.path.commonprefix(entries)
      if prefix:
        for tarinfo in tar.getmembers():
          tail = tarinfo.name.split(prefix, 1)[1]
          tarinfo.name = os.path.join(leading, tail)
    for tarinfo in tar.getmembers():
      tar.extract(tarinfo, path)
    return

  def downloadPackage(self, url, filename, targetDirname):
    '''Download the tarball for a package at url, save it as filename, and untar it into targetDirname'''
    import tarfile, urllib.request, urllib.parse, urllib.error
    filename, headers = urllib.request.urlretrieve(url, filename)
    tar = tarfile.open(filename, 'r:gz')
    self.untar(tar, targetDirname, leading = filename.split('.')[0])
    return

  def validateOptions(self, opts):
    '''This should be overridden for PETSc check_for_option_mistakes()+check_petsc_arch()'''
    return

  def normalizeOptions(self, opts):
    '''Support a few standard configure option types'''
    for l in range(0, len(opts)):
      name = opts[l]
      if name.find('enable-') >= 0:
        if name.find('=') == -1:
          opts[l] = name.replace('enable-','with-')+'=1'
        else:
          head, tail = name.split('=', 1)
          opts[l] = head.replace('enable-','with-')+'='+tail
      if name.find('disable-') >= 0:
        if name.find('=') == -1:
          opts[l] = name.replace('disable-','with-')+'=0'
        else:
          head, tail = name.split('=', 1)
          if tail == '1': tail = '0'
          opts[l] = head.replace('disable-','with-')+'='+tail
      if name.find('without-') >= 0:
        if name.find('=') == -1:
          opts[l] = name.replace('without-','with-')+'=0'
        else:
          head, tail = name.split('=', 1)
          if tail == '1': tail = '0'
          opts[l] = head.replace('without-','with-')+'='+tail
    return

  def setupLog(self):
    '''Sometime symlinks can get broken if the original files are deleted. Delete such broken links'''
    import os
    for logfile in ['configure.log','configure.log.bkp']:
      if os.path.islink(logfile) and not os.path.isfile(logfile): os.remove(logfile)
    return

  def getArch(self, framework, opts):
    import os
    arch = ''
    if hasattr(framework, 'arch'): return framework.arch
    # Check PROJECT_ARCH.
    # We check to see if it is set on the command line, then we check the environment variables, and, finally, we use the 
    # filename of the configure script if it is not 'configure.py'.
    found = 0
    for name in opts:
      if name.find(self.PROJECT+'_ARCH=') >= 0:
        arch  = name.split('=')[1]
        found = 1
        break
    if not found:
      arch = os.environ.get(self.PROJECT+'_ARCH')
      if arch is not None:
        found = 1
      else:
        arch = ''
    # Use the filename of script, if we have been unable to determine this from PROJECT_ARCH.
    if not found:
      filename = os.path.basename(sys.argv[0])
      if not filename.startswith('configure') and not filename.startswith('reconfigure'):
        arch = os.path.splitext(os.path.basename(sys.argv[0]))[0]
        opts.append(self.PROJECT+'_ARCH='+arch)
      else:
        # Unable to determine any suitable value for 'arch' attribute.
        msg ='*******************************************************************************\n'\
            +'UNABLE to DETERMINE AMANZI_ARCH.  Please set AMANZI_ARCH and re-run configure.\n' \
            +'*******************************************************************************\n'
        print(msg)
        if framework is not None:
          framework.logClear()
          import traceback
          if hasattr(framework, 'log'):
            try:
              framework.log.write(msg)
              traceback.print_tb(sys.exc_info()[2], file = framework.log)
              framework.log.close()
              self.cleanupLog(framework)
            except:
              pass
          else:
            traceback.print_tb(sys.exc_info()[2])
        sys.exit(1)
    return arch

  def cleanupLog(self, framework):
    '''Move configure.log to PROJECT_ARCH/conf - and update configure.log.bkp in both locations appropriately'''
    arch    = self.getArch(framework, sys.argv)
    logFile = 'configure.log'
    if hasattr(framework, 'logName'): logFile = framework.logName

    if arch:
      import shutil
      import os

      confDir = os.path.join(arch, 'conf')
      if not os.path.isdir(arch):    os.mkdir(arch)
      if not os.path.isdir(confDir): os.mkdir(confDir)

      logFileBkp        = logFile + '.bkp'
      logFileArchive    = os.path.join(confDir, logFile)
      logFileArchiveBkp = logFileArchive + '.bkp'

      # Keep backup in $PROJECT_ARCH/conf location
      if os.path.isfile(logFileArchiveBkp): os.remove(logFileArchiveBkp)
      if os.path.isfile(logFileArchive):    os.rename(logFileArchive, logFileArchiveBkp)
      if os.path.isfile(logFile):
        shutil.copyfile(logFile, logFileArchive)
        os.remove(logFile)
      if os.path.isfile(logFileArchive):    os.symlink(logFileArchive, logFile)
      # If the old bkp is using the same $PROJECT_ARCH/conf, then update bkp link
      if os.path.realpath(logFileBkp) == os.path.realpath(logFileArchive):
        if os.path.isfile(logFileBkp):        os.remove(logFileBkp)
        if os.path.isfile(logFileArchiveBkp): os.symlink(logFileArchiveBkp, logFileBkp)
    return

  def setupBuildSystem(self):
    import subprocess

    self.configDir = os.path.abspath('config')
    self.bsDir     = os.path.join(self.configDir, 'BuildSystem')
    bsURL = 'http://petsc.cs.iit.edu/petsc/BuildSystem'
    hgCommand = 'hg clone ' + bsURL +' '+ self.bsDir
    if not os.path.isdir(self.configDir):
      raise RuntimeError('Run configure from $'+self.PROJECT+'_DIR, not '+os.path.abspath('.'))
    # Try to clone BuildSystem via Mercurial; failing that, we then download the tarball.
    if not os.path.isdir(self.bsDir):
      print('===============================================================================')
      print('''++ Could not locate BuildSystem in %s.''' % self.configDir)
      print('''++ Downloading it from %s.''' % bsURL)
      print('''++ Attempting: %s''' % hgCommand)
      (status,output) = subprocess.getstatusoutput(hgCommand)
      if status:
        print('++ Unable to clone BuildSystem.  Attempting download from ' + bsURL + '/archive/tip.tar.gz')
        self.downloadPackage('http://petsc.cs.iit.edu/petsc/BuildSystem/archive/tip.tar.gz', 'BuildSystem.tar.gz', self.configDir)
      print('===============================================================================')
    # to load ~/.pythonrc.py before inserting correct BuildSystem to path
    import user
    sys.path.insert(0, self.bsDir)
    sys.path.insert(0, self.configDir)
    return

  def checkCygwin(self, opts):
    '''Cygwin 1.5.11-1 is broken
       You cannot use Windows Python with Cygwin
       Cygwin-Python-2.4/2.5/2.6 does not have working threads
       Cygwin linker messes up Intel F90'''
    logMsg = []
    if os.path.exists('/usr/bin/cygcheck.exe'):
      buf = os.popen('/usr/bin/cygcheck.exe -c cygwin').read()
      if buf.find('1.5.11-1') > -1:
        print('===============================================================================')
        print(' *** cygwin-1.5.11-1 detected. ./configure fails with this version ***')
        print(' *** Please upgrade to cygwin-1.5.12-1 or newer version. This can  ***')
        print(' *** be done by running cygwin-setup, selecting "next" all the way.***')
        print('===============================================================================')
        sys.exit(3)
      if sys.platform != 'cygwin':
        print('===============================================================================')
        print(' *** Non-cygwin python detected. Please rerun ./configure **')
        print(' *** with cygwin-python. ***')
        print('===============================================================================')
        sys.exit(3)
      buf = os.popen('/usr/bin/cygcheck.exe -c python').read()
      if (buf.find('2.4') > -1) or (buf.find('2.5') > -1) or (buf.find('2.6') > -1):
        opts.append('--useThreads=0')
        logMsg.append('''\
===============================================================================
** Cygwin-python-2.4/2.5/2.6 detected. Threads do not work correctly with this
** version. Disabling thread usage for this run of ./configure *******
===============================================================================''')
      if os.path.exists('/usr/bin/link.exe'):
        if '--ignore-cygwin-link' in opts: return 0
        for arg in opts:
          if (arg.find('win32fe') >= 0 and (arg.find('f90') >=0 or arg.find('ifort') >=0)):
            print('===============================================================================')
            print(' *** Cygwin /usr/bin/link detected! Compiles with CVF/Intel f90 can break!  **')
            print(' *** To workarround do: "mv /usr/bin/link.exe /usr/bin/link-cygwin.exe"     **')
            print(' *** Or to ignore this check, use configure option: --ignore-cygwin-link    **')
            print('===============================================================================')
            sys.exit(3)
    return logMsg

  def checkRedhat(self, opts):
    '''Disable threads on RHL9'''
    logMsg = []
    if os.path.exists('/etc/redhat-release'):
      try:
        f = open('/etc/redhat-release','r')
        buf = f.read()
        f.close()
      except:
        # can't read file - assume dangerous RHL9
        buf = 'Shrike'
      if buf.find('Shrike') > -1: 
        opts.append('--useThreads=0')
        logMsg.append('''\
==============================================================================
   *** RHL9 detected. Threads do not work correctly with this distribution ***
   ****** Disabling thread usage for this run of ./configure *********
===============================================================================''')
    return logMsg

  def systemKludges(self, opts):
    logMsg = []
    logMsg.extend(self.checkCygwin(opts))
    logMsg.extend(self.checkRedhat(opts))
    return logMsg

  def configure(self, configure_options):
    print('===============================================================================')
    print('             Configuring '+self.Project+' to compile on your system                       ')
    print('===============================================================================')  

    # Arguments passed in take precedence (but don't destroy argv[0])
    sys.argv = sys.argv[:1] + configure_options + sys.argv[1:]
    self.validateOptions(sys.argv)
    self.normalizeOptions(sys.argv)
    logMessages = self.systemKludges(sys.argv)

    import config.base
    import config.framework
    import pickle

    framework = None
    arch = self.getArch(framework, sys.argv)

    try:
      framework = config.framework.Framework(['--configModules='+self.Project+'.Configure','--optionsModule='+self.Project+'.compilerOptions']+sys.argv[1:], loadArgDB = 0)
      framework.arch = arch
      framework.setup()
      framework.logPrint('\n'.join(logMessages))
      framework.configure(out = sys.stdout)
      framework.storeSubstitutions(framework.argDB)
      framework.argDB['configureCache'] = pickle.dumps(framework)
      for i in framework.packages:
        if hasattr(i,'postProcess'):
          i.postProcess()
      framework.printSummary()
      framework.logClear()
      framework.closeLog()
      try:
        self.cleanupLog(framework)
      except:
        # perhaps print an error about unable to shuffle logs?
        pass
      return 0
    except (RuntimeError, config.base.ConfigureSetupError) as e:
      emsg = str(e)
      if not emsg.endswith('\n'): emsg = emsg+'\n'
      msg ='*******************************************************************************\n'\
          +'         UNABLE to CONFIGURE with GIVEN OPTIONS    (see configure.log for details):\n' \
          +'-------------------------------------------------------------------------------\n'  \
          +emsg+'*******************************************************************************\n'
      se = ''
    except (TypeError, ValueError) as e:
      emsg = str(e)
      if not emsg.endswith('\n'): emsg = emsg+'\n'
      msg ='*******************************************************************************\n'\
          +'                ERROR in COMMAND LINE ARGUMENT to ./configure \n' \
          +'-------------------------------------------------------------------------------\n'  \
          +emsg+'*******************************************************************************\n'
      se = ''
    except ImportError as e :
      emsg = str(e)
      if not emsg.endswith('\n'): emsg = emsg+'\n'
      msg ='*******************************************************************************\n'\
          +'                     UNABLE to FIND MODULE for ./configure \n' \
          +'-------------------------------------------------------------------------------\n'  \
          +emsg+'*******************************************************************************\n'
      se = ''
    except OSError as e :
      emsg = str(e)
      if not emsg.endswith('\n'): emsg = emsg+'\n'
      msg ='*******************************************************************************\n'\
          +'                    UNABLE to EXECUTE BINARIES for ./configure \n' \
          +'-------------------------------------------------------------------------------\n'  \
          +emsg+'*******************************************************************************\n'
      se = ''
    except SystemExit as e:
      if e.code is None or e.code == 0:
        return
      msg ='*******************************************************************************\n'\
          +'         CONFIGURATION FAILURE  (Please send configure.log to the Amanzi developers)\n' \
          +'*******************************************************************************\n'
      se  = str(e)
    except Exception as e:
      msg ='*******************************************************************************\n'\
          +'        CONFIGURATION CRASH  (Please send configure.log to the Amanzi developers)\n' \
          +'*******************************************************************************\n'
      se  = str(e)
    print(msg)
    if framework is not None:
      framework.logClear()
      if hasattr(framework, 'log'):
        import traceback
        try:
          framework.log.write(msg+se)
          traceback.print_tb(sys.exc_info()[2], file = framework.log)
          framework.log.close()
          self.cleanupLog(framework)
        except:
          pass
        sys.exit(1)
    else:
      print(se)
      import traceback
      traceback.print_tb(sys.exc_info()[2])
    return

if __name__ == '__main__':
  ConfigurationManager(project = 'Amanzi').configure([])
