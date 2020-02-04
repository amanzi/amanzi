import config.base

import os
import re

class AmanziDir(object):
  dir = os.getcwd()

class AmanziArch(object):
  arch = ''
  def __init__(self, framework=None, opts=None):
    import os
    import sys
    # Check AMANZI_ARCH.
    # We check to see if it is set on the command line, then we check the environment variables, and, finally, we use the 
    # filename of the configure script if it is not 'configure.py'.
    found = 0
    for name in opts:
      if name.find('AMANZI_ARCH=') >= 0:
        arch  = name.split('=')[1]
        found = 1
        break
    if not found:
      arch = os.environ.get('AMANZI_ARCH')
      if arch is not None:
        found = 1
      else:
        arch = ''
    # Use the filename of script, if we have been unable to determine this from PROJECT_ARCH.
    if not found:
      filename = os.path.basename(sys.argv[0])
      if not filename.startswith('configure') and not filename.startswith('reconfigure'):
        arch = os.path.splitext(os.path.basename(sys.argv[0]))[0]
        opts.append('AMANZI_ARCH='+arch)
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
    self.arch = arch

class Configure(config.base.Configure):
  def __init__(self, framework):
    import sys
    config.base.Configure.__init__(self, framework)
    self.Project      = 'Amanzi'
    self.project      = self.Project.lower()
    self.PROJECT      = self.Project.upper()
    self.headerPrefix = self.PROJECT
    self.substPrefix  = self.PROJECT
    self.defineAutoconfMacros()
    self.arch = AmanziArch(framework, sys.argv)
    # Fix these later
    self.dir  = AmanziDir()
    self.PACKAGES_INCLUDES = ''
    self.PACKAGES_LIBS     = ''
    return

  def __str2__(self):
    desc = []
    desc.append('xxx=========================================================================xxx')
    desc.append('   Configure stage complete. Now build '+self.Project+' libraries with:')
    desc.append('   make all')
    desc.append('xxx=========================================================================xxx')
    return '\n'.join(desc)+'\n'

  def setupHelp(self, help):
    import nargs
    help.addArgument(self.Project, '-prefix=<path>', nargs.Arg(None, '', 'Specifiy location to install '+self.Project+' (eg. /usr/local)'))
    help.addArgument(self.Project, '-with-shared', nargs.ArgBool(None, 0, 'Make libraries shared'))
    help.addArgument(self.Project, '-with-dynamic', nargs.ArgBool(None, 0, 'Make libraries dynamic'))
    return

  def setupDependencies(self, framework):
    config.base.Configure.setupDependencies(self, framework)
    self.setCompilers  = framework.require('config.setCompilers',      self)
    self.compilers     = framework.require('config.compilers',         self)
    self.types         = framework.require('config.types',             self)
    self.headers       = framework.require('config.headers',           self)
    self.functions     = framework.require('config.functions',         self)
    self.libraries     = framework.require('config.libraries',         self)

    # This is where we add all Amanzi TPL dependencies, via framework.require().
    # Note that we do this for all packages upon which Amanzi may depend, regardless of whether they are optional.  For optional 
    # packages, we determine if they actually need to be configured in the 'configure' method for that package.
    self.blaslapack    = framework.require('config.packages.BlasLapack', self)
    self.mpi           = framework.require('config.packages.MPI',        self)
    self.boost         = framework.require('config.packages.boost',      self)
    self.netcdf        = framework.require('config.packages.netcdf',     self)
    self.hdf5          = framework.require('config.packages.hdf5',       self)
    self.exodusii      = framework.require('config.packages.exodusii',   self)
    self.moab          = framework.require('config.packages.MOAB',       self)

    # Set some additional options for some of the packages:

    # Set the archProvider and installDirProvider for all packages.
    # There should probably be a more automatic way of doing this, at some point.

    self.blaslapack.archProvider        = self.arch
    self.blaslapack.installDirProvider  = self.dir
    self.mpi.archProvider               = self.arch
    self.mpi.installDirProvider         = self.dir
    self.boost.archProvider             = self.arch
    self.boost.installDirProvider       = self.dir
    self.netcdf.archProvider            = self.arch
    self.netcdf.installDirProvider      = self.dir
    self.hdf5.archProvider              = self.arch
    self.hdf5.installDirProvider        = self.dir
    self.exodusii.archProvider          = self.arch
    self.exodusii.installDirProvider    = self.dir
    self.moab.archProvider              = self.arch
    self.moab.installDirProvider        = self.dir

    # For BLAS/LAPACK:
    #force blaslapack to depend on scalarType so precision is set before BlasLapack is built
    #framework.require('PETSc.utilities.scalarTypes', self.blaslapack)
    #self.blaslapack.precisionProvider = self.scalartypes
    self.blaslapack.headerPrefix       = self.headerPrefix

    # For MPI:
    #self.mpi.languageProvider   = self.languages  #RTM: Do we need to set this?  Appears to depend on PETSc.utilities.languages
    self.mpi.headerPrefix       = self.headerPrefix

    self.compilers.headerPrefix = self.headerPrefix
    self.types.headerPrefix     = self.headerPrefix
    self.headers.headerPrefix   = self.headerPrefix
    self.functions.headerPrefix = self.headerPrefix
    self.libraries.headerPrefix = self.headerPrefix
    headersC = map(lambda name: name+'.h', ['malloc'])
    functions = ['drand48', 'getcwd']
    libraries1 = [(['socket', 'nsl'], 'socket'), (['fpe'], 'handle_sigfpes')]
    self.headers.headers.extend(headersC)
    self.functions.functions.extend(functions)
    self.libraries.libraries.extend(libraries1)
    return

  def defineAutoconfMacros(self):
    self.hostMacro = 'dnl Version: 2.13\ndnl Variable: host_cpu\ndnl Variable: host_vendor\ndnl Variable: host_os\nAC_CANONICAL_HOST'
    return
    
  def Dump(self):
    ''' Actually put the values into the configuration files '''
    # eventually everything between -- should be gone
#-----------------------------------------------------------------------------------------------------    

    # Sometimes we need C compiler, even if built with C++
    self.setCompilers.pushLanguage('C')
    self.addMakeMacro('CC_FLAGS',self.setCompilers.getCompilerFlags())    
    self.setCompilers.popLanguage()

    # C preprocessor values
    self.addMakeMacro('CPP_FLAGS',self.setCompilers.CPPFLAGS)
    
    # compiler values
    self.setCompilers.pushLanguage('Cxx')
    self.addMakeMacro('PCC',self.setCompilers.getCompiler())
    self.addMakeMacro('PCC_FLAGS',self.setCompilers.getCompilerFlags())
    self.setCompilers.popLanguage()
    # .o or .obj 
    self.addMakeMacro('CC_SUFFIX','o')

    # executable linker values
    self.setCompilers.pushLanguage('Cxx')
    pcc_linker = self.setCompilers.getLinker()
    self.addMakeMacro('PCC_LINKER',pcc_linker)
    self.addMakeMacro('PCC_LINKER_FLAGS',self.setCompilers.getLinkerFlags())
    self.setCompilers.popLanguage()
    # '' for Unix, .exe for Windows
    self.addMakeMacro('CC_LINKER_SUFFIX','')
    self.addMakeMacro('PCC_LINKER_LIBS',self.libraries.toStringNoDupes(self.compilers.flibs+self.compilers.cxxlibs+self.compilers.LIBS.split(' ')))

    if hasattr(self.compilers, 'FC'):
      self.setCompilers.pushLanguage('FC')
      # need FPPFLAGS in config/setCompilers
      self.addDefine('HAVE_FORTRAN','1')
      self.addMakeMacro('FPP_FLAGS',self.setCompilers.CPPFLAGS)
    
      # compiler values
      self.addMakeMacro('FC_FLAGS',self.setCompilers.getCompilerFlags())
      self.setCompilers.popLanguage()
      # .o or .obj 
      self.addMakeMacro('FC_SUFFIX','o')

      # executable linker values
      self.setCompilers.pushLanguage('FC')
      # Cannot have NAG f90 as the linker - so use pcc_linker as fc_linker
      fc_linker = self.setCompilers.getLinker()
      if config.setCompilers.Configure.isNAG(fc_linker):
        self.addMakeMacro('FC_LINKER',pcc_linker)
      else:
        self.addMakeMacro('FC_LINKER',fc_linker)
      self.addMakeMacro('FC_LINKER_FLAGS',self.setCompilers.getLinkerFlags())
      # apple requires this shared library linker flag on SOME versions of the os
      if self.setCompilers.getLinkerFlags().find('-Wl,-commons,use_dylibs') > -1:
        self.addMakeMacro('DARWIN_COMMONS_USE_DYLIBS',' -Wl,-commons,use_dylibs ')
      self.setCompilers.popLanguage()

      # F90 Modules
      if self.setCompilers.fortranModuleIncludeFlag:
        self.addMakeMacro('FC_MODULE_FLAG', self.setCompilers.fortranModuleIncludeFlag)
      else: # for non-f90 compilers like g77
        self.addMakeMacro('FC_MODULE_FLAG', '-I')
      if self.setCompilers.fortranModuleIncludeFlag:
        self.addMakeMacro('FC_MODULE_OUTPUT_FLAG', self.setCompilers.fortranModuleOutputFlag)
    else:
      self.addMakeMacro('FC','')

    # shared library linker values
    self.setCompilers.pushLanguage('Cxx')
    # need to fix BuildSystem to collect these separately
    self.addMakeMacro('SL_LINKER',self.setCompilers.getLinker())
    self.addMakeMacro('SL_LINKER_FLAGS','${PCC_LINKER_FLAGS}')
    self.setCompilers.popLanguage()
    # One of 'a', 'so', 'lib', 'dll', 'dylib' (perhaps others also?) depending on the library generator and architecture
    # Note: . is not included in this macro, consistent with AR_LIB_SUFFIX
    if self.setCompilers.sharedLibraryExt == self.setCompilers.AR_LIB_SUFFIX:
      self.addMakeMacro('SL_LINKER_SUFFIX', '')
      self.addDefine('SLSUFFIX','""')
    else:
      self.addMakeMacro('SL_LINKER_SUFFIX', self.setCompilers.sharedLibraryExt)
      self.addDefine('SLSUFFIX','"'+self.setCompilers.sharedLibraryExt+'"')
      
    #SL_LINKER_LIBS is currently same as PCC_LINKER_LIBS - so simplify
    self.addMakeMacro('SL_LINKER_LIBS','${PCC_LINKER_LIBS}')
    #self.addMakeMacro('SL_LINKER_LIBS',self.libraries.toStringNoDupes(self.compilers.flibs+self.compilers.cxxlibs+self.compilers.LIBS.split(' ')))

    if not os.path.exists(os.path.join(self.dir.dir,self.arch.arch,'lib')):
      os.makedirs(os.path.join(self.dir.dir,self.arch.arch,'lib'))

    # add a makefile entry for configure options
    self.addMakeMacro('CONFIGURE_OPTIONS', self.framework.getOptionsString(['configModules', 'optionsModule']).replace('\"','\\"'))
    return

  def dumpConfigInfo(self):
    import time
    if not os.path.isdir(os.path.join(self.arch.arch, 'include')):
      os.makedirs(os.path.join(self.arch.arch,'include'))
    fd = open(os.path.join(self.arch.arch,'include', self.project+'configinfo.h'), 'w')
    fd.write('static const char *'+self.project+'configureruntime = "'+time.ctime(time.time())+'";\n')
    fd.write('static const char *'+self.project+'configureoptions = "'+self.framework.getOptionsString(['configModules', 'optionsModule']).replace('\"','\\"')+'";\n')
    fd.close()
    return

  def dumpMachineInfo(self):
    import platform
    import time
    import script
    if not os.path.isdir(os.path.join(self.arch.arch, 'include')):
      os.makedirs(os.path.join(self.arch.arch,'include'))
    fd = open(os.path.join(self.arch.arch,'include', self.project+'machineinfo.h'),'w')
    fd.write('static const char *'+self.project+'machineinfo = \"\\n\"\n')
    fd.write('\"-----------------------------------------\\n\"\n')
    fd.write('\"Libraries compiled on %s on %s \\n\"\n' % (time.ctime(time.time()), platform.node()))
    fd.write('\"Machine characteristics: %s\\n\"\n' % (platform.platform()))
    fd.write('\"Using '+self.Project+' directory: %s\\n\"\n' % (self.dir.dir))
    fd.write('\"Using '+self.Project+' arch: %s\\n\"\n' % (self.arch.arch))
    fd.write('\"-----------------------------------------\\n\";\n')
    fd.write('static const char *'+self.project+'compilerinfo = \"\\n\"\n')
    self.setCompilers.pushLanguage('Cxx')
    fd.write('\"Using C compiler: %s %s ${COPTFLAGS} ${CFLAGS}\\n\"\n' % (self.setCompilers.getCompiler(), self.setCompilers.getCompilerFlags()))
    self.setCompilers.popLanguage()
    if hasattr(self.compilers, 'FC'):
      self.setCompilers.pushLanguage('FC')
      fd.write('\"Using Fortran compiler: %s %s ${FOPTFLAGS} ${FFLAGS} %s\\n\"\n' % (self.setCompilers.getCompiler(), self.setCompilers.getCompilerFlags(), self.setCompilers.CPPFLAGS))
      self.setCompilers.popLanguage()
    fd.write('\"-----------------------------------------\\n\";\n')
    fd.write('static const char *'+self.project+'compilerflagsinfo = \"\\n\"\n')
    fd.write('\"Using include paths: %s %s %s\\n\"\n' % ('-I'+os.path.join(self.dir.dir, self.arch.arch, 'include'), '-I'+os.path.join(self.dir.dir, 'include'), self.PACKAGES_INCLUDES))
    fd.write('\"-----------------------------------------\\n\";\n')
    fd.write('static const char *'+self.project+'linkerinfo = \"\\n\"\n')
    self.setCompilers.pushLanguage('Cxx')
    fd.write('\"Using C linker: %s\\n\"\n' % (self.setCompilers.getLinker()))
    self.setCompilers.popLanguage()
    if hasattr(self.compilers, 'FC'):
      self.setCompilers.pushLanguage('FC')
      fd.write('\"Using Fortran linker: %s\\n\"\n' % (self.setCompilers.getLinker()))
      self.setCompilers.popLanguage()
    #fd.write('\"Using libraries: %s%s -L%s %s %s %s\\n\"\n' % (self.setCompilers.CSharedLinkerFlag, os.path.join(self.dir.dir, self.arch.arch, 'lib'), os.path.join(self.dir.dir, self.arch.arch, 'lib'), '-lpetscts -lpetscsnes -lpetscksp -lpetscdm -lpetscmat -lpetscvec -lpetscsys', self.PACKAGES_LIBS, self.libraries.toStringNoDupes(self.compilers.flibs+self.compilers.cxxlibs+self.compilers.LIBS.split(' '))))
    fd.write('\"-----------------------------------------\\n\";\n')
    fd.close()
    return

  def configureScript(self):
    '''Output a script in the conf directory which will reproduce the configuration'''
    import nargs
    import sys
    if not os.path.isdir(os.path.join(self.arch.arch, 'conf')):
      os.makedirs(os.path.join(self.arch.arch, 'conf'))
    scriptName = os.path.join(self.arch.arch, 'conf', 'reconfigure-'+self.arch.arch+'.py')
    args = dict([(nargs.Arg.parseArgument(arg)[0], arg) for arg in self.framework.clArgs])
    if 'configModules' in args:
      if nargs.Arg.parseArgument(args['configModules'])[1] == self.Project+'.Configure':
        del args['configModules']
    if 'optionsModule' in args:
      if nargs.Arg.parseArgument(args['optionsModule'])[1] == self.Project+'.compilerOptions':
        del args['optionsModule']
    if self.PROJECT+'_ARCH' not in args:
      args[self.PROJECT+'_ARCH'] = self.PROJECT+'_ARCH='+str(self.arch.arch)
    f = open(scriptName, 'w')
    f.write('#!'+sys.executable+'\n')
    f.write('if __name__ == \'__main__\':\n')
    f.write('  import sys\n')
    f.write('  import os\n')
    f.write('  sys.path.insert(0, os.path.abspath(\'config\'))\n')
    f.write('  import configure\n')
    # pretty print repr(args.values())
    f.write('  configure_options = [\n')
    for itm in args.values():
      f.write('    \''+str(itm)+'\',\n')
    f.write('  ]\n')
    f.write('  configure.ConfigurationManager('+self.Project+').configure(configure_options)\n')
    f.close()
    try:
      os.chmod(scriptName, 0o775)
    except OSError as e:
      self.framework.logPrint('Unable to make reconfigure script executable:\n'+str(e))
    self.framework.actions.addArgument(self.Project, 'File creation', 'Created '+scriptName+' for automatic reconfiguration')
    return

  def configure(self):
    #if not os.path.samefile(self.dir.dir, os.getcwd()):
    #  raise RuntimeError('Wrong PETSC_DIR option specified: '+str(self.dir.dir) + '\n  Configure invoked in: '+os.path.realpath(os.getcwd()))
    #if self.framework.argDB['prefix'] and os.path.isdir(self.framework.argDB['prefix']) and os.path.samefile(self.framework.argDB['prefix'], self.dir.dir):
    #  raise RuntimeError('Incorrect option --prefix='+self.framework.argDB['prefix']+' specified. It cannot be same as PETSC_DIR!')
    self.framework.header          = self.arch.arch+'/include/'+self.project+'conf.h'
    self.framework.cHeader         = self.arch.arch+'/include/'+self.project+'fix.h'
    self.framework.makeMacroHeader = self.arch.arch+'/conf/'+self.project+'variables'
    self.framework.makeRuleHeader  = self.arch.arch+'/conf/'+self.project+'rules'
    self.executeTest(self.configureScript)
    
    self.Dump()
    self.dumpConfigInfo()
    self.dumpMachineInfo()
    self.framework.log.write('================================================================================\n')
    self.logClear()
    return
