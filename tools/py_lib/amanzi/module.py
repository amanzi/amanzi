################################################################################

import os, sys, string
################################################################################

# For debugging only!
def print_environ():
    for var in list(os.environ.keys()):
        print(var + '=' + os.environ[var])

# Module Gobal
LOADED_MODULES_KEY = 'LOADEDMODULES'
MODULE_PATH_KEY    = 'MODULEPATH' 
MODULES_HOME_KEY   = 'MODULESHOME'
MODULE_VERSION_KEY = 'MODULE_VERSION'


class ModuleInterface:

    def __init__(self,module_cmd=None):

        # Set the module command
        if module_cmd is None:
            module_cmd_search = 'which modulecmd'
        else:
            module_cmd_search = 'which ' + module_cmd

        self._modulecmd = os.popen(module_cmd_search).read().strip()
        assert os.path.exists(self._modulecmd), \
                "Can not locate modulecmd: %s returned %s" % ( module_cmd_search, self._modulecmd)

    def modulecmd(self,command,*arguments):
        commands = os.popen('%s python %s %s' % (self._modulecmd,command,string.join(arguments))).read()
        exec(commands)

        # Catch any changes to PYTHONPATH
        if 'PYTHONPATH' in os.environ:
            pp = ['']
            pythonpath = os.environ['PYTHONPATH'].split(":")
            for p in sys.path:
                if ( p not in pp) and (p):
                    pp.append(p)
            sys.path = pp

    def list(self):
        self.modulecmd('list')  
        if LOADED_MODULES_KEY in os.environ:
            return os.environ[LOADED_MODULES_KEY].rsplit(':')
        else:
            return []

    def load(self,modules):
        self.modulecmd('load',modules)

    def unload(self,module_name):
        if self.isloaded(module_name):
            self.modulecmd('unload',module_name)
        else:
            err_mess = "Will not unload %s. Module is not loaded" % (module_name)
            print(err_mess)

    def use(self,path,append=False):
        assert os.path.exists(path), \
                "Can not add path: %s to MODULEPATH does not exist" % (path)
        if append:
            arguments='--append' + ' ' + path
        else:
            arguments=path
        self.modulecmd('use',arguments)
           
    def unuse(self,path):
        self.modulecmd('unuse',path)

    def swap(self,old_module,new_module):
        self.modulecmd('swap',old_module,new_module)

    def purge(self):
        self.modulecmd('purge')

    def isloaded(self,module_name):
        try:
            idx = self.list().index(module_name)
        except ValueError:
            idx = -1
       
        if idx < 0:
            return False
        else:
            return True

    def available(self,regexp_pattern=None):
        assert True, \
                "method available is not implemented at this time" 

    def version(self):
        if MODULE_VERSION_KEY in os.environ:
            return os.environ[MODULE_VERSION_KEY]
        else:
            return None

    def search_paths(self):
        if MODULE_PATH_KEY in os.environ:
            return os.environ[MODULE_PATH_KEY].rsplit(":")
        else:
            return []

################################################################################
if __name__ == '__main__':

    '''
    Change module_cmd to the full path name
    of modulecmd to avoid which search
    '''
    module_cmd = os.environ['MODULESHOME'] + '/bin/modulecmd'
    if module_cmd is None:
        print('Will define module command through a which search')

    module = ModuleInterface(module_cmd)

    print('Module Version:' + module.version())
    print('Current Search Paths :' + str(module.search_paths()))

    hdf5_module = 'hdf5-serial/1.8.5'
    module.load(hdf5_module)
    print('Loaded modules:' + str(module.list()))

    print('HDF5 Module is loaded:' + str(module.isloaded(hdf5_module)))

    try:
        module.use('/some/path/dne')
    except Exception:
        print('Caught assert error  while trying to add path that did not exist') 

    module.use(os.environ['HOME'])
    print('After calling use Search Paths :' + str(module.search_paths()))

    module.purge()
    print('After purging modules:' + str(module.list()))

    module.available()

