import os,sys
try:
    import xml.etree.cElementTree as ET
except:
    import xml.etree.ElementTree as ET
try:
    from amanzi_xml.utils.search import *
except:
    print "Error: Unable to locate amanzi xml_functions module"
    print "Add $AMANZI_SRC_DIR/tools/amanzi_xml to your PYTHONPATH environment variable"
    sys.exit()
from xml.dom import minidom
import re
import sys
import numpy
import subprocess

class ATSXML(object):
    def __init__(self,filename=None,xml_string=None):
        if filename is not None and xml_string is not None:
            print 'Error: Cannot define both filename and string!\n'
            return 
        elif filename is not None: 
            self.filename = filename
            self.read_xml(filename)
        elif xml_string is not None: 
            self.filename = None
            self.fromstring(xml_string)
    def read_xml(self,filename):
        self.tree = ET.parse(filename)
        self.root = self.tree.getroot()
    def fromstring(self,xml_string):
        self.root = ET.fromstring(xml_string)
    @property
    def regions(self):
        return self.root.findall('ParameterList/[@name="Regions"]/ParameterList')
    @property
    def parent_map(self):
        return create_parent_map(self)
    def replace_regions(self,fromatsxml,mapping=None):
        ''' 
        Replace regions in xml with regions in another xml

        :param fromatsxml: ATSXML object with replacement regions
        :type fromatsxml: ATSXML object
        :param mapping: dictionary of replacement region names keyed by current region names
        :type mapping: dict[str] = str
        '''
        replace_regions(self,fromatsxml,mapping=mapping)
    def replace_by_value(self,oldvalue,newvalue):
        """Replace all matches of a given 'value'='oldvalue' with 'newvalue'"""
        replace_by_value(self.root,oldvalue,newvalue)
    def replace_by_name(self,name,value):
        """Replace all matches of a given 'name' with 'newvalue'"""
        replace_by_name(self.root,name,value)
    def replace_file(self,newfilename,oldfilename=None):
        if oldfilename is not None: replace_by_value(self.root,oldfilename,newfilename)
        else: replace_by_name(self.root,'file',newfilename)
    def replace_mesh_file(self,newfilename):
        for e in self.root.findall('.//ParameterList/[@name="Region: Labeled Set"]'):
            e2 = e.find('.//Parameter/[@name="file"]')
            print e2.attrib
            if e2 is not None: e2.set('value',newfilename)
            e2 = e.find('.//Parameter/[@name="File"]')
            if e2 is not None: e2.set('value',newfilename)
        for e in self.root.findall('.//ParameterList/[@name="Read Mesh File"]'):
            e2 = e.find('.//Parameter/[@name="file"]')
            if e2 is not None: e2.set('value',newfilename)
            e2 = e.find('.//Parameter/[@name="File"]')
            if e2 is not None: e2.set('value',newfilename)
    def replace_restart_file(self,filename):
        replace_by_name(self.root,'restart file', filename)
    def replace_elem(self,elem_sink,elem_src):
        """Replace the element 'sink' with the element 'src' in the hierarchy 'xml'"""
        replace_elem(self.root,elem_sink,elem_src)
    # Potentially not needed anymore, xml's written directly with write look ok now
    def write(self,filename):
        #self.tree.write(filename)
        m_str = ET.tostring(self.root)
        m_str_strip = m_str.replace('\n','')
        m_str_strip = re.sub('>\s*<','><',m_str_strip)
        m_reparsed = minidom.parseString(m_str_strip)
        with open(filename, "w") as f:
            f.write(m_reparsed.toprettyxml(indent="  "))
    #def findall_name(self,name,elem=None):
    #    if elem is None:
    #        return self.root.findall('.//*/[@name="'+str(name)+'"]')
    #    else:
    #        return elem.findall('.//*/[@name="'+str(name)+'"]')
    #def find_name(self,name,elem=None):
    #    if elem is None:
    #        elem = self.root
    #    if isinstance(name,str): name = numpy.array([name])
    #    for n in name:
    #        elem = elem.find('.//*/[@name="'+str(n)+'"]')
    #    return elem
    #def findall_value(self,value,elem=None):
    #    if elem is None:
    #        return self.root.findall('.//*/[@value="'+str(value)+'"]')
    #    else:
    #        return elem.findall('.//*/[@value="'+str(value)+'"]')
    #def find_value(self,value,elem=None):
    #    if elem is None:
    #        return self.root.find('.//*/[@value="'+str(value)+'"]')
    #    else:
    #        return elem.find('.//*/[@value="'+str(value)+'"]')
    #def findall_name_value(self,name,value,elem=None):
    #    if elem is None:
    #        return self.root.findall('.//*/[@name="'+str(name)+'"]/[@value="'+str(value)+'"]')
    #    else:
    #        return elem.findall('.//*/[@name="'+str(name)+'"]/[@value="'+str(value)+'"]')
    #def find_name_value(self,name,value,elem=None):
    #    if elem is None:
    #        return self.root.find('.//*/[@name="'+str(name)+'"]/[@value="'+str(value)+'"]')
    #    else:
    #        return elem.find('.//*/[@name="'+str(name)+'"]/[@value="'+str(value)+'"]')
    def uptree(self,elem,level=1):
        pm = self.parent_map
        elems = [elem]
        for i in range(level): 
            if elem in pm:
                elems.append(pm[elem])
                elem = pm[elem]
            else:
                break
        ind = ''
        s = ''
        for e in reversed(elems[1:]):
            s += "\n"+ind+e.attrib['name']
            ind += '    '
        s += "\n"+ind+elems[0].attrib['name']
        if 'value' in elems[0].attrib:
            s += ": "+elems[0].attrib['value']
        for e in elems[1].getchildren():
            if 'value' in e.attrib and e is not elems[0]:
                s += "\n"+ind+e.attrib['name']+": "+e.attrib['value']
        print s
    def dump(self,elem=None,uplevel=0):
        if elem is None: return ET.dump(self)
        else:
            if uplevel > 0:
                pm = self.parent_map
                for i in range(uplevel): elem = pm[elem]
            ET.dump(elem)
    def remove_elem(self,elem):
        print "remove_elem deprecated, it has been replaced with remove"
    def remove(self,elem):
        """Removes the xml 'elem' wherever it is locaed in 'xml'"""
        remove(self.root,elem)
    def set_verbosity(self,value):
        '''
        Set verbosity level to value in all available Verbose Objects

        :param value: Verbosity level: high,low
        :type value: str
        '''
        for e in findall_name(self.root,'Verbosity Level'):
            e.set('value',value)
    # Now in amanzi_xml.utils.search
    #def add_Parameter(self,name,value,vtype,elem=None):
    #    if elem is None: ET.SubElement(self.root,'Parameter',{'name':name,'type':vtype,'value':str(value)})
    #    else: ET.SubElement(elem,'Parameter',{'name':name,'type':vtype,'value':str(value)})
    #def add_ParameterList(self,name,elem=None):
    #    if elem is None: ET.SubElement(self.root,'ParameterList',{'name':name,'type':'ParameterList'})
    #    else: ET.SubElement(elem,'ParameterList',{'name':name,'type':'ParameterList'})

#class Region(object):
#    def __init__(self,element):
#        self.element = element
#        self.name = element.attrib['name']
#        self.type = element.find('ParameterList').attrib['name']
#        d = {}
#        pl = element.find('ParameterList')
#        for v in pl.findall('Parameter'):
#            d[v.attrib['name']] = v.attrib['value']
#        self.dict = d
#    def __repr__(self):
#        return self.name

def get_root(filename=None,xml_string=None):
    if filename is not None and xml_string is not None:
        print 'Error: Cannot define both filename and string!\n'
        return 
    elif filename is not None: 
        t = ET.parse(filename)
        r = t.getroot()
    elif xml_string is not None: 
        r = ET.fromstring(xml_string)
    return r

def get_regions(xml):
    return xml.findall('ParameterList/[@name="Regions"]/ParameterList')

def write_xml(xml,filename):
    #self.tree.write(filename)
    m_str = ET.tostring(xml)
    m_str_strip = m_str.replace('\n','')
    m_str_strip = re.sub('>\s*<','><',m_str_strip)
    m_reparsed = minidom.parseString(m_str_strip)
    with open(filename, "w") as f:
        f.write(m_reparsed.toprettyxml(indent="  "))

def replace_regions(toatsxml,fromatsxml,mapping=None):
    ''' 
    Replace regions in xml with regions in another xml

    :param toatsxml: ATSXML object to replace regions in
    :type fromatsxml: ATSXML object
    :param fromatsxml: ATSXML object with replacement regions
    :type fromatsxml: ATSXML object
    :param mapping: dictionary of replacement region names keyed by current region names
    :type mapping: dict[str] = str
    '''
    # Create mapping
    m = {}
    for r in get_regions(fromatsxml): 
        m[r.attrib['name']] = r.attrib['name']
    # Override mappings with those provided
    if mapping is not None:
        for k,v in mapping.items(): m[k] = v
    # Change regions
    for rb in toatsxml.findall('.//Parameter/[@name="region"]/../..'):
        for rp in rb.findall('*/Parameter/[@name="region"]/..'):
            r = rp.find('./Parameter/[@name="region"]')
            name = r.attrib['value']
            if name in m:
                r.set('value',m[name])
            else: 
                rb.remove(rp)
    for r in toatsxml.findall('.//Parameter/[@name="regions"]'):
        rs = []
        for rn in r.attrib['value'].strip('{}').split(','):
            if rn in m: rs.append(m[rn])
        #rs = [rn for rn in r.attrib['value'].strip('{}').split(',')]
        rs = '{'+','.join(rs)+'}'
        r.set('value',rs)
    for r in toatsxml.findall('.//Parameter/[@name="domain name"]'):
        r.set('value',m[r.attrib['value']])
    # Replace Mesh and Regions blocks
    replace_elem(toatsxml,toatsxml.find('./ParameterList/[@name="Regions"]'),fromatsxml.find('./ParameterList/[@name="Regions"]'))
    replace_elem(toatsxml,find_name(toatsxml,'Mesh'),find_name(fromatsxml,'Mesh'))
    newinput = get_value(toatsxml,'Native Unstructured Input')
    replace_by_name(toatsxml,'Native Unstructured Input',newinput)
    newgridoption = get_value(toatsxml,'grid_option')
    replace_by_name(toatsxml,'grid_option',newgridoption)

def run(xml, ats_exe=None, nproc=1, mpiexec='mpiexec', run_file_name='run.xml',stdout=None,stderr=None,hostname=None,cpuset=None):
    ''' Run ats model based on tpl_xml_file using parameters defined in pars dictionary '''

    # ensure that ATS's executable exists and that it's module is loaded
    if ats_exe is None:
        try:
            path = os.path.join(os.environ['ATS_DIR'],'bin')
        except KeyError:
            raise RuntimeError("Missing ATS installation, please set the ATS_DIR environmental variable.")
        ats_exe = os.path.join(path, "ats")
        if not os.path.isfile(ats_exe):
            raise RuntimeError("Missing ATS installation, please build and install ATS.")

    write_xml(xml,run_file_name)
    
    try:
        outfile = open(stdout,'w')
        if not stderr is None:
            errfile = open(stderr,'w')
        cmdlst = [mpiexec,'-n',str(nproc)]
        if cpuset is not None:
            if isinstance(cpuset,list): ','.join([str(v) for v in cpuset])
            cmdlst += ['-cpu-set',cpuset]
        cmdlst += [ats_exe,"--xml_file="+run_file_name]
        ierr = subprocess.call(cmdlst,stdout=outfile,stderr=errfile)
        print ierr
    except:
        print "Error: ATS simulation failed"



