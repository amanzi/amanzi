try:
    import xml.etree.cElementTree as ET
except:
    import xml.etree.ElementTree as ET
from copy import deepcopy
from xml.dom import minidom
import re
import sys
import numpy
from glob import glob
import parse_xml

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
        return {c:p for p in self.tree.iter() for c in p}
    @property
    def simulated_times(self):
        base_name = self.find_name('visualization').find('.//*/[@name="file name base"]').attrib['value']
        # Assume base_name+'_data*' will get the times ok, should work unless things change in ATS output
        base_name += '_data'
        ts = []
        for fn in glob(base_name+'.h5.*.xmf'):
            t = ET.parse(fn)
            r = t.getroot()
            ts.append(float(r.find('.//*/Time').attrib['Value']))
        return numpy.array(ts)
    def replace_regions(self,fromatsxml,mapping=None):
        ''' 
        Replace regions in xml with regions in another xml

        :param fromatsxml: ATSXML object with replacement regions
        :type fromatsxml: ATSXML object
        :param mapping: dictionary of replacement region names keyed by current region names
        :type mapping: dict[str] = str
        '''
        # Create mapping
        m = {}
        for r in fromatsxml.regions: 
            m[r.attrib['name']] = r.attrib['name']
        # Override mappings with those provided
        if mapping is not None:
            for k,v in mapping.items(): m[k] = v
        # Change regions
        for rb in self.root.findall('.//Parameter/[@name="region"]/../..'):
            for rp in rb.findall('*/Parameter/[@name="region"]/..'):
                r = rp.find('./Parameter/[@name="region"]')
                name = r.attrib['value']
                if name in m:
                    r.set('value',m[name])
                else: 
                    rb.remove(rp)
        for r in self.root.findall('.//Parameter/[@name="regions"]'):
            rs = []
            for rn in r.attrib['value'].strip('{}').split(','):
                if rn in m: rs.append(m[rn])
            #rs = [rn for rn in r.attrib['value'].strip('{}').split(',')]
            rs = '{'+','.join(rs)+'}'
            r.set('value',rs)
        for r in self.root.findall('.//Parameter/[@name="domain name"]'):
            r.set('value',m[r.attrib['value']])
        # Replace Mesh and Regions blocks
        self.replace_elem(self.find_name('Regions'),fromatsxml.find_name('Regions'))
        self.replace_elem(self.find_name('Mesh'),fromatsxml.find_name('Mesh'))
        newinput = self.getvalue('Native Unstructured Input')
        self.replace_by_name('Native Unstructured Input',newinput)
        newgridoption = self.getvalue('grid_option')
        self.replace_by_name('grid_option',newgridoption)
    def replace_by_value(self,oldvalue,newvalue):
        for r in self.root.findall('.//Parameter/[@value="'+oldvalue+'"]'):
            r.set('value',str(newvalue))
    def replace_by_name(self,name,value):
        for r in self.root.findall('.//Parameter/[@name="'+str(name)+'"]'):
            r.set('value',str(value))
    #def replace_by_name(self,name,value,elem=None):
    #    e = self.find_name(name,elem=None)
    #    e.set('value',str(value))
    def replace_file(self,newfilename,oldfilename=None):
        if oldfilename is not None: self.replace_by_value(oldfilename,newfilename)
        else: self.replace_by_name('file',newfilename)
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
        self.replace_all_by_name('restart file', filename)
    #def replace_elems(self,fromatsxml,element_names):
    #    for nm in element_names:
    #        e = fromatsxml.root.findall(".//ParameterList[@name='"+nm+"']")[0]
    #        self.root.remove(self.root.findall(".//ParameterList[@name='"+nm+"']")[0])
    #        self.root.append(e)        
    def replace_elem(self,elem_sink,elem_src):
        pm = self.parent_map
        p_elem_sink = pm[elem_sink]
        p_elem_sink.remove(elem_sink)
        p_elem_sink.append(elem_src)
    def getvalue(self,name):
        out = []
        for r in self.root.findall('.//Parameter/[@name="'+str(name)+'"]'):
            out.append(r.attrib['value'])
        if len(out) > 1: return out
        else: return out[0]
    def write(self,filename):
        #self.tree.write(filename)
        m_str = ET.tostring(self.root)
        m_str_strip = m_str.replace('\n','')
        m_str_strip = re.sub('>\s*<','><',m_str_strip)
        m_reparsed = minidom.parseString(m_str_strip)
        with open(filename, "w") as f:
            f.write(m_reparsed.toprettyxml(indent="  "))
    def findall_name(self,name,elem=None):
        if elem is None:
            return self.root.findall('.//*/[@name="'+str(name)+'"]')
        else:
            return elem.findall('.//*/[@name="'+str(name)+'"]')
    def find_name(self,name,elem=None):
        if elem is None:
            elem = self.root
        if isinstance(name,str): name = numpy.array([name])
        for n in name:
            elem = elem.find('.//*/[@name="'+str(n)+'"]')
        return elem
    def findall_value(self,value,elem=None):
        if elem is None:
            return self.root.findall('.//*/[@value="'+str(value)+'"]')
        else:
            return elem.findall('.//*/[@value="'+str(value)+'"]')
    def find_value(self,value,elem=None):
        if elem is None:
            return self.root.find('.//*/[@value="'+str(value)+'"]')
        else:
            return elem.find('.//*/[@value="'+str(value)+'"]')
    def findall_name_value(self,name,value,elem=None):
        if elem is None:
            return self.root.findall('.//*/[@name="'+str(name)+'"]/[@value="'+str(value)+'"]')
        else:
            return elem.findall('.//*/[@name="'+str(name)+'"]/[@value="'+str(value)+'"]')
    def find_name_value(self,name,value,elem=None):
        if elem is None:
            return self.root.find('.//*/[@name="'+str(name)+'"]/[@value="'+str(value)+'"]')
        else:
            return elem.find('.//*/[@name="'+str(name)+'"]/[@value="'+str(value)+'"]')
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
        pm = self.parent_map
        pm[elem].remove(elem)
    def set_verbosity(self,value):
        '''
        Set verbosity level to value in all available Verbose Objects

        :param value: Verbosity level: high,low
        :type value: str
        '''
        for e in self.findall_name('Verbosity Level'):
            e.set('value',value)
    def add_Parameter(self,name,value,vtype,elem=None):
        if elem is None: ET.SubElement(self.root,'Parameter',{'name':name,'type':vtype,'value':str(value)})
        else: ET.SubElement(elem,'Parameter',{'name':name,'type':vtype,'value':str(value)})
    def add_ParameterList(self,name,elem=None):
        if elem is None: ET.SubElement(self.root,'ParameterList',{'name':name,'type':'ParameterList'})
        else: ET.SubElement(elem,'ParameterList',{'name':name,'type':'ParameterList'})
class Region(object):
    def __init__(self,element):
        self.element = element
        self.name = element.attrib['name']
        self.type = element.find('ParameterList').attrib['name']
        d = {}
        pl = element.find('ParameterList')
        for v in pl.findall('Parameter'):
            d[v.attrib['name']] = v.attrib['value']
        self.dict = d
    def __repr__(self):
        return self.name




