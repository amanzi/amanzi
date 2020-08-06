import warnings
from . import base
from . import parameter
from amanzi_xml.utils import parser
from amanzi_xml.utils.parser import ValidationLevel
from amanzi_xml.utils import io
from amanzi_xml.utils import errors
from amanzi_xml.utils import search

_doublespace = ["field evaluators", "PKs", "initial conditions"]
_doublespace_two = ["PKs"]
_postspace = ["mesh", "regions", "coordinator", "checkpoint",
              "visualization", "visualization_surface", "PKs"]
_prespace = ["PKs", "field evaluators", "initial conditions"]


class ParameterList(base.TeuchosBaseXML):
    """Base class for all Teuchos::ParameterList objects

    To inherit this class, a derived class should (at a minimum), provide the
    following class data:

      _required_parameters = [ (pname, ptype), ... ]
      _optional_parameters = [ (pname, ptype, default), ... ]

    where:
        pname:   is a string describing the element's name attribute,
        ptype:   is a string describing the element's type attribute (which
                 should also be the element class's _type data)
        default: is the default value used in the code (or None) (not
                 currently used)
    """
    _required_parameters = []
    _optional_parameters = []
    _type = "ParameterList"

    def __init__(self, name):
        """Initialize ParameterList with name"""
        super(ParameterList, self).__init__('ParameterList',{'name': str(name), 'type': self._type})

    @classmethod
    def from_Element(cls, elem, validation=None):
        """Constructor from etree.Element class, with checking for required Parameters"""
        # new the class
        this = super(ParameterList,cls).__new__(cls)
        # init the class with tag
        super(ParameterList,this).__init__('ParameterList')

        # pull in all named items from the Element
        for el in elem.items():
            this.set(el[0],el[1])

        # recursively call the factory, populating the children
        this.extend([parser.fromElement(el) for el in elem])

        # validate
        if validation is None:
            if cls == ParameterList:
                # this simply notes that ParameterList has no required or
                # optional information and can safely bypass validation.
                validation = ValidationLevel(None,None,None,None)
            else:
                validation = ValidationLevel()

        this.validate(validation)
        return this

    def __str__(self):
        return io.toString(self)

    def validate(self, vlevel):
        if not vlevel.none:
            # check required parameters
            for (pname,ptype) in self._required_parameters:
                # check name exists
                if not self.isElement(pname):
                    if vlevel["missing_required"] is not None:
                        if vlevel["missing_required"]:
                            raise RuntimeError("Missing required parameter %s in list %s"%(pname, self.get("name")))
                        else:
                            warnings.warn(RuntimeWarning("Missing required parameter %s in list %s"%(pname, self.get("name"))))

                # check type
                else:
                    if vlevel["type_required"] is not None:
                        p = self.getElement(pname)
                        if not isinstance(p, parser.objects[ptype]):
                            if vlevel["type_required"]:
                                raise RuntimeError("Parameter %s in list %s is of type %s, not required type %s"%(pname, self.get("name"), p.get("type"), ptype))
                            else:
                                warnings.warn(RuntimeWarning("Parameter %s in list %s is of type %s, not required type %s"%(pname, self.get("name"), p.get("type"), ptype)))

            # check optional parameters
            for elem in self:
                # make sure not checking a required parameter
                try:
                    ename,etype = ((ename,etype) for (ename,etype) in self._required_parameters if ename == elem.get("name")).next()
                except StopIteration:
                    pass
                else:
                    continue

                # look for the element in optionals
                try:
                    ename, etype, edef = ((ename,etype,edef) for (ename,etype,edef) in self._optional_parameters if ename == elem.get("name")).next()
                except StopIteration:
                    # element not in optionals
                    if vlevel["existing_not_optional"] is not None:
                        if vlevel["existing_not_optional"]:
                            raise RuntimeError("Parameter %s in list %s is not on the list of optional parameters"%(elem.get("name"), self.get("name")))
                        else:
                            warnings.warn(RuntimeWarning("Parameter %s in list %s is not on the list of optional parameters"%(elem.get("name"), self.get("name"))))
                else:
                    # found element
                    if vlevel["optional_type"] is not None:
                        if elem.get("type") != etype:
                            if vlevel["optional_type"]:
                                raise RuntimeError("Parameter %s in list %s is of type %s, not expected optional type %s"%(pname, self.get("name"), elem.get("type"), etype))
                            else:
                                warnings.warn(RuntimeWarning("Parameter %s in list %s is of type %s, not expected optional type %s"%(pname, self.get("name"), elem.get("type"), etype)))


    def sort(self):
        """Sorts mychildren() by params first, then lists."""
        self.getchildren().sort(key=lambda x: x.tag)

    def indent(self, ntabs, doublespace=False, doublespace_two=False):
        """Properly indent this list (and its Parameters/sublists) with [ntabs] tabs."""
        self.sort()

        doublespace_sublists = self.get("name") in _doublespace
        doublespace_sublists_two_deep = self.get("name") in _doublespace_two
        postspace = doublespace or self.get("name") in _postspace
        prespace = self.get("name") in _prespace
        
        if len(self) == 0:
            self.text = "\n" + " "*(ntabs)*base._tabsize
        else:
            self.text = "\n" + " "*(ntabs+1)*base._tabsize

            for el in self:
                el.indent(ntabs+1, doublespace_sublists or doublespace_two,
                          doublespace_sublists_two_deep)
            self[-1].tail = "\n" + " "*(ntabs)*base._tabsize
            if doublespace_sublists or doublespace_two:
                self[-1].tail = "\n"+self[-1].tail

        if prespace:
            self.text = "\n"+self.text
    
        if postspace:
            self.tail = "\n\n" + " "*(ntabs)*base._tabsize
        else:
            self.tail = "\n" + " "*(ntabs)*base._tabsize

    def isElement(self, name):
        try:
            val = self.getElement(name)
        except errors.MissingXMLError:
            return False
        else:
            return True

    def sublist(self, name):
        """Get a sublist of with the given name, creating if necessary."""
        try:
            val = self.getElement(name)
        except errors.MissingXMLError:
            self.append(ParameterList(name))
            return self.getElement(name)
        else:
            assert isinstance(val, ParameterList)
            return val

    def getName(self):
        return self.get("name")
    
    def setName(self, name):
        self.set("name", name)
        
    def setParameter(self, name, ptype, val):
        """Set parameter name of type ptype to value, creating a new Parameter if necessary."""
        try:
            param = self.getElement(name)
        except errors.MissingXMLError:
            pass
        else:
            assert ptype == param.get('type')
            self.remove(param)

        self.append(parameter.Parameter(name, ptype, val))

    def getElement(self, name):
        """Get Parameter/sublist from the ParameterList"""
        return search.child_by_name(self, name)

    def pop(self, name):
        return search.remove_element(self, name, False, True)
    
# register
parser.objects['ParameterList'] = ParameterList
