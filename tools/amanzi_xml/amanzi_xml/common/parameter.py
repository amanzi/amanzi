from . import base
import amanzi_xml.utils.io as io

_valid_parameter_types = ['double', 'int', 'string', 'bool']
def _valid_parameter_from_string(ptype, value):
    # ensure types
    retval = None

    if ptype == "double":
        try:
            retval = float(value)
        except ValueError:
            
            raise RuntimeError("Parameter of type double with invalid value \"%s\""%str(value))

    elif ptype == "int":
        try:
            retval = int(value)
        except ValueError:
            raise RuntimeError("Parameter of type int with invalid value \"%s\""%str(value))

    elif ptype == "bool":
        if value is True:
            retval = value
        elif value is False:
            retval = value
        elif value == "true":
            retval = True
        elif value == "True":
            retval = True
        elif value == "TRUE":
            retval = True
        elif value == "false":
            retval = False
        elif value == "False":
            retval = False
        elif value == "FALSE":
            retval = False
        else:
            raise RuntimeError("Parameter of type bool with invalid value \"%s\""%str(value))

    elif ptype == "string":
        try:
            assert type(value) is str
        except AssertionError:
            raise RuntimeError("Parameter of type string with invalid value \"%s\""%str(value))
        else:
            retval = str(value)
    return retval


def _valid_parameter_from_type(ptype, value):
    # null string is valid
    if type(value) == str and len(value) == 0:
        return ''
    
    # ensure types
    retval = None
    if ptype == "double":
        retval = str(float(value))
    elif ptype == "int":
        retval = str(int(value))
    elif ptype == "bool":
        if value is True:
            retval = "true"
        elif value is False:
            retval = "false"
    elif ptype == "string":
        retval = str(value)
    return retval


class Parameter(base.TeuchosBaseXML):
    """Base class for all Teuchos::Parameter entries"""

    def __init__(self, name, ptype, value):
        """Initialize Parameter with name, type, and value"""
        super(Parameter, self).__init__('Parameter',{'name': str(name)})
        self.setType(str(ptype))
        self.setValue(value)

    @classmethod
    def from_Element(cls, elem):
        """Constructor from etree.Element class, with checking for valid types and values"""
        # new the class
        this = super(Parameter,cls).__new__(cls)

        # init the class with tag
        super(Parameter,this).__init__('Parameter')

        # add the named items
        for el in elem.items():
            this.set(el[0],el[1])

        mytype = this.get("type")
        myval = this.get("value")
        try:
            this.setType(mytype) #validation
            this.setValue(myval) #validation
        except RuntimeError as err:
            raise RuntimeError(str(err)+"\n from xml entry: \"%s\""%str(elem.attrib))
        return this

    def __str__(self):
        return io.toString(self)

    def indent(self, ntabs, doublespace=False, doublespace_two=False):
        self.text = ""
        self.tail = "\n" + " "*ntabs*base._tabsize

    def setName(self, name):
        """Set the name.  Maintained for interface consistency"""
        self.set('name', name)
        
    def setType(self, ptype):
        """Set the type.  This invalidates any current value!"""
        if ptype in _valid_parameter_types:
            self.set('type', ptype)
            self._basetype = ptype
            self._isarray = False
        elif ptype in ["Array(%s)"%lcv for lcv in _valid_parameter_types]:
            self.set('type', ptype)
            self._basetype = ptype[6:-1]
            self._isarray = True
        elif ptype in ["Array %s"%lcv for lcv in _valid_parameter_types]:
            self._basetype = ptype.split()[1]
            self.set('type', "Array(%s)"%self._basetype)
            self._isarray = True
        else:
            raise RuntimeError("Unknown Parameter type %s"%ptype)

        self.set('value', None)

    def setValue(self, value):
        """Set the value.  This checks it can be cast to the correct type."""
        if self._isarray:
            if type(value) is str:
                vals = value.strip().strip('{').strip('}').split(',')
                if len(vals) == 1 and vals[0] == '':
                    self.value = ['',]
                else:
                    self.value = [self._checkSingleValueFromString(val) for val in vals]
            elif type(value) is list:
                if len(value) == 1 and value[0] == '':
                    self.value = ['',]
                else:
                    self.value = [self._checkSingleValueFromString(val) for val in value]

            self.set("value", "{"+",".join([self._checkSingleValueFromType(val) for val in self.value])+"}")
        else:
            self.value = self._checkSingleValueFromString(value)
            self.set("value", self._checkSingleValueFromType(self.value))

    def getValue(self):
        """Gets the value, cast to the correct type."""
        return self._checkSingleValueFromString(self.get('value'))

    def _checkSingleValueFromString(self, value):
        retval = _valid_parameter_from_string(self._basetype, value)
        assert retval is not None
        return retval

    def _checkSingleValueFromType(self, value):
        retval = _valid_parameter_from_type(self._basetype, value)
        assert retval is not None
        return retval


from amanzi_xml.utils import parser
def make_class(typename, array=False):
    if array:
        class C(Parameter):
            """Parameter of type Array(%s)"""%typename
            _typename = "Array(%s)"%typename
            def __init__(self, name, val):
                super(C,self).__init__(name, self._typename, val)

    else:
        class C(Parameter):
            """Parameter of type %s"""%typename
            _typename = typename
            def __init__(self, name, val):
                super(C,self).__init__(name, self._typename, val)
    return C

# create and register
IntParameter = make_class("int")
parser.objects[IntParameter._typename] = IntParameter

DoubleParameter = make_class("double")
parser.objects[DoubleParameter._typename] = DoubleParameter

BoolParameter = make_class("bool")
parser.objects[BoolParameter._typename] = BoolParameter

StringParameter = make_class("string")
parser.objects[StringParameter._typename] = StringParameter

ArrayIntParameter = make_class("int",True)
parser.objects[ArrayIntParameter._typename] = ArrayIntParameter

ArrayDoubleParameter = make_class("double",True)
parser.objects[ArrayDoubleParameter._typename] = ArrayDoubleParameter

ArrayBoolParameter = make_class("bool",True)
parser.objects[ArrayBoolParameter._typename] = ArrayBoolParameter

ArrayStringParameter = make_class("string",True)
parser.objects[ArrayStringParameter._typename] = ArrayStringParameter



