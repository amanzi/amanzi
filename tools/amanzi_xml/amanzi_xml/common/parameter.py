from . import base
from decimal import Decimal
import amanzi_xml.utils.io as io

_valid_parameter_types = ['double', 'int', 'string', 'bool']
def _validValueFromString(ptype, value):
    """Given a string, return a ptype-typed value."""
    assert isinstance(value, str)
    value = value.strip()

    if ptype == 'double':
        if '.' not in value and 'e' not in value.lower():
            value = value + '.0'

        try:
            retval = Decimal(value)
        except ValueError:
            raise ValueError("Parameter of type double with invalid value \"%s\""%str(value))

    elif ptype == 'int':
        try:
            retval = int(value)
        except ValueError:
            raise ValueError("Parameter of type int with invalid value \"%s\""%str(value))

    elif ptype == 'bool':
        if value.lower() == 'true':
            retval = True
        elif value.lower() == 'false':
            retval = False
        else:
            raise ValueError("Parameter of type bool with invalid value \"%s\""%str(value))

    elif ptype == 'string':
        retval = value

    else:
        raise ValueError(f"Invalid ptype {ptype}")

    return retval


def _validValueFromType(ptype, value):
    """Given a typed value, convert it to ptype."""
    if ptype == 'string':
        if not isinstance(value, str):
            raise ValueError("Parameter of type string given non-string value.")
        return value

    elif ptype == 'bool':
        if not isinstance(value, bool):
            raise ValueError("Parameter of type bool given non-bool value.")
        return value

    elif ptype == 'int':
        if not isinstance(value, int):
            try:
                int_val = np.round(value)
            except ValueError:
                raise ValueError("Parameter of type int given non-int-convertible value.")

            if abs(int_val - val) > 1.e-10:
                raise ValueError("Parameter of type int given non-int-convertible value.")
            value = int(int_val)
        return int(value)

    elif ptype == 'double':
        if not isinstance(value, Decimal):
            try:
                value = Decimal(value)
            except ValueError:
                raise ValueError("Parameter of type float given non-Decimal-convertible value.")

        if '.' not in str(value) and 'e' not in str(value).lower():
            value = Decimal(str(value)+'.0')
        return value

    else:
        raise ValueError(f"Invalid ptype {ptype}")


def _validValue(ptype, value):
    """Valid typed value from string or typed value."""
    if ptype != 'string' and isinstance(value, str):
        return _validValueFromString(ptype, value)
    else:
        return _validValueFromType(ptype, value)
    

def _checkType(ptype, value):
    """Simply asserts that value is the right type."""
    if ptype == 'string':
        assert isinstance(value, str)
    elif ptype == 'bool':
        assert isinstance(value, bool)
    elif ptype == 'int':
        assert isinstance(value, int)
    elif ptype == 'double':
        assert isinstance(value, Decimal)
    else:
        assert(False)

    
def _convertTypeToString(ptype, value):
    """Given a type, write it to a string."""
    _checkType(ptype, value)
    
    retval = None

    if ptype == 'bool':
        if value is True:
            retval = "true"
        elif value is False:
            retval = "false"

    elif ptype == 'double':
        retval = f'{value:g}'
        if '.' not in retval and 'e' not in retval.lower():
            retval = retval + '.0'

    else:
        retval = str(value).strip()
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
        except ValueError as err:
            raise ValueError(str(err)+f"\n from xml entry: \"{str(elem.attrib)}\"")
        return this

    def __str__(self):
        return io.toString(self)

    def indent(self, ntabs, doublespace=False, doublespace_two=False):
        self.text = ""
        self.tail = "\n" + " "*ntabs*base._tabsize

    def getName(self):
        """Get the name.  Maintained for interface consistency"""
        return self.get('name')
        
    def setName(self, name):
        """Set the name.  Maintained for interface consistency"""
        self.set('name', name)
        
    def getType(self):
        """Get the type string.  Maintained for interface consistency"""
        return self.get('type')

    def setType(self, ptype):
        """Set the type.  This invalidates any current value!"""
        if ptype in _valid_parameter_types:
            self.set('type', ptype)
            self._basetype = ptype
            self._isarray = False
        elif ptype in [f'Array({t})' for t in _valid_parameter_types]:
            self.set('type', ptype)
            self._basetype = ptype[6:-1]
            self._isarray = True
        else:
            raise RuntimeError("Unknown Parameter type %s"%ptype)

    def getValue(self):
        """Gets the value, cast to the correct type. Note this is here for API consistency only."""
        return self.value
    
    def setValue(self, value):
        """Set the value given either string or a typed value.

        This checks it can be cast to the correct type. Note that
        setType() must already have been called.

        Note this must handle both strings and typed values.
        """
        if self._isarray:
            if type(value) is str:
                value = value.strip().strip('{').strip('}').split(',')
                if len(value) == 1 and value[0] == '':
                    value = []

            try:
                vals = [i for i in value]
            except TypeError:
                raise ValueError(f'Cannot convert "{value}" to list for an array-typed parameter.')

            self.value = [_validValue(self._basetype, val) for val in vals]
            self.set("value", "{" + ", ".join([_convertTypeToString(self._basetype, val) for val in self.value]) + "}")

        else:
            self.value = _validValue(self.getType(), value)
            self.set("value", _convertTypeToString(self.getType(), self.value))
            
    def __copy__(self):
        return Parameter(self.getName(), self.getType(), self.getValue())
    
    def __deepcopy__(self, memo):
        cp = self.__copy__()
        memo[id(self)] = cp
        return cp


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



