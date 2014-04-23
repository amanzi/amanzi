class XMLError(RuntimeError):
    """Generic error in AmanziXML parsing"""
    pass

class MissingXMLError(XMLError):
    """Missing child in AmanziXML"""
    pass

class NonUniqueXMLError(XMLError):
    """Multiple children of name in AmanziXML"""
    pass

