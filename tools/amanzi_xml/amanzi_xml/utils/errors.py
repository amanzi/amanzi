class XMLError(RuntimeError):
    """Generic error in AmanziXML parsing"""
    pass

class MissingXMLError(XMLError):
    """Missing child in AmanziXML"""
    pass

class NonUniqueXMLError(XMLError):
    """Multiple children of name in AmanziXML"""
    pass

class NotNativeSpecError(RuntimeError):
    """Error for attempting to read a non-native spec file."""
    pass


def deprecated(msg):
    def decorator(func):
        def wrapper(*args, **kwargs):
            warnings.warn(msg, warnings.DeprecatedWarning)
            return func(*args, **kwargs)
        return wrapper
    return decorator

