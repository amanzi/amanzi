from xml.etree.ElementTree import Element as Element

_tabsize = 2

class TeuchosBaseXML(Element):
    def __init__(self, *args, **kwargs):
        super(TeuchosBaseXML,self).__init__(*args, **kwargs)



