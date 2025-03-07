import xml.etree.ElementTree as etree

from . import base
from amanzi_xml.utils import parser

class Comment(etree.Element):
    def __init__(self, text):
        super(Comment, self).__init__(etree.Comment)
        self.text = text
    
    def getName(self):
        return str(self.tag)

    def indent(self, ntabs, doublespace=False, doublespace_two=False):
        self.tail = "\n" + " "*ntabs*base._tabsize

    @classmethod
    def from_Element(cls, elem):
        return Comment(elem.text)

# register
parser.objects['function Comment'] = Comment
