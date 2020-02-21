# -*- coding: utf-8 -*-
"""
    amanzi_xml
    ~~~~~~~~~~

    Original code from sphinx.directives.code - literalinclude directive

    :copyright: Copyright 2007-2011 by the Sphinx team, see Sphinx AUTHORS.
    :license: BSD, see Sphinx LICENSE for details.

    Derivative work - custom amanzi_xml_include directive

    :copyright: Copyright 2013 by the Amanzi team, see AUTHORS.
    :license: Three-clause BSD, see COPYRIGHT for details.

"""

import sys
import codecs

from docutils import nodes
from docutils.parsers.rst import Directive, directives

from sphinx import addnodes
from sphinx.util import parselinenos
from sphinx.util.nodes import set_source_info

class AmanziXMLInclude(Directive):
    """
    Based on the Sphinx ``literalinclude``, this derctive customizes output for 
    showing snippets of Amanzi XML files in out Tutorials.

    Like ``.. include:: :literal:``, but only warns if the include file is
    not found, and does not raise errors.  Also has several options for
    selecting what to include.
    """

    has_content = False
    required_arguments = 1
    optional_arguments = 0
    final_argument_whitespace = False

    option_spec = {
        'linenos': directives.flag,
        'encoding': directives.encoding,
        'lines': directives.unchanged_required,
        'string': directives.unchanged_required,
        'indent': directives.unchanged_required,
        'start': directives.unchanged_required,
        'end': directives.unchanged_required,
        'prepend': directives.unchanged_required,
        'append': directives.unchanged_required,
    }

    def run(self):

        # grab state of the document - not sure all this holds
        document = self.state.document

        # Ensure lines of the XML file can be inserted
        if not document.settings.file_insertion_enabled:
            return [document.reporter.warning('File insertion disabled',
                                              line=self.lineno)]
        env = document.settings.env
        rel_filename, filename = env.relfn2path(self.arguments[0])

        #
        # It looks like lines and start,end,append, and prepend do not
        # work together??  Should probably add more tests here.
        # 

        # extract the encoding 
        encoding = self.options.get('encoding', env.config.source_encoding)
        # set the codec
        codec_info = codecs.lookup(encoding)

        # Read the xml file
        try:
            f = codecs.StreamReaderWriter(open(filename, 'rb'),
                    codec_info[2], codec_info[3], 'strict')
            lines = f.readlines()
            f.close()
        except (IOError, OSError):
            return [document.reporter.warning(
                'Include file %r not found or reading it failed' % filename,
                line=self.lineno)]
        except UnicodeError:
            return [document.reporter.warning(
                'Encoding %r used for reading included file %r seems to '
                'be wrong, try giving an :encoding: option' %
                (encoding, filename))]

        #    
        # Collect lines specified in the "lines" option
        #
        linespec = self.options.get('lines')
        if linespec is not None:
            try:
                linelist = parselinenos(linespec, len(lines))
            except ValueError as err:
                return [document.reporter.warning(str(err), line=self.lineno)]
            # just ignore nonexisting lines
            nlines = len(lines)
            lines = [lines[i] for i in linelist if i < nlines]
            if not lines:
                return [document.reporter.warning(
                    'Line spec %r: no lines pulled from include file %r' %
                    (linespec, filename), line=self.lineno)]


        #        
        # Collect lines from section specified with [start,end]
        #
        start    = self.options.get('start')
        end      = self.options.get('end')
        prepend  = self.options.get('prepend')
        append   = self.options.get('append')

        if start is not None or end is not None:
            use = not start
            res = []
            for line in lines:
                if not use and start and start in line:
                    use = True
                    res.append(line)
                elif use and end and end in line:
                    use = False
                    res.append(line)
                    break
                elif use:
                    res.append(line)
            lines = res

        # Prepend and Append if options used
        if prepend:
           lines.insert(0, prepend + '\n')
        if append:
           lines.append(append + '\n')
        
        #
        #  Adjust indent of single line strings
        #
        indent=self.options.get('indent')
        if indent:
            instr=" "*int(indent)
        else:
            instr=""

        #
        #  Collect the first line containing the string
        # 
        line_string = self.options.get('string')
        if line_string:
            result=[]
            for line in lines:
                if line_string in line:
                    if indent:
                        result.append(instr+line.lstrip())
                    else:
                        result.append(line)
                    break
            lines = result

        # Join lines into a single string
        text = ''.join(lines)

        # Insert text in docutils node
        retnode = nodes.literal_block(text, text, source=filename)
        
        # Set some info for Sphinx (not sure what this really does).
        set_source_info(self, retnode)

        # Language: xml
        retnode['language']='xml'

        # Display line numbers
        if 'linenos' in self.options:
            retnode['linenos'] = True

        # Record dependency for Sphinx build 
        # (note: env = document.settings.env from sphinx)
        env.note_dependency(rel_filename)

        # return
        return [retnode]

# 
#   Add extension to Sphinx
# 
def setup(app):
    app.add_directive('amanzi_xml_include', AmanziXMLInclude)
