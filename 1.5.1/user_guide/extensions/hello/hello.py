
from docutils import nodes
from docutils.parsers.rst import Directive

class HelloWorld(Directive):
    """
    This simple directive is intended to help us learn about creating
    sphinx extensions.  To try it out! In your rst document include,

    .. helloworld::

       Your message content here.
       It can have more than one line.

    Then the processed document will include a block containing the text

    Hello World
    Your message content here.
    It can have more than one line.

    It's not all that exciting, but it is a good way to get started.
    """

    # Enable Content in the Directive
    has_content = True

    # No arguments yet
    required_arguments = 0
    optional_arguments = 0

    # No options yet
    option_spec = {
    }

    def run(self):

        message_title = "Hello World"
        # Content is a list of unicode strings, create a body list
        message_body_l = [message_title]
        for str in self.content:
            message_body_l.append(str)
        # Convert list to a tuple so we can join into single string
        message_body_t=tuple(message_body_l)
        message_body = u'\n'.join(message_body_t)

        # Create a docutils node 
        content = nodes.literal_block(message_body, message_body)

        return [content]

def setup(app):
    app.add_directive('helloworld', HelloWorld)
