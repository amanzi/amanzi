#include <stdio.h>
 
#include <xercesc/util/XMLString.hpp>
#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/framework/LocalFileInputSource.hpp>
#include <xercesc/sax/ErrorHandler.hpp>
#include <xercesc/sax/SAXParseException.hpp>
#include <xercesc/validators/common/Grammar.hpp>

class ErinErrorHandler : public ErrorHandler
{
    private:
        void reportParseException(const SAXParseException& ex)
        {
            char* msg = XMLString::transcode(ex.getMessage());
            fprintf(stderr, "at line %llu column %llu, %s\n",
                    ex.getColumnNumber(), ex.getLineNumber(), msg);
            XMLString::release(&msg);
        }
 
    public:
        void warning(const SAXParseException& ex)
        {
            reportParseException(ex);
        }
 
        void error(const SAXParseException& ex)
        {
            reportParseException(ex);
        }
 
        void fatalError(const SAXParseException& ex)
        {
            reportParseException(ex);
        }
 
        void resetErrors()
        {
        }
};
