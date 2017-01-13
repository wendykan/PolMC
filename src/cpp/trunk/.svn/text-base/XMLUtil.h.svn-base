#ifndef XMLUTIL_H
#define XMLUTIL_H
#if HAVE_CONFIG_H
#include "config.h"
#endif

#include <libxml/parser.h>
#include <libxml/xpath.h>
#include <sstream>
#include <iostream>
#include <string>
#include "fastinterpolate.h"

using namespace std;

enum {
    kNoErr,
    kCantAllocateXMLContext,
    kCantProcessXpath,
    kXPathResultIsEmpty,
    kTooManyNodesReturned,
	kOutOfMemoryErr,
	kIncorrectBase64Length,
	kInvalidBase64Character
};

namespace XMLUtil {
#if HAVE_LIBXML2
    int extractSingleXPath(xmlDocPtr doc, xmlChar* xpath, char* outString);
    int	extractSingleXPathInNewBuffer(xmlDocPtr doc, xmlChar* xpath, char* *outBuffer, long* outLen) ;
    int	extractArrayOfStringsXPath(xmlDocPtr doc, xmlChar* xpath, char* * outArray, long * outHowMany);
    
    template<class T> bool InitVariable(xmlXPathContextPtr inContext, char* inXPath, T& ioVariable);
    
    template<class T> bool InitVariable(xmlXPathContextPtr inContext, char* inXPath, T& ioVariable) 
    {
        xmlXPathObjectPtr pobj = NULL;
        xmlChar* strPtr;
        istringstream input;

        pobj = xmlXPathEvalExpression(BAD_CAST(inXPath), inContext);
        
        if (pobj == NULL)
            return true;

        if(xmlXPathNodeSetIsEmpty(pobj->nodesetval)){
            xmlXPathFreeObject (pobj);
            return true;
        }
        
        strPtr = xmlXPathCastNodeSetToString(pobj->nodesetval);
        
        input.str((char*)strPtr);
        xmlFree(strPtr);
        xmlXPathFreeObject(pobj);
        input >> ioVariable;
        
        return false;
    }

    template<> inline bool InitVariable(xmlXPathContextPtr inContext, char* inXPath, fastinterpolate& ioVariable) 
    {
        xmlXPathObjectPtr pobj = NULL;
        xmlChar* strPtr;
        istringstream input;
        
        pobj = xmlXPathEvalExpression(BAD_CAST(inXPath), inContext);
        
        if (pobj == NULL)
            return true;

        if(xmlXPathNodeSetIsEmpty(pobj->nodesetval)){
            xmlXPathFreeObject (pobj);
            return true;
        }
        
        strPtr = xmlXPathCastNodeSetToString(pobj->nodesetval);
        
        input.str((char*)strPtr);
        xmlFree(strPtr);
        xmlXPathFreeObject(pobj);
        ioVariable.fromstream(input);
        
        if (ioVariable.size() == 0) {
            pobj = xmlXPathEvalExpression(BAD_CAST(".//name/text()"), inContext);
            strPtr = xmlXPathCastNodeSetToString(pobj->nodesetval);

            clog <<  "Warning: table " << string(inXPath) << " for object named " << string( (char*)strPtr) << " is empty" << endl;
            xmlFree(strPtr);
            xmlXPathFreeObject(pobj);

            return true;
        }
        
        return false;
    }

    template<> inline bool InitVariable(xmlXPathContextPtr inContext, char* inXPath, string& ioVariable) 
    {
        xmlXPathObjectPtr pobj = NULL;
        xmlChar* strPtr;
        
        pobj = xmlXPathEvalExpression(BAD_CAST(inXPath), inContext);
        
        if (pobj == NULL)
            return true;
        
        if(xmlXPathNodeSetIsEmpty(pobj->nodesetval)){
            xmlXPathFreeObject (pobj);
            return true;
        }
        
        strPtr = xmlXPathCastNodeSetToString(pobj->nodesetval);
        
        ioVariable = string((char*)strPtr);
        xmlFree(strPtr);
        xmlXPathFreeObject(pobj);
        
        return false;
    }
    
	int decode_base64(unsigned char * buf, long l, char* * outBuffer, long *outLength);
	int encode_base64(unsigned char * buf, long l, char* * outBuffer, long *outLength);

#endif

    
}
#endif
