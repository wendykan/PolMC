
#if HAVE_CONFIG_H
#include <config.h>
#endif

#include "XMLUtil.h"
#include <sstream>

using namespace std;

int
XMLUtil::extractSingleXPath(xmlDocPtr doc, xmlChar* xpath, char* outString) {
    
    strcpy(outString,"");
    
    xmlNodeSetPtr nodeset;
    xmlXPathObjectPtr result;
    xmlChar *keyword;
    
    xmlXPathContextPtr context;
    
    if ( (context = xmlXPathNewContext(doc)) == NULL ) {
        return kCantAllocateXMLContext;
    }
    
    if ( (result = xmlXPathEvalExpression(xpath, context)) == NULL ) {
        return kCantProcessXpath;
    }
    
    if(xmlXPathNodeSetIsEmpty(result->nodesetval)){
        xmlXPathFreeObject (result);
        return kXPathResultIsEmpty;
    }
    xmlXPathFreeContext(context);
    
    if (result) {
        nodeset = result->nodesetval;
        
        if (nodeset->nodeNr == 1) {
            xmlNodePtr theNode = nodeset->nodeTab[0];
            
            if (theNode->children != NULL )
                strcpy(outString,(char*)theNode->children->content);
            else
                strcpy(outString,(char*)theNode->content);
        } else {
            return kTooManyNodesReturned;
        }
        xmlXPathFreeObject (result);
        xmlCleanupParser();
        
    }
    
    return kNoErr;
}

int
XMLUtil::extractSingleXPathInNewBuffer(xmlDocPtr doc, xmlChar* xpath, char* *outBuffer, long* outLen) {
    
    xmlNodeSetPtr nodeset;
    xmlXPathObjectPtr result;
    xmlChar *keyword;
    
    xmlXPathContextPtr context;
    
    if ( (context = xmlXPathNewContext(doc)) == NULL ) {
        return kCantAllocateXMLContext;
    }
    
    if ( (result = xmlXPathEvalExpression(xpath, context)) == NULL ) {
        return kCantProcessXpath;
    }
    
    xmlXPathFreeContext(context);
    
    if(xmlXPathNodeSetIsEmpty(result->nodesetval)){
        xmlXPathFreeObject (result);
        return kXPathResultIsEmpty;
    }
    
    if (result) {
        nodeset = result->nodesetval;
        
        if (nodeset->nodeNr == 1) {
            xmlNodePtr theNode = nodeset->nodeTab[0];
            
            if (theNode->children != NULL ) {
                *outBuffer = (char*)malloc(strlen((char*)theNode->children->content)*sizeof(char)+1);
                strcpy(*outBuffer,(char*)theNode->children->content);
            } else {
                *outBuffer = (char*)malloc(strlen((char*)theNode->content)*sizeof(char)+1);
                strcpy(*outBuffer,(char*)theNode->content);
            }
        } else {
            return kTooManyNodesReturned;
        }
        xmlXPathFreeObject (result);
        xmlCleanupParser();
        
    }
    
    
    return kNoErr;
}


int
XMLUtil::extractArrayOfStringsXPath(xmlDocPtr doc, xmlChar* xpath, char* * outArray, long * outHowMany) {
    
    xmlNodeSetPtr nodeset;
    xmlXPathObjectPtr result;
    xmlChar *keyword;
    long i;
    
    xmlXPathContextPtr context;
    
    if ( (context = xmlXPathNewContext(doc)) == NULL ) {
        return kCantAllocateXMLContext;
    }
    
    if ( (result = xmlXPathEvalExpression(xpath, context)) == NULL ) {
        return kCantProcessXpath;
    }
    
    xmlXPathFreeContext(context);
    
    if(xmlXPathNodeSetIsEmpty(result->nodesetval)){
        outArray = NULL;
        outHowMany = 0;
        xmlXPathFreeObject (result);
        return kXPathResultIsEmpty;
    }
    
    if (result) {
        nodeset = result->nodesetval;
        *outHowMany = nodeset->nodeNr;
        for ( i = 0; i < *outHowMany; ++i ) {
            xmlNodePtr theNode = nodeset->nodeTab[0];
            
            if (theNode->children != NULL ) {
                outArray[i] = (char*)malloc(strlen((char*)theNode->children->content)*sizeof(char)+1);
                strcpy((outArray[i]),(char*)theNode->children->content);
            } else {
                outArray[i] = (char*)malloc(strlen((char*)theNode->content)*sizeof(char)+1);
                strcpy((outArray[i]),(char*)theNode->content);
            }
        }
        
        xmlXPathFreeObject (result);
        xmlCleanupParser();
        
    }
    
    return kNoErr;
}

/* Encodes a buffer (buf) of length l into a newly created buffer outBuffer. Final length of buffer is returned
through outLength.  User is responsible for freeing outBuffer.
*/

int
XMLUtil::encode_base64(unsigned char * buf, long l, char* * outBuffer, long *outLength)
{
	/* Base64 encoding (transform three 8-bit numbers into four 6-bit numbers (from 0 to 63)
	then map this number to a letter A-Z a-z,0-9 + and -.
	When the number of characters to be encoded is not a multiple of three, pad with '='.
	*/
	char *encode = (char*) "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

	long i,j;
	
	short c1,c2,c3;
	long b1,b2,b3,b4;
	
	*outLength = l % 3 == 0 ? l / 3 * 4 : (l/3 + 1) * 4; 
	if ( (*outBuffer = (char *)malloc(*outLength + 1)) == NULL ) {
		return kOutOfMemoryErr;
	}
	
	j = 0;
	for (i = 0; i < l; i+=3) {
		c1 = buf[i];
		if ( i+1 < l)
			c2 = buf[i+1];
		else
			c2 = 512;
		
		if ( i+2 < l)
			c3 = buf[i+2];
		else
			c3 = 512;
		
		b1 = c1 >> 2;
		b2 = ((c1 & 0x03) << 4) | ((c2 & 0xf0) >> 4);
		b3 = ((c2 & 0x0f) << 2) | ((c3 & 0xc0) >> 6);
		b4 = (c3 & 0x3f);
		
		(*outBuffer)[j++] = encode[b1];
		(*outBuffer)[j++] = encode[b2];
		(*outBuffer)[j++] = c2 != 512 ? encode[b3]: '=';
		(*outBuffer)[j++] = c3 != 512 ? encode[b4]: '=';
	}
	(*outBuffer)[j] = 0;
	
	if (j != *outLength) {
		clog << "j is not calculated length. oops in base64 encoding" << endl;
		*outLength=j;
	}
	
	return kNoErr;
}

/* Decodes a buffer (buf) of length l into a newly created buffer outBuffer. Final length of buffer is returned
through outLength.  User is responsible for freeing outBuffer.
*/

int
XMLUtil::decode_base64(unsigned char * buf, long l, char* * outBuffer, long *outLength)
{
	int decode[256] = { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
		-1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
		-1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
		-1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
		-1, -1, -1, 62, -1, -1, -1, 63, 52, 53,
		54, 55, 56, 57, 58, 59, 60, 61, -1, -1,
		-1, 512,-1, -1, -1, 0,  1,  2,  3,  4,
		5,  6,  7,  8,  9, 10, 11, 12, 13, 14,
		15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 
		25, -1, -1, -1, -1, -1, -1, 26, 27, 28,
		29, 30, 31, 32, 33, 34, 35, 36, 37,38,
		39, 40, 41, 42, 43, 44, 45, 46, 47,48,
		49, 50, 51, -1, -1, -1, -1, -1, -1,-1,
		-1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
		-1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
		-1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
		-1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
		-1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
		-1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
		-1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
		-1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
		-1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
		-1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
		-1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
		-1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
		-1, -1, -1, -1, -1, -1
	};
	
	long i,j;
	
	unsigned short c1,c2,c3;
	long b1,b2,b3,b4;
	
	if (l % 4 == 0) {
		*outLength = l / 4 * 3;
	} else {
		return kIncorrectBase64Length;
	}
	
	if ( (*outBuffer = (char *)malloc(*outLength)) == NULL ) {
		return kOutOfMemoryErr;
	}
	
	j = 0;
	for (i = 0; i < l; i+=4) {
		if ( (b1 = decode[buf[i]]) == -1)
			return kInvalidBase64Character;
		if ( (b2 = decode[buf[i+1]]) == -1)
			return kInvalidBase64Character;
		if ( (b3 = decode[buf[i+2]]) == -1)
			return kInvalidBase64Character;
		if ( (b4 = decode[buf[i+3]]) == -1)
			return kInvalidBase64Character;
		
		c1 = b1 << 2 | ((b2 & 0x30) >> 4);
		c2 = ((b2 & 0x0f) << 4) | ((b3 & 0x3c) >> 2);
		c3 = ((b3 & 0x03) << 6) | (b4 & 0x3f);
		
		if (b4 == 512) {
			if (i+3 == l-1) {
				(*outLength)--;
			} else {
				return kInvalidBase64Character; 
			}
		}
		if (b3 == 512) {
			if (i+3 == l-1 && b4 == 512) {
				(*outLength)--;
			} else {
				return kInvalidBase64Character; 
			}
		}
		
		(*outBuffer)[j++] = c1;
		(*outBuffer)[j++] = c2;
		(*outBuffer)[j++] = c3;
	}
	(*outBuffer)[*outLength]=0;
	
	return kNoErr;
}

