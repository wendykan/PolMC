#ifndef __MCGENERICOBJECT
#define __MCGENERICOBJECT

#include "MCObject.h"

#if HAVE_CONFIG_H
#include <config.h>
#endif

#if HAVE_LIBXML2
#include <libxml/parser.h>
#include <libxml/xpath.h>
#endif


class MCObject;

class MCGenericObject : public MCObject {

public:
    MCGenericObject() {}
    MCGenericObject(map<string,string>& inDict):MCObject(inDict) { init(inDict); }

    void init(map<string,string>& inDict);
    virtual bool	IsOutsideObject(RealV& inLocalPosition); 


    MCObject* BuildObjectFromFile(string inFilename);
    MCObject* BuildObjectFromXMLFile(string inFilename);
    long BuildObjectFromOBJStream(istream& input);
    void      split(const string& str,
                    vector<string>& tokens,
                    const string& delimiters);
        
#if HAVE_LIBXML2
    bool	ReadRealVFromString(string inString,  RealV& outVector);
#endif

protected:
    
};

#endif

