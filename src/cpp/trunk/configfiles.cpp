
#if HAVE_CONFIG_H
    #include <config.h>
#endif

#include "configfiles.h"

void
ReadParametersFromStream(istream& inStream, map<string, string>& varDict)
{   
    string s;
    string::size_type s_equal;

    while ( inStream ) {
        getline(inStream, s);

        if (inStream.fail() && ! inStream.eof())
            break;

        if (s[0] == '#' || s[0] == '/') {
            // Skip comments, loosely defined as any line beginning by # or /
            continue;
        } else if ( (s_equal = s.find("=")) != string::npos) {
            string var(s.begin(),s.begin()+s_equal);
            string value(s.begin()+s_equal+1,s.end());

            if (VariableExists(varDict, var)) {
                cerr << "Warning: variable " << var << " defined twice (second definition:" << value <<" skipped)\n";
            } else {
				varDict[var] = value;
			}
            
        } else if (!s.empty()) {
            cerr << "Warning: line invalid " << s << "\n";
        }

    }
    
    if (inStream.fail() && ! inStream.eof())
        throw runtime_error("Processing of file failed.");

    return;
}

void
WriteParametersToStream(ostream& inStream, map<string, string >& varDict)
{   
 
    map<string, string, less<string> >::iterator pos;

    for (pos =  varDict.begin(); pos != varDict.end(); ++pos) {
        inStream << pos->first << "=" << pos->second << endl;
    }
    
}

void
WriteParametersToXMLStream(ostream& inStream, map<string, string >& varDict, char * inXMLElement)
{   
	
	string xmlelement;
	
	if (inXMLElement == 0) {
		xmlelement.assign("variable");
	} else {
		xmlelement.assign(inXMLElement);
	}
	
    map<string, string, less<string> >::iterator pos;
	
    for (pos =  varDict.begin(); pos != varDict.end(); ++pos) {
        inStream << "<" << xmlelement << " name=\""<< pos->first << "\">" << pos->second << "</" << xmlelement << ">" << endl;
    }
    
}

bool
VariableExists(map<string, string>& varDict, string inName)
{
    map<string,string>::iterator p = varDict.find(inName);
    if (p != varDict.end()) {
        return true;
    } else {
        return false;
    }

}
