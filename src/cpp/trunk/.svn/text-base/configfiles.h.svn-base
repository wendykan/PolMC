#include <string>
#include <sstream>
#include <iostream>
#include <map>

#include <stdexcept>

using namespace std;

void ReadParametersFromStream(istream& inStream, map<string, string>& varDict);

void WriteParametersToStream(ostream& inStream, map<string, string>& varDict);
void WriteParametersToXMLStream(ostream& inStream, map<string, string>& varDict, char* inXMLElement = 0);

bool VariableExists(map<string, string>& varDict, string inName);

template<class T> void InitVariable(map<string, string>& varDict, string inName, T& ioVariable);
 
template<class T> void InitVariable(map<string, string>& varDict, string inName, T& ioVariable) 
{
    map<string,string>::iterator p = varDict.find(inName);
    
    if (p != varDict.end()) {
        istringstream s(varDict[inName]);
        s >> ioVariable ;
    } else {
        throw runtime_error("variable \""+inName+"\" undefined");
    }
}
