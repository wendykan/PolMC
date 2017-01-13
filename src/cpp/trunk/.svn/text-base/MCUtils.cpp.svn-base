#include <cmath>
#include "MCUtils.h"
#include "mydebug.h"

using namespace std;


double
MinBinBoundary(long inIndex, double inMin, double inMax, long inN)
{
    if (inIndex == -1)
        return inMin;

    return inMin + inIndex * BinWidth_(inMin, inMax, inN);
}

double
MaxBinBoundary(long inIndex, double inMin, double inMax, long inN)
{
    if (inIndex == -1)
        return inMax;
    
    return inMin + (inIndex+1) * BinWidth_(inMin, inMax, inN);
    
}

double
BinCenter(long inIndex, double inMin, double inMax, long inN)
{
  return (MinBinBoundary(inIndex, inMin, inMax, inN) + MaxBinBoundary(inIndex, inMin, inMax, inN))/2.;
}

void 
split(const string& str,
   vector<string>& tokens,
   const string& delimiters)
{
    // Skip delimiters at beginning.
    string::size_type start = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type current     = str.find_first_of(delimiters, start);
    
    if (string::npos == start) 
        return;
    
    while (string::npos != current || string::npos != start)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(start, current - start));
        start = current != string::npos ? current+1 : string::npos;
        // Find next "non-delimiter"
        current = str.find_first_of(delimiters, start);
    }
}

void 
trim(string& str)
{
    long first = 0;
    long last = str.size()-1;

    while (str[first] == ' ' || str[first] == '\t' || str[first] == '\n') {
        first++;
    }
    
    while (str[last] == ' ' || str[last] == '\t' || str[last] == '\n') {
        last--;
    }

    string temp;
    for (long i = first; i <= last; i++ ) {
        temp.push_back(str[i]);
    }

    str = temp;
}

string
removeBracketedProperty(string& str, char inLeft, char inRight)
{
    string::size_type first     = str.find_first_of(inLeft);
    string::size_type last     = str.find_first_of(inRight, first + 1);

	if (first == string::npos || last == string::npos) 
		return string("");
	
    string out, main;
	
    for (long i = first + 1; i < last; i++ ) {
        out.push_back(str[i]);
    }

    for (long i = 0; i < first; i++ ) {
        main.push_back(str[i]);
    }

    for (long i = last + 1; i < str.size(); i++ ) {
        main.push_back(str[i]);
    }
	
	str = main;
	
	return out;
}
