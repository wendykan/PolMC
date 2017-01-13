#ifndef __UTILS_H
#define __UTILS_H

#include <vector>
#include <string>
//#include <ppc_intrinsics.h> 

using namespace std;

//long    WhichBin(const double& inValue, const double& inMin, const double& inMax, const long inN);
//long	WhichBin(double inValue, double inMin, double inMax, long inN);

#define BinWidth_(inMin, inMax, inN) ((inMax-inMin)/inN)
#define To1DIndexFrom2D_(i, j, NX) (j*(NX) + i )
#define To1DIndexFrom3D_(i, j, k, NX, NY) (k*(NX*NY)+ j*(NX) + i )
#define To1DIndexFrom4D_(i, j, k, l, NX, NY, NZ) (l*(NX*NY*NZ) + k*(NX*NY)+ j*(NX) + i )
#define To3DXIndex_(u, NX, NY, NZ) ( ( u % (NX * NY) ) % NX )
#define To3DYIndex_(u, NX, NY, NZ) ( ( u % (NX * NY) ) / NX )
#define To3DZIndex_(u, NX, NY, NZ) ( u / (NX * NY) )

inline long
WhichBin(const double& inValue, const double& inMin, const double& inMax,const long& inN)
{
    // Return bin index, 0-based (i.e. 0 to inN-1)
    // Returns -1 if out of range
    
	//double r = __fres(inMax-inMin); 
	//r *= (inValue-inMin);
    double r = (inValue-inMin)/(inMax-inMin);
    
	
	// This bizarre !(r<0) is a trick to compare faster on most machines
    if ( ! (r < 0.) && r < 1. ) {
        return long(r *= inN);
    } else if (inValue == inMax) {
        return inN-1;
    } else if (inValue == inMin) {
		return 0;
    }
	
    return -1;
}
/*
inline long
WhichBin(double inValue, double inMin, double inMax, long inN)
{
    // Return bin index, 0-based (i.e. 0 to inN-1)
    // Returns -1 if out of range
    
    if (inValue == inMax) {
        return inN-1;
    } else if (inValue == inMin) {
		return 0;
    }
    
    long r = long((inValue-inMin)/BinWidth_(inMin, inMax, inN));
    
    if ( r < 0 || r >= inN ) {
        return -1;
    } else if ( r == inN ) {
        // There sometimes is a round off error that makes
        // r be inN although we have checked for inValue== inMax above.
        // There might be a better way to do this, but I simply 
        // correct here
        return r-1;
    }
    
    return r;
}

*/
double	BinWidth(double inMin, double inMax , long inN);
double	MinBinBoundary(long inIndex, double inMin, double inMax, long inN);
double	MaxBinBoundary(long inIndex, double inMin, double inMax, long inN);
double  BinCenter(long inIndex, double inMin, double inMax, long inN);
void split(const string& str, vector<string>& tokens, const string& delimiters);
void trim(string& str);
string removeBracketedProperty(string& str, char inLeft, char inRight);

#endif


