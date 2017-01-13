#ifndef __FASTINTERPOLATE_H
#define __FASTINTERPOLATE_H

#include <istream>
#include <vector>
#include <valarray>

using namespace std;


class fastinterpolate {

public:

    fastinterpolate();
    fastinterpolate(string inFilename);
    fastinterpolate(istream& inStream);
    fastinterpolate(double *inX, double* inY, long inN);
    fastinterpolate(const fastinterpolate& inRhs); // Copy
    fastinterpolate& operator=(const fastinterpolate& inRhs); // Assignment
    
    ~fastinterpolate();
    void init();
    void enableCache() ;
    void disableCache() ;
    void enableApproximation();
    void disableApproximation();
    void fromstream(istream& inStream);
    void tostream(ostream& inStream);

    void attachFile(string inFilename);

    void RemoveDuplicates();
    
    void CheckEquallySpaced();
    void            insertpoints(double *inX, double* inY, long inN);
    void            insertpoint(double inX, double inY);
    double           y(const double& inY);
    valarray<double> y(valarray<double>& inY);
    void	NormalizeX(double inNorm);
    void	NormalizeY(double inNorm);
    double * xArray();
    double * yArray();
    size_t  size() ;
    double  xelement(size_t i) ;
    double  yelement(size_t i) ;
    double  yMaximum();
    double  xMaximum();
    void   clearTables();
    
protected:

    double *mX;
    double *mY;
    size_t mSize;
    size_t mCapacity;
    

    bool    mInitialized;
    bool    mCacheEnabled;
    bool    mCacheInitialized;
    bool    mXisEquallySpaced;
    bool    mStartsAtZeroAndEquallySpaced;
    bool    mApproximate;
    double  mDX;
    double  cX;
    double  cY;
    double  cYMaximum;
    double  cXMaximum;
};

#endif

