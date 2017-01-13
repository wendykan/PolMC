#include <istream>
#include <vector>
#include <valarray>

using namespace std;


class cubicspline {

    public:
    
        cubicspline();
        cubicspline(string inFilename);
        cubicspline(string inFilename, double inYp1, double inYpn);
        cubicspline(istream& inStream);

        ~cubicspline();
                void init();
                void enableCache() { mCacheEnabled = true; }
                void disableCache() { mCacheEnabled = false; }
        void fromstream(istream& inStream);

        void attachFile(string inFilename);

        void            insertpoint(double inX, double inY);
        double           y(const double& inY);
        valarray<double> y(valarray<double>& inY);

        double   dydx(double inY);
        double   d2ydx2(double inY);

#ifndef __RANGE_CHECK
        size_t  size()  { return mSize; }
        double  xelement(size_t i)  { return mX[i]; }
        double  yelement(size_t i)  { return mY[i]; }
#else
        size_t  size()  { return mX.size(); }
        double  xelement(size_t i)  { return mX.at(i); }
        double  yelement(size_t i)  { return mY.at(i); }
#endif      
        void    initSplineTables();
        void    clearTables();
        
    protected:

#ifndef __RANGE_CHECK
        double *mX;  
        double *mY;  
        double *mY2; 
                size_t mSize;
#else
        vector<double>   mX; 
        vector<double>   mY; 
        vector<double> mY2;  
#endif

        bool    mInitialized;
        bool    mIsNatural;
        bool    mCacheEnabled;
        bool    mCacheInitialized;
                bool    mXisEquallySpaced;
                double  mDX;
                double  cX;
        double  cY;

                double  mYp1;
                double  mYpn;
};
