
#if HAVE_CONFIG_H
    #include <config.h>
#endif

#include <iostream>
#include <fstream>
#include <string>
#include "cubicspline.h"
#include <stdexcept>

#include "mydebug.h"

using namespace std;

cubicspline::cubicspline()
{
    init();
}

cubicspline::cubicspline(string inFilename)
{
    init();
        attachFile(inFilename);
}

cubicspline::cubicspline(string inFilename, double inYp1, double inYpn)
{
    init();
        mIsNatural = false;
        mYp1 = inYp1;
        mYpn = inYpn;
        
    attachFile(inFilename);
}


cubicspline::cubicspline(istream& inStream)
{
    init();
    fromstream(inStream);
}

void
cubicspline::init()
{
    mInitialized = false;
        mIsNatural = false;
        mCacheEnabled = false;
        mXisEquallySpaced = false;
        mDX = 0.;
        mX = 0;
        mY = 0;
        mY2 = 0;
        mSize = 0;

}

void
cubicspline::fromstream(istream& inStream)
{
    double x,y;
    
    try {
            while (inStream >> x) {
                inStream >> y;
                insertpoint(x,y);
            }
    } catch(...) {
        cerr << "error reading file";
    }

    initSplineTables();
}


cubicspline::~cubicspline()
{

}

void
cubicspline::attachFile(string inFilename)
{

    ifstream stream(inFilename.c_str());

    if (stream.fail() ) {
        throw runtime_error("unable to open file");
    }
    
    
    clearTables();
    
    fromstream(stream);

}

void
cubicspline::clearTables()
{
#ifndef __RANGE_CHECK
    delete mX;
    delete mY;
    delete mY2;
#else
    mX.clear();
    mY.clear();
    mY2.clear();
#endif
}

void
cubicspline::insertpoint(double inX, double inY)
{

#ifndef __RANGE_CHECK   
        double* newX = new double[mSize+1];
        double* newY = new double[mSize+1];
        
        size_t i = 0;
        size_t where = mSize;
        size_t found = 0;
        while (i < mSize) {
            newX[i+found] = mX[i];
            newY[i+found] = mY[i];

            if (mX[i] > inX)  {
                where = i;
                found = 1;
            }
            ++i;
        }
        
        newX[where] = inX;
        newY[where] = inY;

        delete mX;
        delete mY;

//        mX = new double[mSize+1];
//        mY = new double[mSize+1];
        mSize++;
        
        mX = newX;
        mY = newY;
#else
    // the complicated way for the STL to find something.
    vector<double>::iterator wherex = find_if(mX.begin(), mX.end(), bind1st(less<double>(), inX) );

    // We should enforce that no X value is used twice.  Algorithm would fail.
    // I don;'t know how to do this properly here, I'll do it in initSpline
/*  if (mX.at( distance(mX.begin(), wherex)) == inX) {
        cout << mX.at( wherex ) << " " << inX << endl;
        cerr << "Same X used twice.  Aborting.";
        throw;
    }
*/
    // Must insert vy first (par. 16.3.6) because wherex is invalid after insert()
    mY.insert(mY.begin() + distance(mX.begin(), wherex), inY);
    mX.insert(wherex, inX);
#endif




}

void
cubicspline::initSplineTables()
{
    cX = cY = NAN;
    mCacheEnabled = false;
        mCacheInitialized = false;

    // setting up the tables for cubic spline interpolation (second derivative zero at edges)
    // (Numerical recipes, section 3.3)
#ifndef __RANGE_CHECK
        delete mY2;
        mY2 = new double[mSize];
        double* u = new double[mSize];
    size_t n = mSize;
        
    // Check to see if same X used twice.  If so, abort.
    for (long k = 0; k < n - 1; ++k) {
        if (mX[k] == mX[k+1]) {
            cerr << "Same X used twice.  Check element [" << k << "]. Aborting.";
            throw logic_error("same X used twice in cubicspline.");
        }
    }

        if (mIsNatural) {
            mY2[0] = 0;
            u[0] = 0;
        } else {
            mY2[0] = -0.5;
            u[0] = (3./(mX[1]-mX[0])) * ((mY[1]-mY[0])/(mX[1]-mX[0])-mYp1);
        }

    for (size_t i = 1; i <= n - 2; i++) {

        double sig = (mX[i] - mX[i-1]) / (mX[i+1] - mX[i-1]);
        double p = sig * mY2[i-1] + 2;
        mY2[i] = (sig - 1.)/p;
        u[i] =   (mY[i+1] - mY[i]) / (mX[i+1] - mX[i]) - (mY[i] - mY[i-1]) / (mX[i] - mX[i-1]);
        u[i] = (6. * u[i] / (mX[i+1] - mX[i-1]) - sig * u[i-1]) / p;

    }

        double qn,un;
        
        if (mIsNatural) {
            un = 0.;
            qn = 0.;
        } else {
            qn = 0.5;
            un = (3./(mX[n-1]-mX[n-2])) * (mYpn - (mY[n-1]-mY[n-2])/(mX[n-1]-mX[n-2]));
        }

        mY2[n-1] = (un-qn*u[n - 2]) * (qn * mY2[n-2] + 1.);
    
    for (long k = n - 2; k >= 0; --k) {
        mY2[k] = mY2[k] * mY2[k+1] + u[k];
    }
        
        delete u;

#else
    size_t n = mY.size();
    mY2.resize(n);
    vector<double> u(n);
    
    // Check to see if same X used twice.  If so, abort.
    for (long k = 0; k < n - 1; ++k) {
        if (mX.at(k) == mX.at(k+1)) {
            cerr << "Same X used twice.  Check element [" << k << "]. Aborting.";
            throw logic_error("same X used twice in cubicspline.");
        }
    }

        if (mIsNatural) {
            mY2.at(0) = 0;
            u.at(0) = 0;
        } else {
            mY2.at(0) = -0.5;
            u.at(0) = (3./(mX.at(1)-mX.at(0))) * ((mY.at(1)-mY.at(0))/(mX.at(1)-mX.at(0))-mYp1);
        }

    for (size_t i = 1; i <= n - 2; i++) {

        double sig = (mX.at(i) - mX.at(i-1)) / (mX.at(i+1) - mX.at(i-1));
        double p = sig * mY2.at(i-1) + 2;
        mY2.at(i) = (sig - 1.)/p;
        u.at(i) =   (mY.at(i+1) - mY.at(i)) / (mX.at(i+1) - mX.at(i)) - (mY.at(i) - mY.at(i-1)) / (mX.at(i) - mX.at(i-1));
        u.at(i) = (6. * u.at(i) / (mX.at(i+1) - mX.at(i-1)) - sig * u.at(i-1)) / p;

    }

        double qn,un;
        
        if (mIsNatural) {
            un = 0.;
            qn = 0.;
        } else {
            qn = 0.5;
            un = (3./(mX.at(n-1)-mX.at(n-2))) * (mYpn - (mY.at(n-1)-mY.at(n-2))/(mX.at(n-1)-mX.at(n-2)));
        }

        mY2.at(n-1) = (un-qn*u.at(n - 2)) * (qn * mY2.at(n-2) + 1.);
    
    for (long k = n - 2; k >= 0; --k) {
        mY2.at(k) = mY2.at(k) * mY2.at(k+1) + u.at(k);
    }
#endif
    mDX = mX[1]-mX[0];
    mXisEquallySpaced = false;
    // Check to see if X is equally spaced.
    for (long k = 1; k < mSize - 2; ++k) {
        if ( log10(abs(mX[k+1] - mX[k]) - fabs(mDX)) > -5.) {
            mXisEquallySpaced = false;
            clog << "Not equally spaced? : " << mX[k+1] - mX[k] << " != " << mDX << endl;
        }
    }
/*    if (mXisEquallySpaced)
        clog << RightNow() << "Equally spaced X in cubispline function dX = " << mDX <<  endl;
    else
        clog << RightNow() << "Not equally spaced X in cubispline function" << endl;
*/        
    mInitialized = true;    
}

/*
double
cubicspline::y(double inX)
{
    // interpolation routine (Numerical recipes, section 3.3)
    if (!mInitialized) {
        cerr << "Spline table not initialized" << endl;
        throw;
    }
    if (mX.size() == 0) {
        cerr << "Spline table empty" << endl;
        throw;
    }
    
    if (inX == cX && mCacheEnabled == true) {
        return cY;
    }
    
    size_t lo = 0, hi = mX.size() - 1;  
    size_t k;
    
    while (hi - lo > 1) {
        k = (hi + lo) / 2;
        
        if (mX.at(k) > inX)
            hi = k;
        else
            lo = k;
                
    }
        
    double h;
    h = mX.at(hi) - mX.at(lo);

    
    if (h == 0.) {
        cerr << "Error in spline interpolation (bad table: duplicate x values)\n";
        throw;
    } else if (inX > mX.at(mX.size() - 1) || inX < mX.at(0) ) {
        cerr << "Error in cubicspline: out of range " << inX << 
                " not between " << mX.at(0) << " and " << mX.at(mX.size()-1) << endl;
        throw;
    }
    

    double a = (mX.at(hi) - inX) / h;
    double b = (inX - mX.at(lo)) / h;

    cout << "Debug ON\n";
        cout << mX.at(hi) << " " << inX << " " <<  h << endl;;
        cout << inX << " " << mX.at(lo) << " " << b << endl;
        cout << a << " " << b << endl;

        cout << "my2lo" << mY2.at(lo) << endl;
        cout << "my2hi" << mY2.at(hi) << endl;
    cout << "Debug OFF\n";

    cX = inX;
    cY = a * mY.at(lo) + b * mY.at(hi) + ( (a*a*a - a) * mY2.at(lo) + (b*b*b - b) * mY2.at(hi) ) * (h * h)/6.;  
    return cY;
}
*/
double
cubicspline::y(const double& inX)
{
        
        size_t xsize,xsize_1;

        if ( ! mCacheEnabled ) {

            if (! mXisEquallySpaced) {
                // interpolation routine (Numerical recipes, section 3.3)
                xsize_1 = xsize = mSize;
                xsize_1--;
        
                if (inX > mX[xsize_1] || inX < mX[0] ) {
                    cerr << "Error in cubicspline: out of range " << inX << 
                                        " not between " << mX[0] << " and " << mX[xsize_1] << endl;
                        throw runtime_error("out of range value in cubicspline.");
                }
            
                size_t lo = 0;
                size_t hi = xsize_1;
                double a,b;
                double h;
                double xhi,xlo;
                
                size_t k;
 
                while (hi - lo > 1) {
                    k = (hi + lo) >> 1;
                    
                    if (mX[k] > inX)
                            hi = k;
                    else
                            lo = k;
                }

                xhi = mX[hi];
                xlo = mX[lo];
    
                h = xhi-xlo;
                a = (xhi - inX) / h;
                b = (inX - xlo) / h;

                return a * mY[lo] + b * mY[hi] + ( (a*a*a - a) * mY2[lo] + (b*b*b - b) * mY2[hi] ) * (h * h)/6.;
            } else {
                // interpolation routine (Numerical recipes, section 3.3)
                // for equally spaced absissa
                /*
                size_t hi = long((inX-mX[0])/mDX);
                size_t lo = hi++;

                if (hi > mSize-1 || lo < 0 ) {
                    cerr << "Error in cubicspline: out of range " << inX <<
                    " not between " << mX[0] << " and " << mX[mSize-1] << " (indices are " << hi << " and " << lo << ")" << endl;
                    throw runtime_error("out of range value in cubicspline.");
                }
                
                */
                double xnorm = (inX-*mX)/mDX;
                
                size_t hi = size_t(xnorm);
                size_t lo = hi++;

                double a = ( hi - xnorm) ;
                double b = (xnorm -  lo) ;
                
                return a * mY[lo] + b * mY[hi] + ( (a*a*a - a) * mY2[lo] + (b*b*b - b) * mY2[hi] ) * (mDX * mDX)/6.;

            }

        } else {
            // interpolation routine (Numerical recipes, section 3.3)
            if (!mInitialized) {
                    cerr << "Spline table not initialized" << endl;
                    throw runtime_error("Spline table not initialized in cubicspline.");
            }
    
            if (inX == cX && mCacheEnabled && mCacheInitialized) {
                    return cY;
            }
    
            if ( ( xsize_1 = xsize = mSize )== 0) {
                    cerr << "Spline table empty" << endl;
                    throw runtime_error("Spline table empty in cubicspline.");
            }
    
            xsize_1--;
    
            if (inX > mX[xsize_1] || inX < mX[0] ) {
                cerr << "Error in cubicspline: out of range " << inX << 
                                    " not between " << mX[0] << " and " << mX[xsize_1] << endl;
                    throw runtime_error("out of range value in cubicspline.");
            }
        
            size_t lo = 0;
            size_t hi = xsize_1;    
            size_t k;
            double h;
            double xhi,xlo;
    
            while (hi - lo > 1) {
                k = (hi + lo) >> 1;
                
                if (mX[k] > inX)
                        hi = k;
                else
                        lo = k;
                                    
            }
    
            xhi = mX[hi];
            xlo = mX[lo];
            
            h = xhi-xlo;
    
            
            if (h == 0.) {
                    cerr << "Error in spline interpolation (bad table: duplicate x values)\n";
                    throw runtime_error("duplicate x values in cubicspline.");
            }   
    
            double a = (xhi - inX) / h;
            double b = (inX - xlo) / h;
    
            cX = inX;
            cY = a * mY[lo] + b * mY[hi] + ( (a*a*a - a) * mY2[lo] + (b*b*b - b) * mY2[hi] ) * (h * h)/6.;  
            mCacheInitialized = true;
    
            return cY;
        
        }
}

valarray<double>
cubicspline::y(valarray<double>& inX)
{
    valarray<double> inY(inX.size());
    
    for (long i = 0; i < inY.size(); i++) {
        inY[i] = y(inX[i]);
    }
    
    return inY;
}


double
cubicspline::dydx(double inX)
{
    // interpolation routine (Numerical recipes, section 3.3)
    if (!mInitialized) {
        cerr << "Spline table not initialized" << endl;
        throw logic_error("Spline table not initialized in cubicspline.");
    }


#ifndef __RANGE_CHECK
    size_t lo = 0, hi = mSize - 1;  
    size_t k;

    while (hi - lo > 1) {
        k = (hi + lo) / 2;
        
        if (mX[k] > inX)
            hi = k;
        else
            lo = k;
                
    }
    
    double h = mX[hi] - mX[lo];
    
    if (h == 0.) {
        cerr << "Error in spline interpolation (bad table: duplicate x values)\n";
        throw runtime_error("duplicate x values in cubicspline.");
    }
    
    double a = (mX[hi] - inX) / h;
    double b = (inX - mX[lo]) / h;

    return (mY[hi] - mY[lo]) / h - (3*a*a - 1)/6 * h * mY2[lo] + (3*b*b - 1) / 6 * h * mY2[hi];
#else
    size_t lo = 0, hi = mX.size() - 1;  
    size_t k;

    while (hi - lo > 1) {
        k = (hi + lo) / 2;
        
        if (mX.at(k) > inX)
            hi = k;
        else
            lo = k;
                
    }
    
    double h = mX.at(hi) - mX.at(lo);
    
    if (h == 0.) {
        cerr << "Error in spline interpolation (bad table: duplicate x values)\n";
        throw runtime_error("duplicate x values in cubicspline.");
    }
    
    double a = (mX.at(hi) - inX) / h;
    double b = (inX - mX.at(lo)) / h;

    return (mY.at(hi) - mY.at(lo)) / h - (3*a*a - 1)/6 * h * mY2.at(lo) + (3*b*b - 1) / 6 * h * mY2.at(hi);
#endif

}

double
cubicspline::d2ydx2(double inX)
{
    // interpolation routine (Numerical recipes, section 3.3)
    if (!mInitialized) {
        cerr << "Spline table not initialized" << endl;
        throw runtime_error("Spline table not initialized in cubicspline.");

    }

#ifndef __RANGE_CHECK
    size_t lo = 0, hi = mSize - 1;  
    size_t k;

    while (hi - lo > 1) {
        k = (hi + lo) / 2;
        
        if (mX[k] > inX)
            hi = k;
        else
            lo = k;
                
    }
    
    double h = mX[hi] - mX[lo];
    
    if (h == 0.) {
        cerr << "Error in spline interpolation (bad table: duplicate x values)\n";
        throw runtime_error("duplicate x values in cubicspline.");

    }
    
    double a = (mX[hi] - inX) / h;
    double b = (inX - mX[lo]) / h;

    return a * mY2[lo] + b*mY2[hi];
#else   

    size_t lo = 0, hi = mX.size() - 1;  
    size_t k;

    while (hi - lo > 1) {
        k = (hi + lo) / 2;
        
        if (mX.at(k) > inX)
            hi = k;
        else
            lo = k;
                
    }
    
    double h = mX.at(hi) - mX.at(lo);
    
    if (h == 0.) {
        cerr << "Error in spline interpolation (bad table: duplicate x values)\n";
        throw runtime_error("duplicate x values in cubicspline.");

    }
    
    double a = (mX.at(hi) - inX) / h;
    double b = (inX - mX.at(lo)) / h;

    return a * mY2.at(lo) + b*mY2.at(hi);
#endif
}

