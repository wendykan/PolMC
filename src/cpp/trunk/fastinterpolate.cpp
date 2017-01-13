
#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <fstream>
#include <string>
#include "fastinterpolate.h"
#include <stdexcept>

#include "mydebug.h"

using namespace std;

fastinterpolate::fastinterpolate()
{
    init();
}

fastinterpolate::fastinterpolate(string inFilename)
{
    init();
    attachFile(inFilename);
}


fastinterpolate::fastinterpolate(istream& inStream)
{
    init();
    fromstream(inStream);
}

fastinterpolate::fastinterpolate(double *inX, double* inY, long inN)
{
    init();

    insertpoints(inX, inY, inN);

    CheckEquallySpaced();
    
}

fastinterpolate::fastinterpolate(const fastinterpolate& inRhs)
{
    init();

    insertpoints(inRhs.mX, inRhs.mY, inRhs.mSize);

    CheckEquallySpaced();
}

fastinterpolate&
fastinterpolate::operator=(const fastinterpolate& inRhs)
{
    clearTables();
    init();

    insertpoints(inRhs.mX, inRhs.mY, inRhs.mSize);

    CheckEquallySpaced();

    return *this;
}


void
fastinterpolate::CheckEquallySpaced()
{
    if (mSize >= 2) {
        mDX = mX[1]-mX[0];
        mXisEquallySpaced = true;
        // Check to see if X is equally spaced.
        for (long k = 1; k < mSize - 2; ++k) {
            if ( ! areTheSame(abs(mX[k+1] - mX[k]), abs(mDX), 4) ) {
                mXisEquallySpaced = false;
            }
        }

 //       #ifdef __MYDEBUG
        if (mXisEquallySpaced) {
            PrintMessageIfLevel_("Equally spaced X in fastinterpolate function dX = " << mDX, kVerbose );
        } else {
            PrintMessageIfLevel_( "Not equally spaced X in fastinterpolate function", kVerbose );
        }
//        #endif
        
        if (mX[0] == 0 && mXisEquallySpaced)
            mStartsAtZeroAndEquallySpaced = true;

    } else {
        mXisEquallySpaced = false;
        mStartsAtZeroAndEquallySpaced = false;
    }

    mInitialized = true;
}

void
fastinterpolate::init()
{
    mInitialized = false;
    mCacheEnabled = false;
    mXisEquallySpaced = false;
    mStartsAtZeroAndEquallySpaced = false;
    mApproximate = false;
    mDX = 0.;
    mCapacity = 1000;
    
    mX = new double[mCapacity];
    mY = new double[mCapacity];
    mSize = 0;

    cYMaximum = NAN;
    cXMaximum = NAN;
    cY = NAN;
    cX = NAN;
}

void
fastinterpolate::fromstream(istream& inStream)
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

    CheckEquallySpaced();    
}

void
fastinterpolate::tostream(ostream& inStream)
{
    for (long i = 0; i < mSize; i++) {
        inStream << mX[i] << "\t" << mY[i] << endl;
    }
}


fastinterpolate::~fastinterpolate()
{
    clearTables();
}

void
fastinterpolate::enableCache()
{
    mCacheEnabled = true; 
}

void
fastinterpolate::disableCache()
{
    mCacheEnabled = false; 
}

void
fastinterpolate::enableApproximation()
{
    mApproximate = true; 
}

void
fastinterpolate::disableApproximation()
{
    mApproximate = false; 
}

size_t
fastinterpolate::size()
{
    return mSize; 
}

double
fastinterpolate::xelement(size_t i)
{
    return mX[i]; 
}

double
fastinterpolate::yelement(size_t i)
{
    return mY[i]; 
}

void
fastinterpolate::attachFile(string inFilename)
{

    ifstream stream(inFilename.c_str());

    if (stream.fail() ) {
        throw runtime_error(string("unable to open file ")+inFilename);
    }

    clearTables();
    
    fromstream(stream);

}

void
fastinterpolate::clearTables()
{
    mSize = 0;
}

void
fastinterpolate::RemoveDuplicates()
{
    long i = 0;
    while (i < mSize - 1) {
        if ( abs(mX[i+1] - mX[i]) < 1e-6 ) {
            // Overwrite lower value, unless it is the first one
            for (long j = i+(i==0); j < mSize - 1; j++) {
                mX[j] = mX[j+1];
                mY[j] = mY[j+1];
            }
            mSize--;
        } else {
            i++;
        }
    }
    
}

void
fastinterpolate::insertpoints(double* inX, double *inY, long inN)
{
    for (long i = 0; i < inN; i++) {
        insertpoint(inX[i],inY[i]);
    }

    CheckEquallySpaced();    

}

void
fastinterpolate::insertpoint(double inX, double inY)
{
    long i = 0;
    if (mSize+1 > mCapacity) {
        mCapacity *= 2;
        double* newX = new double[mCapacity];
        double* newY = new double[mCapacity];

        while (i < mSize) {
            newX[i] = mX[i];
            newY[i] = mY[i];
            ++i;
        }

        delete mX;
        delete mY;

        mX = newX;
        mY = newY;
    }

    i = 0;
    if (mSize > 0) {
        for (i = mSize-1; i >= 0; i--) {
            if (mX[i] == inX) {
  //              clog << "Warning: duplicate X at point " << i << " with (" << mX[i] << ", " << mY[i] << ") and " << "(" << inX << ", " << inY << ")\n";
                break;
            } else if (mX[i] < inX) {
                // Insert after i (i.e. i+1)
                for (long j = mSize - 1; j >= i + 1; j--) {
                    mX[j+1] = mX[j];
                    mY[j+1] = mY[j];
                }
                i++;
                mSize++;
                break;
            }
    
        }
    } else {
        mSize++;
    }
    
    mX[i] = inX;
    mY[i] = inY;

    cYMaximum = NAN;
    cXMaximum = NAN;
}

double
fastinterpolate::y(const double& inX)
{

    if (mSize == 1 && areTheSame(mX[0], inX,7) )
        return mY[0];
    
    if (mStartsAtZeroAndEquallySpaced && mApproximate) {
        size_t lo = size_t(inX/mDX);

        if (lo < mSize) {
            return  mY[lo];
        } else {
            cerr << "Error in fastinterpolate: out of range " << inX <<
            " not between " << mX[0] << " and " << mX[mSize-1] << endl;
            throw runtime_error("out of range value in fastinterpolate.");
        }
    } else if (mStartsAtZeroAndEquallySpaced) {
        size_t lo = size_t(inX/mDX);

        if (lo < mSize) {
            return  mY[lo]+(inX - mX[lo])*(mY[lo+1]-mY[lo])/mDX;
        } else {
            cerr << "Error in fastinterpolate: out of range " << inX <<
            " not between " << mX[0] << " and " << mX[mSize-1] << endl;
            throw runtime_error("out of range value in fastinterpolate.");
        }
    } else if (mXisEquallySpaced) {
        size_t lo = size_t((inX-*mX)/mDX);

        if (lo < mSize) {
            return  mY[lo]+(inX - mX[lo])*(mY[lo+1]-mY[lo])/mDX;
        } else {
            cerr << "Error in fastinterpolate: out of range " << inX <<
            " not between " << mX[0] << " and " << mX[mSize-1] << endl;
            throw runtime_error("out of range value in fastinterpolate.");
        }
        
    } else {

        size_t xsize,xsize_1;

        xsize_1 = xsize = mSize;
        xsize_1--;

        if (inX > mX[xsize_1] || inX < mX[0] ) {
            cerr << "Error in fastinterpolate: out of range " << inX <<
            " not between " << mX[0] << " and " << mX[xsize_1] << endl;
            throw runtime_error("out of range value in fastinterpolate.");
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
        b = (inX - xlo);

        return mY[lo]+b*(mY[hi]-mY[lo])/h;

    }

}

valarray<double>
fastinterpolate::y(valarray<double>& inX)
{
    valarray<double> inY(inX.size());

    for (long i = 0; i < inY.size(); i++) {
        inY[i] = y(inX[i]);
    }

    return inY;
}

void
fastinterpolate::NormalizeX(double inNorm)
{
    for (long i = 0; i < mSize; i++) {
        mX[i] /= inNorm;
    }

    RemoveDuplicates();
}

void
fastinterpolate::NormalizeY(double inNorm)
{
    for (long i = 0; i < mSize; i++) {
        mY[i] /= inNorm;
    }

    cYMaximum /= inNorm;
}

double *
fastinterpolate::xArray()
{
    return mX;
}

double * 
fastinterpolate::yArray()
{
    return mY;
}

double
fastinterpolate::yMaximum()
{
    cYMaximum = 0;
    for (long i = 0; i < mSize; i++) {
        cYMaximum = mY[i] > cYMaximum? mY[i] : cYMaximum;
    }
    
    return cYMaximum;

}

double
fastinterpolate::xMaximum()
{
    cXMaximum = 0;
    for (long i = 0; i < mSize; i++) {
        cXMaximum = mX[i] > cXMaximum? mX[i] : cXMaximum;
    }
    
    return cXMaximum;
    
}

