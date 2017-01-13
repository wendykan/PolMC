#include "MCEllipsoid.h"

#include "configfiles.h"
#include "MCObject.h"
#include "Photon.h"
#include "StokesV.h"
#include "MuellerM.h"
#include "MCRandomScatterer.h"
#include "MCUtils.h"

MCEllipsoid::MCEllipsoid(double inWidth,
                         double inHeight,
                         double inDepth,
                         long inPtsWidth,
                         long inPtsHeight,
                         double inIndexMedium,
                         double inIndexOutside,
                         double inRotationPerCmInClearSpace,
						 long inDetectionType,
                         long inAcceptanceCosineElements,
                         bool inUseNextEventEstimator)
:MCObject(inIndexMedium, inIndexOutside, inRotationPerCmInClearSpace, inUseNextEventEstimator)
{
    init(inWidth, inHeight, inDepth, inPtsWidth, inPtsHeight, inDetectionType, inAcceptanceCosineElements);
}

void
MCEllipsoid::init(double inWidth,
                  double inHeight,
                  double inDepth,
                  long inPtsWidth,
                  long inPtsHeight,
				  long inDetectionType,
				  long inAcceptanceCosineElements)
{
    mObjectName="ellipsoid";
    mWidth = inWidth;
    mHeight = inHeight;
    mDepth = inDepth;

    /* Three corners of triangle */
    RealV a,b,origin;

    double z = mDepth/2. * cos(PI/4);

    long interfaceIndex = 0;

    mNumOfInterfaces = 36;

    for (long i = 0; i < mNumOfInterfaces; i++) {
        // We have to allocate individually because
        // new[] must be deleted with delete[], and we can't do that
        // because we give ownership to MCWorld (later on)
        SurfaceElement* se = new SurfaceElement;
        mSurfaceElements.push_back(se);
    }
    
    for (long i  = 0; i < 6; i++) {
        a = RealV(mWidth/2. * cos(2*PI/6*double(i))*sin(PI/4), mHeight/2. * sin(2*PI/6*double(i))*sin(PI/4), z);
        b = RealV(mWidth/2. * cos(2*PI/6*double(i+1))*sin(PI/4), mHeight/2. * sin(2*PI/6*double(i+1))*sin(PI/4), z);
        origin = RealV(0,0, mDepth/2.);

        mSurfaceElements[interfaceIndex]->uniqueID = (interfaceIndex);
        mSurfaceElements[interfaceIndex]->origin = origin;
        mSurfaceElements[interfaceIndex]->a = (a-origin);
        mSurfaceElements[interfaceIndex]->b = (b-origin);
        mSurfaceElements[interfaceIndex]->Na = inPtsWidth;
        mSurfaceElements[interfaceIndex]->Nb = inPtsHeight;
        mSurfaceElements[interfaceIndex]->el = RealV(1,0,0);
        mSurfaceElements[interfaceIndex]->er = RealV(0,1,0);
        mSurfaceElements[interfaceIndex]->normal = RealV::CrossProduct(mSurfaceElements[interfaceIndex]->a,mSurfaceElements[interfaceIndex]->b);
        mSurfaceElements[interfaceIndex]->normal.normalize();
        mSurfaceElements[interfaceIndex]->surfaceShape = kTriangle;
		mSurfaceElements[interfaceIndex]->mDetection = inDetectionType;
		mSurfaceElements[interfaceIndex]->mCosineElements = inAcceptanceCosineElements;
		mSurfaceElements[interfaceIndex]->init();
		
        interfaceIndex++;
    }


    for (long i  = 0; i < 6; i++) {
        a = RealV(mWidth/2. * cos(2*PI/6*double(i+1))*sin(PI/4), mHeight/2. * sin(2*PI/6*double(i+1))*sin(PI/4), z);
        b = RealV(mWidth/2. * cos(2*PI/6*double(i))*sin(PI/4), mHeight/2. * sin(2*PI/6*double(i))*sin(PI/4), z);
        origin = RealV(mWidth/2. * cos(2*PI/6*double(i+0.5)), mHeight/2. * sin(2*PI/6*double(i+0.5)), 0);

        mSurfaceElements[interfaceIndex]->uniqueID=(interfaceIndex);
        mSurfaceElements[interfaceIndex]->origin = origin;
        mSurfaceElements[interfaceIndex]->a = (a-origin);
        mSurfaceElements[interfaceIndex]->b = (b-origin);
        mSurfaceElements[interfaceIndex]->Na = inPtsWidth;
        mSurfaceElements[interfaceIndex]->Nb = inPtsHeight;
        mSurfaceElements[interfaceIndex]->el = RealV(1,0,0);
        mSurfaceElements[interfaceIndex]->er = RealV(0,1,0);
        mSurfaceElements[interfaceIndex]->normal = RealV::CrossProduct(mSurfaceElements[interfaceIndex]->a,mSurfaceElements[interfaceIndex]->b);
        mSurfaceElements[interfaceIndex]->normal.normalize();
        mSurfaceElements[interfaceIndex]->surfaceShape = kTriangle;
		mSurfaceElements[interfaceIndex]->mDetection = inDetectionType;
		mSurfaceElements[interfaceIndex]->mCosineElements = inAcceptanceCosineElements;
		mSurfaceElements[interfaceIndex]->init();
		
        interfaceIndex++;

        b = a;
        a = RealV(mWidth/2. * cos(2*PI/6*double(i+1.5)), mHeight/2. * sin(2*PI/6*double(i+1.5)), 0);

        mSurfaceElements[interfaceIndex]->uniqueID=(interfaceIndex);
        mSurfaceElements[interfaceIndex]->origin = origin;
        mSurfaceElements[interfaceIndex]->a = (a-origin);
        mSurfaceElements[interfaceIndex]->b = (b-origin);
        mSurfaceElements[interfaceIndex]->Na = inPtsWidth;
        mSurfaceElements[interfaceIndex]->Nb = inPtsHeight;
        mSurfaceElements[interfaceIndex]->el = RealV(1,0,0);
        mSurfaceElements[interfaceIndex]->er = RealV(0,1,0);
        mSurfaceElements[interfaceIndex]->normal = RealV::CrossProduct(mSurfaceElements[interfaceIndex]->a,mSurfaceElements[interfaceIndex]->b);
        mSurfaceElements[interfaceIndex]->normal.normalize();
        mSurfaceElements[interfaceIndex]->surfaceShape = kTriangle;
		mSurfaceElements[interfaceIndex]->mDetection = inDetectionType;
		mSurfaceElements[interfaceIndex]->mCosineElements = inAcceptanceCosineElements;
		mSurfaceElements[interfaceIndex]->init();
		
        interfaceIndex++;

    }

    for (long i  = 0; i < 6; i++) {
        a = RealV(mWidth/2. * cos(2*PI/6*double(i))*sin(PI/4), mHeight/2. * sin(2*PI/6*double(i))*sin(PI/4), -z);
        b = RealV(mWidth/2. * cos(2*PI/6*double(i+1))*sin(PI/4), mHeight/2. * sin(2*PI/6*double(i+1))*sin(PI/4), -z);
        origin = RealV(mWidth/2. * cos(2*PI/6*double(i+0.5)), mHeight/2. * sin(2*PI/6*double(i+0.5)), 0);

        mSurfaceElements[interfaceIndex]->uniqueID=(interfaceIndex);
        mSurfaceElements[interfaceIndex]->origin = origin;
        mSurfaceElements[interfaceIndex]->a = (a-origin);
        mSurfaceElements[interfaceIndex]->b = (b-origin);
        mSurfaceElements[interfaceIndex]->Na = inPtsWidth;
        mSurfaceElements[interfaceIndex]->Nb = inPtsHeight;
        mSurfaceElements[interfaceIndex]->el = RealV(1,0,0);
        mSurfaceElements[interfaceIndex]->er = RealV(0,1,0);
        mSurfaceElements[interfaceIndex]->normal = RealV::CrossProduct(mSurfaceElements[interfaceIndex]->a,mSurfaceElements[interfaceIndex]->b);
        mSurfaceElements[interfaceIndex]->normal.normalize();
        mSurfaceElements[interfaceIndex]->surfaceShape = kTriangle;
		mSurfaceElements[interfaceIndex]->mDetection = inDetectionType;
		mSurfaceElements[interfaceIndex]->mCosineElements = inAcceptanceCosineElements;
		mSurfaceElements[interfaceIndex]->init();
		
        interfaceIndex++;

        a = b;
        b = RealV(mWidth/2. * cos(2*PI/6*double(i+1.5)), mHeight/2. * sin(2*PI/6*double(i+1.5)), 0);

        mSurfaceElements[interfaceIndex]->uniqueID=(interfaceIndex);
        mSurfaceElements[interfaceIndex]->origin = origin;
        mSurfaceElements[interfaceIndex]->a = (a-origin);
        mSurfaceElements[interfaceIndex]->b = (b-origin);
        mSurfaceElements[interfaceIndex]->Na = inPtsWidth;
        mSurfaceElements[interfaceIndex]->Nb = inPtsHeight;
        mSurfaceElements[interfaceIndex]->el = RealV(1,0,0);
        mSurfaceElements[interfaceIndex]->er = RealV(0,1,0);
        mSurfaceElements[interfaceIndex]->normal = RealV::CrossProduct(mSurfaceElements[interfaceIndex]->a,mSurfaceElements[interfaceIndex]->b);
        mSurfaceElements[interfaceIndex]->normal.normalize();
        mSurfaceElements[interfaceIndex]->surfaceShape = kTriangle;
		mSurfaceElements[interfaceIndex]->mDetection = inDetectionType;
		mSurfaceElements[interfaceIndex]->mCosineElements = inAcceptanceCosineElements;
		mSurfaceElements[interfaceIndex]->init();
		
        interfaceIndex++;

    }

    for (long i  = 0; i < 6; i++) {
        a = RealV(mWidth/2. * cos(2*PI/6*double(i+1))*sin(PI/4), mHeight/2. * sin(2*PI/6*double(i+1))*sin(PI/4), -z);
        b = RealV(mWidth/2. * cos(2*PI/6*double(i))*sin(PI/4), mHeight/2. * sin(2*PI/6*double(i))*sin(PI/4), -z);

        origin = RealV(0,0,-mDepth/2.);

        mSurfaceElements[interfaceIndex]->uniqueID=(interfaceIndex);
        mSurfaceElements[interfaceIndex]->origin =origin;
        mSurfaceElements[interfaceIndex]->a = (a-origin);
        mSurfaceElements[interfaceIndex]->b = (b-origin);
        mSurfaceElements[interfaceIndex]->Na = inPtsWidth;
        mSurfaceElements[interfaceIndex]->Nb = inPtsHeight;
        mSurfaceElements[interfaceIndex]->el = RealV(1,0,0);
        mSurfaceElements[interfaceIndex]->er = RealV(0,1,0);
        mSurfaceElements[interfaceIndex]->normal = RealV::CrossProduct(mSurfaceElements[interfaceIndex]->a,mSurfaceElements[interfaceIndex]->b);
        mSurfaceElements[interfaceIndex]->normal.normalize();
        mSurfaceElements[interfaceIndex]->surfaceShape = kTriangle;
		mSurfaceElements[interfaceIndex]->mDetection = inDetectionType;
		mSurfaceElements[interfaceIndex]->mCosineElements = inAcceptanceCosineElements;
		mSurfaceElements[interfaceIndex]->init();
		
        interfaceIndex++;
    }
    
    double minNormToOrigin = 1;
    for (long i  = 0; i <= mNumOfInterfaces - 1; i++) {
        minNormToOrigin = RealV(mSurfaceElements[i]->origin + mSurfaceElements[i]->a/2).norm();

        minNormForComplicatedCalculation = minNormToOrigin < minNormForComplicatedCalculation ? minNormToOrigin : minNormForComplicatedCalculation;
    }    
    
    FinishCreate();
    
}

void 
MCEllipsoid::init(map<string,string>& inDict)
{
    dict = inDict;

    long Nx, Ny, acceptanceCosElements;
    double xwidth, yheight, zdepth;
    
    InitVariable(dict, "Nx",Nx);
    InitVariable(dict, "Ny",Ny);
    InitVariable(dict, "zdepth",zdepth);
    InitVariable(dict, "xwidth",xwidth);
    InitVariable(dict, "yheight",yheight);
    InitVariable(dict, "acceptanceCosElements",acceptanceCosElements);

    init(xwidth, yheight, zdepth, Nx, Ny, kFromInside, acceptanceCosElements);

}

bool
MCEllipsoid::IsOutsideObject(RealV& inLocalPosition)
{
 
    if (inLocalPosition.norm() <= minNormForComplicatedCalculation) {
        return false;
    } else {
        for (long i  = 0; i < mSurfaceElements.size(); i++) {
            SurfaceElement* se = mSurfaceElements[i];
            // http://www.thepolygoners.com/tutorials/lineplane/lineplane.html
            RealV& normal = se->normal;

            double nDotd = RealV::DotProduct(normal, inLocalPosition);

            if ( nDotd <= 0 )
                continue;

            RealV vecPosSrc = RealV(0,0,0) - se->origin;

            double nDotPos = RealV::DotProduct(normal, vecPosSrc);

            // If t is between 0 and 1, we are outside
            double t = - nDotPos / nDotd;

            if (t >= 0 && t <= 1.) {
                  return true; 
            }

        }

        return false;
    }
    
}

RealV    
MCEllipsoid::GetRandomPointInsideObjectUniformDistribution()
{

    RealV localPt;
    
    do {
        double x = 2 * (mRandomScatterer->RandomFloat() - 0.5);
        double y = 2 * (mRandomScatterer->RandomFloat() - 0.5);
        double z = 2 * (mRandomScatterer->RandomFloat() - 0.5);
        
        localPt = RealV(x, y, z);
        } while (IsOutsideObject(localPt));
    
    return mGlobalOrigin + localPt;
    
}

