#include "MCCylinder.h"

#include "configfiles.h"
#include "MCObject.h"
#include "Photon.h"
#include "StokesV.h"
#include "MuellerM.h"
#include "MCRandomScatterer.h"
#include "MCUtils.h"

MCCylinder::MCCylinder(double inRadius,
					   double inHeight,
					   long inPtsRadius,
					   long inPtsHeight,
					   double inIndexMedium,
					   double inIndexOutside,
					   double inRotationPerCmInClearSpace,
					   long inDetectionType,
					   long inAcceptanceCosineElements,
					   bool inUseNextEventEstimator)
:MCObject(inIndexMedium, inIndexOutside, inRotationPerCmInClearSpace, inUseNextEventEstimator)
{
    init(inRadius, inHeight, inPtsRadius, inPtsHeight, inDetectionType, inAcceptanceCosineElements);
}

void
MCCylinder::init(double inRadius,
          double inHeight,
          long inPtsRadius,
          long inPtsHeight,
		  long inDetectionType, 
		  long inAcceptanceCosineElements)
{
    mObjectName="cylinder";
    mRadius = inRadius;
    mHeight = inHeight;

    /* Three corners of triangle */
    RealV a,b,origin;

    double z = 0.5 * mHeight;

    long interfaceIndex = 0;

	inPtsRadius = inPtsRadius < 6 ? 6: inPtsRadius;
	inPtsHeight = inPtsHeight < 1 ? 1: inPtsHeight;
	
	mNumOfInterfaces = inPtsRadius + inPtsRadius + inPtsRadius * inPtsHeight;

    for (long i = 0; i < mNumOfInterfaces; i++) {
        // We have to allocate individually because
        // new[] must be deleted with delete[], and we can't do that
        // because we give ownership to MCWorld (later on)
        SurfaceElement* se = new SurfaceElement;
        mSurfaceElements.push_back(se);
    }
    
	long Nsides = inPtsRadius;
	
    for (long i  = 0; i < Nsides; i++) {
        a = RealV(mRadius * cos(2*PI/Nsides*double(i)), mRadius * sin(2*PI/Nsides*double(i)), z);
        b = RealV(mRadius * cos(2*PI/Nsides*double(i+1)), mRadius * sin(2*PI/Nsides*double(i+1)), z);
        origin = RealV(0,0,z);

        mSurfaceElements[interfaceIndex]->uniqueID = (interfaceIndex);
        mSurfaceElements[interfaceIndex]->origin = origin;
        mSurfaceElements[interfaceIndex]->a = (a-origin);
        mSurfaceElements[interfaceIndex]->b = (b-origin);
        mSurfaceElements[interfaceIndex]->Na = inPtsRadius;
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

    for (long i  = 0; i < Nsides; i++) {
        // Faces
        origin = RealV(mRadius * cos(2*PI/Nsides*double(i)), mRadius * sin(2*PI/Nsides*double(i)), z);

        a = RealV(mRadius * cos(2*PI/Nsides*double(i)), mRadius * sin(2*PI/Nsides*double(i)), -z);

        b = RealV(mRadius * cos(2*PI/Nsides*double(i+1)), mRadius * sin(2*PI/Nsides*double(i+1)), z);

        mSurfaceElements[interfaceIndex]->uniqueID = (interfaceIndex);
        mSurfaceElements[interfaceIndex]->origin = origin;
        mSurfaceElements[interfaceIndex]->a = (a-origin);
        mSurfaceElements[interfaceIndex]->b = (b-origin);
        mSurfaceElements[interfaceIndex]->Na = 1;
        mSurfaceElements[interfaceIndex]->Nb = 1;
        mSurfaceElements[interfaceIndex]->el = RealV(1,0,0);
        mSurfaceElements[interfaceIndex]->er = RealV(0,1,0);
        mSurfaceElements[interfaceIndex]->normal = RealV::CrossProduct(mSurfaceElements[interfaceIndex]->a,mSurfaceElements[interfaceIndex]->b);
        mSurfaceElements[interfaceIndex]->normal.normalize();
        mSurfaceElements[interfaceIndex]->surfaceShape = kParallelogram;
		mSurfaceElements[interfaceIndex]->mDetection = inDetectionType;
		mSurfaceElements[interfaceIndex]->mCosineElements = inAcceptanceCosineElements;
		mSurfaceElements[interfaceIndex]->init();
		
        interfaceIndex++;
    }



    for (long i  = 0; i < Nsides; i++) {
        /* Opposite of top face */
        a = RealV(mRadius * cos(2*PI/Nsides*double(i+1)), mRadius * sin(2*PI/Nsides*double(i+1)), -z);
        b = RealV(mRadius * cos(2*PI/Nsides*double(i)), mRadius * sin(2*PI/Nsides*double(i)), -z);
        origin = RealV(0,0,-z);

        mSurfaceElements[interfaceIndex]->uniqueID = (interfaceIndex);
        mSurfaceElements[interfaceIndex]->origin = origin;
        mSurfaceElements[interfaceIndex]->a = (a-origin);
        mSurfaceElements[interfaceIndex]->b = (b-origin);
        mSurfaceElements[interfaceIndex]->Na = 1;
        mSurfaceElements[interfaceIndex]->Nb = 1;
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
    
    FinishCreate();

}

void 
MCCylinder::init(map<string,string>& inDict)
{
    dict = inDict;

    long Nr, Nh, acceptanceCosElements;
    double radius, height;
    
    InitVariable(dict, "Nr",Nr);
    InitVariable(dict, "Nh",Nh);
    InitVariable(dict, "radius",radius);
    InitVariable(dict, "height",height);
    InitVariable(dict, "acceptanceCosElements",acceptanceCosElements);

    init(radius, height, Nr, Nh, kFromInside, acceptanceCosElements);

}

bool
MCCylinder::IsOutsideObject(RealV& inLocalPosition)
{
 
//    throw runtime_error("IsOutsideObject Not implemented in MCCylinder");
	return false;

}

